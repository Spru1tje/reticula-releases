"""Instrument file readers.

Every reader returns pandas DataFrames with clean float columns. Native binary
formats are parsed with established open-source libraries where they exist
(nmrglue for Bruker TopSpin, brukeropusreader for OPUS, tifffile for SEM TIFF);
everything else goes through a robust text parser that copes with comma or
semicolon delimiters, decimal commas (Belgian/European locale), metadata
headers and footers.
"""
from __future__ import annotations
import io, os, re, zipfile, tempfile
import numpy as np
import pandas as pd

NUM = r"[-+]?(?:\d+[.,]?\d*|[.,]\d+)(?:[eE][-+]?\d+)?"


# --------------------------------------------------------------------------- utils
def _to_text(src) -> str:
    if hasattr(src, "read"):
        data = src.read()
        if hasattr(src, "seek"):
            try: src.seek(0)
            except Exception: pass
    elif isinstance(src, (bytes, bytearray)):
        data = bytes(src)
    else:
        with open(src, "rb") as fh:
            data = fh.read()
    if isinstance(data, str):
        return data
    encs = ["utf-8-sig", "cp1252", "latin-1"]
    if data[:2] in (b"\xff\xfe", b"\xfe\xff") or data[:400].count(b"\x00") > 50:
        encs.insert(0, "utf-16")
    for enc in encs:
        try:
            return data.decode(enc)
        except Exception:
            continue
    return data.decode("latin-1", "ignore")


def _to_bytes(src) -> bytes:
    if hasattr(src, "read"):
        b = src.read()
        try: src.seek(0)
        except Exception: pass
        return b if isinstance(b, bytes) else b.encode()
    if isinstance(src, (bytes, bytearray)):
        return bytes(src)
    with open(src, "rb") as fh:
        return fh.read()


def _split(line: str, delim):
    if delim == "ws":
        return line.split()
    return [c.strip() for c in line.split(delim)]


def _guess_delim(lines):
    sample = [l for l in lines if l.strip()][-150:]
    if not sample:
        return "ws"
    scores = {}
    for d in ["\t", ";", ","]:
        counts = [l.count(d) for l in sample]
        vals, freq = np.unique(counts, return_counts=True)
        mode = vals[np.argmax(freq)]
        cons = freq.max() / len(counts)
        if mode > 0 and cons > 0.7:
            scores[d] = cons
    for d in ["\t", ";"]:
        if d in scores and (scores[d] >= scores.get(",", 0) - 0.05 or "," not in scores):
            return d
    if "," in scores:
        return ","
    return "ws"


def _num(cell: str):
    c = cell.strip().strip('"').replace("\u2212", "-")
    if not c:
        return np.nan
    try:
        return float(c)
    except ValueError:
        pass
    try:
        return float(c.replace(",", "."))
    except ValueError:
        return None


def read_table(src, min_numeric_cols=2) -> pd.DataFrame:
    """Generic numeric table reader. Finds the largest contiguous numeric block,
    uses the nearest non-numeric line above it as header."""
    text = _to_text(src)
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    delim = _guess_delim(lines)
    # if delimiter is comma and the data rows contain decimal commas inside quotes, pandas handles quotes poorly -> custom
    rows, numeric = [], []
    for l in lines:
        cells = _split(l, delim) if l.strip() else []
        vals = [_num(c) for c in cells]
        ok = len(cells) >= min_numeric_cols and sum(v is not None and not (isinstance(v, float) and np.isnan(v)) for v in vals) >= min_numeric_cols \
            and all(v is not None for v in vals[:min_numeric_cols])
        rows.append((l, cells, vals))
        numeric.append(ok)
    # longest run of numeric rows
    best, cur, best_start, start = 0, 0, 0, 0
    for i, ok in enumerate(numeric):
        if ok:
            if cur == 0:
                start = i
            cur += 1
            if cur > best:
                best, best_start = cur, start
        else:
            cur = 0
    if best == 0:
        raise ValueError("No numeric data block found")
    block = rows[best_start: best_start + best]
    ncol = int(np.median([len(r[1]) for r in block]))
    header = None
    for j in range(best_start - 1, max(-1, best_start - 15), -1):
        cells = rows[j][1]
        if len(cells) >= min(ncol, 2) and any(re.search("[A-Za-z°%]", c) for c in cells):
            header = [c.strip().strip('"') for c in cells][:ncol]
            # units-only line (e.g. Mettler STARe '[s] [°C] [mg]') -> combine with the names line above
            if all(re.fullmatch(r"\[.*\]|\(.*\)", h) for h in header if h) and j > 0:
                above = [c.strip() for c in rows[j - 1][1]]
                if len(above) >= len(header):
                    k = len(above) - len(header)
                    header = above[:k] + [f"{a} {u}" for a, u in zip(above[k:], header)]
            break
    data = []
    for _, cells, vals in block:
        v = [(x if x is not None else np.nan) for x in vals[:ncol]]
        v += [np.nan] * (ncol - len(v))
        data.append(v)
    if not header or len(header) < ncol:
        header = (header or []) + [f"col{i+1}" for i in range(len(header or []), ncol)]
    # de-duplicate
    seen = {}
    hdr = []
    for h in header:
        h = h or "col"
        if h in seen:
            seen[h] += 1; h = f"{h}_{seen[h]}"
        else:
            seen[h] = 0
        hdr.append(h)
    df = pd.DataFrame(data, columns=hdr).dropna(axis=1, how="all")
    df.attrs["header_text"] = "\n".join(r[0] for r in rows[:best_start])[-4000:]
    return df


def read_csv_table(src, keywords=(r"lane", r"band", r"sample")) -> pd.DataFrame:
    """Delimited table with text columns: header row = first line matching one of the keywords."""
    raw = _to_bytes(src)
    if raw[:4] in (b"\xd0\xcf\x11\xe0", b"PK\x03\x04"):
        return pd.read_excel(io.BytesIO(raw))
    text = _to_text(raw)
    lines = text.replace("\r", "").split("\n")
    d = _guess_delim(lines)
    sep = {"\t": "\t", ";": ";", ",": ",", "ws": r"\s+"}[d]
    hi = 0
    for i, l in enumerate(lines):
        if sum(bool(re.search(k, l, re.I)) for k in keywords) >= min(2, len(keywords)) and l.count(sep if sep != r"\s+" else " ") >= 1:
            hi = i; break
    df = pd.read_csv(io.StringIO("\n".join(lines[hi:])), sep=sep, engine="python", comment=None)
    df.columns = [str(c).strip() for c in df.columns]
    if sep == ";":
        for c in df.columns:
            if df[c].dtype == object:
                conv = pd.to_numeric(df[c].astype(str).str.replace(",", "."), errors="coerce")
                if conv.notna().mean() > 0.8:
                    df[c] = conv
    return df


def find_col(df: pd.DataFrame, *patterns, default=None):
    for p in patterns:
        for c in df.columns:
            if re.search(p, str(c), re.I):
                return c
    return default


def read_xy(src):
    """Two-column spectrum/pattern (x, y[, esd]). Returns DataFrame x, y (, esd)."""
    df = read_table(src)
    out = pd.DataFrame({"x": df.iloc[:, 0].astype(float), "y": df.iloc[:, 1].astype(float)})
    if df.shape[1] >= 3:
        out["esd"] = df.iloc[:, 2].astype(float)
    return out.dropna(subset=["x", "y"]).reset_index(drop=True)


# --------------------------------------------------------------------------- PXRD
def read_pxrd(src, name: str = ""):
    """STOE WinXPOW exports (.xy/.xye/.dat/.txt, 2theta-intensity[-esd]),
    Bruker .xy, PANalytical .xy/.asc and any two-column text. STOE .raw (binary)
    is read only if the optional `xylib` package is installed."""
    ext = os.path.splitext(name or getattr(src, "name", ""))[1].lower()
    if ext in (".rasx", ".ras"):
        from mfs import rigaku as RG
        scans = RG.read_any(_to_bytes(src), name)
        pick = next((s_ for s_ in scans if s_["is_2theta_theta"]), scans[0])
        df = pd.DataFrame({"x": pick["x"], "y": pick["y"]})
        df.attrs["rigaku"] = pick
        df.attrs["rigaku_all"] = scans
        return df
    if ext == ".raw":
        try:
            import xylib  # optional
            with tempfile.NamedTemporaryFile(suffix=".raw", delete=False) as t:
                t.write(_to_bytes(src)); p = t.name
            d = xylib.load_file(p, "")
            b = d.get_block(0)
            x = np.array([b.get_column(1).get_value(i) for i in range(b.get_point_count())])
            y = np.array([b.get_column(2).get_value(i) for i in range(b.get_point_count())])
            return pd.DataFrame({"x": x, "y": y})
        except ImportError:
            raise ValueError("STOE .raw is a binary format: export it from WinXPOW as .xy/.dat "
                             "(File > Export / 'Save as ASCII'), or install the optional xylib package.")
    df = read_xy(src)
    return df


# --------------------------------------------------------------------------- FTIR / Raman
def read_opus(src):
    """Bruker OPUS binary (.0, .1, ...) via brukeropusreader. Returns x (cm-1), y and the block used."""
    from brukeropusreader import read_file
    with tempfile.NamedTemporaryFile(suffix=".0", delete=False) as t:
        t.write(_to_bytes(src)); p = t.name
    d = read_file(p)
    for block in ["AB", "ScSm", "ScRf", "TR", "RAM", "ScSmRaman"]:
        if block in d:
            y = np.asarray(d[block], float)
            par = d.get(f"{block} Data Parameter", {})
            fxv, lxv, npt = par.get("FXV"), par.get("LXV"), par.get("NPT", len(y))
            y = y[:npt]
            x = np.linspace(fxv, lxv, len(y)) if fxv is not None else np.arange(len(y), dtype=float)
            out = pd.DataFrame({"x": x, "y": y})
            out.attrs["block"] = block
            return out
    raise ValueError("No absorbance/transmittance/Raman block found in OPUS file")


def read_ftir(src, name: str = ""):
    ext = os.path.splitext(name or getattr(src, "name", ""))[1].lower()
    if re.fullmatch(r"\.\d+", ext or ""):
        return read_opus(src)
    raw = _to_bytes(src)
    if raw[:4] == b"\n\n\xfe\xfe" or (b"\x00" in raw[:64] and not raw[:2] in (b"\xff\xfe", b"\xfe\xff")):
        return read_opus(io.BytesIO(raw))
    return read_xy(io.BytesIO(raw))  # OPUS 'Save as Data Point Table' (.dpt) / CSV


# --------------------------------------------------------------------------- NMR
def read_nmr(src, name: str = ""):
    """Bruker TopSpin processed spectrum.
    Accepts (1) a .zip of the experiment folder (contains acqus + pdata/1/1r, procs),
    (2) a TopSpin 'totxt' export (# LEFT = .. ppm. RIGHT = .. ppm. header + intensities),
    (3) any two-column ppm/intensity text."""
    raw = _to_bytes(src)
    if raw[:2] == b"PK":
        import nmrglue as ng
        tmp = tempfile.mkdtemp()
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            z.extractall(tmp)
        pdirs = []
        for root, dirs, files in os.walk(tmp):
            if "1r" in files and "procs" in files:
                pdirs.append(root)
        if not pdirs:
            raise ValueError("Zip does not contain pdata/<n>/1r + procs (process the spectrum in TopSpin first)")
        pd_dir = sorted(pdirs)[0]
        dic, data = ng.bruker.read_pdata(pd_dir, scale_data=True)
        pr = dic["procs"]
        n = len(data)
        # TopSpin convention: OFFSET = ppm of the first (left) point, SW_p in Hz, SF in MHz
        ppm = pr["OFFSET"] - np.arange(n) * (pr["SW_p"] / pr["SF"]) / n
        out = pd.DataFrame({"x": ppm, "y": np.asarray(data, float)})
        out.attrs["source"] = "TopSpin pdata " + os.path.relpath(pd_dir, tmp)
        try:
            out.attrs["SF"] = dic["procs"]["SF"]
            out.attrs["NS"] = dic.get("acqus", {}).get("NS"); out.attrs["D1_s"] = dic.get("acqus", {}).get("D", [0, None])[1]
        except Exception:
            pass
        return out
    text = _to_text(raw)
    m = re.search(r"LEFT\s*=\s*(" + NUM + r")\s*ppm.*?RIGHT\s*=\s*(" + NUM + r")\s*ppm", text, re.S | re.I)
    if m:
        left, right = float(m.group(1).replace(",", ".")), float(m.group(2).replace(",", "."))
        vals = []
        for l in text.splitlines():
            s = l.strip()
            if not s or s.startswith("#"):
                continue
            v = _num(s.split()[0])
            if v is not None and not np.isnan(v):
                vals.append(v)
        y = np.array(vals)
        x = np.linspace(left, right, len(y))
        return pd.DataFrame({"x": x, "y": y})
    return read_xy(io.BytesIO(raw))


# --------------------------------------------------------------------------- N2 physisorption
STP_CM3_PER_MMOL = 22.414
_RE_P = re.compile(r"relative\s*pressure|\bp\s*/\s*p\s*[o0°ₒ]|\bp/p0\b|p_rel", re.I)
_RE_Q = re.compile(r"quantity\s*adsorbed|volume\s*adsorbed|amount\s*adsorbed|uptake|n_ads|\bq\b", re.I)
_RE_W = re.compile(r"pore\s*(width|diameter|size)", re.I)
_RE_KIND = [  # (kind, title regex, header regex)
    ("bet", r"\bbet\b", r"1\s*/\s*\[?\s*q|1/\[w"),
    ("langmuir", r"langmuir", r"p\s*/\s*q\b"),
    ("tplot", r"t-?\s*plot", r"thickness"),
    ("psd", r"dft|nldft|bjh|pore\s*size|horvath|hk\b|pore\s*distribution", r"pore\s*(width|diameter|size)"),
    ("other", r"freundlich|temkin|dubinin|astakhov|\bdr\b|\bda\b|\bmp\b|alpha-?s|f-ratio|heat\s*of", r"log\s*\(|ln\s*\(|w/w0"),
    ("iso_tabular", r"tabular", r"elapsed\s*time|absolute\s*pressure|saturation\s*pressure"),
    ("iso_plot", r"isotherm.*plot|plot.*isotherm|linear\s*plot|log\s*plot", None),
]


def _cells(row):
    return [str(c).strip() if str(c).strip() not in ("nan", "None", "<NA>") else "" for c in row]


def _q_factor(header: str):
    """Factor converting the quantity column to cm3 g-1 STP, plus the unit as written."""
    h = header.lower().replace("³", "3").replace(" ", "")
    if "mmol/g" in h or "mol/kg" in h or "mmolg" in h:
        return STP_CM3_PER_MMOL, "mmol/g"
    if "mg/g" in h:
        return STP_CM3_PER_MMOL / 28.0134, "mg/g"
    if "mol/g" in h:
        return STP_CM3_PER_MMOL * 1000, "mol/g"
    return 1.0, "cm3/g STP"


def _find_tables(grid, sheet=""):
    """Locate every table (header row followed by numeric rows) in a sheet grid.
    A MicroActive full-report export stacks many reports on one sheet; each comes back as its own table."""
    tables, r, n = [], 0, len(grid)
    while r < n:
        cells = _cells(grid[r])
        heads = [i for i, c in enumerate(cells) if c and (_RE_P.search(c) or _RE_Q.search(c) or _RE_W.search(c)
                                                        or re.search(r"pressure|thickness|quantity|volume", c, re.I))
                 and _num(c) is None]
        if len(heads) < 2:
            r += 1; continue
        # need numeric data in the header columns in one of the next 3 rows
        first = None
        for k in range(r + 1, min(n, r + 4)):
            nxt = _cells(grid[k])
            if sum(1 for i in heads if i < len(nxt) and nxt[i] and _num(nxt[i]) is not None) >= 2:
                first = k; break
        if first is None:
            r += 1; continue
        # title and group labels (e.g. 'Adsorption | | Desorption') above the header
        title, groups = "", {}
        for k in range(r - 1, max(-1, r - 12), -1):
            prev = _cells(grid[k])
            if not any(prev):
                if title: break
                continue
            if tables and k <= tables[-1]["end"]:
                break
            texts = [(i, c) for i, c in enumerate(prev) if c]
            if k == r - 1 and all(re.fullmatch(r"(?i)\s*(adsorption|desorption)\s*", c) for _, c in texts):
                groups = {i: c.lower() for i, c in texts}; continue
            if len(texts) == 1 and ":" not in texts[0][1] and _num(texts[0][1]) is None:
                title = texts[0][1]; break
        rows, k, branch_rows = [], first, []
        cur_branch = None
        while k < n:
            cc = _cells(grid[k])
            vals = [_num(cc[i]) if i < len(cc) else np.nan for i in range(len(cells))]
            numeric = sum(1 for i in heads if i < len(vals) and vals[i] is not None and not np.isnan(vals[i]))
            if numeric == 0:
                txt = " ".join(c for c in cc if c).strip().lower()
                if re.fullmatch(r"(adsorption|desorption)( branch)?", txt):
                    cur_branch = txt.split()[0]; k += 1; continue
                break
            rows.append(vals); branch_rows.append(cur_branch); k += 1
        tables.append({"sheet": sheet, "row": r, "end": k, "title": title, "header": cells, "cols": heads,
                       "groups": groups, "data": rows, "row_branch": branch_rows})
        r = k
    for t in tables:
        t["kind"] = _classify(t)
    return tables


def _classify(t):
    title = (t["title"] + " " + t["sheet"]).lower()
    hdr = " ".join(t["header"][i] for i in t["cols"]).lower()
    for kind, tre, hre in _RE_KIND:
        if (tre and re.search(tre, t["title"].lower() or t["sheet"].lower())) or (hre and re.search(hre, hdr)):
            if kind.startswith("iso") and not (_RE_P.search(hdr) and _RE_Q.search(hdr)):
                continue
            return kind
    if _RE_P.search(hdr) and _RE_Q.search(hdr):
        return "isotherm"
    if _RE_W.search(hdr):
        return "psd"
    return "other"


_ISO_PRIORITY = {"iso_tabular": 0, "isotherm": 1, "iso_plot": 2}


def _iso_from_table(t):
    """Return list of (p_rel, q_cm3STP, branch, order) from an isotherm table; keeps measurement order."""
    H = t["header"]
    ip = [i for i in t["cols"] if _RE_P.search(H[i])]
    iq = [i for i in t["cols"] if _RE_Q.search(H[i]) and not _RE_P.search(H[i])]
    pairs = []
    for a in ip:
        qs = [q for q in iq if q > a]
        if qs:
            pairs.append((a, qs[0]))
    if not pairs:
        return [], None
    pts = []
    for k, (a, q) in enumerate(pairs):
        fac, unit = _q_factor(H[q])
        label = t["groups"].get(a) or t["groups"].get(q)
        branch_fixed = label if label in ("adsorption", "desorption") else (("adsorption" if k == 0 else "desorption") if len(pairs) == 2 else None)
        seq = [(row[a], row[q], rb) for row, rb in zip(t["data"], t["row_branch"])
               if a < len(row) and q < len(row) and row[a] is not None and row[q] is not None
               and not np.isnan(row[a]) and not np.isnan(row[q])]
        if not seq:
            continue
        p = np.array([s[0] for s in seq], float); y = np.array([s[1] for s in seq], float) * fac
        if p.max() > 1.5:  # an absolute-pressure column slipped through
            continue
        if branch_fixed:
            br = [branch_fixed[:3]] * len(p)
        elif any(s[2] for s in seq):
            br = [(s[2] or "adsorption")[:3] for s in seq]
        else:
            imax = int(np.argmax(p))   # ads up to the highest pressure, then desorption in measurement order
            br = ["ads"] * (imax + 1) + ["des"] * (len(p) - imax - 1)
        pts += [(float(pp), float(qq), b, i) for i, (pp, qq, b) in enumerate(zip(p, y, br))]
        if len(pairs) > 2 and k == 1:
            break  # several sample overlays side by side: keep the first ads/des pair
    return pts, unit


def _split_bands(grid):
    """Micromeritics TriStar/ASAP 'Full Report Set' exports put every report SIDE BY SIDE on one sheet, in column
    bands separated by columns of '|'. Split such a grid into one grid per band (left to right); other grids are
    returned unchanged. Without this, header/meta lines of different reports are read as one row."""
    if not grid:
        return [grid]
    ncol = max(len(r) for r in grid)
    seps = []
    for c in range(ncol):
        vals = [str(r[c]).strip() for r in grid if c < len(r) and str(r[c]).strip() not in ("", "nan", "None")]
        if len(vals) >= 5 and sum(v == "|" for v in vals) >= 0.6 * len(vals):
            seps.append(c)
    if not seps:
        return [grid]
    edges = [-1] + seps + [ncol]
    bands = []
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a <= 1:
            continue
        band = [[(row[c] if c < len(row) else "") for c in range(a + 1, b)] for row in grid]
        if any(str(v).strip() not in ("", "nan", "None") for row in band for v in row):
            bands.append(band)
    return bands


def _decimal_point(v: str) -> str:
    """'0,0425' / '77,350' (decimal comma, as in Belgian-locale MicroActive exports) -> '0.0425'."""
    return v.replace(",", ".") if re.fullmatch(r"[-+]?\d+,\d+", v) else v


def _sorption_meta(grids):
    meta = {}
    pats = [("sample", r"^\s*Sample:\s*(.+)$"), ("mass_g", r"Sample\s*mass:\s*(" + NUM + r")\s*g"),
            ("bet_reported", r"BET\s+Surface\s+Area:\s*(" + NUM + ")"),
            ("bet_reported_err", r"BET\s+Surface\s+Area:\s*" + NUM + r"\s*±\s*(" + NUM + ")"),
            ("bet_C_reported", r"^\s*C:\s*(" + NUM + ")"), ("bet_Qm_reported", r"^\s*Qm:\s*(" + NUM + ")"),
            ("bet_R_reported", r"^\s*Correlation\s*Coefficient:\s*(" + NUM + ")"),
            ("langmuir_reported", r"Langmuir\s+Surface\s+Area:\s*(" + NUM + ")"),
            ("single_point_reported", r"Single\s+point\s+surface\s+area.*?:\s*(" + NUM + r")\s*m"),
            ("vtotal_reported", r"total\s+pore\s+volume.*?:\s*(" + NUM + r")\s*cm"),
            ("tplot_vmicro_reported", r"t-?Plot\s+micropore\s+volume:\s*(" + NUM + ")"),
            ("tplot_sext_reported", r"t-?Plot\s+External\s+Surface\s+Area:\s*(" + NUM + ")"),
            ("adsorptive", r"Analysis\s+adsorptive:\s*(\S+)"), ("bath_K", r"bath\s+temp\.?:\s*(" + NUM + r")\s*K"),
            ("instrument", r"^\s*(3Flex|ASAP|TriStar|Gemini|MicroActive)[^:]*$")]
    for grid in grids:
        for row in grid:
            cells = [c for c in _cells(row) if c]
            for cell in cells + [" ".join(cells)]:
                for key, pat in pats:
                    if key in meta:
                        continue
                    mm = re.search(pat, cell, re.I)
                    if mm:
                        v = mm.group(1).strip()
                        if key == "sample":
                            v = re.split(r"\s{2,}|\t", v)[0].strip()
                            if not v or re.match(r"(?i)operator|submitter", v):
                                continue
                        meta[key] = _decimal_point(v)
    return meta


def read_sorption(src, name: str = ""):
    """Micromeritics MicroActive / 3Flex / ASAP export (Excel .xls/.xlsx 'Export' of the report, or CSV/TXT).

    v2: table-aware. The file is split into its individual reports (Isotherm Tabular, BET, Langmuir, t-Plot,
    Isotherm Linear/Log Plot overlay data, DFT/BJH...). Only ONE isotherm table is used, in this priority:
    Isotherm Tabular Report > untitled isotherm table > plot-overlay columns. Adsorption and desorption are split
    in measurement order (ads up to the highest p/p0, then des) or by the 'Adsorption'/'Desorption' labels,
    the quantity is converted to cm3 g-1 STP (mmol/g x 22.414), and the BET report's own points/range are kept
    for comparison. Returns dict(iso=DataFrame[p_rel, q, branch, order], psd=DataFrame|None, meta=dict,
    tables=list of {sheet,title,kind,n,used})."""
    ext = os.path.splitext(name or getattr(src, "name", ""))[1].lower()
    raw = _to_bytes(src)
    grids = []
    if ext in (".xls", ".xlsx", ".xlsm") or raw[:4] in (b"\xd0\xcf\x11\xe0", b"PK\x03\x04"):
        engine = "xlrd" if raw[:4] == b"\xd0\xcf\x11\xe0" else None
        sheets = pd.read_excel(io.BytesIO(raw), sheet_name=None, header=None, engine=engine)
        for sn, sh in sheets.items():
            g = sh.astype(object).where(sh.notna(), "").astype(str).values.tolist()
            bands = _split_bands(g)
            for k, band in enumerate(bands):
                grids.append((str(sn) if len(bands) == 1 else f"{sn} [report {k + 1}]", band))
    else:
        text = _to_text(raw)
        lines = text.replace("\r", "").split("\n")
        d = _guess_delim(lines)
        grids.append(("", [_split(l, d) for l in lines]))
    meta = _sorption_meta([g for _, g in grids])
    tables = []
    for sn, g in grids:
        tables += _find_tables(g, sn)
    iso_t = sorted([t for t in tables if t["kind"] in _ISO_PRIORITY], key=lambda t: (_ISO_PRIORITY[t["kind"]], -len(t["data"])))
    pts, unit, used = [], None, None
    for t in iso_t:
        pts, unit = _iso_from_table(t)
        if len(pts) >= 3:
            used = t; break
    if used is None:
        # plain 2-3 column text: p/p0, q  (ads then des in measurement order)
        df = read_table(io.BytesIO(raw))
        p, q = df.iloc[:, 0].values.astype(float), df.iloc[:, 1].values.astype(float)
        imax = int(np.argmax(p))
        pts = [(float(a), float(b), "ads" if i <= imax else "des", i) for i, (a, b) in enumerate(zip(p, q))]
        unit = "cm3/g STP"
    iso = pd.DataFrame(pts, columns=["p_rel", "q", "branch", "order"])
    iso = iso[(iso.p_rel > 0) & (iso.p_rel <= 1.05)].reset_index(drop=True)
    meta["q_unit_in_file"] = unit
    if used is not None:
        meta["isotherm_table"] = f"{used['title'] or 'untitled table'} (sheet '{used['sheet']}', row {used['row'] + 1})"
    # BET report points -> the range Micromeritics used
    bt = next((t for t in tables if t["kind"] == "bet"), None)
    if bt is not None:
        ip = [i for i in bt["cols"] if _RE_P.search(bt["header"][i])]
        if ip:
            pp = [row[ip[0]] for row in bt["data"] if row[ip[0]] is not None and not np.isnan(row[ip[0]])]
            if pp:
                meta["bet_reported_range"] = (float(min(pp)), float(max(pp)))
                meta["bet_reported_npts"] = len(pp)
    # PSD
    psd = None
    for t in tables:
        if t["kind"] != "psd":
            continue
        H = t["header"]
        iw = [i for i in t["cols"] if _RE_W.search(H[i])]
        cols = {"w": iw[0]}
        for i in t["cols"]:
            c = H[i].lower()
            if i == iw[0]: continue
            if "cumul" in c and "volume" in c: cols.setdefault("cum", i)
            elif "dv/dlog" in c.replace(" ", ""): cols.setdefault("dvdlogw", i)
            elif "dv/d" in c.replace(" ", ""): cols.setdefault("dvdw", i)
            elif "increment" in c and "volume" in c: cols.setdefault("inc", i)
        if len(cols) < 2:
            continue
        psd = pd.DataFrame([{k: row[i] if i < len(row) else np.nan for k, i in cols.items()} for row in t["data"]]).astype(float)
        if "dvdw" not in psd and "inc" in psd:
            psd["dvdw"] = psd["inc"] / np.gradient(psd["w"])
        unit_w = H[iw[0]]
        psd.attrs["width_unit"] = "nm" if "nm" in unit_w else "Å"
        psd.attrs["model"] = t["title"]
        break
    meta.setdefault("sample", os.path.splitext(os.path.basename(name or getattr(src, "name", "") or "sample"))[0])
    summary = [{"sheet": t["sheet"], "row": t["row"] + 1, "title": t["title"] or "(untitled)", "kind": t["kind"],
                "points": len(t["data"]), "columns": " | ".join(t["header"][i] for i in t["cols"]),
                "used": t is used or (t["kind"] == "psd" and psd is not None and t["title"] == psd.attrs.get("model"))}
               for t in tables]
    return {"iso": iso, "psd": psd, "meta": meta, "tables": summary}


# --------------------------------------------------------------------------- TGA
def read_tga(src, name: str = ""):
    """TA Instruments (Universal Analysis / TRIOS text export), Netzsch Proteus ASCII
    (##Temp./°C;Time/min;Mass/%), Mettler-Toledo STARe export, Setaram/PerkinElmer Pyris
    text and any CSV with temperature + mass columns.
    Returns DataFrame T (°C), w (% of initial mass), plus t (min) if present."""
    text = _to_text(src)
    # TA Universal Analysis: data after 'StartOfData'
    if "StartOfData" in text:
        text = text.split("StartOfData", 1)[1]
        df = read_table(io.StringIO(text))
        df.columns = [f"c{i}" for i in range(df.shape[1])]
        # UA default signal order: Time (min), Temperature (°C), Weight (mg) [, Weight (%)]
        sig = re.findall(r"Sig\d+\s+(.*)", _to_text(src))
        cols = {}
        for i, s in enumerate(sig[: df.shape[1]]):
            cols[f"c{i}"] = s
        df = df.rename(columns=cols)
    else:
        # Netzsch headers start with '##'
        text = re.sub(r"^##", "", text, flags=re.M)
        df = read_table(io.StringIO(text))
    tcol = find_col(df, r"^ts\b|sample\s*temp|temp", r"°c|deg")
    mcol = find_col(df, r"mass\s*/?\s*%|weight\s*\(%\)|weight\s*%|mass.*%|tg\s*/?\s*%|%")
    mabs = find_col(df, r"weight|mass|value|tg")
    timec = find_col(df, r"time|^t\b|min")
    if tcol is None:
        tcol = df.columns[1] if df.shape[1] > 2 else df.columns[0]
    out = pd.DataFrame({"T": df[tcol].astype(float)})
    if mcol is not None and mcol != tcol:
        w = df[mcol].astype(float)
    else:
        mabs = mabs if mabs not in (None, tcol) else df.columns[-1]
        w = df[mabs].astype(float)
        w = 100 * w / w.iloc[0]
    if w.iloc[: max(3, len(w)//50)].mean() < 5:   # fraction 0..1
        w = 100 * w
    out["w"] = w.values
    if timec is not None and timec not in (tcol, mcol):
        out["t"] = df[timec].astype(float).values
    return out.dropna().reset_index(drop=True)


# --------------------------------------------------------------------------- ICP-OES (Syngistix)
def read_icp(src, name: str = ""):
    """PerkinElmer Syngistix 'Export > Report/ CSV' (one row per sample x analyte)
    or a wide table (rows = samples, columns = elements). Returns long DataFrame
    sample, analyte, conc, sd, unit."""
    raw = _to_bytes(src)
    if raw[:4] in (b"\xd0\xcf\x11\xe0", b"PK\x03\x04"):
        df = pd.read_excel(io.BytesIO(raw))
    else:
        text = _to_text(raw)
        lines = text.replace("\r", "").split("\n")
        d = _guess_delim(lines)
        sep = "\t" if d == "\t" else (";" if d == ";" else ",")
        # find header line
        hi = next((i for i, l in enumerate(lines) if re.search(r"sample\s*id|analyte", l, re.I)), 0)
        df = pd.read_csv(io.StringIO("\n".join(lines[hi:])), sep=sep, engine="python")
    df.columns = [str(c).strip() for c in df.columns]
    sc = find_col(df, r"sample\s*id|sample\s*name|^sample")
    ac = find_col(df, r"analyte")
    if ac is not None:
        cc = find_col(df, r"conc.*samp|sample\s*conc|conc.*calib|calib.*conc|^conc|concentration")
        sdc = find_col(df, r"^sd|std|stdev|sd\s*\(")
        uc = find_col(df, r"unit")
        out = pd.DataFrame({"sample": df[sc].astype(str), "analyte": df[ac].astype(str),
                            "conc": pd.to_numeric(df[cc].astype(str).str.replace(",", "."), errors="coerce")})
        out["sd"] = pd.to_numeric(df[sdc].astype(str).str.replace(",", "."), errors="coerce") if sdc else np.nan
        out["unit"] = df[uc].astype(str) if uc else ""
        out["element"] = out["analyte"].str.extract(r"([A-Z][a-z]?)")[0]
        return out.dropna(subset=["conc"])
    # wide
    long = df.melt(id_vars=[sc], var_name="analyte", value_name="conc").rename(columns={sc: "sample"})
    long["conc"] = pd.to_numeric(long["conc"], errors="coerce")
    long["sd"] = np.nan; long["unit"] = ""
    long["element"] = long["analyte"].str.extract(r"([A-Z][a-z]?)")[0]
    return long.dropna(subset=["conc"])


# --------------------------------------------------------------------------- Image Lab (SDS-PAGE)
def read_imagelab(src, name: str = ""):
    """Bio-Rad Image Lab 'Lane and Bands' table export (Analysis Table > Export > to CSV/Excel).
    Returns DataFrame lane, band, mw, band_pct, lane_pct, volume."""
    raw = _to_bytes(src)
    if raw[:4] in (b"\xd0\xcf\x11\xe0", b"PK\x03\x04"):
        df = pd.read_excel(io.BytesIO(raw))
    else:
        df = read_csv_table(io.BytesIO(raw), (r"lane", r"band"))
    lane = find_col(df, r"^lane$|lane\s*(no|number|#)?$")
    band = find_col(df, r"band\s*(no|number|#)|^band$")
    mw = find_col(df, r"mol.*w|kda|mw")
    bp = find_col(df, r"band\s*%")
    lp = find_col(df, r"lane\s*%")
    vol = find_col(df, r"adj.*vol|^volume")
    out = pd.DataFrame({"lane": df[lane] if lane else np.arange(len(df)) + 1})
    out["band"] = df[band] if band else 1
    out["mw"] = pd.to_numeric(df[mw], errors="coerce") if mw else np.nan
    out["band_pct"] = pd.to_numeric(df[bp], errors="coerce") if bp else np.nan
    out["lane_pct"] = pd.to_numeric(df[lp], errors="coerce") if lp else np.nan
    out["volume"] = pd.to_numeric(df[vol], errors="coerce") if vol else np.nan
    lab = find_col(df, r"band\s*label|label|name")
    out["band_label"] = df[lab].astype(str) if lab else out["band"].astype(str)
    tc = find_col(df, r"time")
    if tc:
        out["time"] = pd.to_numeric(df[tc], errors="coerce")
    return out


# --------------------------------------------------------------------------- Zetasizer
def read_zetasizer(src, name: str = ""):
    """Malvern Zetasizer software 'Export' (File > Export > CSV with selected parameters).
    Returns dict(summary=DataFrame, dist=DataFrame|None) with columns sample, z_avg, pdi, zeta, zeta_sd."""
    text = _to_text(src)
    lines = text.replace("\r", "").split("\n")
    d = _guess_delim(lines)
    sep = "\t" if d == "\t" else (";" if d == ";" else ",")
    df = pd.read_csv(io.StringIO(text), sep=sep, engine="python")
    df.columns = [str(c).strip() for c in df.columns]
    sn = find_col(df, r"sample\s*name|^sample")
    za = find_col(df, r"z-?ave")
    pdi = find_col(df, r"pdi")
    zp = find_col(df, r"zeta\s*potential\s*\(mv\)$|zeta\s*potential(?!.*(sd|dev))")
    zsd = find_col(df, r"zeta.*(sd|dev)")
    summ = pd.DataFrame({"sample": df[sn].astype(str) if sn else [f"S{i+1}" for i in range(len(df))]})
    for key, col in [("z_avg", za), ("pdi", pdi), ("zeta", zp), ("zeta_sd", zsd)]:
        summ[key] = pd.to_numeric(df[col], errors="coerce") if col else np.nan
    sizes = [c for c in df.columns if re.match(r"sizes?\s*\[\d+\]", c, re.I)]
    ints = [c for c in df.columns if re.match(r"intensit(y|ies)\s*\[\d+\]", c, re.I)]
    dist = None
    if sizes and ints:
        recs = []
        for i, row in df.iterrows():
            for s, v in zip(sizes, ints):
                recs.append((summ["sample"].iloc[i], float(row[s]), float(row[v])))
        dist = pd.DataFrame(recs, columns=["sample", "size_nm", "intensity"])
    return {"summary": summ, "dist": dist}


# --------------------------------------------------------------------------- XPS
def read_xps(src, name: str = ""):
    """CasaXPS / Thermo Avantage / generic exported spectrum: binding energy (eV) + counts.
    VAMAS (.vms) single-block files are parsed natively (regular-scan blocks)."""
    text = _to_text(src)
    if text.startswith("VAMAS Surface Chemical Analysis Standard"):
        return _read_vamas(text)
    df = read_xy(io.StringIO(text))
    return df


def _read_vamas(text):
    L = [l.strip() for l in text.split("\n")]
    # Minimal regular-scan VAMAS reader: find 'abscissa start/increment' via known layout
    nums = []
    try:
        idx = max(i for i, l in enumerate(L) if l.lower() in ("kinetic energy", "binding energy"))
        start, step = float(L[idx + 2]), float(L[idx + 3])
        npts_idx = idx + 4
        # skip corresponding variables until data
        n = None
        for j in range(npts_idx, npts_idx + 40):
            if L[j].isdigit() and int(L[j]) > 50:
                n = int(L[j]); data_start = j + 3; break
        vals = [float(v) for v in L[data_start:data_start + n]]
        ke = start + step * np.arange(n)
        x = ke if L[idx].lower() == "binding energy" else 1486.6 - ke
        return pd.DataFrame({"x": x, "y": vals})
    except Exception as e:
        raise ValueError("VAMAS file layout not recognised - export the region as CSV/TXT from CasaXPS") from e


# --------------------------------------------------------------------------- SEM
_UNIT_NM = {"nm": 1.0, "a": 0.1, "å": 0.1, "ang": 0.1, "angstrom": 0.1, "pm": 1e-3, "um": 1e3, "µm": 1e3, "μm": 1e3,
            "micron": 1e3, "microns": 1e3, "mm": 1e6, "m": 1e9, "cm": 1e7}
EM_EXT = (".dm3", ".dm4", ".ser", ".emd", ".mrc", ".emi")


def unit_to_nm(unit) -> float | None:
    u = (unit or "").strip().lower().replace("\\u00b5", "µ")
    return _UNIT_NM.get(u)


def read_em(src, name: str = ""):
    """TEM/STEM files via ncempy: Gatan DigitalMicrograph .dm3/.dm4, FEI/TIA .ser, Velox/Berkeley .emd, .mrc.
    Returns the same dict as read_sem (image 8-bit, pixel_nm, meta). Stacks -> first frame."""
    import tempfile
    import ncempy.io as nio
    raw = _to_bytes(src)
    ext = os.path.splitext(name or getattr(src, "name", ""))[1].lower() or ".dm4"
    with tempfile.TemporaryDirectory() as td:
        fn = os.path.join(td, "img" + ext)
        open(fn, "wb").write(raw)
        d = nio.read(fn)
    data = np.asarray(d["data"])
    while data.ndim > 2:
        data = data[0]
    if data.ndim < 2 or min(data.shape) < 8:
        raise ValueError("This file holds a spectrum or line profile, not an image")
    ps = d.get("pixelSize"); pu = d.get("pixelUnit")
    pixel_nm = None
    if ps is not None and len(ps):
        f = unit_to_nm(pu[-1] if isinstance(pu, (list, tuple)) else pu)
        if f:
            pixel_nm = float(ps[-1]) * f
    img = np.asarray(data, float)
    lo, hi = np.percentile(img, [0.1, 99.9])
    img8 = np.clip((img - lo) / max(hi - lo, 1e-9) * 255, 0, 255).astype(np.uint8)
    meta = {"source": f"ncempy ({ext})", "pixelSize": str(ps), "pixelUnit": str(pu)}
    return {"image": img8, "pixel_nm": pixel_nm, "meta": meta, "databar_row": None, "kind": "TEM"}


def parse_sem_sidecar(text: str, width_px=None):
    """Instrument text sidecars next to the image:
    JEOL .txt  ($CM_MAG, $$SM_MICRON_BAR <px>, $$SM_MICRON_MARKER <1um>, $CM_FULL_SIZE <w h>),
    Hitachi .txt (PixelSize=<nm>, Magnification=), TESCAN .hdr (PixelSizeX=<m>, Magnification=),
    Zeiss/other 'Image Pixel Size = 4.96 nm'. Returns dict(pixel_nm, magnification, source)."""
    out = {"pixel_nm": None, "magnification": None, "source": None}
    if not text:
        return out
    t = text.replace("\r", "")
    g = lambda rx: (re.search(rx, t, re.I | re.M) or [None, None])[1]
    if "$CM_" in t or "$$SM_" in t:
        out["source"] = "JEOL sidecar"
        bar, mark = g(r"^\$\$SM_MICRON_BAR\s+(" + NUM + ")"), g(r"^\$\$SM_MICRON_MARKER\s+(\S+)")
        full = re.search(r"^\$CM_FULL_SIZE\s+(\d+)\s+(\d+)", t, re.M)
        if bar and mark:
            mm = re.fullmatch(r"(" + NUM + r")\s*(nm|um|µm|mm)", mark.strip(), re.I)
            if mm:
                f = {"nm": 1, "um": 1e3, "µm": 1e3, "mm": 1e6}[mm.group(2).lower()]
                px = float(bar)
                if full and width_px:            # bar length is given in CM_FULL_SIZE pixels
                    px *= width_px / float(full.group(1))
                out["pixel_nm"] = float(mm.group(1)) * f / px
        m = g(r"^\$CM_MAG\s+(" + NUM + ")")
        out["magnification"] = float(m) if m else None
        return out
    m = g(r"^\s*PixelSizeX?\s*[=:]\s*(" + NUM + ")")
    if m:
        v = float(m.replace(",", "."))
        if re.search(r"^\s*\[MAIN\]|TESCAN|PixelSizeX", t, re.I | re.M) and v < 1e-3:
            out["pixel_nm"], out["source"] = v * 1e9, "TESCAN .hdr"          # metres
        else:
            out["pixel_nm"], out["source"] = v, "Hitachi sidecar"           # nm
    if out["pixel_nm"] is None:
        m = re.search(r"Image Pixel Size\s*=\s*(" + NUM + r")\s*(nm|µm|um|pm)", t)
        if m:
            out["pixel_nm"] = float(m.group(1)) * {"nm": 1, "µm": 1e3, "um": 1e3, "pm": 1e-3}[m.group(2)]
            out["source"] = "Zeiss text"
    m = re.search(r"^\s*Mag(?:nification)?\s*[=:]\s*([^\n;]{1,30})", t, re.I | re.M)
    if m:
        from .semcal import parse_magnification
        out["magnification"] = parse_magnification("Mag = " + m.group(1) + ("" if re.search(r"[xX×]", m.group(1)) else " x"))
    return out


def read_sem(src, name: str = "", sidecar: str | None = None):
    """SEM/TEM image: TIFF/PNG/JPG (+ .dm3/.dm4/.ser/.emd/.mrc through read_em). Returns dict(image=2D uint8, pixel_nm=float|None, meta=dict, databar_row=int|None).
    Pixel size is read from FEI/Thermo (Philips XL30 successor) [Scan] PixelWidth,
    Zeiss SmartSEM (tag 34118 'Image Pixel Size'), Hitachi/JEOL text tags or
    the XL30 magnification + 127 mm reference width when only magnification is stored."""
    import tifffile
    from PIL import Image
    if os.path.splitext(name or getattr(src, "name", ""))[1].lower() in EM_EXT:
        return read_em(src, name or getattr(src, "name", ""))
    raw = _to_bytes(src)
    meta, pixel_nm = {}, None
    img = None
    try:
        with tifffile.TiffFile(io.BytesIO(raw)) as tf:
            page = tf.pages[0]
            img = page.asarray()
            if page.photometric == 3 and page.colormap is not None and img.ndim == 2:     # palette TIFF -> grey levels
                cm = np.asarray(page.colormap, float)
                img = (cm[:, img].transpose(1, 2, 0) @ [0.299, 0.587, 0.114]) / 257.0
            sis = getattr(tf, "sis_metadata", None)                                       # Olympus analySIS / SIS
            if sis and sis.get("pixelsizex"):
                pixel_nm = float(sis["pixelsizex"]) * 1e9
                meta["calibration"] = f"Olympus SIS ({sis.get('cameraname', 'analySIS')})"
                if sis.get("magnification"):
                    meta["magnification"] = float(sis["magnification"])
                meta["instrument"] = str(page.tags["Artist"].value) if "Artist" in page.tags else "Olympus SIS"
            for tag in page.tags.values():
                v = tag.value
                s = v if isinstance(v, str) else (str(v) if isinstance(v, dict) else None)
                if s:
                    meta[tag.name] = s[:20000]
            ij = getattr(tf, "imagej_metadata", None) or {}
            if ij and ij.get("unit") and "XResolution" in page.tags:
                num, den = page.tags["XResolution"].value
                f = unit_to_nm(str(ij.get("unit")))
                if f and num:
                    pixel_nm = den / num * f          # ImageJ: pixels per unit
                    meta["calibration"] = f"ImageJ ({ij.get('unit')})"
            fei = getattr(tf, "fei_metadata", None) or getattr(tf, "sem_metadata", None)
            if fei and isinstance(fei, dict):
                for sec in ("Scan", "EScan", "Image"):
                    if sec in fei and "PixelWidth" in fei[sec]:
                        pixel_nm = float(fei[sec]["PixelWidth"]) * 1e9; break
    except Exception:
        img = np.array(Image.open(io.BytesIO(raw)))
    blob = "\n".join(v for v in meta.values() if isinstance(v, str))
    if pixel_nm is None:
        m = re.search(r"PixelWidth\s*=\s*(" + NUM + ")", blob)
        if m:
            pixel_nm = float(m.group(1).rstrip(",.").replace(",", ".")) * 1e9
    if pixel_nm is None:
        m = re.search(r"Image Pixel Size\s*=\s*(" + NUM + r")\s*(nm|µm|um|pm)", blob)
        if m:
            f = {"nm": 1, "µm": 1e3, "um": 1e3, "pm": 1e-3}[m.group(2)]
            pixel_nm = float(m.group(1)) * f
    if pixel_nm is None:
        m = re.search(r"PixelSize\s*[=:]\s*(" + NUM + ")", blob)
        if m:
            pixel_nm = float(m.group(1))
    if pixel_nm is not None and "calibration" not in meta:
        meta["calibration"] = "ImageJ" if "calibration" in meta else ("FEI/Thermo PixelWidth" if "PixelWidth" in blob or "FEI" in blob
                                                                       else ("Zeiss Image Pixel Size" if "Image Pixel Size" in blob else "file tags"))
    mag = re.search(r"Mag(?:nification)?'?\s*[=:]\s*(" + NUM + r")\s*([kK](?=\s*[xX]))?", blob)
    if mag:
        meta["magnification"] = float(mag.group(1).rstrip(",.").replace(",", ".")) * (1e3 if mag.group(2) else 1)
    if sidecar:
        sc = parse_sem_sidecar(sidecar, img.shape[1] if img is not None and img.ndim >= 2 else None)
        if pixel_nm is None and sc["pixel_nm"]:
            pixel_nm = sc["pixel_nm"]; meta["calibration"] = sc["source"]
        if sc["magnification"] and "magnification" not in meta:
            meta["magnification"] = sc["magnification"]
    if img.ndim == 3:
        img = (img[..., :3].astype(float) @ [0.299, 0.587, 0.114])
    img = np.asarray(img, float)
    lo, hi = np.percentile(img, [0.1, 99.9])
    img8 = np.clip((img - lo) / max(hi - lo, 1e-9) * 255, 0, 255).astype(np.uint8)
    return {"image": img8, "pixel_nm": pixel_nm, "meta": meta, "databar_row": detect_databar(img8)}


def pixel_from_magnification(mag: float, width_px: int, ref_width_mm: float = 127.0) -> float:
    """Pixel size (nm) for Philips/FEI XL-series images where magnification is quoted
    relative to a 127 mm (Polaroid 4x5) display width. Always verify against the databar."""
    return ref_width_mm * 1e6 / (mag * width_px)


def detect_databar(img8: np.ndarray):
    """First row of a burned-in databar (uniform black or white band with sparse text) at the bottom, or None."""
    h = img8.shape[0]
    med = np.median(img8, axis=1)
    dark = (med <= 8) | (med >= 247)
    if not dark[-1]:
        return None
    for r in range(int(h * 0.6), h - 4):
        if dark[r:r + 4].all() and dark[r:].mean() > 0.6:
            return r
    return None


# --------------------------------------------------------------------------- generic kinetics / tables
def read_generic(src, name: str = ""):
    return read_table(src)
