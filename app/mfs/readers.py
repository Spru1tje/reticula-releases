"""Instrument file readers.

Every reader returns pandas DataFrames with clean float columns. Native binary
formats are parsed with established open-source libraries where they exist
(nmrglue for Bruker TopSpin, brukeropusreader for OPUS, tifffile for SEM TIFF);
everything else goes through a robust text parser that copes with comma or
semicolon delimiters, decimal commas (Belgian/European locale), metadata
headers and footers.
"""
from __future__ import annotations
import io, os, re, struct, zipfile, tempfile
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


def read_spa(src):
    """v2.18.2: Thermo OMNIC .spa (one spectrum per file). The directory of 16-byte entries starts at byte 304 (type, offset, size);
    entry 2 holds the point count, first and last x (cm-1), entry 3 the float32 intensities. Verified against OMNIC CSV exports
    (x within 0.001 cm-1, y identical). Returns x ascending and y."""
    b = _to_bytes(src)
    if not b.startswith(b"Spectral Data File") and b[:5] != b"Spect":
        raise ValueError("Not an OMNIC .spa file")
    ent = {}
    for i in range(64):
        o = 304 + 16 * i
        if o + 16 > len(b) or (b[o] == 0 and i > 2):
            break
        ent.setdefault(b[o], struct.unpack_from("<II", b, o + 2))
    if 2 not in ent or 3 not in ent:
        raise ValueError("OMNIC .spa: no data block found")
    p2, p3 = ent[2][0], ent[3][0]
    nx = struct.unpack_from("<i", b, p2 + 4)[0]
    fx, lx = struct.unpack_from("<ff", b, p2 + 16)
    if nx < 2 or p3 + 4 * nx > len(b):
        raise ValueError("OMNIC .spa: unexpected data layout")
    y = np.frombuffer(b, "<f4", count=nx, offset=p3).astype(float)
    x = np.linspace(fx, lx, nx)
    if x[0] > x[-1]:
        x, y = x[::-1], y[::-1]
    out = pd.DataFrame({"x": x, "y": y})
    out.attrs["source"] = "OMNIC .spa"
    return out


def read_ftir(src, name: str = ""):
    ext = os.path.splitext(name or getattr(src, "name", ""))[1].lower()
    if ext == ".spa" or _to_bytes(src)[:18] == b"Spectral Data File":
        return read_spa(src)
    if re.fullmatch(r"\.\d+", ext or ""):
        return read_opus(src)
    raw = _to_bytes(src)
    if raw[:4] == b"\n\n\xfe\xfe" or (b"\x00" in raw[:64] and not raw[:2] in (b"\xff\xfe", b"\xfe\xff")):
        return read_opus(io.BytesIO(raw))
    return read_xy(io.BytesIO(raw))  # OPUS 'Save as Data Point Table' (.dpt) / CSV


def read_ftir_multi(src, name: str = ""):
    """v2.18.2: every spectrum in a text file with one column per spectrum (first column = x, later columns = y, the header names them).
    Returns a list of (label, DataFrame x, y); a single-spectrum file gives one entry. OPUS binaries give one entry."""
    ext = os.path.splitext(name or getattr(src, "name", ""))[1].lower()
    if ext == ".spa" or _to_bytes(src)[:18] == b"Spectral Data File":
        return [(None, read_spa(src))]                  # an .spa holds one spectrum (series are .spg)
    if re.fullmatch(r"\.\d+", ext or "") or _to_bytes(src)[:4] == b"\n\n\xfe\xfe":
        return [(None, read_ftir(src, name))]
    df = read_table(src)
    x = pd.to_numeric(df.iloc[:, 0], errors="coerce")
    out = []
    for j in range(1, df.shape[1]):
        y = pd.to_numeric(df.iloc[:, j], errors="coerce")
        if y.notna().mean() < 0.5:                  # text column (sample name, comment): not a spectrum
            continue
        label = str(df.columns[j]).strip()
        if not label or label.lower().startswith("unnamed") or label == "nan" or re.fullmatch(r"(col|c)\d+", label.lower()):
            label = f"spectrum {len(out) + 1}"
        d = pd.DataFrame({"x": x, "y": y}).dropna(subset=["x", "y"]).reset_index(drop=True)
        out.append((label, d))
    if not out:
        raise ValueError("No spectrum column found (expected: first column = wavenumber, next columns = intensities)")
    if len(out) == 1 and out[0][0] == "spectrum 1":
        out = [(None, out[0][1])]                       # one unnamed spectrum: keep the sample name as it is
    return out


# --------------------------------------------------------------------------- Edinburgh Instruments .FS (Fluoracle / F980 / FLS)
def read_edinburgh_fs(src, name: str = ""):
    """v2.18.4: Edinburgh Instruments Fluoracle binary `.FS` file (one scan, or 'EAI Multiple Scans' with several). Each scan block holds the file type, the axis names and units,
    the program version, the sample name, the instrument settings and the data as float32 (counts) on an even wavelength grid between the start and end wavelength that are
    stored in the block. Verified on 9 files of an FLS (NIR PMT) emission series: 351 points 900-1250 nm per scan.
    Returns a list of dicts {name, x, y, x_label, y_label, kind, detector, excitation, x_start, x_end, step} (one per scan; wavelength ascending)."""
    b = _to_bytes(src)
    if not b.startswith(b"EAI "):
        raise ValueError("Not an Edinburgh Instruments .FS file (it should start with 'EAI ...')")
    head = re.match(rb"EAI [A-Za-z ]+?\x00", b)
    kind_file = head.group().rstrip(b"\x00").decode("latin1") if head else "EAI"
    marks = [m.start() for m in re.finditer(rb"EAI (?:Emission|Excitation|Synchronous|Kinetics|Lifetime|Decay|Absorption|Scan)[A-Za-z ]*\x00", b)]
    if b.startswith(b"EAI Emission Scan") or kind_file.endswith("Scan"):
        marks = [m for m in marks if m > 0] or marks                   # a single-scan file repeats the type string once at the start
    scans = []
    for k, m0 in enumerate(marks):
        m1 = marks[k + 1] if k + 1 < len(marks) else len(b)
        blk = b[m0:m1]
        kind = re.match(rb"EAI [A-Za-z ]+", blk).group().decode("latin1")
        jw = blk.find(b"Wavelength")
        if jw < 8:
            continue
        x0, x1 = struct.unpack_from("<ff", blk, jw - 8)
        strs = [(mm.start(), mm.end(), mm.group().decode("latin1")) for mm in re.finditer(rb"[\x20-\x7e]{2,}", blk[:700])]
        sample, name_end = "", jw + 20
        ver = next(((s0, s1) for s0, s1, s in strs if re.fullmatch(r"\d+(\.\d+)+", s)), None)          # the program version string ends the axis block
        for s0, s1, s in strs:
            if (ver is None or s0 >= ver[1]) and s0 > jw + 20 and not re.fullmatch(r"\d+(\.\d+)+", s) and s.strip() and len(s) >= 2:
                sample, name_end = s, s1
                break
        tail = [s for s0, s1, s in strs if s0 >= name_end and len(s) > 3]               # settings strings (detector, excitation source) between the name and the data
        # data: the last N float32 of the block; N from the wavelength grid ('nice' steps) - the first position after the name that gives all-plausible numbers
        start = None
        for s in range(name_end, len(blk) - 12):
            if (len(blk) - s) % 4:
                continue
            n = (len(blk) - s) // 4
            step = (x1 - x0) / (n - 1)
            if not (abs(step * 1000 - round(step * 1000)) < 1e-6 * max(1.0, abs(step) * 1000) and abs(step) >= 0.001):
                continue
            a = np.frombuffer(blk, "<f4", n, s)
            if np.all(np.isfinite(a)) and np.abs(a).max() < 1e9:
                start = s
                break
        if start is None:
            raise ValueError("Could not locate the data of a scan in this .FS file")
        y = np.frombuffer(blk, "<f4", (len(blk) - start) // 4, start).astype(float)
        x = np.linspace(x0, x1, len(y))
        scans.append({"name": sample, "x": x, "y": y, "x_label": "Wavelength (nm)", "y_label": "Counts", "kind": kind, "x_start": float(x0), "x_end": float(x1),
                      "step": float(x[1] - x[0]) if len(x) > 1 else float("nan"), "detector": next((t for t in tail if "PMT" in t or "InGaAs" in t or "CCD" in t), ""),
                      "excitation": next((t for t in tail if "aser" in t or "amp" in t or "Xe" in t), ""), "n": len(y)})
    if not scans:
        raise ValueError("No scan found in this .FS file")
    return scans

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
def _zeta_excel_table(raw: bytes) -> pd.DataFrame:
    """The sheet of a Zetasizer Excel export that holds the results table: the header row is the first row (within 40) that names
    a sample column and at least one of Z-Average / PDI / Zeta potential."""
    sheets = pd.read_excel(io.BytesIO(raw), sheet_name=None, header=None)
    key = re.compile(r"z-?ave|pd?i\b|zeta|polydispersity", re.I)
    for nm, raw_df in sheets.items():
        for i in range(min(40, len(raw_df))):
            cells = [str(c).strip() for c in raw_df.iloc[i].tolist() if str(c) != "nan"]
            if any(re.search(r"sample|^name$", c, re.I) for c in cells) and sum(bool(key.search(c)) for c in cells) >= 1:
                df = raw_df.iloc[i + 1:].copy()
                df.columns = [str(c).strip() if str(c) != "nan" else f"col{j}" for j, c in enumerate(raw_df.iloc[i].tolist())]
                return df.dropna(how="all").reset_index(drop=True)
    raise ValueError("No Zetasizer results table found in this Excel file (needs a header row with 'Sample Name' and Z-Average / PDI / "
                     "Zeta Potential). Export the records from ZS Xplorer (or the Zetasizer software) as Excel/CSV with those parameters.")


DTS_VERIFIED_VERSION = "7.11"      # Zetasizer software version of the .dts files for which Z-average / PDI were checked against the instrument


def _dts_cumulants(d: bytes):
    """Z-average (nm) and PDI of ONE .dts record from its stored mean correlation function (cumulant analysis, ISO 22412):
    ln g1 = a - G t + (mu2/2) t^2 with g1 = sqrt((g2 - 1)/beta), Z-average = kB T / (3 pi eta D), D = G / q^2, PDI = mu2 / G^2,
    q = 4 pi n / lambda sin(theta/2). Angle and wavelength are two consecutive doubles in the record, n and eta are float32 values,
    T a double (25 degrees C when it is not found). Confirmed against the values the Zetasizer software shows (Nano ZS, software 7.11).
    Returns (z_avg, pdi) or (nan, nan)."""
    nan = (float("nan"), float("nan"))
    try:
        a8 = [(al + 8 * int(i), float(v)) for al in range(8) for i, v in enumerate(np.frombuffer(d[al:al + ((len(d) - al) // 8) * 8], "<f8")) if 632.0 <= v <= 633.5]
        lam = ang = None
        for off, v in sorted(a8):
            prev = np.frombuffer(d[off - 8:off], "<f8")[0] if off >= 8 else np.nan
            if 1.0 <= prev <= 180.0:
                ang, lam = float(prev), v
                break
        if lam is None:
            return nan
        f4 = np.frombuffer(d[:(len(d) // 4) * 4], "<f4")
        ri = next((float(v) for v in f4[100:200] if 1.30 <= v <= 1.60), 1.33)                     # refractive index of the dispersant (offset 572 B)
        eta = None                                                                                 # viscosity in mPa s: a float32 in the first bytes (not 4-aligned)
        for al in range(4):
            for v in np.frombuffer(d[al:al + 400], "<f4"):
                if 0.2 <= v <= 3.5 and abs(v - round(float(v))) > 1e-3 and 4 * abs(v * 1e5 - round(float(v) * 1e5)) < 1.0:
                    eta = float(v)
                    break
            if eta is not None:
                break
        if eta is None:
            return nan
        tc = next((float(v) for v in np.frombuffer(d[800:3000 - (3000 - 800) % 8], "<f8") if 4.0 <= v <= 95.0 and abs(v - round(v)) < 1e-9), 25.0)
        # mean correlation function: the first long decaying float32 array; lag times: the long increasing array (microseconds)
        corr = None
        for al in range(4):
            n = (len(d) - al) // 4
            b = np.frombuffer(d[al:al + n * 4], "<f4")
            for i in range(0, n - 150):
                s = b[i:i + 150]
                if np.all(np.isfinite(s)) and 0.3 < s[0] < 1.8 and s[0] > s[75] > s[149] and np.all(np.diff(s) < 0.03) and s[149] < 0.3 * s[0]:
                    j = i
                    while j + 1 < n and np.isfinite(b[j + 1]) and b[j + 1] <= b[j] + 0.03:
                        j += 1
                    corr = (al + 4 * i, j - i + 1)
                    break
            if corr:
                break
        lag = None
        for al in range(4):
            n = (len(d) - al) // 4
            b = np.frombuffer(d[al:al + n * 4], "<f4")
            i = 0
            while i < n - 100:
                j = i
                while j + 1 < n and np.isfinite(b[j + 1]) and b[j + 1] > b[j] > 0:
                    j += 1
                if j - i + 1 >= 150 and b[j] > 1e6 and (corr is None or al + 4 * i > corr[0]):
                    lag = (al + 4 * i, j - i + 1)
                    break
                i = max(j, i) + 1
            if lag:
                break
        if corr is None or lag is None:
            return nan
        m = min(corr[1], lag[1])
        g = np.frombuffer(d[corr[0]:corr[0] + m * 4], "<f4").astype(float)
        t = np.frombuffer(d[lag[0]:lag[0] + m * 4], "<f4").astype(float) * 1e-6
        beta = float(g[:3].mean())
        sel = g > 0.15 * beta
        if sel.sum() < 8 or beta <= 0:
            return nan
        y = 0.5 * np.log(np.clip(g / beta, 1e-9, None))
        A = np.vstack([np.ones(int(sel.sum())), -t[sel], 0.5 * t[sel] ** 2]).T
        c, *_ = np.linalg.lstsq(A, y[sel], rcond=None)
        gam, mu2 = float(c[1]), float(c[2])
        if gam <= 0:
            return nan
        q = 4 * np.pi * ri / (lam * 1e-9) * np.sin(np.radians(ang) / 2)
        dh = 1.380649e-23 * (tc + 273.15) / (3 * np.pi * eta * 1e-3 * (gam / q ** 2)) * 1e9
        return float(dh), float(mu2 / gam ** 2)
    except Exception:
        if os.environ.get("MFS_DEBUG"):
            raise
        return nan


def _dts_stored_summary(d: bytes):
    """v2.18.2: Z-average (nm) and PDI as the software stored them in ONE .dts record (software 8.x and later). The two float32 values
    stand directly before the start of the correlogram: a Z-average between 10 and 10 000 nm, then a PDI between 0.005 and 0.9, then ten
    values that fall from about 0.9 (the correlation function). Checked on 6 records of 2 files (software 8.02): all six equal the Z-average
    and PDI the Zetasizer software shows. Returns (z_avg, pdi), or None when no such pair is found (exactly one pair is expected per record)."""
    try:
        hits = []
        for al in range(4):
            m = (len(d) - al) // 4
            a = np.frombuffer(d[al:al + m * 4], "<f4")
            for k in range(m - 14):
                z, p = a[k], a[k + 1]
                if not (np.isfinite(z) and np.isfinite(p)) or not (10.0 <= z <= 1e4 and 0.005 < p < 0.9):
                    continue
                seg = a[k + 3:k + 13]
                if np.all(np.isfinite(seg)) and np.all((seg > 0.5) & (seg < 1.05)) and np.all(np.diff(seg) <= 0.0005):
                    hits.append((float(z), float(p)))
        return hits[0] if hits else None
    except Exception:
        if os.environ.get("MFS_DEBUG"):
            raise
        return None


def _dts_zeta_record(d: bytes):
    """v2.18.2: a ZETA record of a Zetasizer .dts (software 8.x). The result block lies at fixed distances from the END of the record
    (checked on 11 records against the Zeta Potential Reports of the software: zeta, deviation, conductivity, count rate, the three
    peaks with mean / area / SD all equal). Returns a dict or None when this is not a plausible zeta record."""
    try:
        n = len(d)
        if n < 1000:
            return None
        f = lambda r: struct.unpack_from("<f", d, n + r)[0]
        zp, zsd = f(-552), f(-548)
        if not (np.isfinite(zp) and np.isfinite(zsd)) or abs(zp) > 200 or not (0.0 < zsd < 200):
            return None
        peaks = [(f(-584 + 4 * i), f(-488 + 4 * i), f(-512 + 4 * i)) for i in range(3)]       # mean (mV), area (%), SD (mV)
        area = sum(p[1] for p in peaks)
        if not (60.0 <= area <= 100.5):                                      # only three peaks are stored: a fourth small one may be missing
            return None
        cond = struct.unpack_from("<d", d, n - 564)[0]
        cr = struct.unpack_from("<d", d, n - 200)[0]
        out = {"zeta": float(zp), "zeta_sd": float(zsd), "cond_mS_cm": float(cond) if np.isfinite(cond) else float("nan"),
               "count_kcps": float(cr) if np.isfinite(cr) and 0 < cr < 1e5 else float("nan"),
               "peaks": [(float(m), float(a), float(s)) for m, a, s in peaks if a > 0]}
        ax = None                                                          # apparent zeta potential axis (89 points, about -150 … +150 mV)
        for al in range(4):
            b = np.frombuffer(d, "<f4", (n - al) // 4, al)
            for i in np.where((b > -160) & (b < -140))[0]:
                s = b[i:i + 89]
                if len(s) == 89 and np.all(np.isfinite(s)) and np.all(np.diff(s) > 0) and 140 < s[-1] < 160:
                    ax = al + 4 * int(i); break
            if ax is not None:
                break
        if ax is not None and ax + 1336 + 356 <= n:
            x = np.frombuffer(d, "<f4", 89, ax).astype(float)
            c = np.frombuffer(d, "<f4", 89, ax + 1336).astype(float)
            if np.all(np.isfinite(c)) and c.min() >= 0 and c.sum() > 0 and abs((x * c).sum() / c.sum() - zp) < 0.5 and c.max() > 1e3:                  # strict: records whose block does not reproduce zeta get no distribution
                out["dist"] = (x, c)
        return out
    except Exception:
        if os.environ.get("MFS_DEBUG"):
            raise
        return None


def _zeta_dts(raw: bytes):
    """v2.14: Zetasizer Nano / Classic software `.dts` file (an OLE container with one stream REC<n> per record). Verified on real files:
    the sample name and the 70-bin size axis (0.4-10 000 nm) with the INTENSITY size distribution (sums to 100 %). Z-average and PDI are CALCULATED from the stored
    correlation function (see _dts_cumulants); the zeta potential is not decoded and stays empty. Returns the same dict as read_zetasizer or None when this is not a .dts (a .xls is an OLE file too)."""
    try:
        import olefile
    except ImportError:
        return None
    try:
        ole = olefile.OleFileIO(io.BytesIO(raw))
    except Exception:
        return None
    names = [s[0] for s in ole.listdir(streams=True, storages=False) if len(s) == 1]
    recs = sorted((n for n in names if re.fullmatch(r"REC\d+", n)), key=lambda n: int(n[3:]))
    if not recs or "Header" not in names:
        return None
    sig = np.array([0.4, 0.46323168, 0.53645897], "<f4").tobytes()
    rows, summ, sw_ver, zrows = [], [], "", []
    for n in recs:
        d = ole.openstream(n).read()
        i = d.find(sig)
        if i < 0 or i + 560 > len(d):
            zr = _dts_zeta_record(d)                                        # v2.18.2: a zeta potential record
            if zr is None:
                continue
            strs = [x.decode("utf-16le") for x in re.findall(rb"(?:[\x20-\x7e]\x00){3,}", d)]
            nm = next((s_.strip() for s_ in strs if re.fullmatch(r".+ \d+\s*", s_) and not s_.startswith(("C:", "Z:"))), f"record {int(n[3:])}")
            if strs and not sw_ver:
                sw_ver = strs[0].strip()
            mm = re.fullmatch(r"(.*\S)\s+(\d+)", nm)
            base, run = (mm.group(1), int(mm.group(2))) if mm else (nm, len(summ) + 1)
            k_rec = len(summ)
            pk = zr["peaks"] + [(float("nan"), float("nan"), float("nan"))] * (3 - len(zr["peaks"]))
            row = {"sample": base, "run": run, "rec": k_rec, "z_avg": np.nan, "pdi": np.nan, "zeta": zr["zeta"], "zeta_sd": zr["zeta_sd"],
                   "cond_mS_cm": zr["cond_mS_cm"], "count_kcps": zr["count_kcps"]}
            for j, (m_, a_, s_) in enumerate(pk, 1):
                row[f"peak{j}_mV"], row[f"peak{j}_area_pct"], row[f"peak{j}_sd_mV"] = m_, a_, s_
            summ.append(row)
            if "dist" in zr:
                zrows += [(base, float(a), float(b), k_rec) for a, b in zip(*zr["dist"])]
            continue
        sizes = np.frombuffer(d[i:i + 280], "<f4").astype(float)
        inten = np.frombuffer(d[i + 280:i + 560], "<f4").astype(float)
        if abs(inten.sum() - 100.0) > 1.0:
            continue
        strs = [x.decode("utf-16le") for x in re.findall(rb"(?:[\x20-\x7e]\x00){3,}", d)]
        nm = next((s_.strip() for s_ in strs if re.fullmatch(r".+ \d+\s*", s_) and not s_.startswith("C:")), f"record {int(n[3:])}")
        stored = _dts_stored_summary(d)                                                     # v2.18.2: Z-ave and PDI as stored by the software (8.x)
        verified = bool(strs) and strs[0].strip().startswith(DTS_VERIFIED_VERSION)        # 7.11: the cumulant route was checked against the instrument
        if stored is not None:
            za, pd_ = stored
        elif verified:
            za, pd_ = _dts_cumulants(d)
        else:
            za, pd_ = float("nan"), float("nan")
        if strs and not sw_ver:
            sw_ver = strs[0].strip()
        mm = re.fullmatch(r"(.*\S)\s+(\d+)", nm)                      # 'ZrO2-aliquot2 3' = sample 'ZrO2-aliquot2', run 3
        base, run = (mm.group(1), int(mm.group(2))) if mm else (nm, len(summ) + 1)
        k_rec = len(summ)
        summ.append({"sample": base, "run": run, "rec": k_rec, "z_avg": za, "pdi": pd_, "zeta": np.nan, "zeta_sd": np.nan})
        rows += [(base, float(a), float(b), np.nan, np.nan, k_rec) for a, b in zip(sizes, inten)]
    if not summ:
        return None
    sm = pd.DataFrame(summ)
    dist = pd.DataFrame(rows, columns=["sample", "size_nm", "intensity", "volume", "number", "rec"])
    if zrows:                                                          # v2.18.2: zeta potential records
        zd = pd.DataFrame(zrows, columns=["sample", "zeta_mV", "counts", "rec"])
        n_z = int(sm["zeta"].notna().sum())
        note = (f"From this .dts file (software {sw_ver or 'unknown'}): {n_z} zeta potential record(s) – zeta potential, its deviation (SD of the distribution), "
                "conductivity, count rate, the three zeta peaks (mean, area %, SD) and the zeta distribution, as stored by the software "
                "(each value equals the Zeta Potential Report).")
        if sm["z_avg"].notna().any():
            note += " Size records give the intensity size distribution, Z-average and PDI."
        return {"summary": sm, "dist": dist if len(dist) else None, "zdist": zd, "note": note}
    if sm["z_avg"].isna().all():
        return {"summary": sm, "dist": dist,
                "note": f"From this .dts file (software {sw_ver or 'unknown'}): sample names and the intensity size distribution. Z-average and PDI could not be read from these records "
                        f"(they are read for software 8.x, and for {DTS_VERIFIED_VERSION} by calculation); for ζ, use the Excel / CSV export."}
    return {"summary": sm, "dist": dist,
            "note": f"From this .dts file (software {sw_ver or 'unknown'}): sample names, the intensity size distribution, and Z-average and PDI as stored in each record "
                    "by the Zetasizer software (for software 7.11 calculated from the stored correlation function). "
                    "The ζ-potential is not in size records – export it from the Zetasizer software (File → Export → Excel / CSV)."}


def select_records(z: dict, keep) -> dict:
    """The same result with only the records whose number (column 'rec') is in `keep` - the means over a sample then follow."""
    keep = set(int(k) for k in keep)
    out = dict(z)
    out["summary"] = z["summary"][z["summary"]["rec"].isin(keep)].reset_index(drop=True)
    if z.get("dist") is not None:
        out["dist"] = z["dist"][z["dist"]["rec"].isin(keep)].reset_index(drop=True)
    if z.get("zdist") is not None:
        out["zdist"] = z["zdist"][z["zdist"]["rec"].isin(keep)].reset_index(drop=True)
    return out


def read_zetasizer(src, name: str = ""):
    """Malvern Zetasizer software 'Export' (CSV or Excel with the selected parameters; Zetasizer Nano / Ultra).
    Returns dict(summary=DataFrame, dist=DataFrame|None) with columns sample, z_avg, pdi, zeta, zeta_sd."""
    raw = _to_bytes(src)
    if raw[:4] == b"\xd0\xcf\x11\xe0":
        dts = _zeta_dts(raw)                                            # v2.14: Zetasizer .dts (size distribution + names)
        if dts is not None:
            return dts
    if raw[:2] == b"PK" or raw[:4] == b"\xd0\xcf\x11\xe0":            # v2.14: Excel export (.xlsx / .xls) as well as CSV
        df = _zeta_excel_table(raw)
    else:
        text = _to_text(raw)
        lines = text.replace("\r", "").split("\n")
        d = _guess_delim(lines)
        sep = "\t" if d == "\t" else (";" if d == ";" else ",")
        df = pd.read_csv(io.StringIO(text), sep=sep, engine="python")
    df.columns = [str(c).strip() for c in df.columns]
    sn = find_col(df, r"sample\s*name|^sample|^name$")
    za = find_col(df, r"z-?ave|z-?average")
    pdi = find_col(df, r"\bpd?i\b|polydispersity")
    zp = find_col(df, r"zeta\s*potential\s*\(mv\)$|zeta\s*potential(?!.*(sd|dev))|^zeta\s*\(mv\)")
    zsd = find_col(df, r"zeta.*(sd|dev)")
    summ = pd.DataFrame({"sample": df[sn].astype(str) if sn else [f"S{i+1}" for i in range(len(df))]})
    summ["rec"] = np.arange(len(summ))                      # v2.14: record number, to select / deselect single runs
    for key, col in [("z_avg", za), ("pdi", pdi), ("zeta", zp), ("zeta_sd", zsd)]:
        summ[key] = pd.to_numeric(df[col], errors="coerce") if col else np.nan
    def _bins(pattern):                           # {bin number: column} for Sizes[k] / Intensities[k] / Volumes[k] / Numbers[k]
        out = {}
        for c in df.columns:
            m = re.match(pattern + r"\s*\[(\d+)\]", c, re.I)
            if m:
                out[int(m.group(m.lastindex))] = c
        return out
    sizes, ints = _bins(r"sizes?"), _bins(r"intensit(?:y|ies)")
    vols, nums = _bins(r"volumes?"), _bins(r"numbers?")     # v2.14: volume and number distributions too
    dist = None
    if sizes and ints:
        recs = []
        num = lambda row, cols, k: float(row[cols[k]]) if k in cols and pd.notna(row[cols[k]]) else np.nan
        for i, row in df.iterrows():
            for k in sorted(set(sizes) & set(ints)):
                recs.append((summ["sample"].iloc[i], float(row[sizes[k]]), float(row[ints[k]]), num(row, vols, k), num(row, nums, k), i))
        dist = pd.DataFrame(recs, columns=["sample", "size_nm", "intensity", "volume", "number", "rec"])
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
    return {"image": img8, "raw": np.clip(img, 0, 255).astype(np.uint8), "pixel_nm": pixel_nm, "meta": meta, "databar_row": None, "kind": "TEM"}


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
    # v2.22.2: the contrast stretch above can erase a thin scale bar; calibration reads the unstretched values
    raw8 = np.clip(img, 0, 255).astype(np.uint8)
    return {"image": img8, "raw": raw8, "pixel_nm": pixel_nm, "meta": meta, "databar_row": detect_databar(img8)}


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
