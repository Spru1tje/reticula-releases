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
def read_sorption(src, name: str = ""):
    """Micromeritics MicroActive / 3Flex export (Excel .xls/.xlsx or CSV/TXT report).
    Looks for 'Relative Pressure (P/Po)' and 'Quantity Adsorbed (cm³/g STP)' columns
    (Isotherm Tabular Report) and, if present, 'Pore Width' + 'dV/dw' / 'dV/dlog(w)' / 'Cumulative Pore Volume' (DFT report).
    Returns dict(iso=DataFrame[p_rel, q, branch], psd=DataFrame|None, meta=dict)."""
    ext = os.path.splitext(name or getattr(src, "name", ""))[1].lower()
    raw = _to_bytes(src)
    frames = []
    if ext in (".xls", ".xlsx") or raw[:4] in (b"\xd0\xcf\x11\xe0", b"PK\x03\x04"):
        sheets = pd.read_excel(io.BytesIO(raw), sheet_name=None, header=None)
        for sn, sh in sheets.items():
            frames.append(sh.astype(str).replace("nan", "").values.tolist())
    else:
        text = _to_text(raw)
        lines = text.replace("\r", "").split("\n")
        d = _guess_delim(lines)
        frames.append([_split(l, d) for l in lines])
    iso_ads, iso_des, psd, meta = [], [], None, {}
    for grid in frames:
        # metadata
        for row in grid:
            line = " ".join(str(c) for c in row if str(c).strip() not in ("", "nan", "None", "<NA>"))
            for key, pat in [("sample", r"Sample:\s*(.+)"), ("mass_g", r"Sample Mass:\s*(" + NUM + ")"),
                             ("bet_reported", r"BET Surface Area:\s*(" + NUM + ")")]:
                mm = re.search(pat, line)
                if mm and key not in meta:
                    meta[key] = mm.group(1).strip()
        # locate header rows
        for r, row in enumerate(grid):
            cells = [str(c) for c in row]
            ip = [i for i, c in enumerate(cells) if re.search(r"relative\s*pressure|p\s*/\s*p\s*[o0°]", c, re.I)]
            iq = [i for i, c in enumerate(cells) if re.search(r"quantity\s*adsorbed|volume\s*adsorbed|cm.?3.?/?\s*g", c, re.I)]
            if ip and iq:
                # there may be ads + des pairs side by side
                pairs = []
                for a in ip:
                    q = min([q for q in iq if q > a] or [None], key=lambda z: z if z is not None else 1e9)
                    if q is not None:
                        pairs.append((a, q))
                for k, (a, q) in enumerate(pairs):
                    pts = []
                    for row2 in grid[r + 1:]:
                        if a >= len(row2) or q >= len(row2):
                            break
                        x, y = _num(str(row2[a])), _num(str(row2[q]))
                        if x is None or y is None or np.isnan(x) or np.isnan(y):
                            if pts:
                                break
                            continue
                        pts.append((x, y))
                    hdr = " ".join(cells).lower()
                    if pts:
                        if len(pairs) > 1 and k == 1:
                            iso_des += pts
                        else:
                            iso_ads += pts
            iw = [i for i, c in enumerate(cells) if re.search(r"pore\s*(width|diameter|size)", c, re.I)]
            iv = [i for i, c in enumerate(cells) if re.search(r"dV/d|incremental|differential|cumulative", c, re.I)]
            if iw and iv and psd is None:
                cols = {"w": iw[0]}
                for i in iv:
                    c = cells[i].lower()
                    if "cumul" in c: cols.setdefault("cum", i)
                    elif "log" in c: cols.setdefault("dvdlogw", i)
                    elif "dv/d" in c: cols.setdefault("dvdw", i)
                    elif "increment" in c: cols.setdefault("inc", i)
                pts = []
                for row2 in grid[r + 1:]:
                    vals = {k: _num(str(row2[i])) if i < len(row2) else None for k, i in cols.items()}
                    if vals["w"] is None or np.isnan(vals["w"]):
                        if pts: break
                        continue
                    pts.append(vals)
                if pts:
                    psd = pd.DataFrame(pts)
                    unit = cells[iw[0]]
                    psd.attrs["width_unit"] = "Å" if re.search("Å|ang|A\\)", unit) else ("nm" if "nm" in unit else "Å")
    if not iso_ads:
        # fallback: plain 2- or 3-column text p/p0, q
        df = read_table(io.BytesIO(raw))
        p, q = df.iloc[:, 0].values, df.iloc[:, 1].values
        imax = int(np.argmax(p))
        iso_ads = list(zip(p[: imax + 1], q[: imax + 1]))
        iso_des = list(zip(p[imax + 1:], q[imax + 1:]))
    elif not iso_des:
        # single column pair holding ads then des
        arr = np.array(iso_ads)
        imax = int(np.argmax(arr[:, 0]))
        if imax < len(arr) - 2:
            iso_des = [tuple(r) for r in arr[imax + 1:]]
            iso_ads = [tuple(r) for r in arr[: imax + 1]]
    iso = pd.DataFrame([(a, b, "ads") for a, b in iso_ads] + [(a, b, "des") for a, b in iso_des],
                       columns=["p_rel", "q", "branch"])
    return {"iso": iso, "psd": psd, "meta": meta}


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
def read_sem(src, name: str = ""):
    """SEM TIFF/PNG/JPG. Returns dict(image=2D uint8, pixel_nm=float|None, meta=dict, databar_row=int|None).
    Pixel size is read from FEI/Thermo (Philips XL30 successor) [Scan] PixelWidth,
    Zeiss SmartSEM (tag 34118 'Image Pixel Size'), Hitachi/JEOL text tags or
    the XL30 magnification + 127 mm reference width when only magnification is stored."""
    import tifffile
    from PIL import Image
    raw = _to_bytes(src)
    meta, pixel_nm = {}, None
    img = None
    try:
        with tifffile.TiffFile(io.BytesIO(raw)) as tf:
            page = tf.pages[0]
            img = page.asarray()
            for tag in page.tags.values():
                v = tag.value
                s = v if isinstance(v, str) else (str(v) if isinstance(v, dict) else None)
                if s:
                    meta[tag.name] = s[:20000]
            fei = getattr(tf, "fei_metadata", None) or getattr(tf, "sem_metadata", None)
            if fei and isinstance(fei, dict):
                for sec in ("Scan", "EScan", "Image"):
                    if sec in fei and "PixelWidth" in fei[sec]:
                        pixel_nm = float(fei[sec]["PixelWidth"]) * 1e9; break
    except Exception:
        img = np.array(Image.open(io.BytesIO(raw)))
    blob = "\n".join(meta.values())
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
    mag = re.search(r"Mag(?:nification)?'?\s*[=:]\s*(" + NUM + ")", blob)
    if mag:
        meta["magnification"] = float(mag.group(1).rstrip(",.").replace(",", "."))
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
