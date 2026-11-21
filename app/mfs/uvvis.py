"""UV/Vis spectroscopy (v2.9): readers + analysis.

Readers
-------
* Implen NanoPhotometer / NPOS **Excel export** (.xlsx/.xls): a results/summary sheet + one tab per sample with the
  full 200–900 nm absorbance spectrum.  The parser is layout-agnostic: on every sheet it looks for a wavelength axis
  (≥ 15 monotonic values between 150 and 3300 nm) in a column *or* a row and takes every parallel numeric series as one
  spectrum; header cells / the first text cell give the names.  Key–value lines (method, path length, dilution …)
  and a results table (a header row with 'Sample'…) are returned as metadata.
* Implen **.ids** (Implen Document Source): an *encrypted* archive that only NPOS can open.  We sniff it (zip / XML /
  text containers are read when present); a truly encrypted file raises :class:`IDSEncrypted` with export advice.
* Generic CSV / TXT / DAT (any delimiter, decimal comma ok), Agilent Cary CSV (paired columns), JCAMP-DX (XYDATA=(X++(Y..Y))
  and XYPOINTS).

Analysis
--------
A ⇄ %T ⇄ %R conversions, blank subtraction, Savitzky–Golay smoothing / derivatives, normalisation, peak picking (λmax),
Beer–Lambert calibration (linear / through origin / quadratic; R², s_y/x, LOD = 3.3 s/m, LOQ = 10 s/m, unknowns with
SE), single-λ kinetics with a first-order fit, Kubelka–Munk F(R), Tauc plots with automatic or manual linear window
(direct n = 2 / indirect n = ½), and protein / peptide assay presets (BCA, Bradford, Lowry, Biuret, ninhydrin, OPA,
A280, A205).
"""
from __future__ import annotations

import io
import re
import zipfile
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

HC = 1239.841984       # eV·nm


class IDSEncrypted(ValueError):
    pass


@dataclass
class Spectrum:
    name: str
    x: np.ndarray                 # wavelength (nm), ascending
    y: np.ndarray
    ytype: str = "A"              # "A" absorbance, "T" %T, "R" %R
    meta: dict = field(default_factory=dict)
    source: str = ""


# ----------------------------------------------------------------------------- generic helpers
_YTYPE = [(re.compile(r"%\s*t\b|transmi", re.I), "T"), (re.compile(r"%\s*r\b|reflect", re.I), "R"),
          (re.compile(r"\babs|absorb|\bod\b|\ba\b", re.I), "A")]


def guess_ytype(*texts, y=None):
    for t in texts:
        for rx, k in _YTYPE:
            if t and rx.search(str(t)):
                return k
    if y is not None and np.nanmax(y) > 10 and np.nanmin(y) >= -1:
        return "T"
    return "A"


def _num(v):
    if v is None:
        return np.nan
    if isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool):
        return float(v)
    s = str(v).strip().replace("\u2212", "-")
    if re.fullmatch(r"-?\d+,\d+(e-?\d+)?", s, re.I):
        s = s.replace(",", ".")
    try:
        return float(s)
    except ValueError:
        return np.nan


def _is_axis(vals, lo=150.0, hi=3300.0, min_n=15):
    v = np.asarray(vals, float)
    ok = np.isfinite(v)
    if ok.sum() < min_n:
        return False
    w = v[ok]
    if w.min() < lo or w.max() > hi:
        return False
    d = np.diff(w)
    return bool((np.all(d > 0) or np.all(d < 0)) and np.ptp(w) >= 20)


def _longest_run(mask):
    best, cur, start, bs = 0, 0, 0, 0
    for i, m in enumerate(list(mask) + [False]):
        if m:
            if cur == 0:
                start = i
            cur += 1
        else:
            if cur > best:
                best, bs = cur, start
            cur = 0
    return bs, best


def _label_above(grid, r0, c, n_up=6):
    for r in range(r0 - 1, max(-1, r0 - 1 - n_up), -1):
        v = grid[r][c]
        if isinstance(v, str) and v.strip() and not np.isfinite(_num(v)):
            return v.strip()
    return None


def _label_left(row, c0):
    for c in range(c0 - 1, -1, -1):
        v = row[c]
        if isinstance(v, str) and v.strip() and not np.isfinite(_num(v)):
            return v.strip()
    return None


def spectra_from_grid(grid, sheet="", source=""):
    """Find wavelength axes in a 2-D list of cells (columns first, then rows) and return the spectra next to them."""
    R_ = len(grid)
    C_ = max((len(r) for r in grid), default=0)
    g = [list(r) + [None] * (C_ - len(r)) for r in grid]
    num = np.array([[_num(v) for v in r] for r in g], float) if R_ and C_ else np.zeros((0, 0))
    out, used_cols = [], set()
    # --- column-wise
    for c in range(C_):
        if c in used_cols:
            continue
        r0, n = _longest_run(np.isfinite(num[:, c]))
        if n < 15 or not _is_axis(num[r0:r0 + n, c]):
            continue
        x = num[r0:r0 + n, c]
        xlab = _label_above(g, r0, c) or ""
        if not re.search(r"nm|wave|λ|lambda|wl", xlab, re.I) and xlab and re.search(r"time|min|sec|conc", xlab, re.I):
            continue
        for c2 in range(c + 1, C_):
            y = num[r0:r0 + n, c2]
            if np.isfinite(y).sum() < 0.8 * n:
                if np.isfinite(y).sum() == 0 and c2 > c + 1:
                    break
                continue
            if _is_axis(y) and c2 + 1 < C_ and np.isfinite(num[r0:r0 + n, c2 + 1]).sum() >= 0.8 * n and \
                    np.allclose(np.sort(y[np.isfinite(y)]), np.sort(x[np.isfinite(x)]), atol=0.51):
                break                                     # next x/y pair (Cary-style paired columns)
            lab = _label_above(g, r0, c2) or (_label_above(g, r0, c) if c2 == c + 1 else None) or f"{sheet or 'spectrum'} {c2 - c}"
            if C_ == 2 or (c2 == c + 1 and (not lab or re.search(r"^(abs|absorbance|a|od|%t|%r|y)\b", lab, re.I))):
                lab = sheet or lab
            used_cols.add(c2)
            out.append(_mk(lab, x, y, guess_ytype(_label_above(g, r0, c2), xlab), source, sheet))
        used_cols.add(c)
    if out:
        return out
    # --- row-wise (wavelengths across the top, one sample per row)
    for r in range(R_):
        c0, n = _longest_run(np.isfinite(num[r]))
        if n < 15 or not _is_axis(num[r, c0:c0 + n]):
            continue
        x = num[r, c0:c0 + n]
        for r2 in range(r + 1, R_):
            y = num[r2, c0:c0 + n]
            if np.isfinite(y).sum() < 0.8 * n:
                if np.isfinite(y).sum() == 0:
                    break
                continue
            lab = _label_left(g[r2], c0) or f"{sheet or 'spectrum'} {r2 - r}"
            out.append(_mk(lab, x, y, guess_ytype(lab, _label_left(g[r], c0)), source, sheet))
        if out:
            break
    return out


def _mk(name, x, y, ytype, source, sheet):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    o = np.argsort(x)
    return Spectrum(str(name), x[o], y[o], ytype, {"sheet": sheet}, source)


_KV = re.compile(r"^\s*([A-Za-z][\w \-/().%°µ]*?)\s*[:=]\s*(.+?)\s*$")


def meta_from_grid(grid):
    """Key–value pairs (two adjacent cells or 'key: value' in one cell) from the top of a sheet."""
    meta = {}
    for r in grid[:80]:
        cells = [v for v in r if v is not None and str(v).strip() not in ("", "nan")]
        if not cells:
            continue
        if len(cells) == 1 and isinstance(cells[0], str):
            m = _KV.match(cells[0])
            if m:
                meta.setdefault(m.group(1).strip(), m.group(2).strip())
        elif len(cells) == 2 and isinstance(cells[0], str) and not np.isfinite(_num(cells[0])):
            if re.search(r"wavelength|\bnm\b|\[nm\]", cells[0], re.I) and not np.isfinite(_num(cells[1])):
                continue                                   # header of a spectrum table, not a setting
            meta.setdefault(cells[0].strip().rstrip(":"), cells[1])
    return meta


def table_from_grid(grid, min_rows=1):
    """A results table: the first row that contains a 'sample' header cell and ≥ 2 more header cells."""
    for i, r in enumerate(grid[:200]):
        hdr = [str(v).strip() if v is not None else "" for v in r]
        if not any(re.fullmatch(r"(sample\s*(name|id)?|name|#|no\.?)", h, re.I) for h in hdr):
            continue
        if sum(1 for h in hdr if h and h.lower() != "nan") < 3:
            continue
        rows = []
        for r2 in grid[i + 1:]:
            vals = list(r2) + [None] * (len(hdr) - len(r2))
            if all(v is None or str(v).strip() in ("", "nan") for v in vals):
                break
            rows.append(vals[:len(hdr)])
        if len(rows) >= min_rows:
            cols, seen = [], {}
            for h in hdr:
                h = h or "col"
                seen[h] = seen.get(h, 0) + 1
                cols.append(h if seen[h] == 1 else f"{h} ({seen[h]})")
            df = pd.DataFrame(rows, columns=cols)
            keep = [c for c in df.columns if not df[c].isna().all() and c != "col"]
            return df[keep]
    return None


# ----------------------------------------------------------------------------- readers
@dataclass
class UVFile:
    spectra: list
    meta: dict = field(default_factory=dict)
    results: pd.DataFrame | None = None
    instrument: str = ""
    software: str = ""
    sheets: list = field(default_factory=list)


def _detect_instrument(text):
    t = text.lower()
    if "nanophotometer" in t or "implen" in t or "npos" in t:
        return "Implen NanoPhotometer", "NPOS (NanoPhotometer OS)"
    if "cary" in t:
        return "Agilent Cary", "Cary WinUV"
    if "nanodrop" in t:
        return "Thermo NanoDrop", "NanoDrop"
    if "shimadzu" in t or "uvprobe" in t:
        return "Shimadzu", "UVProbe / LabSolutions"
    if "lambda" in t and "perkin" in t:
        return "PerkinElmer Lambda", "UV WinLab"
    return "", ""


def read_excel(data: bytes, name=""):
    try:
        xl = pd.read_excel(io.BytesIO(data), sheet_name=None, header=None, engine=None)
    except Exception:
        xl = pd.read_excel(io.BytesIO(data), sheet_name=None, header=None, engine="xlrd")
    spectra, meta, results, alltext = [], {}, None, []
    for sh, df in xl.items():
        grid = df.where(pd.notna(df), None).values.tolist()
        alltext.extend(str(v) for r in grid[:60] for v in r if isinstance(v, str))
        sp = spectra_from_grid(grid, sheet=str(sh), source=name)
        spectra.extend(sp)
        m = meta_from_grid(grid)
        for k, v in m.items():
            meta.setdefault(k, v)
        if results is None:
            t = table_from_grid(grid)
            if t is not None and not sp:
                results = t
    inst, sw = _detect_instrument(" ".join(alltext) + " " + " ".join(f"{k} {v}" for k, v in meta.items()))
    _dedupe(spectra)
    return UVFile(spectra, meta, results, inst, sw, list(xl))


def _dedupe(spectra):
    """Drop exact duplicates (an overlay tab repeating the per-sample tabs) and make names unique."""
    seen, keep = {}, []
    for s in spectra:
        h = (len(s.x), round(float(np.nansum(s.y)), 6), round(float(s.x[0]), 3) if len(s.x) else 0)
        if h in seen:
            continue
        seen[h] = s
        keep.append(s)
    names = {}
    for s in keep:
        names[s.name] = names.get(s.name, 0) + 1
        if names[s.name] > 1:
            s.name = f"{s.name} ({names[s.name]})"
    spectra[:] = keep


def _grid_from_text(text):
    lines = [l for l in text.replace("\r", "").split("\n")]
    sample = "\n".join(lines[:200])
    delim = max(["\t", ";", ","], key=lambda d: sample.count(d))
    if sample.count(delim) < 5:
        delim = None
    grid = []
    for l in lines:
        if delim:
            parts = l.split(delim)
            if delim != "," and re.search(r"\d,\d", l):
                parts = [p.replace(",", ".") if re.fullmatch(r"\s*-?\d+,\d+\s*", p) else p for p in parts]
        else:
            parts = l.split()
        grid.append([p.strip() if p.strip() != "" else None for p in parts])
    return grid


def read_jcamp(text, name=""):
    hdr = dict((m.group(1).upper(), m.group(2).strip()) for m in re.finditer(r"^##([^=]+)=(.*)$", text, re.M))
    fx = float(hdr.get("XFACTOR", 1)); fy = float(hdr.get("YFACTOR", 1))
    xs, ys = [], []
    body = text.split("##XYDATA=")[1] if "##XYDATA=" in text else text.split("##XYPOINTS=")[1]
    body = body.split("\n", 1)[1]
    xy_pairs = "##XYPOINTS=" in text and "##XYDATA=" not in text
    dx = None
    if not xy_pairs:
        fx_, lx = float(hdr.get("FIRSTX", 0)), float(hdr.get("LASTX", 0))
        npts = int(float(hdr.get("NPOINTS", 0)) or 0)
        dx = (lx - fx_) / (npts - 1) if npts > 1 else None
    for line in body.splitlines():
        if line.startswith("##"):
            break
        nums = [float(v) for v in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)]
        if not nums:
            continue
        if xy_pairs:
            for i in range(0, len(nums) - 1, 2):
                xs.append(nums[i] * fx); ys.append(nums[i + 1] * fy)
        else:
            x0 = nums[0] * fx
            for j, v in enumerate(nums[1:]):
                xs.append(x0 + j * (dx if dx else 1.0)); ys.append(v * fy)
    unit = (hdr.get("XUNITS", "") or "").lower()
    x = np.asarray(xs); y = np.asarray(ys)
    if "1/cm" in unit or "wavenumber" in unit:
        x = 1e7 / x
    yt = guess_ytype(hdr.get("YUNITS", ""), y=y)
    return UVFile([_mk(hdr.get("TITLE", name) or name, x, y, yt, name, "")], {k.title(): v for k, v in hdr.items() if len(v) < 80})


def read_ids(data: bytes, name=""):
    if data[:2] == b"PK":
        z = zipfile.ZipFile(io.BytesIO(data))
        sp, meta = [], {"container": "zip", "members": ", ".join(z.namelist()[:20])}
        for n in z.namelist():
            raw = z.read(n)
            try:
                sub = read_any(raw, n)
                sp.extend(sub.spectra); meta.update({k: v for k, v in sub.meta.items() if k not in meta})
            except Exception:
                continue
        if sp:
            return UVFile(sp, meta, None, "Implen NanoPhotometer", "NPOS (.ids)")
    head = data[:4096]
    printable = sum(32 <= b < 127 or b in (9, 10, 13) for b in head) / max(1, len(head))
    if printable > 0.95:
        txt = data.decode("utf-8", "replace")
        if txt.lstrip().startswith("<"):
            nums = re.findall(r"<(?:\w*:)?(?:Wavelength|WL|X)[^>]*>\s*([-\d.eE+]+)\s*<", txt)
            vals = re.findall(r"<(?:\w*:)?(?:Absorbance|Abs|Y|Value)[^>]*>\s*([-\d.eE+]+)\s*<", txt)
            if len(nums) >= 15 and len(nums) == len(vals):
                return UVFile([_mk(name, [float(v) for v in nums], [float(v) for v in vals], "A", name, "")],
                              {"container": "xml"}, None, "Implen NanoPhotometer", "NPOS (.ids)")
        return read_text(txt, name)
    ent = _entropy(data[:65536])
    raise IDSEncrypted(f"{name or 'file'}: Implen Document Source (.ids) – an encrypted NPOS archive (byte entropy "
                       f"{ent:.2f} bit/byte) that only the NanoPhotometer software can open. In NPOS open the .ids file and "
                       "use Export → Excel (all samples, spectra included) or PDF, then load the .xlsx here.")


def _entropy(b):
    if not b:
        return 0.0
    c = np.bincount(np.frombuffer(b, np.uint8), minlength=256) / len(b)
    c = c[c > 0]
    return float(-(c * np.log2(c)).sum())


def read_text(text, name=""):
    if "##JCAMP" in text.upper() or "##XYDATA" in text.upper():
        return read_jcamp(text, name)
    grid = _grid_from_text(text)
    sp = spectra_from_grid(grid, sheet="", source=name)
    stem = re.sub(r"\.[^.]+$", "", name.split("/")[-1])
    for s in sp:
        if len(sp) == 1 and re.fullmatch(r"(spectrum \d+|abs\w*|a|y|od|%\s*[tr]|absorbance|transmittance|reflectance)?", s.name.strip(), re.I):
            s.name = stem
    if not sp:
        raise ValueError(f"{name}: no wavelength column (≥ 15 monotonic values in 150–3300 nm) found")
    meta = meta_from_grid(grid)
    inst, sw = _detect_instrument(text[:5000])
    _dedupe(sp)
    return UVFile(sp, meta, None, inst, sw)


def read_any(data, name=""):
    if hasattr(data, "read"):
        name = name or getattr(data, "name", "")
        data = data.getvalue() if hasattr(data, "getvalue") else data.read()
    low = name.lower()
    if low.endswith(".ids"):
        return read_ids(data, name)
    if data[:4] == b"PK\x03\x04" and (low.endswith((".xlsx", ".xlsm")) or b"xl/" in data[:2000]):
        return read_excel(data, name)
    if data[:8] == b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1":
        return read_excel(data, name)
    for enc in ("utf-8-sig", "utf-16", "cp1252"):
        try:
            txt = data.decode(enc)
            if "\x00" in txt[:200]:
                continue
            return read_text(txt, name)
        except UnicodeDecodeError:
            continue
    raise ValueError(f"{name}: unknown UV/Vis format")


# ----------------------------------------------------------------------------- conversions / processing
def to_absorbance(s: Spectrum):
    if s.ytype == "A":
        return s.y.copy()
    T = np.clip(s.y / 100.0, 1e-6, None)
    return -np.log10(T)


def to_percent_T(s: Spectrum):
    if s.ytype == "T":
        return s.y.copy()
    if s.ytype == "R":
        return s.y.copy()
    return 100.0 * 10 ** (-s.y)


def kubelka_munk(R_percent):
    R = np.clip(np.asarray(R_percent, float) / 100.0, 1e-6, 0.999999)
    return (1 - R) ** 2 / (2 * R)


def subtract_blank(x, y, xb, yb, factor=1.0):
    return y - factor * np.interp(x, xb, yb, left=np.nan, right=np.nan)


def smooth(y, window=11, poly=3, deriv=0, dx=1.0):
    from scipy.signal import savgol_filter
    y = np.asarray(y, float)
    w = int(window) | 1
    if w < poly + 2 or len(y) < w:
        return np.gradient(y, dx) if deriv == 1 else (np.gradient(np.gradient(y, dx), dx) if deriv == 2 else y)
    return savgol_filter(y, w, int(poly), deriv=int(deriv), delta=dx)


def derivative(x, y, order=1, window=11, poly=3):
    x = np.asarray(x, float)
    dx = float(np.median(np.diff(x))) if len(x) > 1 else 1.0
    return smooth(y, window, max(poly, order + 1), deriv=order, dx=dx)


def maxima_from_derivative(x, y, d1, frac=0.05):
    """Wavelengths where dy/dλ crosses + → − (spectral maxima), ignoring crossings where y < frac of its range."""
    x = np.asarray(x, float); y = np.asarray(y, float); d1 = np.asarray(d1, float)
    lo, hi = np.nanmin(y), np.nanmax(y)
    idx = np.where((np.sign(d1[:-1]) > 0) & (np.sign(d1[1:]) <= 0))[0]
    out = []
    for i in idx:
        j = i if abs(d1[i]) < abs(d1[i + 1]) else i + 1
        if y[j] - lo >= frac * (hi - lo):
            xi = x[i] + (x[i + 1] - x[i]) * d1[i] / (d1[i] - d1[i + 1]) if d1[i] != d1[i + 1] else x[i]
            out.append(float(xi))
    return out


def normalise(x, y, mode="none", at=None, rng=None):
    y = np.asarray(y, float)
    m = np.ones_like(y, bool) if rng is None else (x >= rng[0]) & (x <= rng[1])
    if mode == "max":
        d = np.nanmax(y[m]); return y / d if d else y
    if mode == "min-max":
        lo, hi = np.nanmin(y[m]), np.nanmax(y[m]); return (y - lo) / (hi - lo) if hi > lo else y
    if mode == "at λ" and at is not None:
        d = float(np.interp(at, x, y)); return y / d if d else y
    if mode == "area":
        d = np.trapezoid(np.abs(y[m]), x[m]) if hasattr(np, "trapezoid") else np.trapz(np.abs(y[m]), x[m])
        return y / d if d else y
    return y


def peaks(x, y, rng=None, prominence=0.05, max_n=6, minima=False):
    """[(λ, y)] of the most prominent maxima (minima for %T) inside rng; prominence is a fraction of the range."""
    from scipy.signal import find_peaks
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(y) & (np.ones_like(x, bool) if rng is None else (x >= rng[0]) & (x <= rng[1]))
    xs, ys = x[m], y[m]
    if len(xs) < 5:
        return []
    yy = -ys if minima else ys
    pr = prominence * (np.nanmax(yy) - np.nanmin(yy) or 1)
    idx, props = find_peaks(yy, prominence=pr)
    order = np.argsort(props["prominences"])[::-1][:max_n]
    return sorted([(float(xs[i]), float(ys[i])) for i in idx[order]])


def lambda_max(x, y, rng=None):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(y) & (np.ones_like(x, bool) if rng is None else (x >= rng[0]) & (x <= rng[1]))
    if not m.any():
        return np.nan, np.nan
    i = int(np.nanargmax(np.where(m, y, -np.inf)))
    return float(x[i]), float(y[i])


def value_at(x, y, lam, half_width=0.0):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if half_width > 0:
        m = (x >= lam - half_width) & (x <= lam + half_width)
        if m.any():
            return float(np.nanmean(y[m]))
    return float(np.interp(lam, x, y))


# ----------------------------------------------------------------------------- Beer–Lambert / calibration
def calibration(c, A, model="linear", weights=None):
    """Fit A = m·c + b (or m·c, or a·c² + m·c + b). Returns a dict with slope, intercept, SEs, R², s_y/x, LOD, LOQ."""
    c = np.asarray(c, float); A = np.asarray(A, float)
    ok = np.isfinite(c) & np.isfinite(A)
    c, A = c[ok], A[ok]
    n = len(c)
    if model == "through origin":
        X = c[:, None]
    elif model == "quadratic":
        X = np.column_stack([c ** 2, c, np.ones_like(c)])
    else:
        X = np.column_stack([c, np.ones_like(c)])
    p = X.shape[1]
    if n <= p:
        raise ValueError(f"need more than {p} standards for a {model} fit")
    W = np.ones(n) if weights is None else np.asarray(weights, float)[ok]
    Xw = X * np.sqrt(W)[:, None]; Aw = A * np.sqrt(W)
    coef, *_ = np.linalg.lstsq(Xw, Aw, rcond=None)
    fit = X @ coef
    res = A - fit
    dof = n - p
    s = float(np.sqrt(np.sum(W * res ** 2) / dof))
    cov = s ** 2 * np.linalg.inv(Xw.T @ Xw)
    se = np.sqrt(np.diag(cov))
    ss_tot = float(np.sum((A - A.mean()) ** 2)) if model != "through origin" else float(np.sum(A ** 2))
    r2 = 1 - float(np.sum(res ** 2)) / ss_tot if ss_tot else np.nan
    out = {"model": model, "n": n, "coef": coef.tolist(), "cov": cov.tolist(), "se": se.tolist(), "R2": r2, "s_yx": s,
           "c": c.tolist(), "A": A.tolist(), "residuals": res.tolist()}
    if model == "quadratic":
        a, m, b = coef
        out.update(slope=float(m), intercept=float(b), quad=float(a), slope_se=float(se[1]), intercept_se=float(se[2]))
    elif model == "through origin":
        out.update(slope=float(coef[0]), intercept=0.0, slope_se=float(se[0]), intercept_se=0.0)
    else:
        out.update(slope=float(coef[0]), intercept=float(coef[1]), slope_se=float(se[0]), intercept_se=float(se[1]))
    m = out["slope"]
    out["LOD"] = 3.3 * s / m if m else np.nan
    out["LOQ"] = 10.0 * s / m if m else np.nan
    return out


def predict_conc(cal, A, n_rep=1):
    """Concentration (and its SE, ISO 8466 / Miller & Miller) for measured absorbance(s) A."""
    A = np.atleast_1d(np.asarray(A, float))
    m, b = cal["slope"], cal["intercept"]
    if cal["model"] == "quadratic":
        a = cal["quad"]
        out = []
        for v in A:
            disc = m * m - 4 * a * (b - v)
            if disc < 0:
                out.append(np.nan); continue
            r = [(-m + s_ * np.sqrt(disc)) / (2 * a) for s_ in (1, -1)] if a else [(v - b) / m]
            r = [x for x in r if x >= -1e-9] or r
            out.append(min(r, key=lambda x: abs(x - np.mean(cal["c"]))))
        c0 = np.array(out)
        return c0, np.full_like(c0, np.nan)
    c0 = (A - b) / m
    cs = np.asarray(cal["c"]); As = np.asarray(cal["A"])
    n = len(cs)
    sxx = float(np.sum((cs - cs.mean()) ** 2))
    se = cal["s_yx"] / abs(m) * np.sqrt(1.0 / n_rep + 1.0 / n + (A - As.mean()) ** 2 / (m ** 2 * sxx)) if sxx else np.full_like(A, np.nan)
    return c0, se


def epsilon(cal, path_cm=1.0, conc_unit_to_M=1.0):
    """Molar absorptivity ε (M⁻¹ cm⁻¹) from the slope (per concentration unit) and the path length."""
    return cal["slope"] / (path_cm * conc_unit_to_M)


ASSAYS = {
    # name: (λ / nm, default standard, standard concentrations, unit, suggested model, note)
    "BCA (bicinchoninic acid)": (562, "BSA", [0, 25, 125, 250, 500, 750, 1000, 1500, 2000], "µg/mL", "quadratic",
                                 "Cu⁺–BCA₂ complex; slightly curved above ~1000 µg/mL (Pierce: 4-PL or quadratic). Peptides with < 3 residues respond weakly."),
    "Bradford (Coomassie G-250)": (595, "BSA", [0, 125, 250, 500, 750, 1000, 1500], "µg/mL", "quadratic",
                                   "Non-linear; dye binds Arg/Lys/aromatics – small peptides (< 3 kDa) are not detected."),
    "Lowry (Folin–Ciocalteu)": (750, "BSA", [0, 25, 50, 100, 200, 400, 600, 800, 1000], "µg/mL", "quadratic",
                                "Tyr/Trp + Cu; many interferents (Tris, EDTA, reducing agents)."),
    "Biuret": (540, "BSA", [0, 1, 2.5, 5, 7.5, 10], "mg/mL", "linear", "Peptide bonds + Cu²⁺; low sensitivity, linear."),
    "Ninhydrin (free amines)": (570, "Gly / Leu", [0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0], "mM", "linear",
                                "Ruhemann's purple; Pro/Hyp give 440 nm. Counts free α-amino groups → follows peptide-bond hydrolysis."),
    "OPA (o-phthalaldehyde)": (340, "Gly / Ser", [0, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0], "mM", "linear",
                               "Primary amines + thiol; degree of hydrolysis DH = h/h_tot (Nielsen 2001). Fast, no heating."),
    "A280 (direct protein)": (280, "ε from sequence", [], "mg/mL", "through origin",
                              "c = A280 / (ε·l); ε from Trp/Tyr/Cys content (Pace 1995)."),
    "A205 (peptide bond)": (205, "ε205 ≈ 31 mL mg⁻¹ cm⁻¹", [], "mg/mL", "through origin",
                            "Scopes 1974; works for peptides without aromatics; buffer must be transparent."),
}


def degree_of_hydrolysis(h_serine_eq, h_tot):
    """OPA/ninhydrin: DH (%) = 100 · h / h_tot (mmol amino groups per g protein)."""
    return 100.0 * np.asarray(h_serine_eq, float) / h_tot


# ----------------------------------------------------------------------------- kinetics at one wavelength
_TIME = re.compile(r"(?<![\d.])(\d+(?:[.,]\d+)?)\s*(h|hr|hrs|hours?|min|mins|minutes?|s|sec|seconds?|d|days?)(?![A-Za-z])", re.I)


def time_from_name(name, to="min"):
    m = _TIME.search(name or "")
    if not m:
        return np.nan
    v = float(m.group(1).replace(",", ".")); u = m.group(2).lower()
    sec = v * (3600 if u.startswith("h") else 60 if u.startswith("m") else 86400 if u.startswith("d") else 1)
    return sec / {"s": 1, "min": 60, "h": 3600, "d": 86400}[to]


def first_order_fit(t, A):
    """A(t) = A∞ + (A0 − A∞)·exp(−k t). Returns k, SE(k), A0, A∞, half-life, R²."""
    from scipy.optimize import curve_fit
    t = np.asarray(t, float); A = np.asarray(A, float)
    ok = np.isfinite(t) & np.isfinite(A)
    t, A = t[ok], A[ok]
    if len(t) < 4:
        raise ValueError("need ≥ 4 time points")
    f = lambda tt, a0, ai, k: ai + (a0 - ai) * np.exp(-k * tt)
    k0 = 3.0 / (np.ptp(t) or 1.0)
    p, cov = curve_fit(f, t, A, p0=[A[0], A[-1], k0], maxfev=20000)
    res = A - f(t, *p)
    r2 = 1 - np.sum(res ** 2) / np.sum((A - A.mean()) ** 2) if np.ptp(A) else np.nan
    se = np.sqrt(np.diag(cov)) if np.all(np.isfinite(cov)) else [np.nan] * 3
    return {"A0": float(p[0]), "Ainf": float(p[1]), "k": float(p[2]), "k_se": float(se[2]), "t_half": float(np.log(2) / p[2]),
            "R2": float(r2), "fn": lambda tt: f(np.asarray(tt, float), *p)}


# ----------------------------------------------------------------------------- Tauc / band gap
def tauc(x_nm, y, kind="direct", from_reflectance=False):
    """Return (hν / eV, (F·hν)^n) with n = 2 (direct allowed) or ½ (indirect allowed). y is F(R) or absorbance."""
    x = np.asarray(x_nm, float); y = np.asarray(y, float)
    E = HC / x
    a = kubelka_munk(y) if from_reflectance else y
    n = 2.0 if kind.startswith("direct") else 0.5
    T = np.clip(a * E, 0, None) ** n
    o = np.argsort(E)
    return E[o], T[o]


def band_gap(E, T, window=None, min_pts=6, min_frac=0.1):
    """Linear fit of the absorption edge, Eg = x-intercept.  Automatic: among windows of ≥ min_pts points whose
    signal lies between min_frac and 0.9 of the maximum, take the steepest one with R² ≥ 0.98."""
    E = np.asarray(E, float); T = np.asarray(T, float)
    ok = np.isfinite(E) & np.isfinite(T)
    E, T = E[ok], T[ok]
    if window:
        m = (E >= window[0]) & (E <= window[1])
        if m.sum() < 3:
            raise ValueError("fewer than 3 points in the Tauc window")
        s, b = np.polyfit(E[m], T[m], 1)
        r2 = np.corrcoef(E[m], T[m])[0, 1] ** 2
        return {"Eg": float(-b / s), "slope": float(s), "intercept": float(b), "R2": float(r2), "window": (float(E[m].min()), float(E[m].max()))}
    # v2.18.12: the window search is vectorised (closed-form least squares from running sums and sliding min/max). It gives the
    # same windows and the same result as the old per-window polyfit loop (test_v21812_known_answers), but is much faster on long scans.
    Tm = np.nanmax(T)
    best = None
    n = len(E)
    for w in sorted({min_pts, max(min_pts, n // 60), max(min_pts, n // 30), max(min_pts, n // 15)}):
        if n - w <= 0:
            continue
        ncount = n - w                                   # window starts 0 .. n-w-1, as in the original loop
        Ex = np.lib.stride_tricks.sliding_window_view(E, w)[:ncount]
        Tx = np.lib.stride_tricks.sliding_window_view(T, w)[:ncount]
        Sx = Ex.sum(1); Sy = Tx.sum(1)
        Sxx = (Ex * Ex).sum(1); Syy = (Tx * Tx).sum(1); Sxy = (Ex * Tx).sum(1)
        den_x = w * Sxx - Sx * Sx
        den_y = w * Syy - Sy * Sy
        with np.errstate(divide="ignore", invalid="ignore"):
            s = (w * Sxy - Sx * Sy) / den_x
            r2 = (w * Sxy - Sx * Sy) ** 2 / (den_x * den_y)
            b = (Sy - s * Sx) / w
        ymin = Tx.min(1); ymax = Tx.max(1)
        valid = (ymin >= min_frac * Tm) & (ymax <= 0.9 * Tm) & (s > 0) & (r2 >= 0.98)
        valid &= np.isfinite(s) & np.isfinite(r2)
        if not valid.any():
            continue
        cand = np.where(valid, s, -np.inf)
        i = int(np.argmax(cand))                         # first maximum, as the original loop did
        if best is None or cand[i] > best[0]:
            best = (float(s[i]), float(b[i]), float(r2[i]), (float(Ex[i].min()), float(Ex[i].max())))
    if best is None:
        raise ValueError("no linear absorption edge found – pick the window manually")
    s, b, r2, win = best
    return {"Eg": float(-b / s), "slope": float(s), "intercept": float(b), "R2": float(r2), "window": win}


# ----------------------------------------------------------------------------- tables
def spectra_table(spectra, ycol="y"):
    rows = []
    for s in spectra:
        rows.append(pd.DataFrame({"sample": s.name, "wavelength_nm": s.x, ycol: s.y}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["sample", "wavelength_nm", ycol])


def wide_table(spectra):
    """One wavelength column + one column per spectrum (union grid, NaN where a spectrum has no point)."""
    xs = np.unique(np.concatenate([s.x for s in spectra])) if spectra else np.array([])
    df = pd.DataFrame({"wavelength_nm": xs})
    for s in spectra:
        df[s.name] = np.interp(xs, s.x, s.y, left=np.nan, right=np.nan)
    return df


def tauc_auto_kind(x_nm, y, from_reflectance=False):
    """v2.18.11: chooses the direct (n = 2) or indirect (n = 1/2) Tauc form by the better automatic linear edge (higher R² of the
    steepest linear window, direct on a tie). Returns {'kind', 'Eg', 'R2', 'direct': {...}, 'indirect': {...}, 'note'}.
    General practice: the rule is a heuristic; check the Tauc plot of the chosen form and the other one before you report the band gap."""
    out = {}
    for kind in ("direct", "indirect"):
        E, T = tauc(x_nm, y, kind=kind, from_reflectance=from_reflectance)
        try:
            out[kind] = band_gap(E, T)
        except ValueError as exc:
            out[kind] = {"status": str(exc)}
    ok = {k: v for k, v in out.items() if "Eg" in v}
    if not ok:
        raise ValueError("no linear absorption edge found for either direct or indirect form – pick the window manually")
    best = max(ok, key=lambda k: (ok[k]["R2"], k == "direct"))
    return {"kind": best, "Eg": ok[best]["Eg"], "R2": ok[best]["R2"], "direct": out["direct"], "indirect": out["indirect"],
            "note": "chosen by the higher R² of the automatic linear edge (direct on a tie)"}
