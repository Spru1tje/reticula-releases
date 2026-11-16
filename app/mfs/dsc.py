"""Differential scanning calorimetry: readers (Mettler STARe, TA TRIOS / Universal Analysis, Netzsch Proteus,
PerkinElmer Pyris, generic CSV/XLSX), segment splitting, baselines, peak evaluation (onset / peak / endset / ΔH),
glass transition (onset / midpoint / inflection / Δcp) and plotting (exo arrow, offsets, TGA–DSC overlay).

Internal sign convention: `q` is stored EXO-POSITIVE (exothermic up). Each reader converts from its vendor's
convention; the user can flip it. Heat flow is kept in W/g (= mW/mg) when the file gives it normalised or when
the sample mass is known; otherwise in mW."""
from __future__ import annotations
import io, re
import numpy as np
import pandas as pd

from .readers import _to_text, _num, find_col

NUM = r"[-+]?\d+(?:[.,]\d+)?(?:[eE][-+]?\d+)?"


def _f(s):
    v = _num(s)
    return None if v is None or (isinstance(v, float) and np.isnan(v)) else float(v)


# ============================================================================ readers
def _block(lines, start, ncol_min=2):
    """Numeric rows from `start` until the first non-numeric line. Splits on tab/;/whitespace, decimal commas ok."""
    rows = []
    for l in lines[start:]:
        s = l.strip()
        if not s:
            if rows:
                break
            continue
        cells = re.split(r"\t|;", s) if ("\t" in s or ";" in s) else s.split()
        if "," in s and not ("\t" in s or ";" in s) and len(s.split()) == 1:
            cells = s.split(",")
        vals = [_f(c) for c in cells]
        if len(vals) >= ncol_min and all(v is not None for v in vals[:ncol_min]):
            rows.append([v if v is not None else np.nan for v in vals])
        elif rows:
            break
    if not rows:
        return None
    n = max(len(r) for r in rows)
    return np.array([r + [np.nan] * (n - len(r)) for r in rows], float)


def _mass(text):
    for rx in (r"(?:Sample\s*(?:Weight|Mass|size)|SAMPLE\s*MASS\s*/mg|Size)\s*[:=;\t ]*\s*(" + NUM + r")\s*(mg)?",
               r",\s*(" + NUM + r")\s*mg\b", r"(" + NUM + r")\s*mg\b"):
        m = re.search(rx, text, re.I)
        if m:
            v = _f(m.group(1))
            if v and 0.01 < v < 1000:
                return v
    return None


def _unit_kind(u: str):
    u = (u or "").lower().replace(" ", "").replace("^-1", "/g").replace("g-1", "/g").replace("wg", "w/g")
    if "mw/mg" in u or "w/g" in u or "/g" in u:
        return "W/g"
    if "mw" in u:
        return "mW"
    if "w" in u:
        return "W"
    return "?"


def read_mettler(text):
    """Mettler-Toledo STARe 'Curve Values' export: Index  t [s]  Ts [°C]  Tr [°C]  Value [Wg^-1 | mW].
    Several 'Curve Name:' blocks = several curves/segments. Footer: Sample '<name>, 5,6400 mg', Method, Results
    (Integral, normalized, Onset, Peak) as reported by STARe. STARe convention: endothermic = negative."""
    L = text.replace("\r", "").split("\n")
    curves, i = [], 0
    while i < len(L):
        if L[i].strip().startswith("Curve Values"):
            hdr = L[i + 1].split(); units = re.findall(r"\[([^\]]*)\]", L[i + 2])
            a = _block(L, i + 3, 3)
            name = ""
            for k in range(i - 1, max(-1, i - 6), -1):
                if L[k].strip().startswith("Curve Name"):
                    name = L[k + 1].strip().lstrip("!]0123456789[").strip() or L[k + 1].strip(); break
            if a is not None:
                cols = {h.lower(): j for j, h in enumerate(hdr)}
                t = a[:, cols.get("t", 1)] / 60.0 if units and units[cols.get("t", 1) - (1 if len(units) < len(hdr) else 0)] == "s" else a[:, cols.get("t", 1)]
                ts = a[:, cols["ts"]] if "ts" in cols else a[:, 2]
                tr = a[:, cols["tr"]] if "tr" in cols else None
                v = a[:, cols.get("value", a.shape[1] - 1)]
                unit = _unit_kind(units[-1] if units else "")
                curves.append({"name": name, "t": t, "T": ts, "Tr": tr, "q": -v * 0 + v, "unit": unit})
                i += 3 + len(a)
                continue
        i += 1
    meta = {"vendor": "Mettler-Toledo STARe", "mass_mg": _mass("\n".join(l for l in L if "mg" in l and "," in l) or text)}
    m = re.search(r"Method:\s*\n\s*(.+)", text)
    if m:
        meta["method"] = m.group(1).strip()
    res = {}
    for key, rx in (("integral_mJ", r"Integral\s+(" + NUM + r")\s*mJ"), ("dH_Jg", r"normalized\s+(" + NUM + r")\s*Jg"),
                    ("onset_C", r"Onset\s+(" + NUM + r")\s*°C"), ("peak_C", r"Peak\s+(" + NUM + r")\s*°C"),
                    ("endset_C", r"Endset\s+(" + NUM + r")\s*°C"), ("tg_mid_C", r"Midpoint\s*(?:ISO|DIN|ASTM)?\s+(" + NUM + r")\s*°C")):
        mm = re.search(rx, text)
        if mm:
            res[key] = _f(mm.group(1))
    if res:
        meta["reported"] = res
    for c in curves:   # STARe: endothermic negative -> already exo-positive
        c["sign"] = +1
    return curves, meta


def read_trios(text):
    """TA TRIOS text export: [step] blocks, each with a header row, a units row and data. TA: exo up (exo +)."""
    L = text.replace("\r", "").split("\n")
    curves, i = [], 0
    while i < len(L):
        if L[i].strip().lower() == "[step]":
            name = L[i + 1].strip()
            j = i + 2
            while j < len(L) and not re.search(r"temp", L[j], re.I):
                j += 1
            hdr = [h.strip() for h in L[j].split("\t")]; units = [u.strip() for u in L[j + 1].split("\t")]
            a = _block(L, j + 2, 2)
            if a is not None:
                d = {h.lower(): k for k, h in enumerate(hdr)}
                kT = next(k for h, k in d.items() if "temp" in h)
                kq = next((k for h, k in d.items() if "heat flow" in h and "normal" in h), None)
                if kq is None:
                    kq = next(k for h, k in d.items() if "heat flow" in h)
                kt = next((k for h, k in d.items() if h.startswith("time")), None)
                km = next((k for h, k in d.items() if h.startswith("weight") and "%" not in units[k]), None)
                kmp = next((k for h, k in d.items() if h.startswith("weight") and "%" in units[k]), None)
                c = {"name": name, "T": a[:, kT], "q": a[:, kq], "unit": _unit_kind(units[kq]),
                     "t": a[:, kt] if kt is not None else None, "sign": +1}
                if kmp is not None:
                    c["w"] = a[:, kmp]
                elif km is not None:
                    c["w"] = 100 * a[:, km] / a[0, km]
                curves.append(c)
                i = j + 2 + len(a); continue
        i += 1
    return curves, {"vendor": "TA Instruments TRIOS", "mass_mg": _mass(text)}


def read_universal(text):
    """TA Universal Analysis export: 'SigN <name> (<unit>)' header lines, Size, StartOfData, then rows."""
    sigs = re.findall(r"^Sig\d+\s+(.+)$", text, re.M)
    body = text.split("StartOfData", 1)[1]
    a = _block(body.replace("\r", "").split("\n"), 0, 2)
    names = [s.lower() for s in sigs][: a.shape[1]]
    k = lambda pat: next((i for i, s in enumerate(names) if re.search(pat, s)), None)
    kT, kq, kt, kw = k(r"temp"), k(r"heat flow"), k(r"^time"), k(r"weight.*%")
    c = {"name": "UA", "T": a[:, kT], "q": a[:, kq], "unit": _unit_kind(sigs[kq]), "t": a[:, kt] if kt is not None else None, "sign": +1}
    if kw is not None:
        c["w"] = a[:, kw]
    return [c], {"vendor": "TA Universal Analysis", "mass_mg": _mass(text)}


def read_netzsch(text):
    """Netzsch Proteus ASCII: '#KEY:;value' header (SAMPLE MASS /mg, SEGMENTS), '##Temp./°C;Time/min;DSC/(mW/mg)[;Mass/%][;Segment]'.
    Netzsch DSC signal: 'exo down' instruments store endothermic positive -> flipped."""
    L = text.replace("\r", "").split("\n")
    j = next(i for i, l in enumerate(L) if l.startswith("##"))
    hdr = [h.strip() for h in L[j][2:].split(";")]
    a = _block(L, j + 1, 2)
    d = {h.lower(): k for k, h in enumerate(hdr)}
    g = lambda pat: next((k for h, k in d.items() if re.search(pat, h)), None)
    kT, kt, kq, km, ks = g(r"^temp"), g(r"^time"), g(r"^dsc|heat"), g(r"^mass|^tg"), g(r"segment")
    unit = _unit_kind(re.sub(r".*?/", "", hdr[kq], count=1).strip("()"))
    exo_down = not re.search(r"#?\s*EXO\s*[:;]\s*1|exo\s*up", text, re.I)
    base = {"q_all": a[:, kq], "T": a[:, kT], "t": a[:, kt] if kt is not None else None}
    seg = a[:, ks] if ks is not None else np.ones(len(a))
    curves = []
    for s in pd.unique(seg):
        m = seg == s
        c = {"name": f"segment {int(s)}", "T": a[m, kT], "q": a[m, kq], "unit": unit,
             "t": a[m, kt] if kt is not None else None, "sign": -1 if exo_down else +1}
        if km is not None:
            c["w"] = a[m, km]
        curves.append(c)
    return curves, {"vendor": "Netzsch Proteus", "mass_mg": _mass(text)}


def read_pyris(text):
    """PerkinElmer Pyris text: 'Sample Weight: 5.000 mg', header 'Time  Unsubtracted Heat Flow  ...  Sample Temperature'
    with units row '(min) (mW) ... (°C)'. Pyris: endothermic up (positive) -> flipped."""
    L = text.replace("\r", "").split("\n")
    j = next(i for i, l in enumerate(L) if re.search(r"heat\s*flow", l, re.I) and re.search(r"temp", l, re.I))
    hdr = [h.strip() for h in re.split(r"\t|\s{2,}", L[j].strip())]
    units = re.findall(r"\(([^)]*)\)", L[j + 1])
    a = _block(L, j + 2, 2)
    d = {h.lower(): k for k, h in enumerate(hdr)}
    g = lambda pat: next((k for h, k in d.items() if re.search(pat, h)), None)
    kT = g(r"sample\s*temp") if g(r"sample\s*temp") is not None else g(r"temp")
    kq, kt = g(r"heat\s*flow"), g(r"^time")
    unit = _unit_kind(units[kq] if kq < len(units) else "mW")
    return [{"name": "Pyris", "T": a[:, kT], "q": a[:, kq], "unit": unit, "t": a[:, kt] if kt is not None else None,
             "sign": -1}], {"vendor": "PerkinElmer Pyris", "mass_mg": _mass(text)}


def read_generic(src, text=None):
    """Any CSV/TXT/XLSX with temperature + heat-flow columns (+ time, mass %, segment/cycle). Sign assumed exo-up."""
    name = getattr(src, "name", "")
    if str(name).lower().endswith((".xlsx", ".xls")):
        raw = src.read() if hasattr(src, "read") else open(src, "rb").read()
        df = pd.read_excel(io.BytesIO(raw))
    else:
        from .readers import read_table
        df = read_table(io.StringIO(text))
    cT = find_col(df, r"^ts\b|sample\s*temp|temp")
    cq = find_col(df, r"heat\s*flow|dsc|^hf|w/g|mw")
    ct = find_col(df, r"^time|^t\b")
    cw = find_col(df, r"mass|weight|tg")
    cs = find_col(df, r"segment|cycle|step")
    if cT is None or cq is None:
        cT, cq = df.columns[0], df.columns[1]
    seg = df[cs].values if cs is not None else np.ones(len(df))
    curves = []
    for s in pd.unique(seg):
        m = seg == s
        c = {"name": f"segment {int(s) if float(s).is_integer() else s}" if cs is not None else "curve", "T": df.loc[m, cT].astype(float).values,
             "q": df.loc[m, cq].astype(float).values, "unit": _unit_kind(str(cq)),
             "t": df.loc[m, ct].astype(float).values if ct is not None else None, "sign": +1}
        if cw is not None and cw not in (cT, cq):
            c["w"] = df.loc[m, cw].astype(float).values
        curves.append(c)
    return curves, {"vendor": "generic table", "mass_mg": _mass(text or "")}


def read_dsc(src, name: str = "", split_segments=True, min_seg_points=20):
    """Detect the format, return dict(curves=[dict(name, T, q (exo +), unit, t (min) | None, w (%) | None, kind)], meta)."""
    name = name or getattr(src, "name", "")
    text = None if str(name).lower().endswith((".xlsx", ".xls")) else _to_text(src)
    if text is None:
        curves, meta = read_generic(src)
    elif "Curve Values" in text and "STARe" in text or re.search(r"^\s*Index\s+t\s+Ts", text, re.M):
        curves, meta = read_mettler(text)
    elif re.search(r"^\[step\]", text, re.M | re.I):
        curves, meta = read_trios(text)
    elif "StartOfData" in text:
        curves, meta = read_universal(text)
    elif re.search(r"^##", text, re.M):
        curves, meta = read_netzsch(text)
    elif re.search(r"Pyris|Sample Weight", text, re.I) and re.search(r"heat\s*flow", text, re.I):
        curves, meta = read_pyris(text)
    else:
        curves, meta = read_generic(src, text)
    out = []
    curves = make_continuous(curves)          # exports that restart t for every curve block
    for c in curves:
        q = np.asarray(c["q"], float) * c.pop("sign", 1)
        c["q"] = q
        c["T"] = np.asarray(c["T"], float)
        ok = np.isfinite(c["T"]) & np.isfinite(q)
        for k in ("T", "q", "t", "w", "Tr"):
            if c.get(k) is not None:
                c[k] = np.asarray(c[k], float)[ok]
        parts = split_by_rate(c, min_seg_points) if split_segments else [c]
        out.extend(parts)
    for k, c in enumerate(out):
        c["id"] = k
    meta["n_curves_in_file"] = len(curves)
    if text is not None and meta.get("timestamp") is None:
        meta["timestamp"] = file_timestamp(text)
    return {"curves": out, "meta": meta}


def heating_rate(c):
    """K/min from T(t), or None."""
    if c.get("t") is None or len(c["t"]) < 3 or np.ptp(c["t"]) <= 0:
        return None
    return float(np.polyfit(c["t"], c["T"], 1)[0])


def split_by_rate(c, min_pts=20):
    """Split one curve into heating / cooling / isothermal segments from dT (per point or per time)."""
    T = c["T"]; n = len(T)
    if n < 2 * min_pts:
        c["kind"] = _kind(c); return [c]
    x = c["t"] if c.get("t") is not None else np.arange(n, dtype=float)
    k = max(5, n // 200)
    dT = np.gradient(pd.Series(T).rolling(k, center=True, min_periods=1).mean().values, x)
    thr = 0.15 * np.nanmedian(np.abs(dT[np.abs(dT) > 0])) if np.any(np.abs(dT) > 0) else 0
    state = np.where(dT > thr, 1, np.where(dT < -thr, -1, 0))
    # merge short runs into neighbours
    cuts = [0] + [i for i in range(1, n) if state[i] != state[i - 1]] + [n]
    runs = [(cuts[i], cuts[i + 1], state[cuts[i]]) for i in range(len(cuts) - 1)]
    merged = []
    for a, b, s in runs:
        if merged and (b - a < min_pts or merged[-1][2] == s):
            merged[-1] = (merged[-1][0], b, merged[-1][2] if b - a < min_pts else s)
        else:
            merged.append((a, b, s))
    if len(merged) == 1:
        c["kind"] = _kind(c); return [c]
    # sharpen ramp <-> isothermal boundaries (the smoothing above smears them by ~k points)
    dTr = np.gradient(pd.Series(T).rolling(3, center=True, min_periods=1).mean().values, x)
    merged = [list(r) for r in merged]
    for j in range(len(merged) - 1):
        A, B = merged[j], merged[j + 1]
        if A[2] != 0 and B[2] == 0:          # ramp then hold: give the slow tail of the ramp to the hold
            rate = np.nanmedian(np.abs(dTr[A[0]:A[1]]))
            while A[1] - A[0] > min_pts and abs(dTr[A[1] - 1]) < 0.5 * rate:
                A[1] -= 1; B[0] -= 1
        elif A[2] == 0 and B[2] != 0:        # hold then ramp: give the slow start of the ramp to the hold
            rate = np.nanmedian(np.abs(dTr[B[0]:B[1]]))
            while B[1] - B[0] > min_pts and abs(dTr[B[0]]) < 0.5 * rate:
                B[0] += 1; A[1] += 1
    merged = [tuple(r) for r in merged]
    out = []
    for i, (a, b, s) in enumerate(merged):
        d = {kk: (v[a:b] if isinstance(v, np.ndarray) and len(v) == n else v) for kk, v in c.items()}
        d["kind"] = {1: "heating", -1: "cooling", 0: "isothermal"}[s]
        d["name"] = f"{c.get('name', 'curve')} · {i + 1} {d['kind']}"
        out.append(d)
    return out


def _kind(c):
    T = c["T"]
    if np.ptp(T) < 2:
        return "isothermal"
    return "heating" if T[-1] > T[0] else "cooling"


def to_specific(c, mass_mg=None):
    """Heat flow in W/g (= mW/mg). mW / W are divided by the sample mass."""
    q, u = c["q"], c["unit"]
    if u == "W/g":
        return q
    if u in ("mW", "W", "?") and mass_mg:
        return q * (1e3 if u == "W" else 1.0) / mass_mg
    return None


# ============================================================================ analysis
def _smooth(y, k=7):
    return pd.Series(y).rolling(k, center=True, min_periods=1).mean().values if k > 1 else y


def baseline_linear(T, q, i0, i1):
    """Straight line through the (locally averaged) signal at the two integration limits."""
    w = max(1, (i1 - i0) // 50)
    y0, y1 = np.mean(q[max(0, i0 - w):i0 + w + 1]), np.mean(q[max(0, i1 - w):i1 + w + 1])
    return y0 + (y1 - y0) * (T - T[i0]) / (T[i1] - T[i0] if T[i1] != T[i0] else 1)


def baseline_sigmoidal(T, q, i0, i1, iters=6):
    """Sigmoidal ('integral tangential') baseline: switches from the pre-peak tangent to the post-peak tangent in
    proportion to the reacted fraction α(T) (cumulative partial area), iterated to self-consistency."""
    n = len(T); w = max(5, (i1 - i0) // 5)
    a0, b0 = np.polyfit(T[max(0, i0 - w):i0 + 1], q[max(0, i0 - w):i0 + 1], 1) if i0 >= 3 else (0, q[i0])
    a1, b1 = np.polyfit(T[i1:min(n, i1 + w)], q[i1:min(n, i1 + w)], 1) if n - i1 >= 3 else (0, q[i1])
    pre, post = a0 * T + b0, a1 * T + b1
    base = baseline_linear(T, q, i0, i1)
    for _ in range(iters):
        s = np.clip(np.abs(q - base), 0, None)
        s[:i0] = 0; s[i1 + 1:] = 0
        cum = np.concatenate([[0], np.cumsum(0.5 * (s[1:] + s[:-1]) * np.diff(T))])
        alpha = cum / cum[-1] if cum[-1] > 0 else np.linspace(0, 1, n)
        base = (1 - alpha) * pre + alpha * post
    return base


def _tangent_cross(T, y, base, i, j):
    """Intersection of the inflectional tangent of y between indices i..j (steepest point) with the baseline."""
    if j - i < 3:
        return None, None
    d = np.gradient(_smooth(y, 3), T)
    k = i + int(np.argmax(np.abs(d[i:j])))
    m, c = d[k], y[k] - d[k] * T[k]
    # baseline locally linear around the crossing
    bm = np.gradient(base, T)[k]; bc = base[k] - bm * T[k]
    if abs(m - bm) < 1e-12:
        return None, k
    return float((bc - c) / (m - bm)), k


def find_peaks_dsc(c, prominence=None, max_peaks=6, smooth=5):
    """Candidate peaks on the (smoothed) signal, both endothermic (down, exo-up convention) and exothermic."""
    from scipy.signal import find_peaks
    q = _smooth(c["q"], smooth)
    span = np.nanpercentile(q, 99.5) - np.nanpercentile(q, 0.5)
    prom = prominence or max(0.05 * span, 1e-9)
    out = []
    for sgn, kind in ((-1, "endo"), (+1, "exo")):
        idx, pr = find_peaks(sgn * q, prominence=prom, distance=max(3, len(q) // 100))
        for i, p in zip(idx, pr["prominences"]):
            out.append((float(p), int(i), kind))
    out.sort(reverse=True)
    return [(i, k) for _, i, k in out[:max_peaks]]


def peak_limits(T, q, ip, frac=0.02, smooth=5):
    """Integration limits: walk out from the peak until the signal is back within `frac` of the peak height of the
    local linear background (estimated from points beyond the peak's flanks)."""
    qs = _smooth(q, smooth); n = len(q)
    h = qs[ip]
    e = max(3, n // 20)                               # robust end levels (ignore switching transients at the ends)
    q0, q1 = float(np.median(qs[:e])), float(np.median(qs[-e:]))
    # rough background from the nearest local level on both sides
    L = ip
    R = ip
    d = np.gradient(qs, T)
    ref = np.median(np.abs(d)) * 3 + 1e-15
    # leave the peak top first (a broad top is flat too): go below half height on each side
    while L > 1 and abs(qs[L] - q0) > 0.5 * abs(h - q0):
        L -= 1
    while R < n - 2 and abs(qs[R] - q1) > 0.5 * abs(h - q1):
        R += 1
    while L > 1 and not (abs(d[L]) < ref and abs(qs[L] - qs[0 if L < 5 else L - 5]) < frac * abs(h - q0) + 1e-15):
        L -= 1
    while R < n - 2 and not (abs(d[R]) < ref and abs(qs[R] - qs[min(n - 1, R + 5)]) < frac * abs(h - q1) + 1e-15):
        R += 1
    return L, R


def evaluate_peak(c, ip=None, Tlim=None, baseline="linear", mass_mg=None, rate=None):
    """Onset / peak / endset (°C) and ΔH (J/g) for one peak. ΔH is negative for endothermic events
    (STARe / ICTAC sign: heat absorbed by the sample) and integrated over TIME (s); if no time column, dt = dT/β."""
    T, q = c["T"], c["q"]
    spec = to_specific(c, mass_mg)
    y = spec if spec is not None else q
    if Tlim:
        i0, i1 = int(np.argmin(np.abs(T - Tlim[0]))), int(np.argmin(np.abs(T - Tlim[1])))
        i0, i1 = min(i0, i1), max(i0, i1)
        if ip is None:
            seg = y[i0:i1 + 1]; b = baseline_linear(T, y, i0, i1)[i0:i1 + 1]
            ip = i0 + int(np.argmax(np.abs(seg - b)))
    else:
        i0, i1 = peak_limits(T, y, ip)
    base = baseline_sigmoidal(T, y, i0, i1) if baseline.startswith("sigm") else baseline_linear(T, y, i0, i1)
    s = y - base
    if c.get("t") is not None:
        x = c["t"] * 60.0                      # s
    else:
        beta = rate or 10.0
        x = T / beta * 60.0
    sl = slice(i0, i1 + 1)
    area = float(np.trapezoid(s[sl], x[sl])) if hasattr(np, "trapezoid") else float(np.trapz(s[sl], x[sl]))
    kind = "endo" if s[ip] < 0 else "exo"
    onset, _ = _tangent_cross(T, y, base, i0, ip)
    endset, _ = _tangent_cross(T, y, base, ip, i1)
    return {"kind": kind, "i0": i0, "i1": i1, "ip": int(ip), "T_peak": float(T[ip]), "T_onset": onset, "T_endset": endset,
            "T_start": float(T[i0]), "T_end": float(T[i1]), "height": float(s[ip]),
            "dH": area if spec is not None else None, "area_mJ": area if spec is None else (area * mass_mg if mass_mg else None),
            "unit": "J/g" if spec is not None else "mJ", "baseline": base, "baseline_kind": baseline}


def evaluate_all(c, baseline="linear", mass_mg=None, prominence=None, max_peaks=6):
    pk = find_peaks_dsc(c, prominence, max_peaks)
    res = []
    for ip, _ in pk:
        try:
            r = evaluate_peak(c, ip, None, baseline, mass_mg)
            if r["i1"] - r["i0"] >= 4:
                res.append(r)
        except Exception:
            continue
    res.sort(key=lambda r: r["T_peak"])
    return res


def glass_transition(c, T1, T2, mass_mg=None, rate=None):
    """Tg in [T1, T2]: tangents on the first/last 25 % of the window; inflection = steepest point; midpoint = half
    step height between the extrapolated tangents (ISO 11357-2 / ASTM E1356); onset/endset = inflection tangent ∩
    tangents; Δcp = step / β."""
    T, q = c["T"], c["q"]
    spec = to_specific(c, mass_mg); y = spec if spec is not None else q
    m = (T >= min(T1, T2)) & (T <= max(T1, T2))
    Tw, yw = T[m], _smooth(y[m], 9)
    if len(Tw) < 20:
        return None
    n = len(Tw); a = n // 4
    p0 = np.polyfit(Tw[:a], yw[:a], 1); p1 = np.polyfit(Tw[-a:], yw[-a:], 1)
    d = np.gradient(_smooth(yw, 9), Tw)
    core = slice(a, n - a)
    k = a + int(np.argmax(np.abs(d[core])))
    Ti = float(Tw[k]); mi, ci = d[k], yw[k] - d[k] * Tw[k]
    l0, l1 = np.polyval(p0, Tw), np.polyval(p1, Tw)
    half = 0.5 * (l0 + l1)
    j = a + int(np.argmin(np.abs(yw[core] - half[core])))
    Tmid = float(Tw[j])
    onset = float((p0[1] - ci) / (mi - p0[0])) if mi != p0[0] else None
    endset = float((p1[1] - ci) / (mi - p1[0])) if mi != p1[0] else None
    step = float(np.polyval(p1, Tmid) - np.polyval(p0, Tmid))
    beta = rate or heating_rate(c) or 10.0
    dcp = abs(step) / (abs(beta) / 60.0) if spec is not None else None
    return {"Tg_mid": Tmid, "Tg_inflection": Ti, "Tg_onset": onset, "Tg_endset": endset, "step": step,
            "dcp_J_gK": dcp, "tangent0": p0, "tangent1": p1, "window": (float(Tw[0]), float(Tw[-1]))}


# ============================================================================ plotting
def draw_dsc(ax, items, x="T", exo_up=True, offset=0.0, colors=None, peaks=None, tgs=None, annotate=True,
             show_baseline=True, ylabel=None, arrow=True, lw=1.0, linestyles=None):
    """items: list of (label, curve, y) with y the heat-flow array to plot (exo-positive). Draws stacked curves,
    the 'exo ↑/↓' arrow, baselines and peak / Tg annotations."""
    import matplotlib.pyplot as plt
    from . import style
    sgn = 1 if exo_up else -1
    cols = colors or list(style.colours_for([lab for lab, _, _ in items]))
    for k, (lab, c, y) in enumerate(items):
        xx = c["T"] if x == "T" else (c["t"] if c.get("t") is not None else c["T"])
        off = k * offset
        col = cols[k] if cols else None
        ax.plot(xx, sgn * y + off, color=col, lw=lw, label=lab, ls=(linestyles[k] if linestyles else "-"))
        for r in (peaks or {}).get(c["id"], []):
            if show_baseline:
                ax.plot(xx[r["i0"]:r["i1"] + 1], sgn * r["baseline"][r["i0"]:r["i1"] + 1] + off, color=col, lw=0.6, ls="--")
            if annotate:
                yp = sgn * y[r["ip"]] + off
                yb = sgn * r["baseline"][r["ip"]] + off
                txt = f"{r['T_peak']:.1f} °C"
                if r.get("T_onset") is not None:
                    txt = f"onset {r['T_onset']:.1f} °C\npeak " + txt
                if r.get("dH") is not None:
                    txt += f"\nΔH = {r['dH']:.1f} J g$^{{-1}}$"
                elif r.get("area_mJ") is not None:
                    txt += f"\n{r['area_mJ']:.3g} mJ"
                if r.get("T_onset") is not None and x == "T":
                    ax.plot([r["T_onset"]], [sgn * np.interp(r["T_onset"], c["T"], r["baseline"]) + off], marker="|", ms=6,
                            color=col, mew=0.8)
                xr = np.nanmax(xx) - np.nanmin(xx)
                xa, xb = xx[r["i0"]], xx[min(r["i1"], len(xx) - 1)]
                xt = max(xa, xb) + 0.04 * xr                   # right of the whole peak (heating or cooling)
                if xt > np.nanmax(xx) - 0.25 * xr:            # no room on the right -> left side
                    xt, ha = min(xa, xb) - 0.04 * xr, "right"
                else:
                    ha = "left"
                t = ax.annotate(txt, (xx[r["ip"]], yp), xytext=(xt, 0.45 * yp + 0.55 * yb), textcoords="data",
                                fontsize=plt.rcParams["font.size"] - 1, va="center", ha=ha,
                                arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4", shrinkA=1, shrinkB=1))
                t.set_gid("mfs_note")
        for g in (tgs or {}).get(c["id"], []):
            yi = sgn * np.interp(g["Tg_mid"], c["T"], y) + off
            ax.plot([g["Tg_mid"]], [yi], marker="o", ms=3, color=col)
            if annotate:
                t = ax.annotate(f"T$_g$ = {g['Tg_mid']:.1f} °C", (g["Tg_mid"], yi), xytext=(6, 6), textcoords="offset points",
                                fontsize=plt.rcParams["font.size"] - 1)
                t.set_gid("mfs_note")
    ax.set_xlabel("Temperature (°C)" if x == "T" else "Time (min)")
    ax.set_ylabel(ylabel or "Heat flow (W g$^{-1}$)")
    if arrow:
        exo_arrow(ax, exo_up)
    if offset:
        ax.set_yticks([])
    return ax


def exo_arrow(ax, exo_up=True, loc=(0.04, 0.80)):
    """'exo ↑' or 'exo ↓' arrow in axes coordinates."""
    import matplotlib.pyplot as plt
    x, y = loc
    dy = 0.12 if exo_up else -0.12
    y0 = y - (0.06 if exo_up else -0.06)
    ax.annotate("", xy=(x, y0 + dy), xytext=(x, y0), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color="black", mutation_scale=7))
    t = ax.text(x + 0.02, y0 + dy / 2, "exo", transform=ax.transAxes, va="center", ha="left", fontsize=plt.rcParams["font.size"] - 1)
    t.set_gid("mfs_note")


def draw_sta(ax, c, y, exo_up=True, color_dsc=None, color_tga="#555555"):
    """TGA–DSC overlay (STA/SDT): heat flow on the left axis, mass % on a twin right axis."""
    ax.plot(c["T"], (1 if exo_up else -1) * y, color=color_dsc or "#0F766E", lw=1.0, label="DSC")
    ax.set_xlabel("Temperature (°C)"); ax.set_ylabel("Heat flow (W g$^{-1}$)")
    tw = ax.twinx()
    tw.plot(c["T"], c["w"], color=color_tga, lw=1.0, ls="-", label="TGA")
    tw.set_ylabel("Mass (%)")
    exo_arrow(ax, exo_up)
    return tw


def auto_tg_window(c, peaks=None, half=20.0):
    """Suggest a Tg window: the steepest step of the smoothed signal outside the integrated peaks, ± `half` K."""
    T, q = c["T"], _smooth(c["q"], 15)
    d = np.abs(_smooth(np.gradient(q, T), 15))
    m = np.isfinite(d)
    for r in peaks or []:
        m &= ~((T >= r["T_start"] - 5) & (T <= r["T_end"] + 5))
    e = int(len(T) * 0.05)
    m[:e] = False; m[len(T) - e:] = False
    if not m.any():
        return None
    k = int(np.argmax(np.where(m, d, -np.inf)))
    lo, hi = float(T[k] - half), float(T[k] + half)
    for r in peaks or []:                          # do not run into a peak
        if r["T_start"] > T[k]:
            hi = min(hi, r["T_start"] - 2)
        if r["T_end"] < T[k]:
            lo = max(lo, r["T_end"] + 2)
    return (max(lo, float(np.nanmin(T))), min(hi, float(np.nanmax(T))))


# ============================================================================ full heating / cooling cycle
SEG_PREFIX = {"heating": "Heat", "cooling": "Cool", "isothermal": "Iso"}
_DT_RX = [(r"(\d{1,2})[./](\d{1,2})[./](\d{4})[ T,]+(\d{1,2}):(\d{2})(?::(\d{2}))?", "dmy"),
          (r"(\d{4})-(\d{2})-(\d{2})[ T]+(\d{1,2}):(\d{2})(?::(\d{2}))?", "ymd")]


def file_timestamp(text):
    """First date + time found in the file header (export / run time), as a datetime, or None."""
    import datetime as _dt
    head = "\n".join((text or "").splitlines()[:40])
    for rx, kind in _DT_RX:
        m = re.search(rx, head)
        if m:
            g = [int(x) if x else 0 for x in m.groups()]
            try:
                if kind == "dmy":
                    return _dt.datetime(g[2], g[1], g[0], g[3], g[4], g[5])
                return _dt.datetime(g[0], g[1], g[2], g[3], g[4], g[5])
            except ValueError:
                continue
    return None


def segment_kind(c, iso_rate=0.05):
    """'heating' / 'cooling' / 'isothermal' from the sign of dT/dt (K/min) – or from T(end) - T(start) without time."""
    T = np.asarray(c["T"], float)
    if c.get("t") is not None and len(T) > 2 and np.ptp(c["t"]) > 0:
        rate = np.polyfit(c["t"], T, 1)[0]
        if abs(rate) < iso_rate or np.ptp(T) < 1.0:
            return "isothermal"
        return "heating" if rate > 0 else "cooling"
    return _kind(c)


def make_continuous(curves):
    """Make t (min) run on from one curve to the next when an export restarts the clock for every curve/file."""
    out, t_end, dt = [], None, None
    for c in curves:
        c = dict(c)
        t = c.get("t")
        if t is not None and len(t):
            t = np.asarray(t, float)
            step = float(np.nanmedian(np.diff(t))) if len(t) > 1 else 0.01
            if t_end is not None and t[0] <= t_end:
                t = t - t[0] + t_end + (dt or step)
            c["t"] = t; t_end = float(t[-1]); dt = step
        out.append(c)
    return out


def label_segments(curves):
    """Number the segments per kind in order: Heat 1, Cool 1, Heat 2, Iso 1 … (sets c['kind'], c['seg'])."""
    n = {}
    for c in curves:
        c["kind"] = segment_kind(c)
        n[c["kind"]] = n.get(c["kind"], 0) + 1
        c["seg"] = f"{SEG_PREFIX[c['kind']]} {n[c['kind']]}"
        c["seg_no"] = n[c["kind"]]
    return curves


def combine_cycle(datasets, order="time"):
    """Merge several reads of one sample (e.g. separate heating and cooling exports) into one cycle.
    datasets: list of dict(curves, meta[, timestamp]) as returned by read_dsc. order: 'time' (file time stamp, then
    the given order) or 'user' (as given). Returns the ordered, time-continuous, relabelled segment list."""
    idx = list(range(len(datasets)))
    if order == "time":
        def key(i):
            ts = datasets[i]["meta"].get("timestamp")
            return (0, ts, i) if ts is not None else (1, None, i)
        if all(datasets[i]["meta"].get("timestamp") is not None for i in idx):
            idx.sort(key=lambda i: (datasets[i]["meta"]["timestamp"], i))
    segs = []
    for i in idx:
        for c in datasets[i]["curves"]:
            c = dict(c); c["file_index"] = i
            segs.append(c)
    segs = make_continuous(segs)
    return label_segments(segs)


_GRP_WORDS = (r"heat(?:ing)?|cool(?:ing)?|full[ _-]?cycle|cycle|[hc]\d|seg(?:ment)?[ _-]?\d+|\d+(?:st|nd|rd|th)|run[ _-]?\d+|"
              r"\d+-curves|real|example|mettler|stare|trios|netzsch|proteus|pyris|ta|export")


def cycle_group_name(name):
    """Sample name without heating/cooling/segment/vendor words: 'SG_DSC_heating_Mettler_REAL' -> 'SG_DSC'."""
    g = re.sub(r"\.[A-Za-z0-9]{2,4}$", "", str(name))
    prev = None
    while prev != g:
        prev = g
        g = re.sub(r"(?i)(^|[ _\-.])(?:" + _GRP_WORDS + r")(?=$|[ _\-.])", r"\1", g)
        g = re.sub(r"[ _\-.]{2,}", "_", g).strip(" _-.")
    return g or str(name)


def main_peak(rows, kind):
    """Largest |ΔH| peak of the given kind ('endo' / 'exo') in a list of evaluate_all results."""
    rs = [r for r in rows or [] if r["kind"] == kind and r.get("dH") is not None]
    if not rs:
        rs = [r for r in rows or [] if r["kind"] == kind]
        return max(rs, key=lambda r: abs(r.get("height") or 0)) if rs else None
    return max(rs, key=lambda r: abs(r["dH"]))


def cycle_summary(segments, peaks):
    """One row per heating segment paired with the cooling segment that follows it (Heat n -> Cool n):
    melting (onset, peak, ΔHm) on heating, crystallisation (onset, peak, ΔHc) on cooling, supercooling ΔT = Tm − Tc
    (peak and onset based) and |ΔHc / ΔHm|. peaks: {segment uid: evaluate_all rows}."""
    rows = []
    groups = {}
    for c in segments:
        groups.setdefault(c.get("sample", ""), []).append(c)
    for smp, segs in groups.items():
        heats = [c for c in segs if c["kind"] == "heating"]
        cools = [c for c in segs if c["kind"] == "cooling"]
        for h in heats or [None]:
            # the first cooling after this heating (by position), else the cooling with the same number
            after = None
            if h is not None:
                pos = next(i for i, x in enumerate(segs) if x is h)
                after = next((c for c in segs[pos + 1:] if c["kind"] == "cooling"), None)
                if after is None:
                    after = next((c for c in cools if c["seg_no"] == h["seg_no"]), None)
            else:
                after = cools[0] if cools else None
            if h is None and after is None:
                continue
            m = main_peak(peaks.get(h["uid"]), "endo") if h is not None else None
            x = main_peak(peaks.get(after["uid"]), "exo") if after is not None else None
            r = {"sample": smp, "heating": h["seg"] if h is not None else "–", "cooling": after["seg"] if after is not None else "–",
                 "Tm onset (°C)": m and m["T_onset"], "Tm peak (°C)": m and m["T_peak"], "ΔHm (J/g)": m and m["dH"],
                 "Tc onset (°C)": x and x["T_onset"], "Tc peak (°C)": x and x["T_peak"], "ΔHc (J/g)": x and x["dH"]}
            r["ΔT = Tm − Tc, peak (K)"] = (m["T_peak"] - x["T_peak"]) if m and x else None
            r["ΔT, onset (K)"] = (m["T_onset"] - x["T_onset"]) if m and x and m["T_onset"] is not None and x["T_onset"] is not None else None
            r["|ΔHc/ΔHm|"] = abs(x["dH"] / m["dH"]) if m and x and m.get("dH") and x.get("dH") else None
            rows.append(r)
    return pd.DataFrame(rows)


DEFAULT_SEG_STYLE = {"heating": ("#C0392B", "-"), "cooling": ("#1F5FA8", "-"), "isothermal": ("#7F7F7F", ":")}
LINESTYLES = {"solid": "-", "dashed": "--", "dotted": ":", "dashdot": "-."}


def default_styles(samples):
    """{(sample, kind): (colour, linestyle)}. One sample: red heating / blue cooling. Several samples: the sample
    colour, solid for heating, dashed for cooling."""
    from . import style
    out = {}
    if len(samples) <= 1:
        for s in samples or [""]:
            for k, v in DEFAULT_SEG_STYLE.items():
                out[(s, k)] = v
        return out
    cols = style.colours_for(samples)
    for s, col in zip(samples, cols):
        out[(s, "heating")] = (col, "-"); out[(s, "cooling")] = (col, "--"); out[(s, "isothermal")] = ("#7F7F7F", ":")
    return out


def direction_arrows(ax, x, y, color, n=2, ls="-", size=8, lw=0.8):
    """Arrow heads along a curve, pointing in the direction of the measurement (index order). They sit on flat
    baseline stretches (no peak within ±3 % of the x range), spread over the segment."""
    from matplotlib.patches import FancyArrowPatch
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y); x, y = x[ok], y[ok]
    if len(x) < 10:
        return []
    xr = np.ptp(x) or 1.0
    yr = (np.nanpercentile(y, 99.5) - np.nanpercentile(y, 0.5)) or 1.0
    med = np.nanmedian(y)
    dev = pd.Series(np.abs(y - med)).rolling(max(3, int(0.06 * len(x))), center=True, min_periods=1).max().values
    flat = np.where(dev < 0.06 * yr)[0]
    out = []
    for f in np.linspace(0, 1, n + 2)[1:-1]:
        i = int(f * (len(x) - 1))
        if len(flat):
            i = int(flat[np.argmin(np.abs(flat - i))])
        j = i
        while j < len(x) - 1 and abs(x[j] - x[i]) < 0.015 * xr:
            j += 1
        if j == i:
            continue
        a = FancyArrowPatch((x[i], y[i]), (x[j], y[j]), arrowstyle="-|>", mutation_scale=size, color=color, lw=lw,
                            shrinkA=0, shrinkB=0, zorder=5)
        a.set_gid("mfs_dir_arrow")
        ax.add_patch(a); out.append(a)
    return out


def draw_cycle_T(ax, items, exo_up=True, styles=None, peaks=None, annotate=True, show_baseline=True, ylabel=None,
                 arrows=True, connect=True, lw=1.0):
    """Full cycle vs temperature: heating and cooling as one closed loop (consecutive segments joined), direction
    arrows on every segment, style per (sample, kind). items: (label, curve, y) in measurement order."""
    st_ = [(styles or {}).get((c.get("sample", ""), c["kind"]), DEFAULT_SEG_STYLE[c["kind"]]) for _, c, _ in items]
    draw_dsc(ax, items, "T", exo_up, 0.0, colors=[s[0] for s in st_], linestyles=[s[1] for s in st_], peaks=peaks,
             annotate=annotate, show_baseline=show_baseline, ylabel=ylabel, lw=lw)
    sgn = 1 if exo_up else -1
    prev = None
    for (lab, c, y), (col, ls) in zip(items, st_):
        if arrows:
            direction_arrows(ax, c["T"], sgn * np.asarray(y), col, n=2 if c["kind"] != "isothermal" else 1)
        if connect and prev is not None and prev[0].get("sample") == c.get("sample"):
            pc, py = prev
            ax.plot([pc["T"][-1], c["T"][0]], [sgn * py[-1], sgn * y[0]], color=col, lw=0.6 * lw, ls=":", zorder=1)
        prev = (c, np.asarray(y))
    return ax


def draw_cycle_time(ax, items, exo_up=True, styles=None, peaks=None, annotate=True, show_baseline=True, ylabel=None,
                    program="Ts", shade=True, lw=1.0):
    """Full cycle vs time: heat flow on the left axis, the temperature programme (Ts, or Tr if chosen and present)
    on a twin right axis, and shaded backgrounds per segment (heating / cooling / isothermal) with their labels."""
    import matplotlib.pyplot as plt
    st_ = [(styles or {}).get((c.get("sample", ""), c["kind"]), DEFAULT_SEG_STYLE[c["kind"]]) for _, c, _ in items]
    draw_dsc(ax, items, "t", exo_up, 0.0, colors=[s[0] for s in st_], linestyles=[s[1] for s in st_], peaks=peaks,
             annotate=annotate, show_baseline=show_baseline, ylabel=ylabel, lw=lw)
    shade_c = {"heating": "#C0392B", "cooling": "#1F5FA8", "isothermal": "#9E9E9E"}
    tw = ax.twinx()
    for k, (lab, c, y) in enumerate(items):
        t = c["t"]; T = c["Tr"] if program == "Tr" and c.get("Tr") is not None else c["T"]
        if shade:
            sp = ax.axvspan(t[0], t[-1], color=shade_c[c["kind"]], alpha=0.07, lw=0, zorder=0)
            sp.set_gid("mfs_shade")
            tx = ax.text(0.5 * (t[0] + t[-1]), 1.0, c.get("seg", c["kind"]), transform=ax.get_xaxis_transform(), ha="center",
                         va="bottom", fontsize=plt.rcParams["font.size"] - 2, color=shade_c[c["kind"]])
            tx.set_gid("mfs_note")
        tw.plot(t, T, color="0.35", lw=0.7 * lw, ls="--", label="temperature" if k == 0 else None)
    tw.set_ylabel("Temperature (°C)")
    return tw


def draw_cycle_stacked(axes, heat_items, cool_items, exo_up=True, styles=None, peaks=None, annotate=True,
                       show_baseline=True, ylabel=None, arrows=True, lw=1.0):
    """Stacked panel: heating segments in the upper axes, cooling segments in the lower one (shared T axis)."""
    for ax, its in zip(axes, (heat_items, cool_items)):
        st_ = [(styles or {}).get((c.get("sample", ""), c["kind"]), DEFAULT_SEG_STYLE[c["kind"]]) for _, c, _ in its]
        if its:
            draw_dsc(ax, its, "T", exo_up, 0.0, colors=[s[0] for s in st_], linestyles=[s[1] for s in st_], peaks=peaks,
                     annotate=annotate, show_baseline=show_baseline, ylabel=ylabel, lw=lw)
            if arrows:
                sgn = 1 if exo_up else -1
                for (lab, c, y), (col, ls) in zip(its, st_):
                    direction_arrows(ax, c["T"], sgn * np.asarray(y), col, n=1)
        else:
            ax.text(0.5, 0.5, "no segment selected", transform=ax.transAxes, ha="center", va="center", color="0.5")
    axes[0].set_xlabel("")
    axes[0].set_title("heating", fontsize=8, loc="left"); axes[1].set_title("cooling", fontsize=8, loc="left")
    return axes


# ------------------------------------------------------------------ cycle-to-cycle comparison
def cycle_compare(segments, peaks):
    """One row per segment of each kind (Heat 1, Heat 2 …; Cool 1, Cool 2 …) with its main peak and the change vs the
    first cycle of that kind: ΔT_peak, ΔT_onset, ΔH ratio (e.g. melting enthalpy lost after the first heating)."""
    rows, first = [], {}
    for c in segments:
        pk = peaks.get(c.get("uid", c.get("id"))) if isinstance(peaks, dict) else None
        kind = c.get("kind")
        if kind not in ("heating", "cooling"):
            continue
        r = main_peak(pk, "endo" if kind == "heating" else "exo")
        row = {"sample": c.get("sample", ""), "segment": c.get("seg", c.get("name")), "kind": kind,
               "rate (K/min)": c.get("rate") if c.get("rate") is not None else heating_rate(c), "T onset (°C)": r and r.get("T_onset"), "T peak (°C)": r and r.get("T_peak"),
               "ΔH (J/g)": r and r.get("dH")}
        key = (row["sample"], kind)
        if key not in first and r:
            first[key] = row
        f = first.get(key)
        if f and r and f is not row:
            row["ΔT peak vs cycle 1 (K)"] = r["T_peak"] - f["T peak (°C)"]
            row["ΔT onset vs cycle 1 (K)"] = (r["T_onset"] - f["T onset (°C)"]) if r.get("T_onset") is not None and f["T onset (°C)"] is not None else None
            row["ΔH / ΔH(cycle 1)"] = (r["dH"] / f["ΔH (J/g)"]) if r.get("dH") and f["ΔH (J/g)"] else None
        rows.append(row)
    return rows


# ------------------------------------------------------------------ Cp with a sapphire reference (ASTM E1269)
# synthetic sapphire α-Al2O3, J g-1 K-1 (ASTM E1269 Table 1 to 600 K; Ginnings & Furukawa / TA TN-8 above)
SAPPHIRE_CP = np.array([
    (150, 0.3133), (170, 0.3912), (190, 0.4659), (200, 0.5014), (220, 0.5684), (240, 0.6294), (250, 0.6577), (260, 0.6846),
    (270, 0.7102), (280, 0.7344), (290, 0.7574), (300, 0.7792), (310, 0.7999), (320, 0.8194), (330, 0.8380), (340, 0.8556),
    (350, 0.8721), (360, 0.8878), (370, 0.9027), (380, 0.9168), (390, 0.9302), (400, 0.9429), (420, 0.9666), (440, 0.9879),
    (460, 1.0074), (480, 1.0250), (500, 1.0411), (520, 1.0559), (540, 1.0694), (560, 1.0819), (580, 1.0934), (600, 1.1040),
    (620, 1.1137), (650, 1.1271), (700, 1.1467), (750, 1.1636), (800, 1.1783), (850, 1.1913), (900, 1.2030), (1000, 1.2237)])


def cp_sapphire(T_C):
    T = np.asarray(T_C, float) + 273.15
    if np.any(T < SAPPHIRE_CP[0, 0]) or np.any(T > SAPPHIRE_CP[-1, 0]):
        raise ValueError("sapphire Cp table covers 150–1000 K (−123 to 727 °C)")
    return np.interp(T, SAPPHIRE_CP[:, 0], SAPPHIRE_CP[:, 1])


def _heating_part(c):
    T = np.asarray(c["T"], float); q = np.asarray(c["q"], float) * (1e3 if c.get("unit") == "W" else 1.0)
    o = np.argsort(T)
    return T[o], q[o]


def cp_from_sapphire(sample, sapphire, blank, m_sample, m_sapphire, T_range=None, step=1.0):
    """Ratio method (ASTM E1269): Cp(T) = Cp_sap(T) · (q_s − q_b)/(q_sap − q_b) · m_sap/m_s.
    sample, sapphire, blank: heating curves (dicts with T, q in mW) at the same rate and purge. Isothermal parts must be
    removed (use the heating segment). Returns T, Cp, and the direct estimate Cp = (q_s − q_b)/(m_s β) for comparison."""
    Ts, qs = _heating_part(sample); Ta, qa = _heating_part(sapphire); Tb, qb = _heating_part(blank)
    lo, hi = max(Ts.min(), Ta.min(), Tb.min()), min(Ts.max(), Ta.max(), Tb.max())
    if T_range:
        lo, hi = max(lo, T_range[0]), min(hi, T_range[1])
    lo += 5; hi -= 2          # skip the start-up transient
    T = np.arange(lo, hi, step)
    if len(T) < 5:
        raise ValueError("no common temperature range")
    s, a, b = np.interp(T, Ts, qs), np.interp(T, Ta, qa), np.interp(T, Tb, qb)
    den = a - b
    cp = cp_sapphire(T) * (s - b) / den * (m_sapphire / m_sample)
    out = {"T": T, "Cp": cp, "Cp_sapphire": cp_sapphire(T), "method": "ratio vs sapphire (ASTM E1269)"}
    rate = heating_rate({"T": Ts, "t": sample.get("t")}) if sample.get("t") is not None else None
    if rate:
        out["Cp_direct"] = np.abs(s - b) / m_sample / (abs(rate) / 60.0)       # mW/mg / (K/s) = J/(g K)
        out["rate"] = rate
        out["Cp_sapphire_check"] = np.abs(a - b) / m_sapphire / (abs(rate) / 60.0)
    return out


def draw_cp(ax, results, colors=None, show_ref=False):
    cols = colors or ["#C0392B", "#1F5FA8", "#009E73", "#CC79A7", "#E69F00"]
    for k, (nm, r) in enumerate(results):
        ax.plot(r["T"], r["Cp"], color=cols[k % len(cols)], lw=1.0, label=nm)
        if show_ref and k == 0:
            ax.plot(r["T"], r["Cp_sapphire"], color="0.5", lw=0.8, ls="--", label="sapphire (ASTM E1269)")
    ax.set_xlabel("Temperature (°C)"); ax.set_ylabel(r"$c_p$ (J g$^{-1}$ K$^{-1}$)")
    ax.legend(frameon=False)
