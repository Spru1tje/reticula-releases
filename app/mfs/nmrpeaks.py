"""v2.13 NMR peak picking, integration, first-order multiplet analysis and the ACS characterisation string."""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field, asdict

import numpy as np

MULT_NAMES = {1: "s", 2: "d", 3: "t", 4: "q", 5: "p", 6: "sext", 7: "sept"}
PATTERNS = [(), (2,), (3,), (4,), (5,), (6,), (7,), (2, 2), (2, 3), (3, 2), (2, 4), (4, 2), (3, 3), (2, 2, 2), (2, 2, 3),
            (2, 3, 2), (3, 2, 2)]


def mad_sigma(y):
    y = np.asarray(y)
    return 1.4826 * np.median(np.abs(y - np.median(y))) + 1e-12


@dataclass
class Peak:
    ppm: float
    height: float
    hz: float = 0.0


@dataclass
class Multiplet:
    lo: float                     # region (ppm), lo < hi
    hi: float
    centre: float
    integral: float = 0.0
    nH: float | None = None
    mult: str = "m"
    J: list = field(default_factory=list)       # Hz, largest first
    lines: list = field(default_factory=list)   # ppm of the lines used
    fit: float = 0.0                            # 0 (perfect) … 1
    label: str = ""                             # assignment label (atom numbers), set by nmrassign
    note: str = ""                              # e.g. 'CHCl3 (solvent)', 'H2O'
    manual: bool = False
    overlap: str = ""                           # compound multiplet overlapping a solvent/H2O line
    nH_assigned: float | None = None            # H count of the atoms assigned here (set by nmrassign)

    def as_dict(self):
        return asdict(self)


def pick_peaks(ppm, y, k=6.0, min_rel=0.004, sf=400.0):
    """Line positions (parabolic interpolation) of all maxima above k·σ and min_rel·max."""
    from scipy.signal import find_peaks
    s = mad_sigma(y)
    hz_pt = abs(ppm[1] - ppm[0]) * sf
    dist = max(1, int(0.25 / max(hz_pt, 1e-6)))
    pk, pr = find_peaks(y, height=max(k * s, min_rel * np.max(y)), distance=dist, prominence=k * s * 0.7)
    out = []
    for i in pk:
        p = float(ppm[i])
        if 0 < i < len(y) - 1:
            a, b, c = y[i - 1], y[i], y[i + 1]
            den = a - 2 * b + c
            if den:
                p = float(ppm[i] + 0.5 * (a - c) / den * (ppm[i + 1] - ppm[i]))
        out.append(Peak(p, float(y[i]), p * sf))
    return out


def regions(ppm, y, peaks, sf, gap_hz=None, nucleus="1H"):
    """Group lines into multiplets: neighbours closer than gap_hz (default 1H 18 Hz) with signal between them
    (the trace does not return to the noise floor). Returns [(lo, hi, [peaks])] in ppm, high δ first."""
    if not peaks:
        return []
    gap_hz = gap_hz or (18.0 if nucleus == "1H" else 3.0)
    s = mad_sigma(y)
    pk = sorted(peaks, key=lambda p: -p.ppm)
    groups = [[pk[0]]]
    for p in pk[1:]:
        q = groups[-1][-1]
        close = (q.ppm - p.ppm) * sf < gap_hz
        if close:
            m = (ppm < q.ppm) & (ppm > p.ppm)
            dip = y[m].min() if m.any() else 0
            close = dip > max(3 * s, 0.02 * min(p.height, q.height))
            close = close or (q.ppm - p.ppm) * sf < 0.6 * gap_hz and dip > 1.5 * s
        if close:
            groups[-1].append(p)
        else:
            groups.append([p])
    out = []
    for g in groups:
        hi, lo = g[0].ppm, g[-1].ppm
        mx = max(p.height for p in g)
        thr = max(2 * s, 0.01 * mx)
        i_hi = int(np.argmin(abs(ppm - hi))); i_lo = int(np.argmin(abs(ppm - lo)))
        stop = 0.08 * mx                                  # also stop in a valley (never climb into a neighbour)
        while i_hi > 0 and y[i_hi] > thr and not (y[i_hi - 1] > y[i_hi] and y[i_hi] < stop):
            i_hi -= 1
        while i_lo < len(y) - 1 and y[i_lo] > thr and not (y[i_lo + 1] > y[i_lo] and y[i_lo] < stop):
            i_lo += 1
        pad = (1.0 if nucleus == "1H" else 2.0) / sf
        out.append((float(ppm[i_lo]) - pad, float(ppm[i_hi]) + pad, g))
    # merge overlapping windows
    out.sort(key=lambda r: -r[1])
    merged = []
    for lo, hi, g in out:
        if merged and hi >= merged[-1][0]:
            plo, phi, pg = merged[-1]
            merged[-1] = (min(lo, plo), max(hi, phi), pg + g)
        else:
            merged.append((lo, hi, g))
    return merged


def integrate(ppm, y, lo, hi):
    """Area between lo and hi (ppm) above a straight line through the region ends."""
    m = (ppm >= lo) & (ppm <= hi)
    if m.sum() < 3:
        return 0.0
    x, v = ppm[m], y[m]
    k = max(1, min(5, len(v) // 30))
    b0, b1 = np.mean(v[:k]), np.mean(v[-k:])
    vmax = np.max(np.abs(v))
    if max(abs(b0), abs(b1)) > 0.03 * vmax:            # an end sits on a line tail: trust the corrected baseline (0)
        b0 = b1 = 0.0
    base = b0 + (b1 - b0) * (np.arange(len(v)) / max(1, len(v) - 1))
    return float(np.abs(np.trapezoid(v - base, x)))


def _sticks(pattern, Js, merge_hz=0.35):
    pos = np.array([0.0]); inten = np.array([1.0])
    for n, J in zip(pattern, Js):
        offs = (np.arange(n) - (n - 1) / 2) * J
        from math import comb
        w = np.array([comb(n - 1, i) for i in range(n)], float)
        pos = (pos[:, None] + offs[None, :]).ravel(); inten = (inten[:, None] * w[None, :]).ravel()
    order = np.argsort(pos); pos, inten = pos[order], inten[order]
    mp, mi = [], []                                           # merge coincident lines
    for p, i in zip(pos, inten):
        if mp and abs(p - mp[-1]) < merge_hz:
            mi[-1] += i; mp[-1] = (mp[-1] * mi[-1] + p * i) / (mi[-1] + i) if merge_hz > 0.35 else (mp[-1] + p) / 2
        else:
            mp.append(p); mi.append(i)
    return np.array(mp), np.array(mi)


def _merge_close(x, h, min_hz=0.55):
    o = np.argsort(x); x, h = list(np.asarray(x)[o]), list(np.asarray(h)[o])
    i = 0
    while i < len(x) - 1:
        if x[i + 1] - x[i] < min_hz:
            w = h[i] + h[i + 1]; x[i] = (x[i] * h[i] + x[i + 1] * h[i + 1]) / w; h[i] = max(h[i], h[i + 1])
            del x[i + 1]; del h[i + 1]
        else:
            i += 1
    return np.array(x), np.array(h)


def _score(x, hr, ps, ins, tol):
    """Match predicted sticks (ps, ins normalised to max 1) to observed lines (x, hr normalised to max 1)."""
    used = np.zeros(len(x), bool); pe = 0.0; ie = 0.0; miss = 0.0
    sc = np.sum(hr[hr > 0.2]) / max(np.sum(ins), 1e-9)
    for p, I in zip(ps, ins):
        d = np.abs(x - p); d[used] = np.inf
        j = int(np.argmin(d))
        if d[j] < 2.5 * tol:
            used[j] = True; pe += (d[j] / tol) ** 2; ie += abs(hr[j] - I * sc)
        else:
            miss += 1.5 * I
    extra = np.sum(np.where(hr[~used] > 0.35, 1.5 * hr[~used], 0.4 * hr[~used]))
    n = len(ps)
    return 0.5 * np.sqrt(pe / n) + 0.8 * ie / n + miss + extra


def analyse_multiplet(lines_hz, heights, tol_hz=0.6, max_J=20.0):
    """First-order analysis of one multiplet. lines_hz, heights: observed lines (Hz). Weak lines that do not belong to
    the pattern (impurities, 13C satellites, noise) cost little; lines closer than 0.55 Hz are merged.
    Returns (name, [J largest first], score 0 good … 1 bad)."""
    x, h = _merge_close(lines_hz, heights)
    if len(x) == 0:
        return "m", [], 1.0
    hr = h / h.max()
    strong = hr > 0.2
    if strong.sum() == 1:
        return "s", [], float(min(1.0, 0.3 * hr[~strong].sum()))
    cen = float(np.sum(x[strong] * hr[strong]) / np.sum(hr[strong]))
    xs = x[strong]
    spc = sorted({round(abs(a - b), 2) for a, b in itertools.combinations(x[hr > 0.08], 2) if 0.5 < abs(a - b) <= max_J})
    cands = []
    for v in spc:
        if not cands or v - cands[-1] > 0.3:
            cands.append(v)
        else:
            cands[-1] = (cands[-1] + v) / 2
    if len(cands) > 10:                                  # keep spacings between strong lines first
        st = sorted({round(abs(a - b), 2) for a, b in itertools.combinations(xs, 2) if 0.5 < abs(a - b) <= max_J})
        cands = sorted(set([c for c in cands if any(abs(c - t) < 0.3 for t in st)][:10]))
    best = ("m", [], 9.0)
    nstrong = int(strong.sum())
    merge_hz = 0.35
    for pat in PATTERNS[1:]:
        k = len(pat)
        nl = int(np.prod(pat))
        if nl < nstrong * 0.7 or nl > len(x) * 1.5 + 2:
            continue
        for Js in (itertools.combinations_with_replacement(cands, k) if k > 1 else ((c,) for c in cands)):
            for perm in set(itertools.permutations(Js)):
                ps, ins = _sticks(pat, perm, merge_hz)
                ins = ins / ins.max()
                for dc in (0.0,):
                    sc = _score(x, hr, ps + cen + dc, ins, tol_hz) + 0.08 * (k - 1)
                    if sc < best[2]:
                        order = np.argsort(perm)[::-1]
                        best = ("".join(MULT_NAMES[pat[i]] for i in order), [round(float(perm[i]), 1) for i in order], sc)
    if best[2] > 1.0:
        return "m", [], 1.0
    return best[0], best[1], float(best[2])


def _split_chain(sp, lines, max_seg=14):
    """Dynamic-programming split of a chain of lines into contiguous first-order multiplets. A cut is only allowed where
    the trace dips below 30 % of the smaller neighbouring line; each group costs 0.35 + its fit score, an unresolved
    group ('m') costs 1.0 + 0.35. Returns list of line lists."""
    n = len(lines)
    if n <= 2:
        return [lines]
    ls = sorted(lines, key=lambda p: -p.ppm)
    cut_ok = [False] * n                                   # cut after line i
    for i in range(n - 1):
        a, b = ls[i], ls[i + 1]
        m = (sp.ppm < a.ppm) & (sp.ppm > b.ppm)
        dip = sp.y[m].min() if m.any() else min(a.height, b.height)
        cut_ok[i] = dip < 0.30 * min(a.height, b.height)
    cache = {}

    def seg_cost(i, j):                                     # lines i..j-1
        if (i, j) in cache:
            return cache[(i, j)]
        g = ls[i:j]
        if j - i > max_seg:
            c = (1.35, "m")
        else:
            nm, J, fit = analyse_multiplet([p.ppm * sp.sf for p in g], [p.height for p in g])
            c = (1.35, nm) if nm == "m" else (0.35 + fit, nm)
        cache[(i, j)] = c
        return c
    best = [0.0] + [np.inf] * n; prev = [0] * (n + 1)
    for j in range(1, n + 1):
        for i in range(j - 1, -1, -1):
            if i > 0 and not cut_ok[i - 1]:
                continue
            if j - i > 3 * max_seg:
                break
            c = best[i] + seg_cost(i, j)[0]
            if c < best[j]:
                best[j] = c; prev[j] = i
    out = []; j = n
    while j > 0:
        i = prev[j]; out.append(ls[i:j]); j = i
    return out[::-1]


def analyse(sp, k=6.0, solvent=None, nucleus=None, split=True):
    """Full 1H (or 13C) peak/multiplet analysis of a processed Spectrum (mfs.nmrproc.Spectrum)."""
    nucleus = nucleus or sp.nucleus
    peaks = pick_peaks(sp.ppm, sp.y, k=k, sf=sp.sf)
    major, minor = split_major_minor(peaks, sp, solvent)
    out = []
    groups = []
    R_major = regions(sp.ppm, sp.y, major, sp.sf, nucleus=nucleus)
    minor = [p for p in minor if not any(lo <= p.ppm <= hi for lo, hi, _ in R_major)]
    R_minor = []
    for lo, hi, g in regions(sp.ppm, sp.y, minor, sp.sf, nucleus=nucleus):
        for mlo, mhi, _ in R_major:                       # never overlap a major multiplet
            if lo < mhi and hi > mlo:
                if g[0].ppm > mhi:
                    lo = max(lo, mhi)
                else:
                    hi = min(hi, mlo)
        if hi > lo:
            R_minor.append((lo, hi, g))
    for lo, hi, g in R_major + R_minor:
        lines = [p for p in g if p.height > 0.04 * max(q.height for q in g)]
        whole = analyse_multiplet([p.ppm * sp.sf for p in lines], [p.height for p in lines]) if len(lines) <= 16 else ("m", [], 1)
        sv = _solvent_1h(sp, solvent) if nucleus == "1H" else None
        if sv is not None and lo - 0.01 <= sv <= hi + 0.01 and len(lines) >= 3 and (whole[0] == "m" or whole[2] > 0.5):
            # a compound multiplet on top of the residual solvent line: damp the line at the solvent shift and retry
            j = int(np.argmin([abs(p.ppm - sv) for p in lines]))
            if abs(lines[j].ppm - sv) < 0.006:
                hh = [p.height for p in lines]
                hh[j] = float(np.median([x for i, x in enumerate(hh) if i != j]))
                alt = analyse_multiplet([p.ppm * sp.sf for p in lines], hh)
                if alt[0] not in ("m", "s") and alt[2] <= 0.5 and alt[1] and max(alt[1]) >= 5.0:
                    groups.append((lo, hi, lines, alt))
                    continue
        if split and nucleus == "1H" and len(lines) > 2 and (whole[0] == "m" or whole[2] > 0.5):
            parts = _split_chain(sp, lines, max_seg=10)
        else:
            parts = [lines]
        if len(parts) == 1:
            groups.append((lo, hi, parts[0]))
            continue
        for gi, part in enumerate(parts):                   # boundaries at the deepest point between groups
            hi_ = hi if gi == 0 else _min_between(sp, parts[gi - 1][-1].ppm, part[0].ppm)
            lo_ = lo if gi == len(parts) - 1 else _min_between(sp, part[-1].ppm, parts[gi + 1][0].ppm)
            groups.append((lo_, hi_, part))
    for grp in groups:
        lo, hi, lines = grp[:3]
        hz = [p.ppm * sp.sf for p in lines]
        if len(grp) == 4:
            name, J, fit = grp[3]
        elif nucleus == "1H":
            name, J, fit = analyse_multiplet(hz, [p.height for p in lines])
        else:
            name, J, fit = "s", [], 0.0
        cen = float(np.sum([p.ppm * p.height for p in lines]) / np.sum([p.height for p in lines]))
        if name == "m":
            cen = (lo + hi) / 2
        width = (max(hz) - min(hz)) if len(hz) > 1 else 0
        if name == "s" and nucleus == "1H":
            hw = _half_width_hz(sp, lines[0].ppm)
            if hw > 6:
                name = "br s"
        mm = Multiplet(lo=lo, hi=hi, centre=cen, integral=integrate(sp.ppm, sp.y, lo, hi), mult=name, J=J,
                       lines=[p.ppm for p in lines], fit=fit)
        if len(grp) == 4:
            mm.overlap = "solvent residual (integral includes it)"
            mm.centre = float(np.mean([p.ppm for p in lines]))
        out.append(mm)
    out.sort(key=lambda m: -m.centre)
    tag_solvent(out, sp, solvent)
    merged = []                                       # neighbouring pieces with the same solvent tag → one region
    for m in out:
        p = merged[-1] if merged else None
        if p is not None and m.note and m.note == p.note and abs(p.lo - m.hi) < 0.01:
            p.lo = min(p.lo, m.lo); p.integral += m.integral; p.lines += m.lines; p.mult = "m"; p.J = []
            p.centre = (p.lo + p.hi) / 2
        else:
            merged.append(m)
    return merged


def split_major_minor(peaks, sp, solvent=None, rel=0.012):
    """Major lines (≥ rel × reference height; reference = strongest non-solvent line) build the multiplets.
    Minor lines: 13C satellites (pairs ±J(CH)/2 = 55–90 Hz around a major line, ≈0.55 % each) are dropped, the rest are
    kept as minor signals (impurity candidates)."""
    from mfs import nmrproc as NP
    if not peaks:
        return [], []
    key = NP.solvent_key(solvent or sp.info.get("solvent"))
    avoid = [0.0]
    if key and NP.SOLVENTS[key].get(sp.nucleus) is not None:
        avoid.append(NP.SOLVENTS[key][sp.nucleus])
    avoid += [{"CDCl3": 1.56, "DMSO-d6": 3.33, "CD3OD": 4.87, "Acetone-d6": 2.84, "D2O": 4.79, "C6D6": 0.40,
               "CD3CN": 2.13, "CD2Cl2": 1.52}.get(key, -99)]
    hs = [p.height for p in peaks if min(abs(p.ppm - a) for a in avoid) > 0.08]
    ref = max(hs) if hs else max(p.height for p in peaks)
    major = [p for p in peaks if p.height >= rel * ref]
    minor = [p for p in peaks if p.height < rel * ref]
    keep = []
    for p in minor:
        sat = False
        for q in major:
            dhz = abs(p.ppm - q.ppm) * sp.sf
            if 55 <= dhz <= 90 and 0.002 < p.height / q.height < 0.012:
                mirror = q.ppm + (q.ppm - p.ppm)
                if any(abs(r.ppm - mirror) * sp.sf < 2.0 for r in minor) or abs(mirror - q.ppm) > 0:
                    sat = True; break
        if not sat:
            keep.append(p)
    return major, keep


def _min_between(sp, a, b):
    m = (sp.ppm < max(a, b)) & (sp.ppm > min(a, b))
    if not m.any():
        return (a + b) / 2
    idx = np.where(m)[0]
    return float(sp.ppm[idx[np.argmin(sp.y[idx])]])


def _half_width_hz(sp, p):
    i = int(np.argmin(abs(sp.ppm - p))); h = sp.y[i]
    a = i
    while a > 0 and sp.y[a] > h / 2:
        a -= 1
    b = i
    while b < len(sp.y) - 1 and sp.y[b] > h / 2:
        b += 1
    return abs(sp.ppm[a] - sp.ppm[b]) * sp.sf


def _solvent_1h(sp, solvent=None):
    from mfs import nmrproc as NP
    key = NP.solvent_key(solvent or sp.info.get("solvent"))
    return NP.SOLVENTS[key].get("1H") if key else None


def tag_solvent(mults, sp, solvent=None):
    from mfs import nmrproc as NP
    key = NP.solvent_key(solvent or sp.info.get("solvent"))
    marks = []
    if key and NP.SOLVENTS[key].get(sp.nucleus) is not None:
        marks.append((NP.SOLVENTS[key][sp.nucleus], f"{key} (solvent)"))
    if sp.nucleus == "1H":
        marks.append((0.0, "TMS"))
        water = {"CDCl3": 1.56, "DMSO-d6": 3.33, "CD3OD": 4.87, "Acetone-d6": 2.84, "C6D6": 0.40, "CD3CN": 2.13,
                 "CD2Cl2": 1.52, "D2O": 4.79, "THF-d8": 2.46, "Toluene-d8": 0.43}.get(key)
        if water:
            marks.append((water, "H2O"))
    for m in mults:
        if m.overlap:
            continue
        for v, nm in marks:
            tol = 0.03 if sp.nucleus == "1H" else 0.3
            if m.lo - tol <= v <= m.hi + tol and abs(m.centre - v) < (0.06 if sp.nucleus == "1H" else 0.6):
                resolved = (sp.nucleus == "1H" and m.mult not in ("s", "br s", "m") and m.J
                            and max(m.J) >= 5.0 and m.fit < 0.5)
                if resolved:                              # a coupled compound multiplet sitting on the solvent line
                    m.overlap = nm
                else:
                    m.note = nm


def normalise(mults, total_h=None, ref_index=None, ref_n=1.0):
    """Integrals → number of H. With a reference multiplet (index, nH), a known total H count, or by the best integer
    fit over the significant multiplets (weighted by integral). Solvent/TMS/H2O multiplets are skipped; multiplets below
    0.5 H are kept with their fraction (minor / impurity candidates)."""
    use = [m for m in mults if not m.note]
    if not use:
        return mults
    ints = np.array([m.integral for m in use])
    sig = ints >= 0.06 * ints.max()
    if ref_index is not None:
        unit = mults[ref_index].integral / ref_n
    elif total_h:
        unit = ints[sig].sum() / total_h
    else:
        best = None
        top = np.sort(ints[sig])[::-1][:10]
        for base in top:
            for k in (1, 2, 3, 4, 6, 9):
                u = base / k
                r = ints[sig] / u
                if r.min() < 0.25 or r.max() > 60:
                    continue
                w = ints[sig] / ints[sig].sum()
                dev = np.where(r < 0.75, 0.35, np.abs(r - np.round(r)))   # < 0.75 H: minor species (cheap but not free)
                err = np.sum(w * dev) + 0.05 * np.log(r.sum())
                if best is None or err < best[0]:
                    best = (err, u)
        unit = best[1] if best else ints[sig].min()
    for m in mults:
        m.nH = round(m.integral / unit, 2) if unit else None
        if (m.note and m.note.split(" ")[0] in ("CDCl3", "C6D6", "D2O", "CD2Cl2") and "solvent" in m.note
                and m.nH and m.nH >= 0.5 and m.mult not in ("s", "br s")):        # residual normally a singlet here
            m.overlap = "coupled lines at the solvent shift: a compound signal may be hidden here"
    return mults


def _fmt_d(v, dec=2):
    return f"{v:.{dec}f}"


def acs_string(mults, nucleus="1H", mhz=400.0, solvent="CDCl3", dec=None, include_notes=False, labels=True, skip=(),
               min_nH=None):
    """ACS-style characterisation string, e.g.
    ¹H NMR (400 MHz, CDCl3) δ 9.09 (d, J = 8.7 Hz, 1H), 7.68 (t, J = 7.6 Hz, 1H), 3.40 – 3.20 (m, 2H)."""
    from mfs import nmrproc as NP
    sol = NP.solvent_key(solvent) or solvent or ""
    sol_fmt = {"CDCl3": "CDCl<sub>3</sub>"}.get(sol, sol)
    sol_fmt = sol
    if nucleus == "1H":
        head = f"¹H NMR ({mhz:.0f} MHz, {sol_fmt}) δ "
        dec = dec or 2
        parts = []
        for i0, m in sorted(enumerate(mults), key=lambda t: -t[1].centre):
            if (m.note and not include_notes and not m.label) or i0 in skip:
                continue
            if min_nH is not None and not m.label and (m.nH or 0) < min_nH:
                continue
            if m.mult == "m":
                d = f"{_fmt_d(m.hi - (1.0 / mhz), dec)} – {_fmt_d(m.lo + (1.0 / mhz), dec)}"
            else:
                d = _fmt_d(m.centre, dec)
            inner = [m.mult]
            if m.J:
                inner.append("J = " + ", ".join(f"{j:.1f}" for j in m.J) + " Hz")
            nh = m.nH_assigned if m.nH_assigned else m.nH
            if nh is not None:
                inner.append(f"{int(round(nh))}H" if nh >= 0.5 else f"{nh:.2f}H")
            if m.label and labels:
                inner.append(m.label)
            parts.append(f"{d} ({', '.join(inner)})")
        return head + ", ".join(parts) + "."
    dec = dec or 1
    head = f"¹³C{{¹H}} NMR ({mhz:.0f} MHz, {sol_fmt}) δ "
    vals = []
    for m in sorted(mults, key=lambda m: -m.centre):
        if m.note and not include_notes:
            continue
        s = _fmt_d(m.centre, dec)
        extra = [x for x in (getattr(m, "mult", "") if m.mult not in ("s", "m") else "", m.label) if x]
        vals.append(s + (f" ({', '.join(extra)})" if extra else ""))
    return head + ", ".join(vals) + "."


def carbon_lines(sp, k=8.0):
    """13C: one Multiplet per line (decoupled). Solvent tagged."""
    peaks = pick_peaks(sp.ppm, sp.y, k=k, sf=sp.sf, min_rel=0.01)
    from mfs import nmrproc as NP
    s = mad_sigma(sp.y)
    out = []
    for p in sorted(peaks, key=lambda p: -p.ppm):
        if out and abs(out[-1].centre - p.ppm) * sp.sf < 1.5:
            continue
        w = 1.5 / sp.sf * 3
        out.append(Multiplet(lo=p.ppm - w, hi=p.ppm + w, centre=p.ppm, integral=p.height, mult="s", lines=[p.ppm]))
    tag_solvent(out, sp)
    # CDCl3 triplet / DMSO septet: tag lines near the solvent value
    key = NP.solvent_key(sp.info.get("solvent"))
    if key and NP.SOLVENTS[key].get("13C"):
        v = NP.SOLVENTS[key]["13C"]
        for m in out:
            if abs(m.centre - v) < 1.0:
                m.note = f"{key} (solvent)"
    return out


def dept_types(c_mults, dept135=None, dept90=None, tol=0.15):
    """Carbon type from DEPT: Cq (absent in DEPT-135), CH (DEPT-90 +), CH2 (DEPT-135 −), CH3 (DEPT-135 +, absent in 90)."""
    from mfs import nmrproc as NP
    l135 = NP.line_list(dept135, k=8) if dept135 is not None else []
    l90 = NP.line_list(dept90, k=8) if dept90 is not None else []
    for m in c_mults:
        if m.note:
            continue
        a = [h for p, h in l135 if abs(p - m.centre) < tol]
        b = [h for p, h in l90 if abs(p - m.centre) < tol]
        if dept135 is None:
            continue
        if not a:
            m.mult = "Cq"
        elif a[0] < 0:
            m.mult = "CH2"
        elif dept90 is not None and b and b[0] > 0:
            m.mult = "CH"
        elif dept90 is not None:
            m.mult = "CH3"
        else:
            m.mult = "CH/CH3"
    return c_mults
