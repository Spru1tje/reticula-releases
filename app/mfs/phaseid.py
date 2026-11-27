"""v2.16: automatic phase identification (search-match) against the reference library (app/references/library: COD structures and standards).

  1. peaks of the measured pattern (background removed, positions refined, intensities relative);
  2. the stick pattern of every library phase for the radiation of the measurement (pymatgen kinematic intensities, cached);
  3. every phase is matched with a small search over an isotropic lattice scale (solid solutions / thermal expansion / a different metal) and a zero shift:
     recall = intensity-weighted share of the strong reflections of the phase that are present in the pattern,
     precision = share of the measured peak intensity that the phase explains;
  4. greedy multi-phase: after the best phase its peaks are removed and the rest is searched again.
The result is a ranked list that the user confirms - it only proposes phases, the refinement then decides (Rwp, weight fractions)."""
from __future__ import annotations

import json
import os

import numpy as np
from scipy.signal import find_peaks, savgol_filter

HERE = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.join(os.path.dirname(HERE), "references", "library")
_STICKS: dict = {}


def library() -> list[dict]:
    """The reference library: [{file, path, name, kind, source, formula, cell}]."""
    try:
        idx = json.load(open(os.path.join(LIB_DIR, "index.json"), encoding="utf-8"))
    except OSError:
        return []
    for e in idx:
        e["path"] = os.path.normpath(os.path.join(LIB_DIR, e["file"]))
    return [e for e in idx if os.path.exists(e["path"])]


def exp_peaks(x, y, min_prom: float = 0.02, max_peaks: int = 60):
    """Peaks of a measured pattern: positions (°2θ), heights above the background (relative to the strongest, 0-1), widths (FWHM estimate, °)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    o = np.argsort(x)
    x, y = x[o], y[o]
    step = float(np.median(np.diff(x)))
    win = max(5, int(round(0.25 / step)) | 1)
    ys = savgol_filter(y, min(win, len(y) - (1 - len(y) % 2)), 2) if len(y) > win else y
    # background: running minimum over ~2.5° then smoothed
    n = max(5, int(round(2.5 / step)))
    from scipy.ndimage import minimum_filter1d, uniform_filter1d
    bg = uniform_filter1d(minimum_filter1d(ys, n), n)
    h = np.clip(ys - bg, 0, None)
    if h.max() <= 0:
        return np.array([]), np.array([]), np.array([])
    pk, pr = find_peaks(h, prominence=min_prom * h.max(), distance=max(2, int(0.15 / step)))
    pos, hh, fw = [], [], []
    for i in pk:
        a, b = max(i - 2, 0), min(i + 3, len(x))
        if b - a >= 3:
            c = np.polyfit(x[a:b], h[a:b], 2)
            xp = -c[1] / (2 * c[0]) if c[0] < 0 else x[i]
            xp = xp if abs(xp - x[i]) < 2 * step else x[i]
        else:
            xp = x[i]
        half = h[i] / 2
        l = i
        while l > 0 and h[l] > half:
            l -= 1
        r = i
        while r < len(h) - 1 and h[r] > half:
            r += 1
        pos.append(float(xp))
        hh.append(float(h[i]))
        fw.append(float(x[r] - x[l]))
    pos, hh, fw = np.array(pos), np.array(hh), np.array(fw)
    order = np.argsort(hh)[::-1][:max_peaks]
    order = np.sort(order)
    return pos[order], hh[order] / h.max(), fw[order]


def sticks(path: str, wavelength: float, tt_range=(2.0, 120.0)):
    """Reflection positions (°2θ) and relative intensities (0-1) of a library CIF for this wavelength (cached)."""
    key = (path, round(wavelength, 4), tt_range)
    if key in _STICKS:
        return _STICKS[key]
    import warnings
    warnings.filterwarnings("ignore")
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    from pymatgen.core import Structure
    s = Structure.from_file(path)
    pat = XRDCalculator(wavelength=wavelength).get_pattern(s, two_theta_range=tt_range)
    tt = np.array(pat.x, float)
    I = np.array(pat.y, float)
    I = I / I.max() if len(I) else I
    _STICKS[key] = (tt, I)
    return tt, I


def _score(tt_sim, I_sim, pk, ph, lam, scale, zero, tol, tt_lo, tt_hi, min_I=0.03):
    """Recall / precision of one phase for a lattice scale and zero shift. Positions: d -> d·scale (cell scale), then the zero shift."""
    s = np.clip(np.sin(np.radians(tt_sim / 2)) / scale, 0, 0.999)
    tt = 2 * np.degrees(np.arcsin(s)) + zero
    ok = (tt >= tt_lo) & (tt <= tt_hi) & (I_sim >= min_I)
    if not ok.any() or len(pk) == 0:
        return 0.0, 0.0, 0, 0, 0.0
    t, I = tt[ok], I_sim[ok]
    j = np.searchsorted(pk, t)
    left = np.clip(j - 1, 0, len(pk) - 1)
    right = np.clip(j, 0, len(pk) - 1)
    dl, dr = np.abs(pk[left] - t), np.abs(pk[right] - t)
    idx = np.where(dl <= dr, left, right)
    dist = np.minimum(dl, dr)
    hit = dist <= tol
    recall = float(I[hit].sum() / I.sum())
    used = np.unique(idx[hit])
    precision = float(ph[used].sum() / max(ph.sum(), 1e-12)) if len(used) else 0.0
    # do the relative intensities of the matched reflections agree? (cosine similarity of simulated and measured heights)
    a, b = I[hit], ph[idx[hit]]
    cosine = float(a @ b / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-12)) if hit.any() else 0.0
    return recall, precision, int(hit.sum()), int(ok.sum()), cosine


def refine_cell_zero(tt_sim, I_sim, pk, ph, scale0, zero0, tol, min_I=0.03):
    """Least-squares refinement of the lattice scale and the zero shift on the matched peaks (intensity-weighted): the grid of the search is too coarse
    and scale / zero are partly confounded, so the final values come from this fit."""
    from scipy.optimize import least_squares
    sel = I_sim >= min_I
    ts, Is = tt_sim[sel], I_sim[sel]

    def pos(par):
        return 2 * np.degrees(np.arcsin(np.clip(np.sin(np.radians(ts / 2)) / par[0], 0, 0.999))) + par[1]

    def match_idx(par):
        t = pos(par)
        j = np.argmin(np.abs(pk[None, :] - t[:, None]), axis=1)
        return j, np.abs(pk[j] - t) <= tol * 1.5

    par = np.array([scale0, zero0])
    for _ in range(3):
        j, ok = match_idx(par)
        if ok.sum() < 3:
            break
        w = np.sqrt(Is[ok] * ph[j[ok]])
        f = lambda q: w * (pos(q)[ok] - pk[j[ok]])
        sol = least_squares(f, par, bounds=([0.95, -0.4], [1.05, 0.4]), x_scale=[0.002, 0.05])
        par = sol.x
    return float(par[0]), float(par[1])


def match(x, y, wavelength=1.540593, scale_range=0.025, tol=0.12, n_phases=3, lib=None, min_prom=0.02, progress=None):
    """Ranks the library phases for a pattern. Returns a list of dicts (greedy, best first):
    name, file, kind, source, formula, score (0-1), recall, precision, n_matched, n_expected, cell_scale (1 = as in the library), cell_change_pct, zero (°),
    explained (share of the pattern intensity explained after the previous phases were removed)."""
    lib = lib if lib is not None else library()
    # v2.24.0: entries flagged exclude_matching (e.g. MIL-53 until its CIF is corrected) are not offered by the matching
    lib = [e for e in lib if not e.get("exclude_matching")]
    pk, ph, fw = exp_peaks(x, y, min_prom)
    if len(pk) < 3:
        return []
    tt_lo, tt_hi = float(np.min(x)), float(np.max(x))
    results = []
    remaining = np.ones(len(pk), bool)
    pk0, ph0 = pk.copy(), ph.copy()
    taken: set[str] = set()
    scales = np.linspace(1 - scale_range, 1 + scale_range, 51)
    zeros = np.linspace(-0.12, 0.12, 7)
    for rank in range(n_phases):
        best = None
        for k, e in enumerate(lib):
            if e["name"] in taken:
                continue
            try:
                tt_sim, I_sim = sticks(e["path"], wavelength, (max(tt_lo - 3, 1.0), min(tt_hi + 3, 150.0)))
            except Exception:
                continue
            if len(tt_sim) == 0:
                continue
            pkr, phr = pk[remaining], ph[remaining]
            top = (0.0, None)
            for sc in scales:
                for z in zeros:
                    rc, pr, nm, ne, cs = _score(tt_sim, I_sim, pkr, phr, wavelength, sc, z, tol, tt_lo, tt_hi)
                    # presence (recall x precision), intensity agreement, and a mild preference for the cell of the library (a solid solution is allowed, not favoured)
                    sco = float(np.sqrt(rc * pr)) * (0.7 + 0.3 * cs) * (1 - 0.12 * abs(sc - 1) / scale_range) if nm >= 3 else 0.0
                    if sco > top[0]:
                        top = (sco, (rc, pr, nm, ne, sc, z))
            if top[1] is None:
                continue
            if best is None or top[0] > best["score"]:
                rc, pr, nm, ne, sc, z = top[1]
                best = {"name": e["name"], "file": e["path"], "kind": e.get("kind", ""), "source": e.get("source", ""), "formula": e.get("formula", ""),
                        "score": top[0], "recall": rc, "precision": pr, "n_matched": nm, "n_expected": ne, "cell_scale": float(sc), "cell_change_pct": 100 * (float(sc) - 1),
                        "zero": float(z)}
            if progress:
                progress((rank + (k + 1) / len(lib)) / n_phases, f"Searching the library: {e['name']}")
        if best is None or best["n_matched"] < 3 or best["score"] < (0.62 if rank == 0 else 0.5) or (rank == 0 and best["precision"] < 0.45 and best["recall"] < 0.9):
            break
        best["confidence"] = "very likely" if best["score"] >= 0.85 else "likely" if best["score"] >= 0.7 else "possible"
        tt_sim, I_sim = sticks(best["file"], wavelength, (max(tt_lo - 3, 1.0), min(tt_hi + 3, 150.0)))
        sc_f, z_f = refine_cell_zero(tt_sim, I_sim, pk[remaining], ph[remaining], best["cell_scale"], best["zero"], tol)
        best["cell_scale"], best["zero"], best["cell_change_pct"] = sc_f, z_f, 100 * (sc_f - 1)
        # peaks explained by the best phase are removed for the next round
        s = np.clip(np.sin(np.radians(tt_sim / 2)) / best["cell_scale"], 0, 0.999)
        tt = 2 * np.degrees(np.arcsin(s)) + best["zero"]
        used = np.zeros(len(pk), bool)
        for t, I in zip(tt, I_sim):
            if I < 0.01:
                continue
            d = np.abs(pk - t)
            j = int(np.argmin(d))
            if d[j] <= tol * 1.2:
                used[j] = True
        before = ph0[remaining].sum()
        best["explained"] = float(ph0[used & remaining].sum() / max(ph0.sum(), 1e-12))
        remaining &= ~used
        taken.add(best["name"])
        results.append(best)
        if remaining.sum() < 3 or ph0[remaining].sum() < 0.08 * ph0.sum():
            break
    return results
