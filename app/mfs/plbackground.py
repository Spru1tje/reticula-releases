r"""v2.18.4: background (pedestal, dark counts, stray light, broad solvent / host emission) subtraction for luminescence spectra, visible and near-infrared.

The problem. A luminescence spectrum is the emission PLUS a background that the experiment adds: the detector pedestal (dark counts; in the NIR several thousand counts under
a lanthanide spectrum is normal), stray light and laser leakage, second-order light of the grating, broad emission of the solvent / cuvette / host. Band ratios, thermometry
(sigma_T from photon counts), integrated intensities and normalisations all change with it, so it must be removed in a standard, documented, reproducible way.

What this module does (numpy / scipy only).
  * finds the SIGNAL-FREE WINDOWS of a whole series automatically: regions in which no spectrum rises above its own robust lower envelope (SNIP) by more than k noise levels;
  * fits the background in those windows with a polynomial of order 0 (pure pedestal) … 4, by iteratively re-weighted (Tukey bisquare) least squares, and chooses the ORDER with
    the Bayesian information criterion summed over the series (a higher order must earn its extra parameter);
  * alternatives for spectra that have no signal-free window: arPLS (asymmetrically reweighted penalised least squares, Baek 2015), ModPoly (iterative polynomial, Zhao & Lui 2007) and
    SNIP (statistics-sensitive non-linear iterative peak clipping), and the subtraction of a measured dark / background spectrum (scaled in the signal-free windows);
  * ranks the candidate methods by objective numbers (bias and noise in the signal-free windows, over-subtraction below zero, complexity) and returns the plan with the reasons;
  * removes single-point spikes (Hampel filter) when asked;
  * the same plan is applied to every spectrum of the series, so the ratios between temperatures are not distorted by a different correction per spectrum.

Manual control: the user can give the windows, the order or the method, or a dark file; the automatic plan is the default."""
from __future__ import annotations

import math

import numpy as np

METHODS = ("none", "offset", "poly", "rubber", "modpoly", "arpls", "snip", "dark")
LABELS = {"none": "none", "masked": "polynomial under a data-derived emission mask", "offset": "constant pedestal (signal-free windows)", "poly": "polynomial in signal-free windows", "rubber": "rubber band (lower convex hull)",
          "modpoly": "ModPoly (iterative polynomial)", "arpls": "arPLS (penalised least squares)", "snip": "SNIP (peak clipping)", "dark": "measured dark / background spectrum"}


# ------------------------------------------------------------------------------------------------------------------------------ helpers
def noise_sigma(y) -> float:
    """Noise level (one standard deviation) from the point-to-point differences: robust against the (smooth) emission bands."""
    y = np.asarray(y, float)
    d = np.diff(y)
    if len(d) < 3:
        return 1.0
    s = 1.4826 * float(np.median(np.abs(d - np.median(d)))) / math.sqrt(2.0)
    return s if s > 0 else float(np.std(d) / math.sqrt(2.0)) or 1.0


def range_kind(x) -> str:
    """'UV/VIS' (< 750 nm), 'VIS-NIR' (up to 1000 nm), 'NIR' (above 850 nm window start) from the wavelength range."""
    x = np.asarray(x, float)
    lo, hi = float(x.min()), float(x.max())
    if lo >= 850:
        return "NIR"
    if hi > 1000:
        return "VIS-NIR"
    return "UV/VIS"


def _cheb(x, order, lo, hi):
    t = 2.0 * (np.asarray(x, float) - lo) / max(hi - lo, 1e-12) - 1.0
    return np.polynomial.chebyshev.chebvander(t, order)


def snip(y, n_iter=40, transform=True):
    """SNIP baseline (Morháč 1997): n_iter clipping passes with a shrinking window; the (log-log-sqrt) transform keeps the clipping fair for large dynamic ranges."""
    y = np.asarray(y, float)
    n = len(y)
    shift = float(min(y.min(), 0.0))
    v = y - shift
    if transform:
        v = np.log(np.log(np.sqrt(np.maximum(v, 0.0) + 1.0) + 1.0) + 1.0)
    w = v.copy()
    for p in range(int(n_iter), 0, -1):
        a = np.empty_like(w)
        a[:] = w
        if n > 2 * p:
            avg = 0.5 * (w[: n - 2 * p] + w[2 * p:])
            a[p: n - p] = np.minimum(w[p: n - p], avg)
        w = a
    if transform:
        w = (np.exp(np.exp(w) - 1.0) - 1.0) ** 2 - 1.0
    return w + shift


def arpls(y, lam=1e6, ratio=1e-6, niter=60):
    """arPLS (Baek et al., Analyst 2015): a smooth baseline under the data; points above it (peaks) get little weight, points below it are kept."""
    from scipy import sparse
    from scipy.sparse.linalg import spsolve
    y = np.asarray(y, float)
    n = len(y)
    D = sparse.diags([np.ones(n - 2), -2 * np.ones(n - 2), np.ones(n - 2)], [0, 1, 2], shape=(n - 2, n))
    H = lam * (D.T @ D)
    w = np.ones(n)
    z = y.copy()
    for _ in range(niter):
        W = sparse.diags(w, 0, shape=(n, n))
        z = spsolve((W + H).tocsc(), w * y)
        d = y - z
        dn = d[d < 0]
        if len(dn) < 2:
            break
        m, s = dn.mean(), dn.std()
        wt = 1.0 / (1.0 + np.exp(np.clip(2.0 * (d - (2 * s - m)) / max(s, 1e-12), -50, 50)))
        if np.linalg.norm(w - wt) / max(np.linalg.norm(w), 1e-12) < ratio:
            w = wt
            break
        w = wt
    return z


def rubberband(x, y, smooth=7):
    """Rubber band: the lower convex hull of the (lightly smoothed) spectrum, i.e. the straight-segment floor that a stretched rubber band under the data would make.
    It never rises above the data, so it removes the floor (pedestal, slope) but cannot swallow a broad emission band that spans the range - the safe choice without signal-free windows."""
    from scipy.signal import savgol_filter
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ys = savgol_filter(y, smooth, 2) if len(y) > smooth + 2 else y
    hull = []
    for i in range(len(x)):
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (x2 - x1) * (ys[i] - y1) - (y2 - y1) * (x[i] - x1) <= 0:       # the last point lies on or above the chord: not on the lower hull
                hull.pop()
            else:
                break
        hull.append((x[i], ys[i]))
    hx, hy = zip(*hull)
    return np.interp(x, hx, hy)


def modpoly(x, y, order=3, niter=100, tol=1e-3):
    """ModPoly (Lieber & Mahadevan-Jansen; improved by Zhao & Lui 2007): fit a polynomial, clip the data to the fit, repeat - the polynomial sinks below the peaks."""
    x = np.asarray(x, float)
    yy = np.asarray(y, float).copy()
    B = _cheb(x, order, x.min(), x.max())
    prev = None
    for _ in range(niter):
        c, *_ = np.linalg.lstsq(B, yy, rcond=None)
        z = B @ c
        yy = np.minimum(yy, z)
        if prev is not None and np.linalg.norm(z - prev) / max(np.linalg.norm(prev), 1e-12) < tol:
            break
        prev = z
    return B @ c


def poly_in_windows(x, y, mask, order, robust=True, niter=15):
    """Polynomial (Chebyshev basis) fitted to the points of `mask`, with Tukey-bisquare re-weighting against stray peaks. Returns (baseline over all x, coefficients, rss, n)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    B = _cheb(x, order, x.min(), x.max())
    Bm, ym = B[mask], y[mask]
    w = np.ones(len(ym))
    c = np.zeros(order + 1)
    for _ in range(niter if robust else 1):
        sw = np.sqrt(w)
        c, *_ = np.linalg.lstsq(Bm * sw[:, None], ym * sw, rcond=None)
        r = ym - Bm @ c
        if not robust:
            break
        s = max(1.4826 * float(np.median(np.abs(r))), 1e-9)
        u = r / (4.685 * s)
        w = np.where(np.abs(u) < 1, (1 - u * u) ** 2, 0.0) + 1e-6
    r = ym - Bm @ c
    return B @ c, c, float(np.sum(w * r * r) / max(np.sum(w), 1e-12) * len(ym)), int(len(ym))


def despike(y, k=8.0, half=3):
    """Hampel filter: points that differ from the median of their neighbours by more than k noise levels are replaced by that median. Returns (y, number of spikes)."""
    y = np.asarray(y, float)
    out = y.copy()
    n = len(y)
    if n < 2 * half + 3:
        return out, 0
    sig = noise_sigma(y)
    cnt = 0
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        nb = np.delete(y[lo:hi], i - lo)
        med = float(np.median(nb))
        if abs(y[i] - med) > k * max(sig, 1e-12) and abs(y[i] - med) > 3 * float(np.std(nb) + 1e-12):
            out[i] = med
            cnt += 1
    return out, cnt


# ------------------------------------------------------------------------------------------------------------------------------ windows
def _runs(mask):
    """[(start, stop)) of the True runs."""
    m = np.concatenate([[False], np.asarray(mask, bool), [False]])
    d = np.diff(m.astype(int))
    return list(zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]))


def find_free_windows(x, ys, k=3.5, min_pts=8, min_frac=0.025, edge=2):
    """Signal-free windows of a whole series. For each spectrum a robust lower envelope (SNIP, window about a quarter of the range) is subtracted from a lightly smoothed copy;
    a point is FREE when no spectrum of the series rises above its envelope by more than k noise levels. Short runs are dropped, every run is trimmed by `edge` points.
    Returns [(x_start, x_stop)] in the units of x (possibly empty)."""
    from scipy.signal import savgol_filter
    x = np.asarray(x, float)
    n = len(x)
    free = np.ones(n, bool)
    nit = max(8, int(0.22 * n))
    for y in ys:
        y = np.interp(x, np.asarray(y[0], float), np.asarray(y[1], float)) if isinstance(y, tuple) else np.asarray(y, float)
        sig = noise_sigma(y)
        ysm = savgol_filter(y, 7, 2) if n > 9 else y
        env = snip(ysm, nit, transform=False)
        free &= (ysm - env) < k * sig * 0.6 + 0.5 * sig
    keep = np.zeros(n, bool)
    for a, b in _runs(free):
        a2, b2 = a + edge, b - edge
        if b2 - a2 >= min_pts:
            keep[a2:b2] = True
    wins = [(float(x[a]), float(x[b - 1])) for a, b in _runs(keep)]
    if sum(b - a for a, b in wins) < min_frac * (x.max() - x.min()):
        return []
    return wins


def mask_from_windows(x, windows):
    x = np.asarray(x, float)
    m = np.zeros(len(x), bool)
    for lo, hi in windows:
        m |= (x >= min(lo, hi)) & (x <= max(lo, hi))
    return m


# ------------------------------------------------------------------------------------------------------------------------------ plan
def _score(x, y, base, free_mask, sig):
    """Objective numbers of a baseline: bias and rms in the signal-free windows (in noise units), and the over-subtraction (points below -3 sigma)."""
    r = y - base
    if free_mask.sum() >= 3:
        bias = float(np.mean(r[free_mask])) / sig
        rms = float(np.sqrt(np.mean(r[free_mask] ** 2))) / sig
    else:
        bias, rms = float("nan"), float("nan")
    neg = float(np.mean(r < -3.0 * sig))
    return {"bias": bias, "rms": rms, "below": neg}


def autoselect(x, ys, dark=None, windows=None, max_order=2, use_dark=None):
    """The automatic plan for a series. x: common wavelength grid; ys: list of intensity arrays (same grid). Returns a dict
    {method, order, windows, sigma, kind, candidates [...], notes [...], dark_scale}. `windows` overrides the automatic windows (manual mode)."""
    x = np.asarray(x, float)
    ys = [np.asarray(y, float) for y in ys]
    sig = float(np.median([noise_sigma(y) for y in ys]))
    kind = range_kind(x)
    notes = []
    wins = list(windows) if windows else find_free_windows(x, ys)
    free = mask_from_windows(x, wins) if wins else np.zeros(len(x), bool)
    cands = []
    if len(wins) and free.sum() >= 8:
        # polynomial order by BIC (summed over the series)
        bic = {}
        for o in range(0, max_order + 1):
            tot = 0.0
            for y in ys:
                _, _, rss, nn = poly_in_windows(x, y, free, o)
                tot += nn * math.log(max(rss / max(nn, 1), 1e-12)) + (o + 1) * math.log(max(nn, 2))
            bic[o] = tot
        best_o = min(bic, key=bic.get)
        order = best_o
        for o in sorted(bic):                                  # the lowest order that is within ΔBIC = 10 of the best wins: a higher order must clearly earn its parameter
            if bic[o] <= bic[best_o] + 10.0:
                order = o
                break
        for o in sorted(bic):
            base = [poly_in_windows(x, y, free, o)[0] for y in ys]
            sc = [_score(x, y, b, free, sig) for y, b in zip(ys, base)]
            cands.append({"method": "offset" if o == 0 else "poly", "order": o, "bic": float(bic[o]),
                          "bias": float(np.nanmean([abs(s["bias"]) for s in sc])), "rms": float(np.nanmean([s["rms"] for s in sc])), "below": float(np.mean([s["below"] for s in sc]))})
    else:
        notes.append("No signal-free window was found (the emission covers the whole range): the rubber band (lower convex hull) is used. It removes the floor but not a broad band. "
                     "Check it, give the windows by hand, or load a measured dark / blank spectrum (General practice: a blank or solvent measurement is the reliable background for broad emission).")
    # flexible alternatives (also scored; used when no windows)
    for meth in ("rubber", "arpls", "snip", "modpoly"):
        try:
            base = [(rubberband(x, y) if meth == "rubber" else arpls(y, lam=1e6) if meth == "arpls" else snip(y, max(8, int(0.22 * len(y)))) if meth == "snip" else modpoly(x, y, 3)) for y in ys]
        except Exception:
            continue
        sc = [_score(x, y, b, free, sig) for y, b in zip(ys, base)]
        cands.append({"method": meth, "order": None, "bic": float("nan"), "bias": float(np.nanmean([abs(s["bias"]) for s in sc])) if free.any() else float("nan"),
                      "rms": float(np.nanmean([s["rms"] for s in sc])) if free.any() else float("nan"), "below": float(np.mean([s["below"] for s in sc]))})
    plan = {"method": "none", "order": None, "windows": wins, "sigma": sig, "kind": kind, "candidates": cands, "notes": notes, "dark_scale": None}
    if dark is not None and use_dark is not False:
        plan["method"] = "dark"
        plan["notes"].append("A measured dark / background spectrum is subtracted (scaled in the signal-free windows).")
        return plan
    if len(wins) and free.sum() >= 8:
        plan["method"] = "offset" if order == 0 else "poly"
        plan["order"] = order
        plan["notes"].append(f"Signal-free windows: " + ", ".join(f"{a:.0f}–{b:.0f} nm" for a, b in wins) + f". Order {order} chosen by the information criterion (BIC).")
    else:
        # v2.18.6: no signal-free window: a polynomial under a data-derived emission mask (not SNIP/arPLS, which eat broad bands)
        # v2.18.6 safety: the masked baseline must never rise above the spectrum by more than three noise units; if it does, the
        # emission covers too much of the range for a polynomial to be anchored, and the rubber band (which never rises above the data) is used
        worst = min(float(np.min(y - masked_poly(x, y)[0])) for y in ys)
        if worst < -3.0 * sig:
            plan["method"] = "rubber"
            plan["notes"].append("The emission covers so much of the range that no polynomial can be anchored outside it (the masked baseline would rise "
                                 "above the spectrum): the rubber band is used. Check it, and give the signal-free windows or a blank spectrum by hand.")
        else:
            plan["method"] = "masked"
        plan["max_order"] = 1
        plan["notes"].append("No signal-free window: the baseline is a polynomial fitted only outside the emission (mask derived from the data, "
                             "order by BIC). The rubber band is the alternative. General practice: check the baseline under the band; a blank or solvent "
                             "spectrum is the reliable background for broad emission.")
    if kind in ("NIR", "VIS-NIR"):
        plan["notes"].append("NIR spectrum: a detector pedestal of thousands of counts is normal; after subtraction check the signal-free windows (they must be flat around zero). "
                             "For samples in water, absorption of water near 975 and 1195 nm changes band ratios (General practice: compare with a measurement in D₂O or a dry sample).")
    return plan


def masked_poly(x, y, max_order=1, width=9, k=3.0, n_iter=8):
    """v2.18.6: baseline under a broad emission band, from the data only. Fit a polynomial to the points outside the emission,
    mark where the smoothed residual rises above k times its noise (plus half a window on each side), refit, repeat. The order is
    the lowest one within ΔBIC = 10 of the best, capped at max_order (default 1: order 2 extrapolates into the band, measured −30 % area error). Returns (baseline, order, mask of the points used for the fit).
    Why: SNIP and arPLS eat the tails of a broad band that spans the range (test: +13 % and −26 % area error); a polynomial that is
    fitted only where there is no emission does not. Its limit: it extrapolates into the band, so a high order is dangerous."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    free = np.ones(len(y), bool)
    base = None
    order = 0
    for _ in range(n_iter):
        cand = {}
        for o in range(0, max_order + 1):
            if free.sum() < o + 3:
                continue
            b_o, _, rss, nn = poly_in_windows(x, y, free, o)
            cand[o] = (b_o, nn * math.log(max(rss / max(nn, 1), 1e-12)) + (o + 1) * math.log(max(nn, 2)))
        if not cand:
            break
        best = min(cand, key=lambda o: cand[o][1])
        order = min(o for o in cand if cand[o][1] <= cand[best][1] + 10.0)
        base = cand[order][0]
        r = y - base
        sm = np.convolve(r, np.ones(width) / width, mode="same")
        sig = np.std(r[free]) / math.sqrt(width) if free.sum() > 3 else 0.0
        emit = sm > k * max(sig, 1e-12)
        new_free = ~np.convolve(emit.astype(float), np.ones(width), mode="same").astype(bool)
        if np.array_equal(new_free, free):
            break
        free = new_free
    if base is None:
        base = np.full_like(y, float(np.percentile(y, 5)))
        order = 0
    return base, order, free


def alternatives(x, plan, max_order=2):
    """Other reasonable background plans for the same series: [(label, plan)]. Used to measure how much the ratio depends on the background CHOICE
    (the spread of the ratio over these plans is added to its systematic error). Polynomial plans: the other orders (0..max_order) and the same order
    with every signal-free window shortened by a quarter on both sides; flexible plans: the other flexible methods; dark-file plans: a polynomial in the windows."""
    x = np.asarray(x, float)
    m = plan.get("method")
    out = []
    if m in ("offset", "poly") and plan.get("windows"):
        cur = int(plan.get("order") or 0)
        for o in range(0, max_order + 1):
            if o != cur:
                out.append((f"order {o}", {**plan, "method": "offset" if o == 0 else "poly", "order": o}))
        short = [(a + 0.25 * (b - a), b - 0.25 * (b - a)) for a, b in plan["windows"] if (b - a) > 6 * float(np.median(np.diff(x)))]
        if short and sum(1 for _ in short) >= 1:
            out.append(("shorter windows", {**plan, "windows": short}))
    elif m == "masked":
        out.append(("rubber band", {**plan, "method": "rubber", "order": None}))
    elif m in ("rubber", "arpls", "snip", "modpoly"):
        for k in ("rubber", "arpls", "snip", "modpoly"):
            if k != m:
                out.append((LABELS[k].split(" (")[0], {**plan, "method": k, "order": 3 if k == "modpoly" else None}))
    elif m == "dark" and plan.get("windows"):
        for o in (0, 1):
            out.append((f"polynomial order {o} in the windows", {**plan, "method": "offset" if o == 0 else "poly", "order": o}))
    return out


def baseline_for(x, y, plan, dark=None):
    """The baseline of ONE spectrum for a plan."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = plan["method"]
    if m == "none":
        return np.zeros_like(y)
    free = mask_from_windows(x, plan["windows"]) if plan.get("windows") else np.zeros(len(x), bool)
    if m in ("offset", "poly"):
        if free.sum() < 3:
            return np.full_like(y, float(np.percentile(y, 5)))
        return poly_in_windows(x, y, free, int(plan.get("order") or 0))[0]
    if m == "masked":
        return masked_poly(x, y, max_order=int(plan.get("max_order") or 1))[0]
    if m == "rubber":
        return rubberband(x, y)
    if m == "arpls":
        return arpls(y, lam=float(plan.get("lam", 1e6)))
    if m == "snip":
        return snip(y, int(plan.get("iterations", max(8, int(0.22 * len(y))))))
    if m == "modpoly":
        return modpoly(x, y, int(plan.get("order") or 3))
    if m == "dark":
        if dark is None:
            return np.zeros_like(y)
        d = np.interp(x, np.asarray(dark[0], float), np.asarray(dark[1], float))
        if free.sum() >= 4:
            A = np.stack([d[free], np.ones(int(free.sum()))], 1)
            (k, c), *_ = np.linalg.lstsq(A, y[free], rcond=None)
            k = float(np.clip(k, 0.0, 3.0)) if np.isfinite(k) else 1.0
            return k * d + (float(np.mean(y[free] - k * d[free])) if True else 0.0)
        return d
    raise ValueError(m)


def apply_plan(spectra, plan, dark=None, despike_on=False):
    """spectra: [{'wl', 'I'}]; returns a list of {'y', 'baseline', 'spikes'} in the same order (the spectra are NOT changed)."""
    out = []
    for s in spectra:
        x = np.asarray(s["wl"], float)
        y = np.asarray(s["I"], float)
        sp = 0
        if despike_on:
            y, sp = despike(y)
        b = baseline_for(x, y, plan, dark)
        out.append({"y": y - b, "baseline": b, "spikes": sp})
    return out


def check_after(x, y_corr, plan):
    """Quality numbers after the subtraction: bias / noise in the signal-free windows (should be 0 ± 1 noise level) and the over-subtraction below zero."""
    x = np.asarray(x, float)
    free = mask_from_windows(x, plan.get("windows") or [])
    sig = noise_sigma(y_corr)
    return _score(x, np.asarray(y_corr, float), np.zeros(len(x)), free, sig)
