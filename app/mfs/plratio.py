"""v2.15: robust band ratios for luminescence thermometry - a GLOBAL fit of the whole temperature series instead of integrating windows.

Why: with overlapping, wavy, dirty or doubled bands an integration window cannot separate the bands (the overlap is cut at an arbitrary
valley, and the result changes with the window). Here every spectrum of the series is described by the SAME bands:
  * every band is a pseudo-Voigt with unit area; its position and width may drift LINEARLY with temperature (shared by all spectra);
  * the areas and the (polynomial) background of every spectrum are linear parameters, solved exactly in every step (variable projection);
  * extra 'nuisance' bands (a shoulder, a third peak) can be added so that they do not contaminate the two bands of the ratio;
  * the simplest model (shared shape) and the drifting model are compared with the BIC and the better one is used.
The ratio is area1 / area2 with a standard error that includes the covariance of the two areas; a robustness check refits with other windows
and the other model and adds the spread as a systematic error. Quality flags per spectrum (R², reduced χ², residual autocorrelation)
tell the user where not to trust a point. Nothing here depends on the window the user would have drawn."""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import least_squares, nnls
from scipy.signal import find_peaks

LN2 = math.log(2.0)


def _noise(y):
    """Noise level of a spectrum from the point-to-point differences (robust: a smooth band hardly contributes)."""
    d = np.diff(np.asarray(y, float))
    if len(d) < 3:
        return 1.0
    return float(np.median(np.abs(d - np.median(d))) / 0.6745 / math.sqrt(2.0)) or float(np.std(d) / math.sqrt(2.0)) or 1.0


def pv_unit(x, c, w, eta):
    """Unit-area pseudo-Voigt (per nm) with FWHM w and Lorentzian share eta."""
    z = (x - c) / w
    g = (2.0 * math.sqrt(LN2 / math.pi) / w) * np.exp(-4.0 * LN2 * z * z)
    l = (2.0 / (math.pi * w)) / (1.0 + 4.0 * z * z)
    return eta * l + (1.0 - eta) * g


# ------------------------------------------------------------------------------------------------------------- start values
def auto_bands(x, y, lo=None, hi=None, n_extra=0):
    """The two strongest separate bands (+ optionally further ones) of a (mean) spectrum: [(centre, fwhm_guess), ...] sorted by centre.
    The first two are the ratio pair (I1 = shorter wavelength). Returns None when fewer than two bands are found."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    lo = x.min() if lo is None else lo
    hi = x.max() if hi is None else hi
    m = (x >= lo) & (x <= hi)
    xs, ys = x[m], y[m]
    if len(xs) < 15:
        return None
    from scipy.signal import savgol_filter
    step = float(np.median(np.diff(xs)))
    win = max(5, int(round(2.0 / step)) | 1)
    ysm = savgol_filter(ys, min(win, len(ys) - (1 - len(ys) % 2)), 2) if len(ys) > win else ys
    rng = float(ysm.max() - ysm.min())
    if rng <= 0:
        return None
    # candidates: maxima of the spectrum AND of the negative second derivative (that reveals shoulders of overlapping bands)
    w2 = max(7, int(round(0.06 * len(ys))) | 1)
    d2 = savgol_filter(ys, min(w2, len(ys) - (1 - len(ys) % 2)), 3, deriv=2)
    cand = list(find_peaks(ysm, prominence=0.03 * rng)[0])
    pk2, pr2 = find_peaks(-d2, prominence=0.04 * float(np.max(-d2) - np.min(-d2)))
    cand += [int(i) for i in pk2]
    cand = sorted(set(int(i) for i in cand))
    merged = []
    for i in cand:                                                # merge candidates closer than ~1.5 nm
        if merged and xs[i] - xs[merged[-1]] < max(1.5, 3 * step):
            if ysm[i] > ysm[merged[-1]]:
                merged[-1] = i
        else:
            merged.append(i)
    base = float(np.percentile(ysm, 5))
    merged = [i for i in merged if ysm[i] - base > 0.06 * rng]
    if len(merged) < 2:
        return None
    strongest = sorted(merged, key=lambda i: ysm[i], reverse=True)[:2 + max(0, n_extra)]
    chosen = sorted(strongest)
    out = []
    for p_ in chosen:
        h = ysm[p_]
        half = h - (h - base) / 2.0
        i0 = p_
        while i0 > 0 and ysm[i0] > half:
            i0 -= 1
        i1 = p_
        while i1 < len(ysm) - 1 and ysm[i1] > half:
            i1 += 1
        w = max(float(xs[i1] - xs[i0]), 3 * step)
        out.append((float(xs[p_]), w))
    # a shoulder has no half-height crossing of its own: never guess a width larger than the distance to the nearest other band x 1.2
    for k, (c, w) in enumerate(out):
        gap = min([abs(c - o) for j, (o, _) in enumerate(out) if j != k] or [w])
        out[k] = (c, float(min(w, 1.2 * gap)))
    return out


def default_window(bands, x):
    """Fit window: the bands plus 3 FWHM on each side (at least as wide as the data allow)."""
    lo = min(c - 3.0 * w for c, w in bands)
    hi = max(c + 3.0 * w for c, w in bands)
    return float(max(lo, np.min(x))), float(min(hi, np.max(x)))


# ------------------------------------------------------------------------------------------------------------------ the fit
class _Model:
    def __init__(self, spectra, bands, window, bg_order, drift_pos, drift_width, counts):
        self.nb = len(bands)
        self.bg_order = int(bg_order)
        self.drift_pos, self.drift_width = bool(drift_pos), bool(drift_width)
        self.T = np.array([s["T"] for s in spectra], float)
        self.T0 = float(np.mean(self.T))
        self.Ts = max(float(np.ptp(self.T)), 1.0)                      # temperatures enter as (T - T0) / range
        lo, hi = window
        self.data = []
        for s in spectra:
            x = np.asarray(s["x"], float)
            y = np.asarray(s["y"], float)
            m = (x >= lo) & (x <= hi) & np.isfinite(y)
            self.data.append((x[m], y[m]))
        self.counts = counts
        self.w = [1.0 / np.maximum(y, 1.0) if counts else np.ones_like(y) / max(_noise(y), 1e-12) ** 2 for x, y in self.data]
        self.guess = list(bands)
        # parameter layout per band: c0, [s], w0, [g], eta
        self.names, p0, lb, ub = [], [], [], []
        allc = sorted(c for c, _ in bands)
        for k, (c, w) in enumerate(bands):
            gap = min([abs(c - o) for o in allc if o != c] or [w * 4])
            dc = max(min(gap * 0.45, 3.0 * w), 0.5 * w)
            self.names += [f"c{k}"]; p0 += [c]; lb += [c - dc]; ub += [c + dc]
            if drift_pos:
                self.names += [f"s{k}"]; p0 += [0.0]; lb += [-0.5 * w]; ub += [0.5 * w]          # nm over the whole temperature range
            self.names += [f"w{k}"]; p0 += [w]; lb += [0.25 * w]; ub += [4.0 * w]
            if drift_width:
                self.names += [f"g{k}"]; p0 += [0.0]; lb += [-0.6 * w]; ub += [0.6 * w]
            self.names += [f"e{k}"]; p0 += [0.5]; lb += [0.0]; ub += [1.0]
        self.p0, self.lb, self.ub = np.array(p0), np.array(lb), np.array(ub)

    def band_params(self, th, ti):
        """(centres, widths, etas) of all bands at the temperature index ti."""
        t = (self.T[ti] - self.T0) / self.Ts
        c, w, e = [], [], []
        i = 0
        for k in range(self.nb):
            c0 = th[i]; i += 1
            s = 0.0
            if self.drift_pos:
                s = th[i]; i += 1
            w0 = th[i]; i += 1
            g = 0.0
            if self.drift_width:
                g = th[i]; i += 1
            e.append(th[i]); i += 1
            c.append(c0 + s * t)
            w.append(max(w0 + g * t, 1e-3))
        return c, w, e

    def design(self, th, ti):
        x = self.data[ti][0]
        c, w, e = self.band_params(th, ti)
        cols = [pv_unit(x, c[k], w[k], e[k]) for k in range(self.nb)]
        xn = (x - x.mean()) / max(np.ptp(x), 1e-9) * 2
        cols_bg = [np.polynomial.chebyshev.chebval(xn, [0] * j + [1]) for j in range(self.bg_order + 1)]
        return np.stack(cols, 1), np.stack(cols_bg, 1)

    def solve(self, th, ti):
        A, B = self.design(th, ti)
        y, w = self.data[ti][1], self.w[ti]
        sw = np.sqrt(w)
        scale = max(float(np.abs(y).max()), 1e-12)
        M = np.hstack([A, B, -B]) * sw[:, None]
        c, _ = nnls(M, y / scale * sw, maxiter=5000)
        amp = c[:self.nb] * scale
        bg = (c[self.nb:self.nb + B.shape[1]] - c[self.nb + B.shape[1]:]) * scale
        return amp, bg, A, B

    def residual(self, th):
        out = []
        for ti in range(len(self.data)):
            amp, bg, A, B = self.solve(th, ti)
            out.append(np.sqrt(self.w[ti]) * (self.data[ti][1] - A @ amp - B @ bg))
        return np.concatenate(out)


def fit_global(spectra, bands, window=None, bg_order=1, drift="auto", counts=False, max_nfev=80):
    """Global fit of the whole series. spectra = [{'T': K, 'x': nm, 'y': counts}], bands = [(centre, fwhm_guess), ...] with the ratio pair FIRST
    (I1 = bands[0], I2 = bands[1]). drift = 'none' | 'drift' | 'auto' (BIC decides). Returns the result dict (see keys below)."""
    window = window or default_window(bands, np.concatenate([np.asarray(s["x"], float) for s in spectra]))
    if drift == "auto":
        a = fit_global(spectra, bands, window, bg_order, "none", counts, max_nfev)
        b = fit_global(spectra, bands, window, bg_order, "drift", counts, max_nfev)
        best = a if a["bic"] <= b["bic"] + 2.0 else b                    # prefer the simpler model unless the drift model is clearly better
        best["model_compare"] = {"shared shape": a["bic"], "drifting shape": b["bic"]}
        return best
    mdl = _Model(spectra, bands, window, bg_order, drift == "drift", drift == "drift", counts)
    th0 = mdl.p0.copy()
    sol = least_squares(mdl.residual, th0, bounds=(mdl.lb, mdl.ub), x_scale=np.maximum(np.abs(th0) * 0.1, 0.05), diff_step=1e-4, max_nfev=max_nfev,
                        ftol=1e-9, xtol=1e-10, loss="soft_l1", f_scale=3.0)       # robust: one dirty spectrum must not bend the shared shape
    th = sol.x
    n_pts = sum(len(d[0]) for d in mdl.data)
    n_lin = len(mdl.data) * (mdl.nb + mdl.bg_order + 1)
    n_par = len(th) + n_lin
    rss = float(np.sum(sol.fun ** 2))
    bic = n_pts * math.log(max(rss, 1e-300) / n_pts) + n_par * math.log(n_pts)
    rows = []
    for ti, (x, y) in enumerate(mdl.data):
        amp, bg, A, B = mdl.solve(th, ti)
        fit = A @ amp + B @ bg
        res = y - fit
        w = mdl.w[ti]
        dof = max(len(y) - mdl.nb - mdl.bg_order - 1 - 1, 1)
        s2 = float(np.sum(w * res ** 2) / dof)
        Mfull = np.hstack([A, B]) * np.sqrt(w)[:, None]
        try:
            cov = np.linalg.inv(Mfull.T @ Mfull) * s2 * 1.0
        except np.linalg.LinAlgError:
            cov = np.linalg.pinv(Mfull.T @ Mfull) * s2
        # amplitudes were solved on a scaled problem; cov above is already in data units because Mfull uses the real columns
        a0, a1 = amp[0], amp[1]
        var = np.diag(cov)[:mdl.nb]
        r = a0 / a1 if a1 > 0 else float("nan")
        if a0 > 0 and a1 > 0:
            rel2 = var[0] / a0 ** 2 + var[1] / a1 ** 2 - 2 * cov[0, 1] / (a0 * a1)
            r_err = abs(r) * math.sqrt(max(rel2, 0.0))
        else:
            r_err = float("nan")
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - float(np.sum(res ** 2)) / ss_tot if ss_tot > 0 else float("nan")
        rr = res - res.mean()
        dw = float(np.sum(np.diff(res) ** 2) / max(np.sum(res ** 2), 1e-300))         # Durbin-Watson: 2 = uncorrelated residuals
        cc, ww, ee = mdl.band_params(th, ti)
        flag = "good"
        if r2 < 0.98 or dw < 0.8 or (counts and s2 > 3.0) or not np.isfinite(r_err):
            flag = "check"
        if r2 < 0.9 or dw < 0.4 or (counts and s2 > 8.0) or a0 <= 0 or a1 <= 0:
            flag = "poor"
        rows.append({"T": float(mdl.T[ti]), "areas": amp.copy(), "area_err": np.sqrt(np.maximum(var, 0)), "ratio": float(r), "ratio_err": float(r_err),
                     "r2": r2, "chi2_red": s2, "dw": dw, "flag": flag, "x": x, "y": y, "fit": fit, "background": B @ bg,
                     "components": [A[:, k] * amp[k] for k in range(mdl.nb)], "centres": cc, "widths": ww, "etas": ee})
    esd = {}
    try:
        J = sol.jac
        cov_th = np.linalg.inv(J.T @ J) * (rss / max(n_pts - n_par, 1))
        esd = {n: float(math.sqrt(max(cov_th[i, i], 0))) for i, n in enumerate(mdl.names)}
    except Exception:
        pass
    return {"rows": rows, "params": dict(zip(mdl.names, map(float, th))), "esd": esd, "model": "drifting shape" if mdl.drift_pos else "shared shape",
            "bic": float(bic), "window": tuple(map(float, window)), "n_bands": mdl.nb, "bg_order": mdl.bg_order, "rss": rss, "T0": mdl.T0, "Ts": mdl.Ts}


def window_integration(spectra, w1, w2, baseline="linear"):
    """The classic method for comparison: area in two fixed windows above the straight line through each window's ends."""
    I1, I2 = [], []
    for s in spectra:
        x = np.asarray(s["x"], float)
        y = np.asarray(s["y"], float)
        out = []
        for lo, hi in (w1, w2):
            m = (x >= min(lo, hi)) & (x <= max(lo, hi))
            if m.sum() < 3:
                out.append(0.0)
                continue
            xs, v = x[m], y[m]
            if baseline == "linear":
                v = v - np.interp(xs, [xs[0], xs[-1]], [v[0], v[-1]])
            out.append(float(np.trapezoid(v, xs)))
        I1.append(out[0]); I2.append(out[1])
    return np.array(I1), np.array(I2)


def robust_ratio(spectra, bands, window=None, bg_order=1, counts=False, rel_windows=(0.8, 1.25), main=None):
    """Main fit + refits with narrower / wider windows and the other shape model. Returns the main result with 'ratio_sys' (spread between the
    alternatives per spectrum) and 'ratio_total' = √(stat² + sys²)."""
    main = main or fit_global(spectra, bands, window, bg_order, "auto", counts)
    bands = list(main.get("bands", bands)) + list(main.get("extra", []))
    lo, hi = main["window"]
    mid, half = (lo + hi) / 2, (hi - lo) / 2
    alts = []
    other = "drift" if main["model"] == "shared shape" else "none"
    allx = np.concatenate([np.asarray(s["x"], float) for s in spectra])
    for f in rel_windows:
        w = (max(mid - half * f, allx.min()), min(mid + half * f, allx.max()))
        try:
            alts.append(fit_global(spectra, bands, w, bg_order, "none" if main["model"] == "shared shape" else "drift", counts))
        except Exception:
            pass
    try:
        alts.append(fit_global(spectra, bands, main["window"], bg_order, other, counts))
    except Exception:
        pass
    R = np.array([r["ratio"] for r in main["rows"]])
    stat = np.array([r["ratio_err"] for r in main["rows"]])
    if alts:
        A = np.array([[r["ratio"] for r in a["rows"]] for a in alts])
        sys_ = np.sqrt(np.nanmean((A - R[None, :]) ** 2, axis=0))
    else:
        sys_ = np.zeros_like(R)
    main["ratio_sys"] = sys_
    main["ratio_total"] = np.sqrt(np.nan_to_num(stat) ** 2 + sys_ ** 2)
    main["n_alternatives"] = len(alts)
    return main


def suggest_extra_band(res, min_sigma=4.0, exclude=()):
    """Looks for a clear structure left in the residual (summed over all spectra) and proposes an extra band: (centre, fwhm) or None."""
    xs = np.unique(np.round(np.concatenate([r["x"] for r in res["rows"]]), 3))
    acc = np.zeros_like(xs)
    cnt = np.zeros_like(xs)
    for r in res["rows"]:
        d = r["y"] - r["fit"]
        s = max(float(np.std(d)), 1e-12)
        acc += np.interp(xs, r["x"], d / s)
        cnt += 1
    z = acc / max(cnt.max(), 1) * math.sqrt(len(res["rows"]))
    step = float(np.median(np.diff(xs))) if len(xs) > 1 else 1.0
    zp = np.concatenate([[z.min() - 1.0], z, [z.min() - 1.0]])              # pad: a structure at the edge of the window is a peak too
    pk, pr = find_peaks(zp, height=min_sigma, distance=max(3, int(1.0 / max(step, 1e-6))))
    pk = pk - 1
    pk = pk[(pk >= 0) & (pk < len(z))]
    if not len(pk):
        return None
    pk = [q for q in pk if not any(abs(xs[q] - c) < 0.7 * w for c, w in exclude)]       # not inside the cores of the bands of the ratio
    if not len(pk):
        return None
    i = pk[int(np.argmax(z[pk]))]
    half = z[i] / 2
    a = i
    while a > 0 and z[a] > half:
        a -= 1
    b = i
    while b < len(z) - 1 and z[b] > half:
        b += 1
    return float(xs[i]), float(max(xs[b] - xs[a], 3 * step))


def auto_fit(spectra, counts=False, lo=None, hi=None, bg_order=1, max_extra=3, bands=None):
    """Complete automatic analysis: find the two bands (also hidden shoulders), fit the series, then add up to `max_extra` nuisance bands
    while a clear structure remains in the residual AND the BIC improves. Returns the fit result with 'bands' (centre, fwhm) and 'extra' (list)."""
    x = np.concatenate([np.asarray(s["x"], float) for s in spectra])
    grid = np.asarray(spectra[0]["x"], float)
    mean = np.mean([np.interp(grid, np.asarray(s["x"], float), np.asarray(s["y"], float)) for s in spectra], axis=0)
    bands = list(bands) if bands else auto_bands(grid, mean, lo, hi)
    if not bands or len(bands) < 2:
        return None
    bands = bands[:2]
    win = default_window(bands, x)
    if lo is not None and hi is not None:                          # the user's / full spectral range is the fit window (tails and neighbours matter)
        win = (float(max(lo, x.min())), float(min(hi, x.max())))
    best = fit_global(spectra, bands, win, bg_order, "auto", counts)
    extra = []
    for _ in range(max_extra):
        sug = suggest_extra_band(best, exclude=bands + extra)
        if sug is None or any(abs(sug[0] - c) < 0.6 * w for c, w in bands + extra):
            break
        trial = fit_global(spectra, bands + extra + [sug], win, bg_order, "auto", counts)
        if trial["bic"] < best["bic"] - 10.0:
            best, extra = trial, extra + [sug]
        else:
            break
    best["bands"], best["extra"] = bands, extra
    return best
