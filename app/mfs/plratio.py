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


def pv_matrix(x, c, w, eta):
    """Same as pv_unit for several bands at once: returns an array (len(x), number of bands)."""
    x = np.asarray(x, float)[:, None]
    c = np.asarray(c, float)[None, :]
    w = np.asarray(w, float)[None, :]
    eta = np.asarray(eta, float)[None, :]
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
        self.var = []
        for s in spectra:
            x = np.asarray(s["x"], float)
            y = np.asarray(s["y"], float)
            m = (x >= lo) & (x <= hi) & np.isfinite(y)
            self.data.append((x[m], y[m]))
            # v2.18.4: after a background subtraction the Poisson variance is that of the RAW counts (signal + background): spectra may carry it as "var"
            self.var.append(np.asarray(s["var"], float)[m] if s.get("var") is not None else y[m])
        self.counts = counts
        self.w = [1.0 / np.maximum(v, 1.0) if counts else np.ones_like(y) / max(_noise(y), 1e-12) ** 2 for (x, y), v in zip(self.data, self.var)]
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
        A = pv_matrix(x, c, w, e)
        xn = (x - x.mean()) / max(np.ptp(x), 1e-9) * 2
        cols_bg = [np.polynomial.chebyshev.chebval(xn, [0] * j + [1]) for j in range(self.bg_order + 1)]
        return A, np.stack(cols_bg, 1)

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


def _group_ratio(amp, cov, g1, g2):
    """Ratio of the summed areas of two groups of bands and its standard error with the full covariance of the areas."""
    s1, s2 = float(np.sum(amp[g1])), float(np.sum(amp[g2]))
    v1 = float(np.sum(cov[np.ix_(g1, g1)]))
    v2 = float(np.sum(cov[np.ix_(g2, g2)]))
    c12 = float(np.sum(cov[np.ix_(g1, g2)]))
    r = s1 / s2 if s2 > 0 else float("nan")
    if s1 > 0 and s2 > 0:
        rel2 = v1 / s1 ** 2 + v2 / s2 ** 2 - 2 * c12 / (s1 * s2)
        return r, abs(r) * math.sqrt(max(rel2, 0.0)), s1, s2, math.sqrt(max(v1, 0.0)), math.sqrt(max(v2, 0.0))
    return r, float("nan"), s1, s2, float("nan"), float("nan")


def fit_global(spectra, bands, window=None, bg_order=1, drift="auto", counts=False, max_nfev=80, groups=None):
    """Global fit of the whole series. spectra = [{'T': K, 'x': nm, 'y': counts}], bands = [(centre, fwhm_guess), ...] with the ratio pair FIRST
    (I1 = bands[0], I2 = bands[1]). drift = 'none' | 'drift' | 'auto' (BIC decides). Returns the result dict (see keys below).
    v2.18.4: groups = ([band indices of I1], [band indices of I2]) - the ratio of SUMMED areas, for bands that are split into several lines (Stark
    components of a lanthanide manifold): every line has its own shape, the ratio uses the sum of each group. Default ([0], [1])."""
    window = window or default_window(bands, np.concatenate([np.asarray(s["x"], float) for s in spectra]))
    if drift == "auto":
        a = fit_global(spectra, bands, window, bg_order, "none", counts, max_nfev, groups)
        b = fit_global(spectra, bands, window, bg_order, "drift", counts, max_nfev, groups)
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
    # v2.18.5: covariance of the NON-LINEAR shape parameters (centres, widths, eta, drifts) from the Jacobian of the projected residual.
    # Their uncertainty propagates into every ratio through its sensitivity to them; without this term the ratio errors were 10-20 % too small.
    try:
        J_ = sol.jac
        cov_th = np.linalg.pinv(J_.T @ J_) * (rss / max(n_pts - n_par, 1))
    except Exception:
        cov_th = None
    g1_, g2_ = (list(groups[0]), list(groups[1])) if groups else ([0], [1])

    def ratio_at(th_, ti_):
        a_, _, _, _ = mdl.solve(th_, ti_)
        s1_, s2_ = float(np.sum(a_[g1_])), float(np.sum(a_[g2_]))
        return s1_ / s2_ if s2_ > 0 else float("nan")
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
        var = np.diag(cov)[:mdl.nb]
        g1, g2 = (list(groups[0]), list(groups[1])) if groups else ([0], [1])
        r, r_err, a0, a1, e0, e1 = _group_ratio(amp, cov[:mdl.nb, :mdl.nb], g1, g2)
        if cov_th is not None and np.isfinite(r_err):                 # add the shape-parameter term: (dR/dtheta)^T Cov(theta) (dR/dtheta)
            grad = np.zeros(len(th))
            for k in range(len(th)):
                h_ = 1e-5 * max(1.0, abs(th[k]))
                tp_, tm_ = th.copy(), th.copy()
                tp_[k] += h_; tm_[k] -= h_
                grad[k] = (ratio_at(tp_, ti) - ratio_at(tm_, ti)) / (2 * h_)
            var_shape = float(grad @ cov_th @ grad)
            r_err = math.sqrt(r_err ** 2 + max(var_shape, 0.0))
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
        areas_out = amp.copy()
        if groups:                                                  # the two ratio 'areas' are the group sums; the single lines stay in 'line_areas'
            areas_out = np.concatenate([[a0, a1], amp])
        rows.append({"T": float(mdl.T[ti]), "areas": areas_out, "line_areas": amp.copy(), "area_err": np.sqrt(np.maximum(var, 0)), "group_err": (e0, e1),
                     "ratio": float(r), "ratio_err": float(r_err),
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
            "bic": float(bic), "window": tuple(map(float, window)), "n_bands": mdl.nb, "bg_order": mdl.bg_order, "rss": rss, "T0": mdl.T0, "Ts": mdl.Ts,
            "groups": (list(groups[0]), list(groups[1])) if groups else None, "drift": "drift" if mdl.drift_pos else "none"}


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


def robust_ratio(spectra, bands, window=None, bg_order=1, counts=False, rel_windows=(0.8, 1.25), main=None, bg_alternatives=None):
    """Main fit + refits with narrower / wider windows and the other shape model. Returns the main result with 'ratio_sys' (spread between the
    alternatives per spectrum) and 'ratio_total' = √(stat² + sys²).
    v2.18.4: bg_alternatives = [(label, spectra)] - the same series after OTHER reasonable background subtractions (see plbackground.alternatives); each is refitted
    with the same bands, window and shape model, and the spread of the ratio over them is the background part of the systematic error ('ratio_sys_bg';
    'ratio_sys_fit' is the part from windows and shape model)."""
    main = main or fit_global(spectra, bands, window, bg_order, "auto", counts)
    bands = list(main.get("bands", bands)) + list(main.get("extra", []))
    groups = main.get("groups")
    lo, hi = main["window"]
    mid, half = (lo + hi) / 2, (hi - lo) / 2
    alts = []
    other = "drift" if main["model"] == "shared shape" else "none"
    allx = np.concatenate([np.asarray(s["x"], float) for s in spectra])
    win_list = [(max(mid - half * f, allx.min()), min(mid + half * f, allx.max())) for f in rel_windows]
    win_labels = [f"fit range x{f}" for f in rel_windows]
    if main.get("ranges"):                                   # manifold: the two ranges must stay whole; only the margin outside them changes (7 or 25 nm instead of 15)
        r1, r2 = main["ranges"]
        a, b = min(r1[0], r2[0]), max(r1[1], r2[1])
        win_list = [(max(a - pad, allx.min()), min(b + pad, allx.max())) for pad in (7.0, 25.0)]
        win_labels = ["margin 7 nm around the ranges", "margin 25 nm around the ranges"]
    used_labels = []
    for w, lab in zip(win_list, win_labels):
        try:
            bb = [bd for bd in bands if w[0] <= bd[0] <= w[1]] if groups else bands
            if groups:
                g1 = [k for k, (c, _) in enumerate(bb) if min(main["ranges"][0]) <= c <= max(main["ranges"][0])] if main.get("ranges") else groups[0]
                g2 = [k for k, (c, _) in enumerate(bb) if min(main["ranges"][1]) <= c <= max(main["ranges"][1])] if main.get("ranges") else groups[1]
                gg = (g1, g2)
            else:
                gg = None
            alts.append(fit_global(spectra, bb, w, bg_order, "none" if main["model"] == "shared shape" else "drift", counts, 80, gg))
            used_labels.append(lab)
        except Exception:
            pass
    try:
        alts.append(fit_global(spectra, bands, main["window"], bg_order, other, counts, 80, groups))
    except Exception:
        pass
    R = np.array([r["ratio"] for r in main["rows"]])
    stat = np.array([r["ratio_err"] for r in main["rows"]])
    if alts:
        A = np.array([[r["ratio"] for r in a["rows"]] for a in alts])
        sys_ = np.sqrt(np.nanmean((A - R[None, :]) ** 2, axis=0))
    else:
        sys_ = np.zeros_like(R)
    sys_bg = np.zeros_like(R)
    bg_rows, bg_labels = [], []
    for lab, sp_alt in (bg_alternatives or []):
        try:
            a = fit_global(sp_alt, bands, main["window"], bg_order, "none" if main["model"] == "shared shape" else "drift", counts, 40, groups)
            if len(a["rows"]) == len(R):
                bg_rows.append([r["ratio"] for r in a["rows"]])
                bg_labels.append(lab)
        except Exception:
            pass
    if bg_rows:
        Ab = np.array(bg_rows)
        sys_bg = np.sqrt(np.nanmean((Ab - R[None, :]) ** 2, axis=0))
    main["ratio_sys_fit"] = sys_
    main["ratio_sys_bg"] = sys_bg
    main["bg_alt_labels"] = bg_labels
    main["bg_alt_ratios"] = np.array(bg_rows) if bg_rows else None
    main["fit_alt_ratios"] = A if alts else None                      # v2.18.4: kept so that the Boltzmann fit can be repeated on every alternative
    main["fit_alt_labels"] = used_labels + (["other peak model"] if len(alts) > len(used_labels) else [])
    main["ratio_sys"] = np.sqrt(sys_ ** 2 + sys_bg ** 2)
    main["ratio_total"] = np.sqrt(np.nan_to_num(stat) ** 2 + main["ratio_sys"] ** 2)
    main["n_alternatives"] = len(alts) + len(bg_rows)
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


# ------------------------------------------------------------------------------------------------------------- v2.18.4: crowded spectra (manifolds)
def auto_lines(x, y, lo=None, hi=None, max_lines=16, min_rel=0.02):
    """ALL lines of a (mean) spectrum between lo and hi: maxima and shoulders (minima of the second derivative), merged when closer than ~1.5 nm,
    weaker than `min_rel` of the strongest dropped. Returns [(centre, fwhm_guess)] sorted by centre (at most max_lines, strongest kept)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    lo = x.min() if lo is None else lo
    hi = x.max() if hi is None else hi
    m = (x >= lo) & (x <= hi)
    xs, ys = x[m], y[m]
    if len(xs) < 15:
        return []
    from scipy.signal import savgol_filter
    step = float(np.median(np.diff(xs)))
    win = max(5, int(round(2.0 / step)) | 1)
    ysm = savgol_filter(ys, min(win, len(ys) - (1 - len(ys) % 2)), 2) if len(ys) > win else ys
    base = float(np.percentile(ysm, 5))
    rng = float(ysm.max() - base)
    if rng <= 0:
        return []
    w2 = max(7, int(round(2.5 / step)) | 1)
    d2 = savgol_filter(ys, min(w2, len(ys) - (1 - len(ys) % 2)), 3, deriv=2)
    cand = list(find_peaks(ysm, prominence=min_rel * rng)[0]) + list(find_peaks(-d2, prominence=0.02 * float(np.ptp(d2)))[0])
    cand = sorted(set(int(i) for i in cand if ysm[int(i)] - base > min_rel * rng))
    merged = []
    for i in cand:
        if merged and xs[i] - xs[merged[-1]] < max(1.5, 2.5 * step):
            if ysm[i] > ysm[merged[-1]]:
                merged[-1] = i
        else:
            merged.append(i)
    merged = sorted(sorted(merged, key=lambda i: -ysm[i])[:max_lines])
    out = []
    for k, p_ in enumerate(merged):
        gaps = [abs(xs[p_] - xs[q]) for q in merged if q != p_]
        h = ysm[p_]
        half = h - (h - base) / 2.0
        i0, i1 = p_, p_
        while i0 > 0 and ysm[i0] > half:
            i0 -= 1
        while i1 < len(ysm) - 1 and ysm[i1] > half:
            i1 += 1
        wd = max(float(xs[i1] - xs[i0]), 2.5 * step)
        if gaps:
            wd = min(wd, 1.1 * min(gaps))
        out.append((float(xs[p_]), float(max(wd, 2.5 * step))))
    return out


def fit_manifolds(spectra, range1, range2, counts=False, lo=None, hi=None, bg_order=1, max_extra=4, lines=None):
    """Ratio of two MANIFOLDS: every line found between lo and hi gets its own pseudo-Voigt (shape shared by the series, optionally drifting); the lines
    whose centre lies in range1 are summed into I1, those in range2 into I2, the others are nuisance lines. Lines left in the residual are added while the
    BIC improves. Returns the fit result with 'bands' (all lines), 'groups', 'extra' = [] and 'mode' = 'manifold', or None."""
    grid = np.asarray(spectra[0]["x"], float)
    mean = np.mean([np.interp(grid, np.asarray(s["x"], float), np.asarray(s["y"], float)) for s in spectra], axis=0)
    xall = np.concatenate([np.asarray(s["x"], float) for s in spectra])
    lo = float(max(lo if lo is not None else xall.min(), xall.min()))
    hi = float(min(hi if hi is not None else xall.max(), xall.max()))
    lines = list(lines) if lines else auto_lines(grid, mean, lo, hi)

    def grp(ls):
        g1 = [k for k, (c, _) in enumerate(ls) if min(range1) <= c <= max(range1)]
        g2 = [k for k, (c, _) in enumerate(ls) if min(range2) <= c <= max(range2)]
        return g1, g2
    g1, g2 = grp(lines)
    if not g1 or not g2:
        return None
    best = fit_global(spectra, lines, (lo, hi), bg_order, "auto", counts, 120, (g1, g2))
    for _ in range(max_extra):
        sug = suggest_extra_band(best, min_sigma=4.0, exclude=[(c, 0.45 * w) for c, w in lines])       # lines can be closer than a width: only the cores are excluded
        if sug is None or any(abs(sug[0] - c) < 0.3 * w for c, w in lines):
            break
        trial_lines = sorted(lines + [sug])
        t1, t2 = grp(trial_lines)
        trial = fit_global(spectra, trial_lines, (lo, hi), bg_order, "auto", counts, 120, (t1, t2))
        if trial["bic"] < best["bic"] - 10.0:
            best, lines, g1, g2 = trial, trial_lines, t1, t2
        else:
            break
    best["bands"], best["extra"], best["groups"], best["mode"] = lines, [], (g1, g2), "manifold"
    best["ranges"] = (tuple(sorted(range1)), tuple(sorted(range2)))
    return best


def flag_share(fit):
    """Fraction of the spectra of a fit result that are not 'good'."""
    rows = fit.get("rows") or []
    return float(np.mean([r["flag"] != "good" for r in rows])) if rows else 1.0
