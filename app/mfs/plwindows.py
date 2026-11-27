"""PL thermometry on user-chosen windows: a JOINT fit of the two windows that uses only the data inside them.

Model (per spectrum): both bands are present in both windows. Band 1 (centre c1, FWHM w1, area A1) and band 2 (c2, w2, A2) are
unit-area pseudo-Voigts with the Lorentzian share fixed (ETA). Each window has its own straight background line. Because the
neighbour band is part of the model, its tail inside a window is fitted, not ignored. The ratio is R = A1 / A2.

Errors: the covariance of the joint least-squares fit (Jacobian; Poisson weights sigma = sqrt(max(y, 1)) when counts are given),
scaled up by the reduced chi-square when it is above 1. The ratio error uses the full covariance of A1 and A2:
  var(R) = R^2 (sA1^2/A1^2 + sA2^2/A2^2 - 2 cov(A1, A2) / (A1 A2)).

Limits, reported as flags per temperature:
  * 'not significant': a band area is below SIGNIFICANCE x its error; the ratio is then undefined (NaN).
  * 'tail outside measured range': a band centre +- 3 FWHM extends beyond the x-range the fit uses (the union of the two windows). The
    neighbour's tail that falls inside the union is modelled and is not flagged; a tail beyond the union is an extrapolation and is.
  * 'poor fit': reduced chi-square above 2.
  * 'edge': a fitted centre lies on its window edge.
  * 'few points': fewer than MIN_POINTS data points in a window.
  * 'no area': a fitted area is not positive.
What the model still assumes: one pseudo-Voigt shape per band with ETA fixed, and a straight background per window. The synthetic tests
generate data with the same shape, so they do NOT test a shape mismatch in real spectra; the residual error of real spectra is not known.
Not checked here: the energy axis, the background outside the windows, the temperature calibration.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import least_squares

from mfs.plratio import pv_unit

ETA = 0.5          # Lorentzian share of the pseudo-Voigt; fixed (see the module docstring)
MIN_POINTS = 12
SIGNIFICANCE = 3.0  # a band area must exceed this many of its own errors to give a ratio
# S-Figuryte-1: minimum mean counts per point in a window. Below it the ratio is not given ("too little signal").
# Measured on our own synthetic Poisson series (same two-band model, 60 seeds per level, window-1 mean counts):
#   mean counts  1.5   2.5   3.5   5.5   8.5   10.5  15.5  30.4  50.6   100.2  200.5
#   bias (%)     +63   +60  +108   +32  +26   +23   +12   +3.9  +2.0   -0.03  -0.56
#   z-std        1.7   0.9   0.9   0.8   1.0   1.0   0.9   1.0   1.1    0.9    1.2   (1.0 = error bars calibrated)
# 100 is the first level where the bias is within its own standard error (+-0.5 %) and the error bars are calibrated.
# This is a property of the synthetic case; real spectra may need a different value. Pass min_mean_counts to change it.
MIN_MEAN_COUNTS = 100.0


def _check_window(w, name):
    if w is None or len(w) != 2:
        raise ValueError(f"{name}: give the window as (low, high) in nm")
    lo, hi = float(min(w)), float(max(w))
    if not (math.isfinite(lo) and math.isfinite(hi)) or hi <= lo:
        raise ValueError(f"{name}: the window {w} is empty or not finite")
    return lo, hi


def _fail(n, msg):
    nan = float("nan")
    return {"A1": nan, "A2": nan, "eA1": nan, "eA2": nan, "cov12": nan, "c1": nan, "c2": nan, "w1": nan, "w2": nan,
            "chi2_red": nan, "n": n, "flags1": [msg], "flags2": [msg], "mean1": nan, "mean2": nan}


def fit_pair(x, y, w1, w2, var=None, counts=True, free_eta=False, weighting="sqrt", warm=None):
    """Joint fit of the two windows of one spectrum. Only points inside w1 or w2 are used.
    free_eta=False (default, eta fixed at ETA = 0.5); free_eta=True (opt-in): the Lorentzian share eta of both bands is one shared fit parameter in [0, 1]."""
    lo1, hi1 = _check_window(w1, "window 1")
    lo2, hi2 = _check_window(w2, "window 2")
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m1 = (x >= lo1) & (x <= hi1) & np.isfinite(y)
    m2 = (x >= lo2) & (x <= hi2) & np.isfinite(y)
    if m1.sum() < MIN_POINTS or m2.sum() < MIN_POINTS:
        return _fail(int(m1.sum() + m2.sum()), "few points")
    m = m1 | m2
    xs, ys = x[m], y[m]
    sel1 = m1[m]
    sel2 = m2[m]
    if var is not None:
        sig = np.sqrt(np.maximum(np.asarray(var, float)[m], 1.0))
    elif counts:
        sig = np.sqrt(np.maximum(ys, 1.0))
    else:
        sig = np.full_like(ys, max(float(np.std(np.diff(ys))) / math.sqrt(2.0), 1e-12))
    step = float(np.median(np.diff(xs)))
    span1, span2 = 0.5 * (hi1 - lo1), 0.5 * (hi2 - lo2)
    xm1, xm2 = 0.5 * (lo1 + hi1), 0.5 * (lo2 + hi2)
    lin1 = np.where(sel1, (xs - xm1) / span1, 0.0)
    lin2 = np.where(sel2, (xs - xm2) / span2, 0.0)

    def model(p):
        A1, c1, u1, A2, c2, u2, o1, s1, o2, s2 = p[:10]
        eta = float(p[10]) if free_eta else ETA
        bg = np.where(sel1, o1 + s1 * lin1, o2 + s2 * lin2)
        return A1 * pv_unit(xs, c1, u1, eta) + A2 * pv_unit(xs, c2, u2, eta) + bg

    def resid(p):
        if weighting == "poisson":
            # S-Figuryte-1 (opt-in): Poisson deviance residuals. Their sum of squares is the Poisson deviance, so the least-squares
            # minimum is the Poisson maximum-likelihood fit of the same model (no sqrt weights, no low-count weighting bias).
            mu = np.maximum(model(p), 1e-9)
            term = np.where(ys > 0, ys * np.log(np.maximum(ys, 1e-300) / mu), 0.0)
            d = np.maximum(2.0 * (mu - ys + term), 0.0)
            return np.sign(ys - mu) * np.sqrt(d)
        return (model(p) - ys) / sig

    def start(sel, lo, hi):
        xi, yi = xs[sel], ys[sel]
        bg = np.interp(xi, [xi[0], xi[-1]], [yi[0], yi[-1]])
        res = yi - bg
        i = int(np.argmax(res))
        return float(xi[i]), max(float(np.sum(res) * step), 1e-9), float(bg[0])

    c1_0, A1_0, o1_0 = start(sel1, lo1, hi1)
    c2_0, A2_0, o2_0 = start(sel2, lo2, hi2)
    wstart1 = max((hi1 - lo1) / 6.0, 3 * step)
    wstart2 = max((hi2 - lo2) / 6.0, 3 * step)
    p0 = [A1_0, c1_0, wstart1, A2_0, c2_0, wstart2, o1_0, 0.0, o2_0, 0.0]
    lb = [0.0, lo1, 1.5 * step, 0.0, lo2, 1.5 * step, -np.inf, -np.inf, -np.inf, -np.inf]
    ub = [np.inf, hi1, hi1 - lo1, np.inf, hi2, hi2 - lo2, np.inf, np.inf, np.inf, np.inf]
    if free_eta:
        p0 += [ETA]
        lb += [0.0]
        ub += [1.0]
    if warm is not None and len(warm) == len(p0):
        # S-Figuryte-2: warm start from the previous temperature of the series; clipped to the bounds. Same model, same minimum.
        p0 = [min(max(float(v), lo_), hi_) for v, lo_, hi_ in zip(warm, lb, ub)]
    try:
        r = least_squares(resid, p0, bounds=(lb, ub), x_scale="jac", max_nfev=4000)
    except Exception as exc:                                 # noqa: BLE001 - reported, not hidden
        return _fail(int(len(xs)), f"fit failed: {exc}")
    p = r.x
    dof = max(len(xs) - len(p), 1)
    chi2_red = float(np.sum(r.fun ** 2) / dof)
    try:
        cov = np.linalg.inv(r.jac.T @ r.jac) * max(chi2_red, 1.0)
    except np.linalg.LinAlgError:
        cov = np.full((len(p), len(p)), np.nan)
    A1, c1, u1, A2, c2, u2 = (float(v) for v in p[:6])
    eA1 = float(math.sqrt(abs(cov[0, 0]))) if np.isfinite(cov[0, 0]) else float("nan")
    eA2 = float(math.sqrt(abs(cov[3, 3]))) if np.isfinite(cov[3, 3]) else float("nan")
    cov12 = float(cov[0, 3]) if np.isfinite(cov[0, 3]) else float("nan")

    # the neighbour's tail is modelled inside the measured range, so the tail warning is about the band extending BEYOND the data
    data_lo, data_hi = float(xs.min()), float(xs.max())

    def flags_for(A, eA, c, u, lo, hi):
        f = []
        if not (A > 0):
            f.append("no area")
        elif not (A > SIGNIFICANCE * eA):
            f.append("not significant")
        if (c - lo) < 1e-6 * (hi - lo) or (hi - c) < 1e-6 * (hi - lo):
            f.append("edge")
        if c - 3 * u < data_lo or c + 3 * u > data_hi:
            f.append("tail outside measured range")
        if chi2_red > 2.0:
            f.append("poor fit")
        return f

    return {"A1": A1, "A2": A2, "eA1": eA1, "eA2": eA2, "cov12": cov12, "c1": c1, "c2": c2, "w1": u1, "w2": u2,
            "eta": float(p[10]) if free_eta else ETA, "chi2_red": chi2_red, "n": int(len(xs)), "p": np.asarray(p, float),
            "mean1": float(np.mean(ys[sel1])) if sel1.any() else float("nan"),
            "mean2": float(np.mean(ys[sel2])) if sel2.any() else float("nan"),
            "flags1": flags_for(A1, eA1, c1, u1, lo1, hi1), "flags2": flags_for(A2, eA2, c2, u2, lo2, hi2)}


def _ratio_and_error(A1, eA1, A2, eA2, cov12, ok):
    """S-Figuryte-2: the one place where the ratio A1/A2 and its error are formed. Returns (ratio, ratio_err, usable)."""
    if ok and np.isfinite(A1) and np.isfinite(A2) and A1 > 0 and A2 > 0 and np.isfinite(eA1) and np.isfinite(eA2):
        R = A1 / A2
        v = (eA1 / A1) ** 2 + (eA2 / A2) ** 2 - 2.0 * cov12 / (A1 * A2)
        return R, (R * math.sqrt(max(v, 0.0)) if np.isfinite(v) else float("nan")), True
    return float("nan"), float("nan"), False


def window_ratio_series(spectra, w1, w2, counts=True, free_eta=False, min_mean_counts=MIN_MEAN_COUNTS, weighting="sqrt", warm_start=False):
    """Ratio I1/I2 per spectrum from the joint two-window fit.

    spectra = [{'T': K, 'x': nm, 'y': counts, optional 'var': array}, ...]. Returns one dict per spectrum (same order):
    T, I1, I1_err, I2, I2_err, ratio, ratio_err, flags. The ratio is NaN when a band area is not significant or not positive.
    """
    lo1, hi1 = _check_window(w1, "window 1")
    lo2, hi2 = _check_window(w2, "window 2")
    if hi1 > lo2 and hi2 > lo1:
        raise ValueError("the two windows overlap; a band must not be counted in both")
    out = []
    prev_p = None
    for s in spectra:
        f = fit_pair(s["x"], s["y"], (lo1, hi1), (lo2, hi2), s.get("var"), counts, free_eta, weighting,
                     warm=prev_p if warm_start else None)
        if warm_start and np.isfinite(f["A1"]):
            prev_p = f.get("p")
        flags = [f"I1: {x}" for x in f["flags1"]] + [f"I2: {x}" for x in f["flags2"]]
        if counts and not (f["mean1"] >= min_mean_counts):
            flags.append("I1: too little signal (mean counts per point below min_mean_counts, default 100)")
        if counts and not (f["mean2"] >= min_mean_counts):
            flags.append("I2: too little signal (mean counts per point below min_mean_counts, default 100)")
        blocking = ("no area", "not significant", "few points", "fit failed", "too little signal")
        ok = not any(any(b in x for b in blocking) for x in flags)
        ratio, ratio_err, usable = _ratio_and_error(f["A1"], f["eA1"], f["A2"], f["eA2"], f["cov12"], ok)
        if not usable:
            flags.append("ratio undefined")
        out.append({"T": float(s["T"]), "I1": f["A1"], "I1_err": f["eA1"], "I2": f["A2"], "I2_err": f["eA2"],
                    "ratio": ratio, "ratio_err": ratio_err, "flags": flags})
    return out


def window_ratio_series_shared(spectra, w1, w2, counts=True, free_eta=False, min_mean_counts=MIN_MEAN_COUNTS):
    """S-PL-2: ONE joint least-squares fit of the whole temperature series. The band widths (FWHM 1, FWHM 2) and eta are shared by all
    spectra; each spectrum keeps its own areas A1, A2, centres c1, c2 and its own two background lines. Start values come from the
    separate per-spectrum fits (fit_pair), so the joint fit starts warm. Ratio error: the joint covariance of each spectrum's own
    (A1, A2) block. Assumes the band widths do NOT change with temperature; if they do, the shared width biases the ratios (see tests).
    Returns the same keys as window_ratio_series."""
    lo1, hi1 = _check_window(w1, "window 1")
    lo2, hi2 = _check_window(w2, "window 2")
    if hi1 > lo2 and hi2 > lo1:
        raise ValueError("the two windows overlap; a band must not be counted in both")
    specs = []
    for s in spectra:
        x = np.asarray(s["x"], float)
        y = np.asarray(s["y"], float)
        m1 = (x >= lo1) & (x <= hi1) & np.isfinite(y)
        m2 = (x >= lo2) & (x <= hi2) & np.isfinite(y)
        m = m1 | m2
        xs, ys = x[m], y[m]
        sig = np.sqrt(np.maximum(np.asarray(s["var"], float)[m], 1.0)) if s.get("var") is not None else (
            np.sqrt(np.maximum(ys, 1.0)) if counts else np.full_like(ys, max(float(np.std(np.diff(ys))) / math.sqrt(2.0), 1e-12)))
        specs.append({"T": float(s["T"]), "xs": xs, "ys": ys, "sig": sig, "sel1": m1[m], "sel2": m2[m], "n": int(m.sum()),
                      "step": float(np.median(np.diff(xs))) if m.sum() > 2 else 1.0})
    starts = [fit_pair(s["x"], s["y"], (lo1, hi1), (lo2, hi2), s.get("var"), counts, free_eta) for s in spectra]
    ok = [f for f in starts if np.isfinite(f["A1"]) and np.isfinite(f["A2"])]
    if not ok:
        return [{"T": sp["T"], "I1": float("nan"), "I1_err": float("nan"), "I2": float("nan"), "I2_err": float("nan"),
                 "ratio": float("nan"), "ratio_err": float("nan"), "flags": ["shared fit: no valid start"]} for sp in specs]
    u1_0 = float(np.median([f["w1"] for f in ok]))
    u2_0 = float(np.median([f["w2"] for f in ok]))
    eta_0 = float(np.median([f["eta"] for f in ok]))
    n_sh = 3 if free_eta else 2
    step = min(sp["step"] for sp in specs)
    p0, lb, ub = [u1_0, u2_0], [1.5 * step, 1.5 * step], [hi1 - lo1, hi2 - lo2]
    if free_eta:
        p0 += [eta_0]
        lb += [0.0]
        ub += [1.0]
    for sp, f in zip(specs, starts):
        A1 = f["A1"] if np.isfinite(f["A1"]) else 1.0
        A2 = f["A2"] if np.isfinite(f["A2"]) else 1.0
        c1 = f["c1"] if np.isfinite(f["c1"]) else 0.5 * (lo1 + hi1)
        c2 = f["c2"] if np.isfinite(f["c2"]) else 0.5 * (lo2 + hi2)
        p0 += [max(A1, 1e-9), c1, max(A2, 1e-9), c2, 0.0, 0.0, 0.0, 0.0]
        lb += [0.0, lo1, 0.0, lo2, -np.inf, -np.inf, -np.inf, -np.inf]
        ub += [np.inf, hi1, np.inf, hi2, np.inf, np.inf, np.inf, np.inf]
    xm1, xm2 = 0.5 * (lo1 + hi1), 0.5 * (lo2 + hi2)
    span1, span2 = 0.5 * (hi1 - lo1), 0.5 * (hi2 - lo2)

    def model_k(p, k, sp):
        base = n_sh + 8 * k
        u1, u2 = p[0], p[1]
        eta = p[2] if free_eta else ETA
        A1, c1, A2, c2, o1, s1, o2, s2 = p[base:base + 8]
        lin1 = np.where(sp["sel1"], (sp["xs"] - xm1) / span1, 0.0)
        lin2 = np.where(sp["sel2"], (sp["xs"] - xm2) / span2, 0.0)
        bg = np.where(sp["sel1"], o1 + s1 * lin1, o2 + s2 * lin2)
        return A1 * pv_unit(sp["xs"], c1, u1, eta) + A2 * pv_unit(sp["xs"], c2, u2, eta) + bg

    def resid(p):
        return np.concatenate([(model_k(p, k, sp) - sp["ys"]) / sp["sig"] for k, sp in enumerate(specs)])

    try:
        r = least_squares(resid, np.array(p0, float), bounds=(np.array(lb, float), np.array(ub, float)), x_scale="jac", max_nfev=4000)
    except Exception as exc:                                  # noqa: BLE001 - reported, not hidden
        return [{"T": sp["T"], "I1": float("nan"), "I1_err": float("nan"), "I2": float("nan"), "I2_err": float("nan"),
                 "ratio": float("nan"), "ratio_err": float("nan"), "flags": [f"shared fit failed: {exc}"]} for sp in specs]
    p = r.x
    n_tot = sum(sp["n"] for sp in specs)
    dof = max(n_tot - len(p), 1)
    chi2_red = float(np.sum(r.fun ** 2) / dof)
    try:
        C = np.linalg.inv(r.jac.T @ r.jac) * max(chi2_red, 1.0)
    except np.linalg.LinAlgError:
        C = np.full((len(p), len(p)), np.nan)
    out = []
    for k, sp in enumerate(specs):
        base = n_sh + 8 * k
        A1, A2 = float(p[base]), float(p[base + 2])
        c1, c2 = float(p[base + 1]), float(p[base + 3])
        eA1 = float(math.sqrt(abs(C[base, base]))) if np.isfinite(C[base, base]) else float("nan")
        eA2 = float(math.sqrt(abs(C[base + 2, base + 2]))) if np.isfinite(C[base + 2, base + 2]) else float("nan")
        cov12 = float(C[base, base + 2]) if np.isfinite(C[base, base + 2]) else float("nan")
        flags = []
        for tag, A, eA in (("I1", A1, eA1), ("I2", A2, eA2)):
            if not (A > 0):
                flags.append(f"{tag}: no area")
            elif not (A > SIGNIFICANCE * eA):
                flags.append(f"{tag}: not significant")
        data_lo, data_hi = float(sp["xs"].min()), float(sp["xs"].max())
        if c1 - 3 * p[0] < data_lo or c1 + 3 * p[0] > data_hi:
            flags.append("I1: tail outside measured range")
        if c2 - 3 * p[1] < data_lo or c2 + 3 * p[1] > data_hi:
            flags.append("I2: tail outside measured range")
        if chi2_red > 2.0:
            flags.append("poor fit (joint series)")
        if counts:
            if not (float(np.mean(sp["ys"][sp["sel1"]])) >= min_mean_counts):
                flags.append("I1: too little signal (mean counts per point below min_mean_counts, default 100)")
            if not (float(np.mean(sp["ys"][sp["sel2"]])) >= min_mean_counts):
                flags.append("I2: too little signal (mean counts per point below min_mean_counts, default 100)")
        ok_row = not any(("no area" in f) or ("not significant" in f) or ("too little signal" in f) for f in flags)
        ratio, ratio_err, usable = _ratio_and_error(A1, eA1, A2, eA2, cov12, ok_row)
        if not usable:
            flags.append("ratio undefined")
        out.append({"T": sp["T"], "I1": A1, "I1_err": eA1, "I2": A2, "I2_err": eA2, "ratio": ratio, "ratio_err": ratio_err,
                    "flags": flags, "w1": float(p[0]), "w2": float(p[1]), "eta": float(p[2]) if free_eta else ETA})
    return out


NEIGHBOUR_SHARE_FLAG = 0.05  # S-Figuryte-2: flag when the neighbour band contributes more than this fraction of a window's band signal


def neighbour_share(spectra, w1, w2, counts=True, free_eta=False, threshold=NEIGHBOUR_SHARE_FLAG):
    """S-Figuryte-2 diagnosis: per temperature, the fraction of each window's band signal that comes from the NEIGHBOUR band.

    Uses the fitted parameters of the joint fit (fit_pair): share1 = neighbour-band signal / (own + neighbour) inside window 1, and the
    same for window 2. A share above `threshold` (default 5 %) is flagged: the neighbour band dominates that window, so the ratio depends
    strongly on the neighbour model. Returns one dict per spectrum: T, share1, share2, flags. The share is a diagnosis, not a correction.
    """
    lo1, hi1 = _check_window(w1, "window 1")
    lo2, hi2 = _check_window(w2, "window 2")
    out = []
    for s in spectra:
        f = fit_pair(s["x"], s["y"], (lo1, hi1), (lo2, hi2), s.get("var"), counts, free_eta)
        x = np.asarray(s["x"], float)
        if not np.isfinite(f["A1"]):
            out.append({"T": float(s["T"]), "share1": float("nan"), "share2": float("nan"), "flags": ["fit failed: no share"]})
            continue
        eta = f["eta"]
        m1 = (x >= lo1) & (x <= hi1)
        m2 = (x >= lo2) & (x <= hi2)
        own1 = f["A1"] * pv_unit(x[m1], f["c1"], f["w1"], eta)
        nb1 = f["A2"] * pv_unit(x[m1], f["c2"], f["w2"], eta)
        own2 = f["A2"] * pv_unit(x[m2], f["c2"], f["w2"], eta)
        nb2 = f["A1"] * pv_unit(x[m2], f["c1"], f["w1"], eta)
        sh1 = float(nb1.sum() / (own1.sum() + nb1.sum())) if (own1.sum() + nb1.sum()) > 0 else float("nan")
        sh2 = float(nb2.sum() / (own2.sum() + nb2.sum())) if (own2.sum() + nb2.sum()) > 0 else float("nan")
        flags = []
        if np.isfinite(sh1) and sh1 > threshold:
            flags.append(f"I1: neighbour band gives {100 * sh1:.1f} % of the window signal")
        if np.isfinite(sh2) and sh2 > threshold:
            flags.append(f"I2: neighbour band gives {100 * sh2:.1f} % of the window signal")
        out.append({"T": float(s["T"]), "share1": sh1, "share2": sh2, "flags": flags})
    return out


def input_problems(spectra, w1, w2):
    """S-Figuryte-3: checks BEFORE the fit that the user can act on. Returns a list of messages (empty when the input is fine):
    a window with fewer than MIN_POINTS points, and a window whose highest point lies on its edge (the band then probably lies outside
    the window). The numbers are not changed; this only explains why a value may be missing or unreliable."""
    lo1, hi1 = _check_window(w1, "window 1")
    lo2, hi2 = _check_window(w2, "window 2")
    msgs = []
    for s_ in spectra:
        x = np.asarray(s_["x"], float)
        y = np.asarray(s_["y"], float)
        T = float(s_["T"])
        for k, (lo, hi) in enumerate(((lo1, hi1), (lo2, hi2)), start=1):
            m = (x >= lo) & (x <= hi) & np.isfinite(y)
            n = int(m.sum())
            if n < MIN_POINTS:
                msgs.append(f"{T:g} K: window {k} holds {n} points; at least {MIN_POINTS} are needed.")
                continue
            ys = y[m]
            i = int(np.argmax(ys))
            if i == 0 or i == n - 1:
                msgs.append(f"{T:g} K: the highest point of window {k} lies on its edge; the band may lie outside the window. Widen the window.")
    return msgs


def window_systematics(spectra, w1, w2, shift=0.10, counts=True):
    """Systematic error from the choice of the windows: each window is shifted by +-shift x its own width (the sign combinations
    of the two windows), the ratio is refitted, and the largest departure from the nominal ratio is the window systematic.

    Returns one dict per spectrum: T, ratio (nominal), ratio_err (statistical), ratio_sys (window systematic, NaN if no variant could be
    fitted), ratio_err_stat_sys = sqrt(ratio_err^2 + ratio_sys^2), n_variants (variants fitted), flags. A variant whose shifted windows
    overlap or are invalid is skipped and counted as not fitted. ratio_err_stat_sys is NOT the full uncertainty: the shape mismatch of
    real spectra (see the module docstring) is not in it."""
    lo1, hi1 = _check_window(w1, "window 1")
    lo2, hi2 = _check_window(w2, "window 2")
    nominal = window_ratio_series(spectra, (lo1, hi1), (lo2, hi2), counts)
    wd1, wd2 = hi1 - lo1, hi2 - lo2
    variants = []
    for s1 in (-shift, 0.0, shift):
        for s2 in (-shift, 0.0, shift):
            if s1 == 0.0 and s2 == 0.0:
                continue
            try:
                a = (lo1 + s1 * wd1, hi1 + s1 * wd1)
                b = (lo2 + s2 * wd2, hi2 + s2 * wd2)
                variants.append(window_ratio_series(spectra, a, b, counts))
            except ValueError:
                continue
    out = []
    for i, nom in enumerate(nominal):
        dev = []
        for v in variants:
            r = v[i]["ratio"]
            if np.isfinite(r) and np.isfinite(nom["ratio"]):
                dev.append(abs(r - nom["ratio"]))
        sys_ = max(dev) if dev else float("nan")
        tot = math.hypot(nom["ratio_err"], sys_) if np.isfinite(sys_) and np.isfinite(nom["ratio_err"]) else float("nan")
        flags = list(nom["flags"])
        if not dev:
            flags.append("no window variant could be fitted")
        out.append({"T": nom["T"], "ratio": nom["ratio"], "ratio_err": nom["ratio_err"], "ratio_sys": sys_,
                    "ratio_err_stat_sys": tot, "n_variants": len(dev), "flags": flags})
    return out


def crosscheck(window_rows, global_T, global_ratio, global_err, k=3.0):
    """Compare the window ratio with the global-fit ratio at the same temperatures.

    window_rows = output of window_ratio_series; global_T, global_ratio, global_err = arrays from the global fit (same temperatures).
    A temperature is compared only when both ratios and both errors are finite. The difference is judged against the combined error
    sqrt(window_err² + global_err²); |difference| > k * combined error marks a disagreement. The method disagreement is a systematic
    effect that the error bars do not contain (see the window-tail note in the module docstring), so a disagreement is a warning, not a
    correction. Returns {'rows': [...], 'n_compared', 'n_disagree', 'median_rel_diff'}."""
    gT = np.asarray(global_T, float)
    gr = np.asarray(global_ratio, float)
    ge = np.asarray(global_err, float)
    rows = []
    for r in window_rows:
        i = int(np.argmin(np.abs(gT - r["T"]))) if len(gT) else -1
        if i < 0 or abs(gT[i] - r["T"]) > 1e-6:
            rows.append({"T": r["T"], "diff": float("nan"), "combined_err": float("nan"), "disagree": False, "note": "no global-fit point at this T"})
            continue
        a, ea, b, eb = r["ratio"], r["ratio_err"], gr[i], ge[i]
        if not (np.isfinite(a) and np.isfinite(ea) and np.isfinite(b) and np.isfinite(eb)):
            rows.append({"T": r["T"], "diff": float("nan"), "combined_err": float("nan"), "disagree": False, "note": "not compared (undefined value)"})
            continue
        comb = math.sqrt(ea * ea + eb * eb)
        d = a - b
        rows.append({"T": r["T"], "diff": d, "combined_err": comb, "disagree": bool(abs(d) > k * comb),
                     "rel_diff": d / b if b else float("nan"), "note": ""})
    cmp_rows = [r for r in rows if np.isfinite(r["diff"])]
    rel = [abs(r["rel_diff"]) for r in cmp_rows if np.isfinite(r.get("rel_diff", float("nan")))]
    return {"rows": rows, "n_compared": len(cmp_rows), "n_disagree": sum(r["disagree"] for r in cmp_rows),
            "median_rel_diff": float(np.median(rel)) if rel else float("nan")}
