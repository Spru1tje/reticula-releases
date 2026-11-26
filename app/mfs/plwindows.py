"""PL thermometry on user-chosen windows: a fit that uses ONLY the data inside the two windows and reports the error of the ratio.

Why: the global fit (mfs/plratio.fit_global) starts from the whole spectrum (band search, shared start). When the user has already chosen
the two windows, this module fits each window on its own: inside one window the data are described by one unit-area pseudo-Voigt band plus a
straight background line. The band AREA of window 1 (I1) and of window 2 (I2) give the ratio R = I1 / I2 per temperature.

Errors: the covariance of each window fit comes from the Jacobian of the least-squares fit (Poisson weights sigma = sqrt(max(y, 1)) when
counts are given, otherwise the noise of the data). If the fit is worse than Poisson (reduced chi-square > 1) the covariance is scaled up by
it. The two windows are disjoint and the ratio error is propagated as R * sqrt((σI1/I1)² + (σI2/I2)²); the two areas are treated as
independent (they come from different x-ranges of the same spectrum, so this is an approximation, stated in the report).

Limits, all reported as flags per temperature:
  * 'tail outside window': the band centre ± 3 FWHM is not inside the window, so the area is an extrapolation of the band tail.
  * 'poor fit': reduced chi-square above 2 (the band + line do not describe the window).
  * 'edge': the fitted centre lies at the window edge (the band sits on the border).
  * 'few points': fewer than 12 data points in a window.
  * 'no area': the fitted area is not positive, the ratio is undefined.
Not checked here: the energy axis, the background outside the windows, the temperature calibration. A flag does not make a result unusable;
it says where to look.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import least_squares

from mfs.plratio import pv_unit

ETA = 0.5          # Lorentzian share of the pseudo-Voigt; fixed, the windows are too short to fit it reliably
MIN_POINTS = 12
SIGNIFICANCE = 3.0  # a band area must exceed this many of its own errors to give a ratio


def _check_window(w, name):
    if w is None or len(w) != 2:
        raise ValueError(f"{name}: give the window as (low, high) in nm")
    lo, hi = float(min(w)), float(max(w))
    if not (math.isfinite(lo) and math.isfinite(hi)) or hi <= lo:
        raise ValueError(f"{name}: the window {w} is empty or not finite")
    return lo, hi


def fit_window(x, y, window, var=None, counts=True):
    """Fit one window of one spectrum. Returns a dict: area, area_err, centre, width, chi2_red, n, flags (list of str).
    Only the points with lo <= x <= hi are used."""
    lo, hi = _check_window(window, "window")
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = (x >= lo) & (x <= hi) & np.isfinite(y)
    xs, ys = x[m], y[m]
    flags = []
    if len(xs) < MIN_POINTS:
        return {"area": float("nan"), "area_err": float("nan"), "centre": float("nan"), "width": float("nan"),
                "chi2_red": float("nan"), "n": int(len(xs)), "flags": ["few points"]}
    step = float(np.median(np.diff(xs)))
    if var is not None:
        v = np.asarray(var, float)[m]
        sig = np.sqrt(np.maximum(v, 1.0))
    elif counts:
        sig = np.sqrt(np.maximum(ys, 1.0))
    else:
        sig = np.full_like(ys, max(float(np.std(np.diff(ys))) / math.sqrt(2.0), 1e-12))
    # start: straight line through the window ends; band at the largest residual
    bg0 = np.interp(xs, [xs[0], xs[-1]], [ys[0], ys[-1]])
    res = ys - bg0
    i0 = int(np.argmax(res))
    c0 = float(xs[i0])
    area0 = max(float(np.sum(res) * step), 1e-9)
    w0 = max((hi - lo) / 6.0, 3.0 * step)
    xm = 0.5 * (lo + hi)
    span = 0.5 * (hi - lo)

    def model(p):
        A, c, w, b0, b1 = p
        return A * pv_unit(xs, c, w, ETA) + b0 + b1 * (xs - xm) / span

    def resid(p):
        return (model(p) - ys) / sig

    p0 = [area0, c0, w0, float(bg0[0]), 0.0]
    lb = [0.0, lo, 1.5 * step, -np.inf, -np.inf]
    ub = [np.inf, hi, hi - lo, np.inf, np.inf]
    try:
        r = least_squares(resid, p0, bounds=(lb, ub), x_scale="jac", max_nfev=2000)
    except Exception as exc:                                 # noqa: BLE001 - reported, not hidden
        return {"area": float("nan"), "area_err": float("nan"), "centre": float("nan"), "width": float("nan"),
                "chi2_red": float("nan"), "n": int(len(xs)), "flags": [f"fit failed: {exc}"]}
    p = r.x
    dof = max(len(xs) - len(p), 1)
    chi2_red = float(np.sum(r.fun ** 2) / dof)
    try:
        cov = np.linalg.inv(r.jac.T @ r.jac)
    except np.linalg.LinAlgError:
        cov = np.full((len(p), len(p)), np.nan)
    scale = max(chi2_red, 1.0)
    area_err = float(math.sqrt(abs(cov[0, 0]) * scale)) if np.isfinite(cov[0, 0]) else float("nan")
    A, c, w = float(p[0]), float(p[1]), float(p[2])
    if not (A > 0):
        flags.append("no area")
    elif not (A > SIGNIFICANCE * area_err):
        flags.append("not significant")          # area below 3 sigma: the band is not seen above the noise and the ratio is not meaningful
    if (c - lo) < 1e-6 * (hi - lo) or (hi - c) < 1e-6 * (hi - lo):
        flags.append("edge")
    if c - 3 * w < lo or c + 3 * w > hi:
        flags.append("tail outside window")
    if chi2_red > 2.0:
        flags.append("poor fit")
    return {"area": A, "area_err": area_err, "centre": c, "width": w, "chi2_red": chi2_red, "n": int(len(xs)), "flags": flags}


def window_ratio_series(spectra, w1, w2, counts=True):
    """Ratio I1/I2 per spectrum from the two windows, each fitted on its own data only.

    spectra = [{'T': K, 'x': nm, 'y': counts, optional 'var': array}, ...]. Returns a list of dicts (one per spectrum, same order):
    T, I1, I1_err, I2, I2_err, ratio, ratio_err, flags. The ratio is nan when an area is not positive.
    """
    lo1, hi1 = _check_window(w1, "window 1")
    lo2, hi2 = _check_window(w2, "window 2")
    if hi1 > lo2 and hi2 > lo1:
        raise ValueError("the two windows overlap; a band must not be counted in both")
    out = []
    for s in spectra:
        a = fit_window(s["x"], s["y"], (lo1, hi1), s.get("var"), counts)
        b = fit_window(s["x"], s["y"], (lo2, hi2), s.get("var"), counts)
        flags = [f"I1: {f}" for f in a["flags"]] + [f"I2: {f}" for f in b["flags"]]
        ratio = ratio_err = float("nan")
        ok = all("no area" not in f and "not significant" not in f and "few points" not in f and "fit failed" not in f for f in a["flags"] + b["flags"])
        if ok and a["area"] > 0 and b["area"] > 0 and np.isfinite(a["area_err"]) and np.isfinite(b["area_err"]):
            ratio = a["area"] / b["area"]
            ratio_err = ratio * math.sqrt((a["area_err"] / a["area"]) ** 2 + (b["area_err"] / b["area"]) ** 2)
        else:
            flags.append("ratio undefined")
        out.append({"T": float(s["T"]), "I1": a["area"], "I1_err": a["area_err"], "I2": b["area"], "I2_err": b["area_err"],
                    "ratio": ratio, "ratio_err": ratio_err, "flags": flags})
    return out
