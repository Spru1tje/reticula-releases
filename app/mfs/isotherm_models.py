"""Sips (Langmuir-Freundlich) isotherm: q = qm (K P)^n / (1 + (K P)^n). Sips (1948). n = 1 gives Langmuir; n < 1 heterogeneous surface.
Fitted with non-linear least squares; the uncertainty is from the covariance matrix. General practice: report n with its error and check
that qm is not extrapolated far beyond the highest measured pressure."""
from __future__ import annotations

import numpy as np


def sips(p, qm, k, n):
    x = (k * np.asarray(p, float)) ** n
    return qm * x / (1.0 + x)


def fit_sips(p, q):
    from scipy.optimize import curve_fit
    p = np.asarray(p, float)
    q = np.asarray(q, float)
    ok = np.isfinite(p) & np.isfinite(q) & (p > 0) & (q >= 0)
    p, q = p[ok], q[ok]
    if len(p) < 4:
        return {"status": "needs at least 4 points"}
    k0 = 1.0 / max(np.median(p), 1e-12)
    try:
        popt, pcov = curve_fit(sips, p, q, p0=[max(q.max(), 1e-12) * 1.2, k0, 1.0], bounds=([0, 0, 0.05], [np.inf, np.inf, 5.0]), maxfev=20000)
    except (RuntimeError, ValueError) as exc:
        return {"status": f"fit failed: {exc}"}
    err = np.sqrt(np.clip(np.diag(pcov), 0, np.inf))
    res = q - sips(p, *popt)
    ss_tot = float(np.sum((q - q.mean()) ** 2))
    r2 = 1 - float(np.sum(res ** 2)) / ss_tot if ss_tot > 0 else float("nan")
    return {"status": "ok", "qm": float(popt[0]), "qm_err": float(err[0]), "K": float(popt[1]), "K_err": float(err[1]),
            "n": float(popt[2]), "n_err": float(err[2]), "R2": r2, "n_points": int(len(p))}


def toth(p, qm, b, t):
    """Toth isotherm (Toth 1971): q = qm P / (b + P^t)^(1/t). t < 1 is heterogeneity; t = 1 is Langmuir (b = 1/K)."""
    p = np.asarray(p, float)
    return qm * p / (b + p ** t) ** (1.0 / t)


def fit_toth(p, q):
    from scipy.optimize import curve_fit
    p = np.asarray(p, float)
    q = np.asarray(q, float)
    ok = np.isfinite(p) & np.isfinite(q) & (p > 0) & (q >= 0)
    p, q = p[ok], q[ok]
    if len(p) < 4:
        return {"status": "needs at least 4 points"}
    try:
        popt, pcov = curve_fit(toth, p, q, p0=[max(q.max(), 1e-12) * 1.2, 1.0, 0.8], bounds=([0, 1e-12, 0.1], [np.inf, np.inf, 3.0]),
                               maxfev=20000)
    except (RuntimeError, ValueError) as exc:
        return {"status": f"fit failed: {exc}"}
    err = np.sqrt(np.clip(np.diag(pcov), 0, np.inf))
    return {"status": "ok", "qm": float(popt[0]), "qm_err": float(err[0]), "b": float(popt[1]), "t": float(popt[2]),
            "t_err": float(err[2]), "rss": float(np.sum((q - toth(p, *popt)) ** 2)), "n_points": int(len(p))}


def dubinin_radushkevich(p, q, p0=1.0, temperature=77.35):
    """Dubinin-Radushkevich (1947) micropore fit: ln q = ln W0 - (R T ln(p0/p))^2 / E^2. E in J/mol (mean adsorption energy).
    temperature 77.35 K is liquid nitrogen. Returns W0 (same unit as q) and E in kJ/mol. General practice: E between 8 and 16 kJ/mol
    is usually read as physisorption in micropores; outside that range the model is not the right description."""
    R = 8.314462618
    p = np.asarray(p, float)
    q = np.asarray(q, float)
    ok = np.isfinite(p) & np.isfinite(q) & (p > 0) & (p < p0) & (q > 0)
    p, q = p[ok], q[ok]
    if len(p) < 3:
        return {"status": "needs at least 3 points with p < p0"}
    x = (R * temperature * np.log(p0 / p)) ** 2
    y = np.log(q)
    A = np.vstack([x, np.ones_like(x)]).T
    (slope, intercept), *_ = np.linalg.lstsq(A, y, rcond=None)
    if slope >= 0:
        return {"status": "no micropore trend (slope not negative)"}
    E = 1.0 / np.sqrt(-slope)
    return {"status": "ok", "W0": float(np.exp(intercept)), "E_kJ": E / 1000.0, "n_points": int(len(p))}


def aic(q, qfit, k):
    """Akaike information criterion for least squares: n ln(RSS/n) + 2k (lower is better). k = number of fitted parameters."""
    q = np.asarray(q, float)
    qfit = np.asarray(qfit, float)
    n = len(q)
    rss = float(np.sum((q - qfit) ** 2))
    if n == 0 or rss <= 0:
        return float("-inf")
    return n * np.log(rss / n) + 2 * k
