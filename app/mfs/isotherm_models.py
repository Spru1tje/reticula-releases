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
