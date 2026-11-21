"""Avrami (JMAK) crystallisation kinetics (v2.18.11): X(t) = 1 - exp(-k t^n).
Linearised: ln(-ln(1 - X)) = ln k + n ln t. n is the Avrami exponent (dimension and mechanism), k the rate constant.
Half-time t_half = (ln 2 / k)^(1/n). Valid for the early and middle part of the transformation (X roughly 0.05 to 0.95);
General practice: check the linear fit (R2) and drop the points where X is near 0 or 1 (noise dominates there)."""
from __future__ import annotations

import numpy as np


def fit_avrami(t, x, xmin=0.05, xmax=0.95):
    t = np.asarray(t, float)
    x = np.asarray(x, float)
    ok = np.isfinite(t) & np.isfinite(x) & (t > 0) & (x > xmin) & (x < xmax)
    t, x = t[ok], x[ok]
    if len(t) < 3:
        return {"status": "needs at least 3 points between xmin and xmax"}
    X = np.log(t)
    Y = np.log(-np.log(1.0 - x))
    A = np.vstack([X, np.ones_like(X)]).T
    (n, lnk), *_ = np.linalg.lstsq(A, Y, rcond=None)
    fit = A @ np.array([n, lnk])
    ss_tot = float(np.sum((Y - Y.mean()) ** 2))
    r2 = 1 - float(np.sum((Y - fit) ** 2)) / ss_tot if ss_tot > 0 else float("nan")
    k = float(np.exp(lnk))
    return {"status": "ok", "n": float(n), "k": k, "t_half": float((np.log(2) / k) ** (1 / n)), "R2": r2, "n_points": int(len(t))}
