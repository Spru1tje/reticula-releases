"""v2.18.11: Avrami (JMAK) crystallisation kinetics, X(t) = 1 - exp(-k t^n) (Avrami 1939, 1940; Johnson and Mehl 1939).

X is the transformed fraction (0 to 1), t the time, k the rate constant (units of t^-n) and n the Avrami exponent. The fit uses the
linear form ln(-ln(1 - X)) = ln k + n ln t, over the fraction window 0.05 to 0.95 where the form is least affected by noise.
General practice: n is only a guide to the mechanism (nucleation and growth dimension); it is not proof of it. Check the fit visually."""
from __future__ import annotations

import math

import numpy as np


def avrami_fit(t, X, lo: float = 0.05, hi: float = 0.95) -> dict:
    t = np.asarray(t, float)
    X = np.asarray(X, float)
    ok = np.isfinite(t) & np.isfinite(X) & (t > 0) & (X >= lo) & (X <= hi)
    t, X = t[ok], X[ok]
    if len(t) < 3:
        return {"status": "needs at least 3 points between the fraction limits"}
    y = np.log(-np.log(1 - X))
    x = np.log(t)
    n, lnk = np.polyfit(x, y, 1)
    yhat = lnk + n * x
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    k = math.exp(lnk)
    t_half = (math.log(2) / k) ** (1.0 / n) if n > 0 else float("nan")
    # v2.18.12: standard errors of ln k and n from the least-squares covariance (two parameters, n - 2 degrees of freedom)
    A = np.vstack([np.ones_like(x), x]).T
    dof = len(x) - 2
    if dof > 0:
        cov = np.linalg.inv(A.T @ A) * (ss_res / dof)
        lnk_se, n_se = float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[1, 1]))
    else:
        lnk_se = n_se = float("nan")
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    out = {"status": "ok", "k": float(k), "n": float(n), "t_half": float(t_half), "R2": r2,
           "n_se": n_se, "k_se_ln": lnk_se, "n_points": int(len(t)), "window": (lo, hi)}
    # v2.18.12: with few points or a poor line the exponent is not reliable; say so next to the result
    if len(t) < 5:
        out["warning"] = f"only {len(t)} points in the fraction window: n is not reliable below 5 points"
    elif r2 < 0.98:
        out["warning"] = f"R² {r2:.3f} is below 0.98: the ln(−ln(1 − X)) plot is not linear, so the Avrami form may not describe this process"
    return out


def avrami_curve(t, k, n):
    """X(t) = 1 - exp(-k t^n)."""
    t = np.asarray(t, float)
    return 1.0 - np.exp(-k * np.power(t, n))
