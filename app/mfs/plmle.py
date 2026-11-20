"""Poisson maximum-likelihood thermometry on raw band counts (v2.18.6).

Model per spectrum k (temperature T_k), integrated counts of the two bands:
    N1_k ~ Poisson(B1 + S1_k),   S1_k = D_k * S2_k
    N2_k ~ Poisson(B2 + S2_k),   D_k = alpha * exp(-dE / (k_B T_k))
B1, B2 are the known background counts of the two windows (dark plus background, measured).
Free parameters: dE (eV), ln(alpha) and ln(S2_k) for every temperature (profiled out as nuisance).

Why this instead of subtract-then-divide: subtraction gives negative or zero band areas
at low signal, which cannot enter a logarithm, and it drops or biases those points.
The likelihood keeps every point, uses the Poisson variance of the raw counts, and takes
the uncertainty from the observed Fisher information (numerical Hessian of the negative
log-likelihood, inverted).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

KB_EV = 8.617333262e-5          # eV / K
EV_CM = 8065.544                # cm-1 / eV


def _mu(p, T, B1, B2):
    dE, lna, u = p[0], p[1], p[2:]
    S2 = np.exp(u)
    S1 = np.exp(lna) * np.exp(-dE / (KB_EV * T)) * S2
    return B1 + S1, B2 + S2


def nll(p, T, N1, N2, B1, B2):
    """Poisson negative log-likelihood (constant terms dropped)."""
    m1, m2 = _mu(p, T, B1, B2)
    return float(np.sum(m1 - N1 * np.log(m1)) + np.sum(m2 - N2 * np.log(m2)))


def _start(T, N1, N2, B1, B2):
    """Starting point from the data only (never from the truth): weighted ln-linear fit of the
    background-subtracted ratio where it is positive, else a neutral guess."""
    S1 = N1 - B1
    S2 = N2 - B2
    ok = (S1 > 0) & (S2 > 0)
    if ok.sum() >= 3:
        lnD = np.log(S1[ok] / S2[ok])
        slope, icpt = np.polyfit(1.0 / T[ok], lnD, 1)
        dE0 = -slope * KB_EV
        lna0 = icpt
    else:
        dE0, lna0 = 0.05, 0.0
    if not np.isfinite(dE0) or dE0 <= 0:
        dE0 = 0.05
    u0 = np.log(np.maximum(N2 - B2, 1.0))
    return np.concatenate([[dE0, lna0], u0])


def hessian(p, T, N1, N2, B1, B2, rel_step=1e-4):
    n = len(p)
    h = np.full(n, rel_step) * np.maximum(np.abs(p), 1.0)
    h[0] = 1e-6                 # eV
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            ei = np.zeros(n); ej = np.zeros(n); ei[i] = h[i]; ej[j] = h[j]
            v = (nll(p + ei + ej, T, N1, N2, B1, B2) - nll(p + ei - ej, T, N1, N2, B1, B2)
                 - nll(p - ei + ej, T, N1, N2, B1, B2) + nll(p - ei - ej, T, N1, N2, B1, B2)) / (4 * h[i] * h[j])
            H[i, j] = H[j, i] = v
    return H


def _grad(p, T, N1, N2, B1, B2, rel_step=1e-6):
    n = len(p)
    h = np.full(n, rel_step) * np.maximum(np.abs(p), 1.0)
    h[0] = 1e-7
    g = np.zeros(n)
    for i in range(n):
        e = np.zeros(n); e[i] = h[i]
        g[i] = (nll(p + e, T, N1, N2, B1, B2) - nll(p - e, T, N1, N2, B1, B2)) / (2 * h[i])
    return g


def fit(T, N1, N2, B1, B2):
    """Maximum-likelihood fit. Returns dict: dE_cm1, sigma_cm1 (from the observed information),
    alpha, converged, nll, and the covariance of (dE, ln alpha).

    Scaling matters: dE is ~1e-1 eV while the nuisance ln S2 is ~1e1 and the objective ~1e6, so the
    optimiser works on dE in units of 100 cm-1 and on the objective divided by the total counts,
    with tight tolerances, and a few Newton steps polish the optimum."""
    T = np.asarray(T, float); N1 = np.asarray(N1, float); N2 = np.asarray(N2, float)
    B1 = np.broadcast_to(np.asarray(B1, float), T.shape).copy()
    B2 = np.broadcast_to(np.asarray(B2, float), T.shape).copy()
    p0 = _start(T, N1, N2, B1, B2)
    scale = max(float(np.sum(N1) + np.sum(N2)), 1.0)
    to_q = lambda p: np.concatenate([[p[0] * EV_CM / 100.0], p[1:]])
    to_p = lambda q: np.concatenate([[q[0] * 100.0 / EV_CM], q[1:]])
    obj = lambda q: nll(to_p(q), T, N1, N2, B1, B2) / scale
    res = minimize(obj, to_q(p0), method="L-BFGS-B",
                   options={"maxiter": 5000, "ftol": 1e-15, "gtol": 1e-10})
    p = to_p(res.x)
    # Newton polish in the physical parameters (accept a step only if the likelihood does not increase)
    rel_step = np.inf
    for _ in range(6):
        g = _grad(p, T, N1, N2, B1, B2)
        Hn = hessian(p, T, N1, N2, B1, B2)
        try:
            step = np.linalg.solve(Hn, g)
        except np.linalg.LinAlgError:
            break
        if not np.all(np.isfinite(step)):
            break
        cand = p - step
        if nll(cand, T, N1, N2, B1, B2) <= nll(p, T, N1, N2, B1, B2):
            p = cand
            rel_step = float(np.max(np.abs(step[:2]) / np.array([max(abs(p[0]), 1e-3), max(abs(p[1]), 1.0)])))
            if rel_step < 1e-8:
                break
        else:
            break
    converged = bool(np.all(np.isfinite(p)) and (res.success or rel_step < 1e-4))
    H = hessian(p, T, N1, N2, B1, B2)
    try:
        cov = np.linalg.inv(H)
        var_dE = cov[0, 0]
        ok_cov = np.isfinite(var_dE) and var_dE > 0
    except np.linalg.LinAlgError:
        cov = np.full_like(H, np.nan); ok_cov = False
    return {
        "dE_cm1": float(p[0] * EV_CM),
        "sigma_cm1": float(np.sqrt(var_dE) * EV_CM) if ok_cov else float("nan"),
        "alpha": float(np.exp(p[1])),
        "converged": converged,
        "nll": float(res.fun),
        "cov_dE_lnalpha": cov[:2, :2] if np.all(np.isfinite(cov[:2, :2])) else None,
    }


def cramer_rao_sigma(T, S2, B1, B2, dE_cm1, alpha):
    """Expected sigma(dE) in cm-1 at the given true values, from the analytic Poisson Fisher information
    (all nuisance parameters profiled). Used as the reference for the Cramér-Rao consistency check."""
    T = np.asarray(T, float); S2 = np.asarray(S2, float)
    dE = dE_cm1 / EV_CM
    S1 = alpha * np.exp(-dE / (KB_EV * T)) * S2
    mu1 = B1 + S1; mu2 = B2 + S2
    K = len(T)
    I = np.zeros((2 + K, 2 + K))
    for k in range(K):
        g1 = np.zeros(2 + K); g1[0] = -S1[k] / (KB_EV * T[k]); g1[1] = S1[k]; g1[2 + k] = S1[k]
        I += np.outer(g1, g1) / mu1[k]
        g2 = np.zeros(2 + K); g2[2 + k] = S2[k]
        I += np.outer(g2, g2) / mu2[k]
    cov = np.linalg.inv(I)
    return float(np.sqrt(cov[0, 0]) * EV_CM)
