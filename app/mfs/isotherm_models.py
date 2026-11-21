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
            "n": float(popt[2]), "n_err": float(err[2]), "R2": r2, "n_points": int(len(p)),
            "RSS": float(np.sum(res ** 2)), "k": 3}


def langmuir(p, qm, k):
    """Langmuir isotherm: q = qm K P / (1 + K P) (Langmuir 1918). Same as Sips with n = 1."""
    x = k * np.asarray(p, float)
    return qm * x / (1.0 + x)


def toth(p, qm, b, t):
    """Toth isotherm: q = qm b P / (1 + (b P)^t)^(1/t) (Toth 1971). t = 1 gives Langmuir; t < 1 a heterogeneous surface."""
    x = b * np.asarray(p, float)
    return qm * x / np.power(1.0 + np.power(x, t), 1.0 / t)


def _clean(p, q):
    p = np.asarray(p, float)
    q = np.asarray(q, float)
    ok = np.isfinite(p) & np.isfinite(q) & (p > 0) & (q >= 0)
    return p[ok], q[ok]


def _fit(func, p, q, p0, lb, ub, names):
    from scipy.optimize import curve_fit
    try:
        popt, pcov = curve_fit(func, p, q, p0=p0, bounds=(lb, ub), maxfev=20000)
    except (RuntimeError, ValueError) as exc:
        return {"status": f"fit failed: {exc}"}
    err = np.sqrt(np.clip(np.diag(pcov), 0, np.inf))
    res = q - func(p, *popt)
    rss = float(np.sum(res ** 2))
    ss_tot = float(np.sum((q - q.mean()) ** 2))
    out = {"status": "ok", "R2": 1 - rss / ss_tot if ss_tot > 0 else float("nan"), "RSS": rss, "n_points": int(len(p)),
           "k": len(popt), "fitted": popt}
    for nm, v, e in zip(names, popt, err):
        out[nm] = float(v)
        out[nm + "_err"] = float(e)
    return out


def fit_langmuir(p, q):
    p, q = _clean(p, q)
    if len(p) < 3:
        return {"status": "needs at least 3 points"}
    k0 = 1.0 / max(np.median(p), 1e-12)
    return _fit(langmuir, p, q, [max(q.max(), 1e-12) * 1.2, k0], [0, 0], [np.inf, np.inf], ["qm", "K"])


def fit_toth(p, q):
    p, q = _clean(p, q)
    if len(p) < 4:
        return {"status": "needs at least 4 points"}
    k0 = 1.0 / max(np.median(p), 1e-12)
    return _fit(toth, p, q, [max(q.max(), 1e-12) * 1.2, k0, 1.0], [0, 0, 0.05], [np.inf, np.inf, 3.0], ["qm", "b", "t"])


def aicc(rss: float, n: int, k: int) -> float:
    """Corrected Akaike information criterion for least squares (Hurvich and Tsai 1989): lower is better. Needs n - k - 1 > 0."""
    if n - k - 1 <= 0 or rss <= 0:
        return float("nan")
    return n * np.log(rss / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1)


def compare_isotherms(p, q):
    """Fits Langmuir (2 parameters), Sips (3) and Toth (3) and ranks them by AICc. The model with the lowest AICc is 'best'.
    Differences below 2 are not a meaningful preference (General practice). Needs at least 6 points for the three-parameter models."""
    p, q = _clean(p, q)
    if len(p) < 6:
        return {"status": "needs at least 6 points"}
    fits = {"Langmuir": fit_langmuir(p, q), "Sips": fit_sips(p, q), "Toth": fit_toth(p, q)}
    rows = []
    for name, r in fits.items():
        if r.get("status") != "ok":
            continue
        score = aicc(r["RSS"], r["n_points"], r["k"])
        rows.append({"model": name, "AICc": float(score), "R2": r["R2"], "RSS": r["RSS"], "k": r["k"]})
    rows = [r for r in rows if np.isfinite(r["AICc"])]
    if not rows:
        return {"status": "no model could be fitted"}
    rows.sort(key=lambda r: r["AICc"])
    best = rows[0]["model"]
    for r in rows:
        r["delta_AICc"] = r["AICc"] - rows[0]["AICc"]
    return {"status": "ok", "best": best, "table": rows, "fits": fits}


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
