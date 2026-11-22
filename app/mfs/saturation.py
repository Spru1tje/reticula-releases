"""v2.19.2: saturation (binding and enzyme-like) curves. Two standard fits, each returns a dict with the fitted values and R².

- Hill (1910): theta = x^n / (K^n + x^n). theta is the fractional occupancy, K the concentration at half occupancy, n the Hill
  coefficient (n > 1 cooperative, n = 1 no cooperativity). Hill plot: log10(theta / (1 - theta)) = n log10(x) - n log10(K).
- Michaelis-Menten (1913): v = Vmax S / (Km + S). Lineweaver-Burk (1934) linear form: 1/v = (Km/Vmax)(1/S) + 1/Vmax.
  The linear form is shown for the check; the reported values come from the non-linear fit, which weights the points correctly.

General practice: a Hill coefficient fitted from fewer than about 6 concentrations, or from data that do not reach both plateaus,
is not reliable. Check the residuals before you report n or Km."""
from __future__ import annotations

import numpy as np


def _r2(y, fit):
    ss_res = float(np.sum((y - fit) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def hill_fit(x, theta):
    """Fit theta(x) = x^n / (K^n + x^n) by least squares. x > 0; theta in 0..1 (values at exactly 0 or 1 are dropped)."""
    from scipy.optimize import curve_fit
    x = np.asarray(x, float)
    th = np.asarray(theta, float)
    ok = np.isfinite(x) & np.isfinite(th) & (x > 0) & (th > 0) & (th < 1)
    x, th = x[ok], th[ok]
    if len(x) < 4:
        return {"status": "needs at least 4 points with 0 < theta < 1"}
    # starting values: K at the median of x, n from the Hill plot slope (least squares)
    lx, ly = np.log10(x), np.log10(th / (1 - th))
    n0 = float(np.polyfit(lx, ly, 1)[0]) if np.ptp(lx) > 0 else 1.0
    K0 = float(np.median(x))

    def model(xx, K, n):
        xn = np.power(xx, n)
        return xn / (np.power(K, n) + xn)

    try:
        p, cov = curve_fit(model, x, th, p0=[K0, max(n0, 0.1)], bounds=([1e-300, 0.05], [np.inf, 20.0]), maxfev=20000)
    except (RuntimeError, ValueError) as exc:
        return {"status": f"fit failed: {exc}"}
    err = np.sqrt(np.diag(cov)) if np.all(np.isfinite(cov)) else np.full(2, np.nan)
    K, n = float(p[0]), float(p[1])
    return {"status": "ok", "K": K, "n": n, "K_err": float(err[0]), "n_err": float(err[1]), "R2": _r2(th, model(x, K, n)),
            "n_points": int(len(x)), "note": "n > 1 suggests positive cooperativity; n = 1 is the non-cooperative case"}


def michaelis_menten_fit(S, v):
    """Fit v = Vmax S / (Km + S) by non-linear least squares; also the Lineweaver-Burk linear estimate for comparison."""
    from scipy.optimize import curve_fit
    S = np.asarray(S, float)
    v = np.asarray(v, float)
    ok = np.isfinite(S) & np.isfinite(v) & (S > 0) & (v > 0)
    S, v = S[ok], v[ok]
    if len(S) < 3:
        return {"status": "needs at least 3 positive substrate concentrations"}

    def model(s, vmax, km):
        return vmax * s / (km + s)

    vmax0 = float(np.max(v)) * 1.1
    km0 = float(np.median(S))
    try:
        p, cov = curve_fit(model, S, v, p0=[vmax0, km0], bounds=([0, 0], [np.inf, np.inf]), maxfev=20000)
    except (RuntimeError, ValueError) as exc:
        return {"status": f"fit failed: {exc}"}
    err = np.sqrt(np.diag(cov)) if np.all(np.isfinite(cov)) else np.full(2, np.nan)
    vmax, km = float(p[0]), float(p[1])
    # Lineweaver-Burk linear estimate: 1/v = (Km/Vmax)(1/S) + 1/Vmax
    slope, intercept = np.polyfit(1.0 / S, 1.0 / v, 1)
    lb_vmax = 1.0 / intercept if intercept > 0 else float("nan")
    lb_km = slope * lb_vmax if np.isfinite(lb_vmax) else float("nan")
    return {"status": "ok", "Vmax": vmax, "Km": km, "Vmax_err": float(err[0]), "Km_err": float(err[1]),
            "R2": _r2(v, model(S, vmax, km)), "n_points": int(len(S)),
            "lineweaver_burk": {"Vmax": float(lb_vmax), "Km": float(lb_km)}}
