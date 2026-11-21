"""Non-isothermal kinetics from heating-rate series (TGA or DSC runs at several heating rates beta).

- Kissinger (1957): ln(beta / Tp^2) = ln(A R / Ea) - Ea / (R Tp). Tp = peak temperature (K).
- Ozawa-Flynn-Wall (1965/1966): log10(beta) = const - 0.4567 Ea / (R T_alpha), at a fixed conversion alpha.

Both give the apparent activation energy Ea (J/mol, returned in kJ/mol). They assume a single, simple reaction; for multi-step
processes use the conversion-dependent (isoconversional) form per alpha. General practice: check that the fit is linear (R2 near 1)."""
from __future__ import annotations

import numpy as np

R = 8.314462618  # J / (mol K)


def _linfit(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    fit = A @ coef
    ss_res = float(np.sum((y - fit) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    dof = max(len(x) - 2, 1)
    s2 = ss_res / dof
    cov = np.linalg.inv(A.T @ A) * s2
    return coef[0], coef[1], float(np.sqrt(cov[0, 0])), r2


def _kelvin_warning(T):
    """v2.18.12: a temperature below 273.15 cannot be a kelvin peak of a thermal run (it would be below 0 °C), so it most likely
    was entered in °C. Peaks above about 300 °C cannot be told apart from kelvin by size; check the units of the table yourself."""
    import numpy as _np
    if len(T) and float(_np.min(T)) < 273.15:
        return ("a temperature is below 273.15 K (0 °C): the values may be °C, not kelvin; enter them in kelvin (K = °C + 273.15)")
    return None

def kissinger(betas, tp_kelvin):
    """Activation energy from peak temperatures at several heating rates (beta in K/min or K/s; the unit cancels in the slope)."""
    b = np.asarray(betas, float)
    T = np.asarray(tp_kelvin, float)
    ok = np.isfinite(b) & np.isfinite(T) & (b > 0) & (T > 0)
    b, T = b[ok], T[ok]
    if len(b) < 3:
        return {"status": "needs at least 3 heating rates"}
    x = 1.0 / T
    y = np.log(b / T ** 2)
    slope, intercept, slope_err, r2 = _linfit(x, y)
    ea = -slope * R
    out = {"status": "ok", "Ea_kJ": ea / 1000.0, "Ea_err_kJ": abs(slope_err * R) / 1000.0,
           "A_per_s": float(np.exp(intercept) * ea / R) if ea > 0 else float("nan"), "R2": r2, "n_points": int(len(b))}
    w = _kelvin_warning(T)
    if w:
        out["warning"] = w
    return out


def ozawa_flynn_wall(betas, t_alpha_kelvin):
    """Activation energy at ONE conversion alpha from the temperatures T_alpha reached at each heating rate."""
    b = np.asarray(betas, float)
    T = np.asarray(t_alpha_kelvin, float)
    ok = np.isfinite(b) & np.isfinite(T) & (b > 0) & (T > 0)
    b, T = b[ok], T[ok]
    if len(b) < 3:
        return {"status": "needs at least 3 heating rates"}
    slope, intercept, slope_err, r2 = _linfit(1.0 / T, np.log10(b))
    ea = -slope * R / 0.4567
    out = {"status": "ok", "Ea_kJ": ea / 1000.0, "Ea_err_kJ": abs(slope_err * R / 0.4567) / 1000.0, "R2": r2, "n_points": int(len(b))}
    w = _kelvin_warning(T)
    if w:
        out["warning"] = w
    return out
