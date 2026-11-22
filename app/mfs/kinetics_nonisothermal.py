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


def kas(betas, t_alpha_kelvin):
    """v2.19.1: Kissinger-Akahira-Sunose (Kissinger 1957; Akahira and Sunose 1971) isoconversional activation energy at ONE conversion.

    ln(beta / T_alpha²) = const - Ea / (R T_alpha)  ->  slope = -Ea / R of ln(beta / T²) against 1 / T.
    Unlike Ozawa-Flynn-Wall it uses no fixed constant (0.4567), so it is the second check next to OFW at the same conversion.
    General practice: the same reliability rules as for Kissinger apply (linear fit, R2 near 1, temperatures in kelvin)."""
    b = np.asarray(betas, float)
    T = np.asarray(t_alpha_kelvin, float)
    ok = np.isfinite(b) & np.isfinite(T) & (b > 0) & (T > 0)
    b, T = b[ok], T[ok]
    if len(b) < 3:
        return {"status": "needs at least 3 heating rates"}
    slope, intercept, slope_err, r2 = _linfit(1.0 / T, np.log(b / T ** 2))
    ea = -slope * R
    out = {"status": "ok", "Ea_kJ": ea / 1000.0, "Ea_err_kJ": abs(slope_err * R) / 1000.0, "R2": r2, "n_points": int(len(b))}
    w = _kelvin_warning(T)
    if w:
        out["warning"] = w
    return out


# ---- v2.19.2: differential (Friedman) method, beta-spread check and Kissinger-Friedman cross-check

def friedman(curves, alpha_target):
    """Friedman (1964) differential isoconversional method at one conversion alpha.

    curves: list of (beta, T, alpha) with T in kelvin and alpha in 0..1 along the run, one tuple per heating rate beta.
    dalpha/dt = beta · dalpha/dT (dalpha/dT from the curve). Then ln(dalpha/dt) = const - Ea / (R T_alpha): the slope of
    ln(dalpha/dt) against 1/T_alpha gives Ea. No model is assumed; it is sensitive to noise, so smooth the curve first."""
    xs, ys, used = [], [], []
    for beta, T, a in curves:
        T = np.asarray(T, float)
        a = np.asarray(a, float)
        ok = np.isfinite(T) & np.isfinite(a)
        T, a = T[ok], a[ok]
        order = np.argsort(a)
        T, a = T[order], a[order]
        if len(a) < 5 or not (a[0] <= alpha_target <= a[-1]):
            continue
        # keep only the monotone part of the curve (alpha rising with T)
        keep = np.r_[True, np.diff(a) > 0]
        T, a = T[keep], a[keep]
        if len(a) < 5 or not (a[0] <= alpha_target <= a[-1]):
            continue
        dadT = np.gradient(a, T)
        T_at = float(np.interp(alpha_target, a, T))
        dadT_at = float(np.interp(alpha_target, a, dadT))
        rate = abs(float(beta)) * dadT_at            # dalpha/dt in 1/s (beta in K/s) or 1/min (beta in K/min)
        if rate > 0 and np.isfinite(rate):
            xs.append(1.0 / T_at)
            ys.append(np.log(rate))
            used.append(float(beta))
    if len(xs) < 3:
        return {"status": "needs at least 3 heating rates with the conversion inside each curve"}
    slope, intercept, slope_err, r2 = _linfit(xs, ys)
    return {"status": "ok", "Ea_kJ": -slope * R / 1000.0, "Ea_err_kJ": abs(slope_err * R) / 1000.0, "R2": r2,
            "alpha": float(alpha_target), "n_points": len(xs), "betas_used": used,
            "note": "differential method (no model assumed); noise in dalpha/dT strongly affects Ea"}


def beta_spread_warning(betas, peak_temperatures=None) -> dict:
    """Checks a heating-rate series before Kissinger or KAS: the rates must spread (max/min >= 2 and neighbours >= 10 % apart),
    and when peak temperatures are given they must rise with the rate. Returns {'ok': bool, 'warnings': [...]}."""
    b = np.asarray(betas, float)
    b = b[np.isfinite(b) & (b > 0)]
    warnings = []
    if len(b) < 3:
        warnings.append("fewer than 3 heating rates: no activation energy can be fitted")
    else:
        s = np.sort(b)
        if s[-1] / s[0] < 2.0:
            warnings.append(f"heating rates span only {s[-1] / s[0]:.2f}x (need at least 2x for a reliable slope)")
        close = np.sum(s[1:] / s[:-1] < 1.1)
        if close:
            warnings.append(f"{int(close)} pair(s) of heating rates are within 10 % of each other")
    if peak_temperatures is not None and len(b) >= 2:
        T = np.asarray(peak_temperatures, float)
        bb = np.asarray(betas, float)
        ok = np.isfinite(T) & np.isfinite(bb) & (bb > 0)
        order = np.argsort(bb[ok])
        Ts = T[ok][order]
        if len(Ts) >= 2 and np.any(np.diff(Ts) <= 0):
            warnings.append("peak temperature does not rise with the heating rate (check the rates and the peak picks)")
    return {"ok": not warnings, "warnings": warnings}


def kissinger_friedman_check(ea_kissinger_kJ, ea_friedman_kJ, tolerance=0.15):
    """Agreement between two activation energies (kJ/mol): relative difference |a - b| / mean(a, b) must not exceed the tolerance."""
    a, b = float(ea_kissinger_kJ), float(ea_friedman_kJ)
    mean = (a + b) / 2.0
    rel = abs(a - b) / abs(mean) if mean else float("inf")
    return {"agree": bool(rel <= tolerance), "relative_difference": rel, "tolerance": tolerance,
            "note": "the two methods disagree: check the peak picks (Kissinger) and the smoothing (Friedman)" if rel > tolerance else "ok"}
