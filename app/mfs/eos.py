"""v2.19.1: third-order Birch-Murnaghan equation of state (Birch 1947, Murnaghan 1944) fitted to pressure-volume data.

P(V) = (3 B0 / 2) [(V0/V)^(7/3) - (V0/V)^(5/3)] {1 + (3/4) (B0' - 4) [(V0/V)^(2/3) - 1]}
with V0 the zero-pressure volume, B0 the bulk modulus and B0' its pressure derivative. Units: V in the input unit, P and B0 in the
same pressure unit (for example GPa with V in Å³ per formula unit).
General practice: the fit needs data on both sides of V0 (compression and, ideally, tension). Check the residuals and the
uncertainties; B0' is strongly correlated with B0 for narrow pressure ranges."""
from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit


def pressure_bm3(V, V0, B0, B0p):
    """Pressure from the third-order Birch-Murnaghan equation of state."""
    V = np.asarray(V, float)
    x = (V0 / V) ** (2.0 / 3.0)
    return 1.5 * B0 * (x ** 3.5 - x ** 2.5) * (1.0 + 0.75 * (B0p - 4.0) * (x - 1.0))


def fit_birch_murnaghan(V, P, p0=None) -> dict:
    """Fits V0, B0 and B0' to P(V) data. Returns dict(status, V0, B0, B0p, errors, rms_residual, n) or dict(status, reason)."""
    V = np.asarray(V, float)
    P = np.asarray(P, float)
    ok = np.isfinite(V) & np.isfinite(P) & (V > 0)
    V, P = V[ok], P[ok]
    if len(V) < 4:
        return {"status": "needs at least 4 points"}
    if p0 is None:
        p0 = [float(np.max(V)), max(float(np.ptp(P)) / max(float(np.ptp(V)) / float(np.mean(V)), 1e-9), 1e-3), 4.0]
    try:
        popt, pcov = curve_fit(pressure_bm3, V, P, p0=p0, maxfev=20000)
    except (RuntimeError, ValueError) as exc:
        return {"status": "fit failed", "reason": str(exc)}
    err = np.sqrt(np.diag(pcov)) if np.all(np.isfinite(pcov)) else np.full(3, np.nan)
    resid = P - pressure_bm3(V, *popt)
    return {"status": "ok", "V0": float(popt[0]), "B0": float(popt[1]), "B0p": float(popt[2]),
            "errors": {"V0": float(err[0]), "B0": float(err[1]), "B0p": float(err[2])},
            "rms_residual": float(np.sqrt(np.mean(resid ** 2))), "n": int(len(V))}
