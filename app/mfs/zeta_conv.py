"""v2.18.11: electrophoretic mobility -> zeta potential, with the validity check that the Smoluchowski form needs.

Henry (1931): mobility mu = 2 eps zeta f(ka) / (3 eta), with f(ka) between 1 (Hueckel limit, small particles or thick double layer,
small ka) and 1.5 (Smoluchowski limit, large ka). The approximation of Ohshima (1994, J. Colloid Interface Sci. 168, 269) is used for f.
Inverting: zeta = 3 eta mu / (2 eps f(ka)). Smoluchowski (f = 1.5) gives zeta = eta mu / eps.

Units: mobility in micrometre·cm / (V·s) as the Malvern Zetasizer reports it (1 µm·cm/(V·s) = 1e-8 m²/(V·s)); zeta returned in mV.
ka is the Debye parameter times the particle radius; the Debye length follows from the ionic strength of a 1:1 electrolyte (mM).
General practice: Smoluchowski is accurate for ka above about 300 (f within about 1 % of 1.5); for ka below about 1 the Hueckel form applies.
The ionic strength and the radius are not read from the file; the user enters them."""
from __future__ import annotations

import math

E_CHARGE = 1.602176634e-19      # C
K_B = 1.380649e-23              # J/K
N_A = 6.02214076e23             # 1/mol
EPS0 = 8.8541878128e-12         # F/m
MOBILITY_UNIT = 1e-8            # m²/(V·s) per µm·cm/(V·s)


def water_props(temp_C: float = 25.0):
    """Viscosity (Pa s) and relative permittivity of water at a temperature (°C), from standard tables (Vogel-type fit for the
    viscosity, linear fit for the permittivity between 0 and 50 °C; accurate to about 1 % in that range)."""
    T = float(temp_C) + 273.15
    eta = 2.414e-5 * 10 ** (247.8 / (T - 140.0))
    eps_r = 87.74 - 0.40008 * float(temp_C) + 9.398e-4 * float(temp_C) ** 2 - 1.41e-6 * float(temp_C) ** 3
    return eta, eps_r


def debye_length_m(ionic_strength_mM: float, temp_C: float = 25.0) -> float:
    """Debye length (m) of a 1:1 electrolyte: kappa² = 2 N_A e² I / (eps0 eps_r k_B T), I in mol/m³ (1 mM = 1 mol/m³)."""
    if ionic_strength_mM <= 0:
        return float("inf")
    _, eps_r = water_props(temp_C)
    T = float(temp_C) + 273.15
    kappa2 = 2 * N_A * E_CHARGE ** 2 * (ionic_strength_mM) / (EPS0 * eps_r * K_B * T)
    return 1.0 / math.sqrt(kappa2)


def henry_f(ka: float) -> float:
    """Henry's function f(ka) (Ohshima 1994 approximation): 1 for ka -> 0, 1.5 for ka -> infinity."""
    if ka <= 0:
        return 1.0
    return 1.0 + 1.0 / (2.0 * (1.0 + 2.5 / (ka * (1.0 + 2.0 * math.exp(-ka)))) ** 3)


def zeta_from_mobility(mobility_um_cm_Vs: float, ionic_strength_mM: float, radius_nm: float | None = None,
                       temp_C: float = 25.0, model: str = "auto") -> dict:
    """Zeta potential (mV) from the mobility.

    model: 'smoluchowski' (f = 1.5), 'henry' (f from ka) or 'auto' (Henry with ka; the Smoluchowski value is returned for comparison).
    radius_nm is needed for ka: if it is not given, only the Smoluchowski value and the Debye length are returned.
    Returns a dict with zeta_mV, the model used, ka, f, the Smoluchowski value, the difference, and a validity note."""
    if mobility_um_cm_Vs is None or not math.isfinite(float(mobility_um_cm_Vs)):
        return {"status": "no mobility"}
    eta, eps_r = water_props(temp_C)
    eps = eps_r * EPS0
    mu = float(mobility_um_cm_Vs) * MOBILITY_UNIT
    zeta_smol = eta * mu / eps * 1e3
    lam_d = debye_length_m(ionic_strength_mM, temp_C)
    ka = None if (radius_nm is None or not math.isfinite(lam_d)) else float(radius_nm) * 1e-9 / lam_d
    if model == "smoluchowski" or ka is None:
        f = 1.5 if (model != "henry" or ka is None) else henry_f(ka)
        note = ("Smoluchowski form (f = 1.5); give the radius to check its validity" if model != "henry" and ka is None
                else "Smoluchowski form (f = 1.5)")
        zeta = zeta_smol if f == 1.5 else 3 * eta * mu / (2 * eps * f) * 1e3
        used = "smoluchowski" if f == 1.5 else "henry"
        return {"status": "ok", "zeta_mV": zeta, "model": used, "ka": ka, "f": f, "zeta_smoluchowski_mV": zeta_smol,
                "debye_length_nm": lam_d * 1e9 if math.isfinite(lam_d) else None, "note": note}
    f = henry_f(ka)
    zeta = 3 * eta * mu / (2 * eps * f) * 1e3
    if ka >= 300:
        note = "ka >= 300: Smoluchowski valid (f within about 1 % of 1.5)"
    elif ka >= 100:
        note = "ka between 100 and 300: Smoluchowski is about 2 % off; Henry form used"
    elif ka >= 1:
        note = "ka between 1 and 100: intermediate regime; Henry form used"
    else:
        note = "ka < 1: Hueckel regime (small particles); Henry form used"
    if model == "henry" or model == "auto":
        used = "henry"
    return {"status": "ok", "zeta_mV": zeta, "model": used, "ka": ka, "f": f, "zeta_smoluchowski_mV": zeta_smol,
            "debye_length_nm": lam_d * 1e9, "note": note}
