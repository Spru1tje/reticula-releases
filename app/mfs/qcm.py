"""v2.19.1: quartz crystal microbalance (QCM) mass from the frequency shift, Sauerbrey (1959) relation.

For a thin, rigid, uniform film on an AT-cut quartz crystal of fundamental frequency f0:
    delta_f = - 2 f0^2 delta_m / (A sqrt(mu_q rho_q)),  mu_q = 2.947e11 g/(cm s^2),  rho_q = 2.648 g/cm^3.
For 5 MHz quartz the mass sensitivity C_f = sqrt(mu_q rho_q) / (2 f0^2), the mass per area per hertz of shift, is about
17.7 ng/(Hz cm^2) (Sauerbrey constant).
General practice: the relation holds only for thin, rigid films with a mass well below the crystal mass. Viscoelastic or
thick films (for example polymers in liquid) need a model such as Kanazawa-Gordon; the function warns when the shift is large."""
from __future__ import annotations

import math

MU_Q = 2.947e11      # g / (cm s^2), shear modulus of AT-cut quartz
RHO_Q = 2.648        # g / cm^3, density of quartz


def sensitivity_ng_per_hz_cm2(f0_hz: float) -> float:
    """Sauerbrey mass sensitivity C_f (ng per (Hz cm²)) for a crystal of fundamental frequency f0 (Hz)."""
    cf = math.sqrt(MU_Q * RHO_Q) / (2.0 * float(f0_hz) ** 2)   # g / (Hz cm^2): mass per area per Hz of shift
    return cf * 1e9


def mass_per_area_from_shift(delta_f_hz: float, f0_hz: float = 5e6) -> dict:
    """Mass per unit area (µg/cm²) from the frequency shift (Hz). A fall in frequency (negative shift) is a mass gain.

    Returns dict(status, mass_ug_per_cm2, sensitivity_ng_per_hz_cm2, warning). The warning is set when |delta_f| is larger than
    2 % of f0, where the thin-film relation starts to fail."""
    df = float(delta_f_hz)
    f0 = float(f0_hz)
    if not math.isfinite(df) or f0 <= 0:
        return {"status": "invalid input"}
    dm_per_area = -df * math.sqrt(MU_Q * RHO_Q) / (2.0 * f0 ** 2)   # g / cm^2
    out = {"status": "ok", "mass_ug_per_cm2": dm_per_area * 1e6, "sensitivity_ng_per_hz_cm2": sensitivity_ng_per_hz_cm2(f0)}
    if abs(df) > 0.02 * f0:
        out["warning"] = "shift larger than 2 % of f0: the thin, rigid film relation is not reliable"
    return out
