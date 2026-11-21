"""v2.19.1: two electrochemistry helpers.

1. Randles-Sevcik (Randles 1948, Sevcik 1948): cyclic-voltammetry peak current of a reversible, diffusion-controlled couple
   ip = 0.4463 n F A C sqrt(n F v D / (R T)), at 298.15 K unless another temperature is given.
   n: electrons, F: Faraday constant (C/mol), A: electrode area (cm²), C: concentration (mol/cm³), v: scan rate (V/s),
   D: diffusion coefficient (cm²/s). The inverse gives D from a measured peak current.
2. Battery cycling: specific capacity Q = I t / m (mAh/g) and coulombic efficiency = discharge capacity / charge capacity.

General practice: Randles-Sevcik assumes a reversible, planar, diffusion-controlled process. Check that ip grows with the square
root of the scan rate before you use it. Capacities depend on the active mass; state which mass was used."""
from __future__ import annotations

import math

F = 96485.33212      # C / mol
R_GAS = 8.314462618  # J / (mol K)


def randles_sevcik_peak_current(n: int, area_cm2: float, conc_mol_cm3: float, scan_v_s: float, D_cm2_s: float,
                                temp_K: float = 298.15) -> float:
    """Peak current (A) of a reversible, diffusion-controlled couple (Randles-Sevcik)."""
    root = math.sqrt(n * F * scan_v_s * D_cm2_s / (R_GAS * temp_K))
    return 0.4463 * n * F * area_cm2 * conc_mol_cm3 * root


def diffusion_from_peak_current(ip_A: float, n: int, area_cm2: float, conc_mol_cm3: float, scan_v_s: float,
                                temp_K: float = 298.15) -> float:
    """Diffusion coefficient (cm²/s) from a measured peak current (A), inverting Randles-Sevcik."""
    per_sqrtD = 0.4463 * n * F * area_cm2 * conc_mol_cm3
    root = ip_A / per_sqrtD                          # sqrt(n F v D / RT)
    return root * root * R_GAS * temp_K / (n * F * scan_v_s)


def specific_capacity_mAh_g(current_mA: float, time_h: float, mass_g: float) -> float:
    """Capacity (mAh/g) of a constant current (mA) applied for time_h (hours) to an active mass (g)."""
    if mass_g <= 0:
        raise ValueError("active mass must be positive")
    return float(current_mA) * float(time_h) / float(mass_g)


def coulombic_efficiency_percent(discharge_mAh_g: float, charge_mAh_g: float) -> float:
    """Coulombic efficiency (%) = discharge capacity / charge capacity x 100."""
    if charge_mAh_g <= 0:
        raise ValueError("charge capacity must be positive")
    return 100.0 * float(discharge_mAh_g) / float(charge_mAh_g)
