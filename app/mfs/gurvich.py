"""v2.19.2: total pore volume from the Gurvich rule (Gurvich 1915; as used in BET-type porosity reports).

The nitrogen uptake at a relative pressure near 0.95 is taken as filling all pores. The gas volume at standard temperature and
pressure (STP, 0 °C, 1 atm) is converted to the volume of liquid nitrogen: V_liq = n × V_m,liquid, with n = V_STP / V_mol,gas.
Constants: V_mol,gas = 22414 cm³ STP per mol (ideal gas at STP); V_m,liquid of N2 at 77 K = 34.67 cm³/mol.

General practice: the Gurvich volume depends on the chosen relative pressure (0.95 is conventional) and on the pores being
filled by liquid-like nitrogen; compare it with the BET or density-functional values before you report a pore volume."""
from __future__ import annotations

V_MOL_GAS_STP = 22414.0      # cm³ STP per mol
V_M_LIQ_N2 = 34.67           # cm³ per mol, liquid nitrogen at 77 K


def gurvich_volume(v_stp_cm3_per_g: float) -> dict:
    """Liquid-nitrogen pore volume (cm³/g) from the STP gas volume adsorbed at the reference relative pressure (near 0.95)."""
    if v_stp_cm3_per_g is None or v_stp_cm3_per_g != v_stp_cm3_per_g or v_stp_cm3_per_g < 0:
        return {"status": "needs a non-negative uptake in cm³ STP/g"}
    n_mol_per_g = v_stp_cm3_per_g / V_MOL_GAS_STP
    return {"status": "ok", "n_mmol_per_g": n_mol_per_g * 1000.0, "V_pore_cm3_per_g": n_mol_per_g * V_M_LIQ_N2,
            "note": "Gurvich rule at the reference relative pressure (conventionally 0.95)"}
