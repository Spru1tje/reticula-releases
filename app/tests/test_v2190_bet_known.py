"""v2.19.0: BET known-answer test. A noise-free isotherm is generated from the BET equation with a known monolayer
capacity n_m and constant C; the BET functions must return exactly those values and the matching surface area
(area = n_m x N_A x sigma(N2) / V_mol, with sigma = 0.162 nm^2 and V_mol = 22414 cm^3 STP per mol)."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from mfs import analysis as A  # noqa: E402

NM_TRUE = 12.5   # cm3 STP per g
C_TRUE = 80.0


def _bet_isotherm(p):
    """Adsorbed amount from the BET equation: q = n_m C x / ((1 - x)(1 - x + C x))."""
    return NM_TRUE * C_TRUE * p / ((1 - p) * (1 - p + C_TRUE * p))


def _area_from_nm(nm):
    return nm * A.CM3STP_TO_MOL * 6.02214076e23 * A.N2_CROSS


def test_bet_rouquerol_recovers_known_monolayer_and_area():
    p = np.linspace(0.02, 0.30, 60)
    q = _bet_isotherm(p)
    res = A.bet_rouquerol(p, q)
    assert res is not None
    assert abs(res["n_m"] - NM_TRUE) / NM_TRUE < 1e-6
    assert abs(res["C"] - C_TRUE) / C_TRUE < 1e-6
    assert abs(res["area"] - _area_from_nm(NM_TRUE)) / _area_from_nm(NM_TRUE) < 1e-6
    assert res["R2"] > 0.999999


def test_bet_fixed_range_recovers_known_monolayer_and_area():
    p = np.linspace(0.05, 0.30, 50)
    q = _bet_isotherm(p)
    res = A.bet_fixed_range(p, q, p_lo=0.05, p_hi=0.30)
    assert res is not None
    assert abs(res["n_m"] - NM_TRUE) / NM_TRUE < 1e-6
    assert abs(res["C"] - C_TRUE) / C_TRUE < 1e-6
    assert abs(res["area"] - _area_from_nm(NM_TRUE)) / _area_from_nm(NM_TRUE) < 1e-6
    assert not res["intercept_negative"]


def test_known_area_value_for_nm_one_mmol_per_gram():
    """Cross-check the unit conversion on its own: 1 cm3 STP per g of N2 is 0.162e-18 m2 x 6.022e23 / 22414
    = 4.35 m2 per g (the textbook value for 0.162 nm2 molecules)."""
    area = _area_from_nm(1.0)
    assert abs(area - 4.3525) < 1e-3
