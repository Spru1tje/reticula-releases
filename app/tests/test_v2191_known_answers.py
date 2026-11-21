"""v2.19.1: known-answer tests for the round's additions (KAS activation energy, Vegard composition, Sauerbrey QCM mass,
Birch-Murnaghan equation of state, Randles-Sevcik peak current, capacity and coulombic efficiency, unit and format detection)."""
import math

import numpy as np
import pytest

from mfs import autoread as AR
from mfs import electrochem as EC
from mfs import eos as EOS
from mfs import kinetics_nonisothermal as KN
from mfs import qcm as QCM
from mfs import vegard as VG


# ---- KAS (Kissinger-Akahira-Sunose): exact data from a known Ea
def test_kas_recovers_known_activation_energy():
    Ea, c = 120_000.0, 20.0                       # J/mol, constant of ln(beta / T²) = c - Ea / (R T)
    T = np.array([600.0, 610.0, 620.0, 630.0])
    beta = T ** 2 * np.exp(c - Ea / (KN.R * T))   # heating rates that give exactly this T at each rate
    r = KN.kas(beta, T)
    assert r["status"] == "ok"
    assert r["Ea_kJ"] == pytest.approx(120.0, rel=1e-9)
    assert r["R2"] == pytest.approx(1.0, abs=1e-12)


def test_kas_needs_three_rates():
    assert KN.kas([5, 10], [600, 610])["status"] == "needs at least 3 heating rates"


# ---- Vegard: composition from lattice parameter
def test_vegard_linear_round_trip_gives_known_fraction():
    aA, aB = 5.0, 5.4
    a = VG.lattice_from_composition(0.25, aA, aB)
    assert a == pytest.approx(5.3, abs=1e-12)     # 0.25 * 5.0 + 0.75 * 5.4
    r = VG.composition_from_lattice(5.3, aA, aB)
    assert r["status"] == "ok" and r["x"] == pytest.approx(0.25, abs=1e-12)


def test_vegard_bowing_known_value_and_inverse():
    aA, aB, b = 5.0, 5.4, 0.1
    a = VG.lattice_from_composition(0.25, aA, aB, bowing=b)
    assert a == pytest.approx(5.3 - 0.1 * 0.25 * 0.75, abs=1e-12)       # 5.28125
    r = VG.composition_from_lattice(a, aA, aB, bowing=b)
    assert r["x"] == pytest.approx(0.25, abs=1e-9)


def test_vegard_outside_end_members_is_reported_not_guessed():
    assert VG.composition_from_lattice(5.6, 5.0, 5.4)["status"] == "outside range"
    with pytest.raises(ValueError):
        VG.lattice_from_composition(1.5, 5.0, 5.4)


# ---- Sauerbrey QCM: 5 MHz quartz
def test_sauerbrey_known_mass_for_100_hz_shift():
    r = QCM.mass_per_area_from_shift(-100.0, 5e6)
    assert r["status"] == "ok"
    assert r["mass_ug_per_cm2"] == pytest.approx(1.767, abs=0.005)      # 100 Hz on 5 MHz quartz
    assert r["sensitivity_ng_per_hz_cm2"] == pytest.approx(17.7, abs=0.05)  # Sauerbrey constant


def test_sauerbrey_large_shift_warns():
    assert "warning" in QCM.mass_per_area_from_shift(-200_000.0, 5e6)   # 4 % of f0


# ---- Birch-Murnaghan EOS: synthetic P(V) with known V0, B0, B0'
def test_birch_murnaghan_recovers_known_parameters():
    V0, B0, B0p = 100.0, 150.0, 4.5
    V = np.linspace(85.0, 105.0, 21)
    P = EOS.pressure_bm3(V, V0, B0, B0p)
    r = EOS.fit_birch_murnaghan(V, P)
    assert r["status"] == "ok"
    assert r["V0"] == pytest.approx(V0, rel=1e-6)
    assert r["B0"] == pytest.approx(B0, rel=1e-6)
    assert r["B0p"] == pytest.approx(B0p, rel=1e-6)
    assert r["rms_residual"] < 1e-6


def test_birch_murnaghan_zero_pressure_at_v0():
    assert EOS.pressure_bm3(100.0, 100.0, 150.0, 4.5) == pytest.approx(0.0, abs=1e-12)


# ---- Randles-Sevcik: n = 1, A = 1 cm², C = 1 mM, v = 0.1 V/s, D = 1e-5 cm²/s, 298.15 K
def test_randles_sevcik_known_peak_current():
    ip = EC.randles_sevcik_peak_current(1, 1.0, 1e-6, 0.1, 1e-5)
    assert ip == pytest.approx(2.687e-4, rel=2e-3)        # 0.268 mA
    assert EC.diffusion_from_peak_current(ip, 1, 1.0, 1e-6, 0.1) == pytest.approx(1e-5, rel=1e-9)


# ---- Battery cycling
def test_specific_capacity_and_coulombic_efficiency_known_values():
    assert EC.specific_capacity_mAh_g(0.1, 10.0, 0.005) == pytest.approx(200.0)   # 0.1 mA x 10 h / 5 mg
    assert EC.coulombic_efficiency_percent(150.0, 180.0) == pytest.approx(83.3333, rel=1e-5)
    with pytest.raises(ValueError):
        EC.specific_capacity_mAh_g(0.1, 10.0, 0.0)


# ---- Units: heat flow (mW) as power, wavelength in nm
def test_heat_flow_in_mW_converts_to_W():
    info = AR.detect_unit("Heat flow (mW)")
    assert info["quantity"] == "power" and info["base"] == "W"
    assert info["factor"] == pytest.approx(1e-3)


def test_wavelength_in_micrometre_converts_to_nm():
    text = "Wavelength (µm)\tAbsorbance\n0.5\t0.1\n0.6\t0.2\n0.7\t0.3\n"
    r = AR.read_text(text)
    assert r["data"]["Wavelength (µm)"] == pytest.approx([500.0, 600.0, 700.0])
    assert r["data"]["Absorbance"] == pytest.approx([0.1, 0.2, 0.3])   # unknown unit is kept unchanged


def test_heat_flow_table_end_to_end():
    text = "Heat flow (mW)\tTemperature (°C)\n2.5\t25\n3.0\t30\n4.0\t35\n"
    r = AR.read_text(text)
    assert r["data"]["Heat flow (mW)"] == pytest.approx([0.0025, 0.003, 0.004])
    assert r["data"]["Temperature (°C)"] == pytest.approx([298.15, 303.15, 308.15])


# ---- Decimal comma in a tab-separated export
def test_tab_export_with_decimal_comma_reads_numbers_correctly():
    text = "Time (min)\tSignal\n1,5\t2,25\n2,5\t3,5\n"
    r = AR.read_text(text)
    assert r["data"]["Time (min)"] == pytest.approx([90.0, 150.0])     # 1.5 min and 2.5 min in seconds
    assert r["data"]["Signal"] == pytest.approx([2.25, 3.5])


def test_tab_export_with_decimal_point_is_unchanged():
    r = AR.read_text("Signal\tTemp\n1.5\t2.25\n2.5\t3.5\n")
    assert r["data"]["Signal"] == pytest.approx([1.5, 2.5])


# ---- Format detection
def test_detect_format_jcamp_table_and_unknown():
    assert AR.detect_format("##TITLE= IR spectrum\n##JCAMP-DX= 5.00\n##XYDATA=(X++(Y..Y))\n") == "jcamp-dx"
    assert AR.detect_format("T;R\n1;2\n3;4\n") == "delimited table"
    assert AR.detect_format("just some words\nno numbers here\n") == "unknown"
