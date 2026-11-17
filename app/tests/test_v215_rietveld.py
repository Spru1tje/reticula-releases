"""v2.15: Rietveld / Le Bail validated on synthetic patterns with known answers (own forward model + the independent pymatgen calculation)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest
from mfs import rietveld as RV, analysis as A

REF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "references")
UIO = os.path.join(REF, "UiO-66_Zr_COD4512074.cif")
M801 = os.path.join(REF, "MOF-801_Zr_COD4121459.cif")


def synth(cifs, true, scales, seed=3, step=0.02, rng_=(3, 60), peak=20000, rad="Cu Kα1"):
    x = np.arange(rng_[0], rng_[1], step)
    phs = [RV.Phase.from_cif(c) for c in cifs]
    R = RV.Refinement(x, np.ones_like(x) * 1000.0, phs, rad, 6, "rietveld", counts=True)
    for k, v in true.items():
        R.P[k].value = v
    y = R._design(R._values()) @ np.array(scales)
    y = y / y.max() * peak + 300 + 1.2 * x
    return x, np.random.default_rng(seed).poisson(y).astype(float), R


def test_structure_factors_agree_with_pymatgen():
    ph = RV.Phase.from_cif(UIO).build(1.2)
    sim, sticks, _ = A.simulate_pxrd(UIO, (3, 45))
    d = 1 / np.sqrt(np.einsum("ki,ij,kj->k", ph.hkl, np.linalg.inv(RV._cell_metric(ph.cell0)), ph.hkl))
    th = np.arcsin(1.540593 / (2 * d))
    I = ph.I0 * (1 + np.cos(2 * th) ** 2) / (np.sin(th) ** 2 * np.cos(th))
    tt = 2 * np.degrees(th)
    m = (tt > 3) & (tt < 45)
    I = I[m] / I[m].max() * 100
    tt = tt[m]
    for k in np.argsort(I)[::-1][:6]:
        j = np.argmin(abs(sticks.two_theta - tt[k]))
        assert abs(sticks.two_theta[j] - tt[k]) < 0.01 and abs(sticks.I[j] - I[k]) < 0.05 * max(I[k], 1)


def test_single_phase_recovers_cell_B_and_width_fast():
    true = {"c0.a": 20.80, "zero": 0.05, "W": 0.004, "X": 0.12, "Y": 0.05, "B0": 1.5}
    x, y, _ = synth([UIO], true, [1.0])
    import time
    t = time.time()
    r = RV.quick_refine(x, y, [UIO])
    assert time.time() - t < 15                                   # (≈ 1 s on a normal PC; CI machines are slower)
    p = r["phases"][0]
    assert abs(p["cell"][0] - 20.80) < 4 * max(p["cell_esd"]["a"], 2e-4) and abs(p["cell"][0] - 20.80) < 0.003
    assert abs(p["B"] - 1.5) < 0.4 and abs(r["global"]["zero"]["value"] - 0.05) < 0.03
    assert r["stats"]["chi2"] < 1.5 and r["stats"]["Rwp"] < 0.07 and p["Rb"] < 0.08
    assert abs(r["global"]["X"]["value"] - 0.12) < 0.02
    assert p["cell_esd"]["a"] > 0


def test_two_phases_weight_fractions_and_cells():
    true = {"c0.a": 20.78, "c1.a": 17.88, "zero": 0.04, "W": 0.004, "X": 0.10, "B0": 1.0, "B1": 1.0}
    x, y, R = synth([UIO, M801], true, [1.0, 0.35])
    zm = [1.0 * R.phases[0].mass_cell * RV.cell_volume(R.cell(0)), 0.35 * R.phases[1].mass_cell * RV.cell_volume(R.cell(1))]
    wt = zm[0] / sum(zm)
    r = RV.quick_refine(x, y, [UIO, M801])
    assert abs(r["phases"][0]["weight_fraction"] - wt) < 0.02
    assert abs(r["phases"][0]["cell"][0] - 20.78) < 0.004 and abs(r["phases"][1]["cell"][0] - 17.88) < 0.004
    assert r["stats"]["chi2"] < 1.6
    assert abs(sum(p["weight_fraction"] for p in r["phases"]) - 1) < 1e-9


def test_nanocrystals_broad_peaks_and_cell_start_far_off():
    """Broad peaks (≈ 1.5 nm crystals look like X ≈ 0.9°), cell 1.3 % away from the CIF: the start-value search must find it."""
    true = {"c0.a": 21.02, "zero": -0.08, "W": 0.0, "X": 0.5, "Y": 0.0, "B0": 2.0}
    x, y, _ = synth([UIO], true, [1.0], peak=6000)
    r = RV.quick_refine(x, y, [UIO])
    assert abs(r["phases"][0]["cell"][0] - 21.02) < 0.05
    assert r["stats"]["Rwp"] < 0.12


def test_kalpha2_radiation_and_lebail():
    true = {"c0.a": 20.76, "zero": 0.0, "W": 0.003, "X": 0.08, "B0": 1.0}
    x, y, _ = synth([UIO], true, [1.0], rad="Cu Kα1 + Kα2", rng_=(4, 50))
    r = RV.quick_refine(x, y, [UIO], radiation="Cu Kα1 + Kα2")
    assert abs(r["phases"][0]["cell"][0] - 20.76) < 0.003 and r["stats"]["chi2"] < 1.6
    lb = RV.quick_refine(x, y, [UIO], radiation="Cu Kα1 + Kα2", mode="lebail")
    assert abs(lb["phases"][0]["cell"][0] - 20.76) < 0.005
    assert lb["stats"]["Rwp"] < 0.08 and lb["lebail_I"] is not None
    assert lb["phases"][0].get("Rb") is None                       # no structure model: no Bragg R


def test_arbitrary_units_and_cell_systems():
    x, y, _ = synth([UIO], {"c0.a": 20.80, "W": 0.004, "X": 0.1}, [1.0])
    r = RV.quick_refine(x, y / 1000.0 + 0.123, [UIO])                  # not counts
    assert not r["counts"] and abs(r["phases"][0]["cell"][0] - 20.80) < 0.01 and any("not in counts" in n for n in r["notes"])
    assert RV.detect_system((5, 5, 5, 90, 90, 90)) == "cubic" and RV.detect_system((5, 5, 8, 90, 90, 120)) == "hexagonal"
    assert RV.detect_system((5, 6, 7, 90, 100, 90)) == "monoclinic" and RV.detect_system((5, 5, 5, 60, 60, 60)) == "rhombohedral"
