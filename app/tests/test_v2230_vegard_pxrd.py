"""v2.23.0: Vegard shift of reflections and composition estimate (cubic end members), pure functions in mfs/vegard.py.
Numbers are synthetic (Cu and Au as cubic end members for arithmetic only, not a claim about a Cu–Au alloy).
Sources of the formula: Vegard (1921), Z. Phys. 5, 17–26, doi:10.1007/BF01349680; bowing form as in mfs/vegard.py
(Denton & Ashcroft 1991, PRA 43, 3161: listed, formula content not checked against the paper)."""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pytest
from mfs import vegard as VG

WL = 1.540593   # Cu Kα1, Å
A_CU, A_AU = 3.615, 4.0782   # synthetic end members (fcc metals, cubic cells)


def test_shift_to_the_same_cell_is_the_identity():
    th = VG.two_theta_for_lattice(A_CU, 3, WL)
    assert VG.shifted_two_theta(th, A_CU, A_CU) == pytest.approx(th, abs=1e-9)


def test_larger_cell_moves_the_peak_to_lower_angle():
    th = VG.two_theta_for_lattice(A_CU, 3, WL)
    assert VG.shifted_two_theta(th, A_CU, A_AU) < th


def test_shift_agrees_with_direct_bragg_for_the_new_cell():
    th_ref = VG.two_theta_for_lattice(A_CU, 4, WL)
    th_direct = VG.two_theta_for_lattice(A_AU, 4, WL)
    assert VG.shifted_two_theta(th_ref, A_CU, A_AU) == pytest.approx(th_direct, abs=1e-9)


@pytest.mark.parametrize("x", [0.0, 0.25, 0.5, 0.8, 1.0])
def test_linear_vegard_round_trip_from_one_peak(x):
    a_x = VG.lattice_from_composition(x, A_CU, A_AU)
    th = VG.two_theta_for_lattice(a_x, 3, WL)
    r = VG.x_from_two_theta(th, 3, WL, A_CU, A_AU)
    assert r["status"] == "ok"
    assert r["x"] == pytest.approx(x, abs=1e-6)


def test_bowed_round_trip_from_one_peak():
    b = 0.05
    x = 0.4
    a_x = VG.lattice_from_composition(x, A_CU, A_AU, b)
    th = VG.two_theta_for_lattice(a_x, 3, WL)
    r = VG.x_from_two_theta(th, 3, WL, A_CU, A_AU, b)
    assert r["status"] == "ok"
    assert r["x"] == pytest.approx(x, abs=1e-6)


def test_peak_outside_the_end_members_is_not_estimated():
    # a cell larger than both end members: no composition may be returned
    th = VG.two_theta_for_lattice(4.5, 3, WL)
    r = VG.x_from_two_theta(th, 3, WL, A_CU, A_AU)
    assert r["status"] == "outside range"
    assert "x" not in r


def test_non_positive_angle_is_rejected():
    r = VG.x_from_two_theta(0.0, 3, WL, A_CU, A_AU)
    assert r["status"] == "no peak"


def test_reflection_beyond_the_wavelength_limit_is_nan():
    assert math.isnan(VG.two_theta_for_lattice(0.5, 3, WL))   # d < λ/2: no reflection


def test_vegard_panel_runs_and_plots_the_shifted_figure(monkeypatch):
    """Smoke test of the PXRD panel (bare mode, plotting helpers stubbed as the app injects them).
    gemmi and pymatgen are declared requirements: no importorskip, a missing package must fail."""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mfs.pages import spectra as SP
    calls = {}

    def newfig(aspect=1.0, **k):
        return plt.subplots()

    def show(fig, name, df, draw, title, units=None):
        calls["name"] = name
        calls["cols"] = list(df.columns)

    monkeypatch.setattr(SP, "newfig", newfig, raising=False)
    monkeypatch.setattr(SP, "show", show, raising=False)
    x = np.linspace(5, 60, 2000)
    y = 100 + 1000 * np.exp(-((x - 43.4) / 0.1) ** 2)
    SP._vegard_ui({"S1": (x, y)}, 5.0, 60.0)
    assert calls.get("name") == "PXRD_Vegard_S1"
    assert calls.get("cols") == ["two_theta", "observed"]


def test_multi_peak_estimate_gives_the_mean_and_zero_spread_for_exact_peaks():
    """v2.23.1 (S-Structuur-2): several exact peaks of the same composition give that x with zero spread."""
    x = 0.3
    a_x = VG.lattice_from_composition(x, A_CU, A_AU)
    peaks = [(VG.two_theta_for_lattice(a_x, n, WL), n) for n in (3, 4, 8, 11)]
    r = VG.x_from_peaks(peaks, WL, A_CU, A_AU)
    assert r["status"] == "ok" and r["n"] == 4
    assert r["x"] == pytest.approx(x, abs=1e-6)
    assert r["sd"] == pytest.approx(0.0, abs=1e-6)


def test_one_bad_peak_makes_the_whole_estimate_fail():
    a_x = VG.lattice_from_composition(0.5, A_CU, A_AU)
    good = (VG.two_theta_for_lattice(a_x, 3, WL), 3)
    bad = (VG.two_theta_for_lattice(4.5, 4, WL), 4)      # outside the two end members
    r = VG.x_from_peaks([good, bad], WL, A_CU, A_AU)
    assert r["status"] == "outside range"
    assert "x" not in r and len(r["bad_peaks"]) == 1


def test_a_single_peak_is_too_few_for_a_spread():
    a_x = VG.lattice_from_composition(0.5, A_CU, A_AU)
    r = VG.x_from_peaks([(VG.two_theta_for_lattice(a_x, 3, WL), 3)], WL, A_CU, A_AU)
    assert r["status"] == "too few peaks"


# ---------------------------------------------------------------- v2.23.1 (S-Structuur-3): non-cubic systems
# d-spacings are checked against pymatgen's Lattice.d_hkl (an independent implementation), on repo reference cells
CORUNDUM = (4.7606, 4.7606, 12.994, 90.0, 90.0, 120.0)   # hexagonal, repo reference
RUTILE = (4.5937, 4.5937, 2.9587, 90.0, 90.0, 90.0)      # tetragonal, repo reference


@pytest.mark.parametrize("cell,hkl", [(CORUNDUM, (1, 0, 4)), (CORUNDUM, (0, 1, 2)), (CORUNDUM, (1, 1, 0)),
                                      (RUTILE, (1, 1, 0)), (RUTILE, (1, 0, 1)), (RUTILE, (2, 1, 3))])
def test_d_spacing_matches_pymatgen(cell, hkl):
    # no importorskip: pymatgen is a declared requirement; if it is missing this test must FAIL, not pass silently
    from pymatgen.core import Lattice
    ref = Lattice.from_parameters(*cell).d_hkl(hkl)
    assert VG.d_spacing(cell, hkl) == pytest.approx(ref, rel=1e-9)


def test_system_detection_on_repo_cells():
    assert VG.system_of(CORUNDUM) == "hexagonal"
    assert VG.system_of(RUTILE) == "tetragonal"
    assert VG.system_of((5.4113, 5.4113, 5.4113, 90, 90, 90)) == "cubic"
    assert VG.system_of((4.0, 5.0, 6.0, 90, 90, 90)) == ""


def test_hexagonal_vegard_endpoints_and_midpoint():
    B = (4.9, 4.9, 13.5, 90.0, 90.0, 120.0)
    assert VG.cell_from_composition(1.0, CORUNDUM, B) == pytest.approx(CORUNDUM, abs=1e-9)
    assert VG.cell_from_composition(0.0, CORUNDUM, B) == pytest.approx(B, abs=1e-9)
    mid = VG.cell_from_composition(0.5, CORUNDUM, B)
    assert mid[0] == pytest.approx((CORUNDUM[0] + B[0]) / 2, abs=1e-9)
    assert mid[2] == pytest.approx((CORUNDUM[2] + B[2]) / 2, abs=1e-9)


def test_different_systems_are_refused():
    with pytest.raises(ValueError):
        VG.cell_from_composition(0.5, CORUNDUM, (5.4113, 5.4113, 5.4113, 90, 90, 90))
