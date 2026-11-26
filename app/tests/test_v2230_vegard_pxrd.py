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
    """Smoke test of the PXRD panel (bare mode, plotting helpers stubbed as the app injects them)."""
    pytest.importorskip("gemmi")
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
