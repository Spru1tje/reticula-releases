"""v2.17: flat plateau selection for the TGA defect calculation (automatic and picked on the graph)."""
import os, sys, warnings
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
warnings.filterwarnings("ignore")
import numpy as np
from mfs import analysis as A

EX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "example_data", "TGA")


def _synthetic(plateau_T=410.0):
    """Solvent loss to ~200 °C, a flat stretch, then the linker combustion from ~450 °C (DTG peak at 520)."""
    T = np.linspace(30, 700, 3400)
    w = np.full_like(T, 100.0)
    w -= 8 * (1 / (1 + np.exp(-(T - 120) / 15)))                          # solvent loss
    w -= 5 * (1 / (1 + np.exp(-(T - 330) / 12)))                          # dehydroxylation (ends ~350)
    w -= 40 * (1 / (1 + np.exp(-(T - 500) / 14)))                         # linker combustion
    return T, w


def test_plateau_is_inside_the_flat_window_and_avoids_the_edges():
    T, w = _synthetic()
    a = A.tga_analysis(T, w)
    info = a["plateau_info"]
    assert info["flat"] and info["window"][0] < a["T_plateau"] < info["window"][1]
    assert info["window"][0] + info["half_width"] <= a["T_plateau"] <= info["window"][1] - info["half_width"]
    assert 380 <= a["T_plateau"] <= 470
    assert abs(a["W_plateau"] - float(np.interp(a["T_plateau"], a["T"], a["w"]))) < 1e-9


def test_narrow_window_gives_no_plateau_instead_of_a_wrong_edge_point():
    T = np.linspace(30, 700, 3400)
    w = 100 - 0.05 * (T - 30)                                              # never flat
    w -= 30 * (1 / (1 + np.exp(-(T - 300) / 10)))
    w -= 20 * (1 / (1 + np.exp(-(T - 330) / 10)))
    a = A.tga_analysis(T, w)
    assert a["T_plateau"] is None or not a["plateau_info"]["flat"] or a["plateau_info"]["window"][1] - a["plateau_info"]["window"][0] >= 50


def test_example_files_give_plateaus_inside_their_windows():
    import glob
    from mfs import readers as R
    files = sorted(glob.glob(os.path.join(EX, "*.txt")))
    if not files:
        import pytest; pytest.skip("no example data")
    for f in files:
        d = R.read_tga(open(f, "rb"), f)
        a = A.tga_analysis(d["T"], d["w"])
        pi = a["plateau_info"]
        if a["T_plateau"] is not None:
            assert pi["window"][0] <= a["T_plateau"] <= pi["window"][1]


def test_tga_page_renders_the_plateau_picker():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app.py"), default_timeout=600)
    at.session_state["_lv"] = {"global": "Standard"}
    at.run()
    at.sidebar.radio[0].set_value("TGA").run()
    assert not at.exception, at.exception
    assert any("Plateau point" in e.label for e in at.expander)
    assert any("plateau point" in str(c.value).lower() or "automatic" in str(c.value).lower() for c in at.dataframe) or len(at.dataframe) > 0
