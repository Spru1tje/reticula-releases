"""S-PL-UI (window-fit panel): the pure table function behind the experimental window-fit expander in the PL thermometry tab.
Synthetic spectra with a known ratio 3:1 only. No Streamlit is needed; the rendering is not tested here."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest

from mfs import pl_thermo as PT
from mfs.plratio import pv_unit

W1, W2 = (512.0, 540.0), (552.0, 585.0)


def _spec(seeds=(1, 2, 3)):
    x = np.arange(500.0, 600.0, 0.4)
    out = []
    for k, s in enumerate(seeds):
        y = 30.0 * (3000.0 * pv_unit(x, 524.0, 9.0, 0.5) + 1000.0 * pv_unit(x, 566.0, 10.0, 0.5)) + 40.0 + 0.2 * (x - 500.0)
        y = np.random.default_rng(s).poisson(y).astype(float)
        out.append({"T": 290.0 + 10.0 * k, "x": x, "y": y, "var": None})
    return out


def test_table_has_one_row_per_temperature_and_the_expected_columns():
    spec = _spec()
    t = PT.window_fit_table(spec, W1, W2, True, np.full(3, 3.0), np.full(3, 0.1))
    assert len(t["df"]) == 3
    assert list(t["df"].columns) == ["T (K)", "window ratio", "window error", "global ratio", "global error",
                                     "difference", "disagree (> 3 σ)", "flags"]
    assert list(t["df"]["T (K)"]) == [290.0, 300.0, 310.0]


def test_agreeing_global_fit_gives_no_disagreement():
    spec = _spec()
    t = PT.window_fit_table(spec, W1, W2, True, np.full(3, 3.0), np.full(3, 0.5))
    assert t["n_compared"] == 3
    assert t["n_disagree"] == 0


def test_a_clearly_wrong_global_ratio_is_flagged_as_disagreement():
    spec = _spec()
    t = PT.window_fit_table(spec, W1, W2, True, np.full(3, 2.0), np.full(3, 0.02))   # global says 2.0 with a tiny error; window says ~3.0
    assert t["n_disagree"] == 3
    assert bool(t["df"]["disagree (> 3 σ)"].all())


def test_window_errors_are_passed_through_not_invented():
    spec = _spec()
    t = PT.window_fit_table(spec, W1, W2, True, np.full(3, 3.0), np.full(3, 0.1))
    assert np.all(np.isfinite(t["df"]["window error"]))
    assert np.all(t["df"]["window error"] > 0)


def test_invalid_windows_raise_so_the_panel_can_show_a_warning():
    spec = _spec()
    with pytest.raises(ValueError):
        PT.window_fit_table(spec, (520.0, 560.0), (550.0, 580.0), True, np.full(3, 3.0), np.full(3, 0.1))


def test_weighting_choice_changes_the_table_and_sqrt_is_the_default():
    # the panel's weighting choice reaches the table: Poisson gives other numbers, and the default equals the explicit sqrt
    spec = _spec()
    default = PT.window_fit_table(spec, W1, W2, True, np.full(3, 3.0), np.full(3, 0.1))
    sqrt = PT.window_fit_table(spec, W1, W2, True, np.full(3, 3.0), np.full(3, 0.1), "sqrt")
    poisson = PT.window_fit_table(spec, W1, W2, True, np.full(3, 3.0), np.full(3, 0.1), "poisson")
    assert list(default["df"]["window ratio"]) == list(sqrt["df"]["window ratio"])
    assert list(poisson["df"]["window ratio"]) != list(sqrt["df"]["window ratio"])

