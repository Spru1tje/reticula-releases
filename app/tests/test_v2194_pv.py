"""v2.19.4: the vectorised pseudo-Voigt for several bands gives the same numbers as the single-band function (known answer: unit area)."""
import numpy as np

from mfs import plratio as PR


def test_pv_matrix_matches_pv_unit_for_each_band():
    x = np.linspace(530.0, 560.0, 301)
    c = [541.5, 549.5]
    w = [3.0, 3.2]
    e = [0.3, 0.7]
    M = PR.pv_matrix(x, c, w, e)
    for k in range(2):
        assert np.allclose(M[:, k], PR.pv_unit(x, c[k], w[k], e[k]))


def test_pv_matrix_has_unit_area_per_band():
    x = np.linspace(400.0, 700.0, 6001)
    M = PR.pv_matrix(x, [550.0], [4.0], [0.5])
    area = np.trapezoid(M[:, 0], x)
    assert abs(area - 1.0) < 1e-3
