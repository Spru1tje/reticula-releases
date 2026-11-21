"""v2.18.7 hotfix: two crashes reported in the photoluminescence tab (Windows build).

1. fit_power with equal powers raised numpy.linalg.LinAlgError (singular matrix).
2. The panel composer called fig.canvas.get_renderer(), which a plain FigureCanvasBase does not have.
"""
import numpy as np
from matplotlib.figure import Figure

from mfs import plextras as PX
from mfs import composer
from mfs.render import renderer_for


def test_fit_power_known_slope_still_works():
    P = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    I = 3.0 * P ** 2.0                                  # exact quadratic: n = 2
    f = PX.fit_power(P, I)
    assert f["status"] == "ok"
    assert abs(f["n"] - 2.0) < 1e-9


def test_fit_power_equal_powers_does_not_raise():
    f = PX.fit_power(np.array([0.25, 0.25, 0.25]), np.array([1.0, 2.0, 3.0]))
    assert f["status"] == "needs at least two different powers"


def test_renderer_for_plain_canvas_does_not_raise():
    fig = Figure()                                       # FigureCanvasBase: no get_renderer()
    assert not hasattr(fig.canvas, "get_renderer")
    r = renderer_for(fig)
    assert r is not None


def test_group_boxes_on_plain_canvas_figure():
    fig = Figure()
    ax = fig.add_subplot(111)
    ax.plot([0, 1], [0, 1])
    boxes = composer._group_boxes(fig)
    assert boxes, "one panel should give one box"
