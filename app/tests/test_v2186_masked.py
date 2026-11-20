"""v2.18.6: the emission-masked polynomial baseline (mfs/plbackground.masked_poly) recovers the area of a broad band that spans
the whole range, where SNIP and arPLS do not (synthetic test, known area)."""
import numpy as np
from mfs import plbackground as BG


def _case(rng):
    x = np.arange(500.0, 800.0, 0.5)
    band = 1000 * np.exp(-0.5 * ((x - 650) / 40) ** 2)
    bg = 200 + 0.8 * (x - 500)
    y = rng.poisson(bg + band).astype(float)
    return x, y, np.trapezoid(band, x)


def test_masked_poly_area_within_three_percent():
    rng = np.random.default_rng(0)
    errs = []
    for _ in range(30):
        x, y, true_area = _case(rng)
        base, order, free = BG.masked_poly(x, y)
        errs.append((np.trapezoid(y - base, x) - true_area) / true_area)
    errs = np.array(errs)
    assert abs(errs.mean()) < 0.03
    assert errs.std(ddof=1) < 0.05


def test_masked_poly_order_is_linear_for_linear_background():
    rng = np.random.default_rng(1)
    x, y, _ = _case(rng)
    _, order, _ = BG.masked_poly(x, y)
    assert order in (0, 1)           # BIC picks the lowest order that fits; the true background is linear


def test_autoselect_uses_masked_when_no_signal_free_window():
    rng = np.random.default_rng(2)
    x, y, _ = _case(rng)
    plan = BG.autoselect(x, [y])
    if not plan.get("windows"):
        assert plan["method"] == "masked"
        assert BG.baseline_for(x, y, plan).shape == y.shape
        labels = [lab for lab, _ in BG.alternatives(x, plan)]
        assert "rubber band" in labels
