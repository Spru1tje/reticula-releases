"""PL thermometry on user windows (mfs.plwindows): synthetic spectra with a KNOWN ratio. Only the data inside the windows are used."""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest

from mfs import plwindows as PW
from mfs.plratio import pv_unit

W1, W2 = (512.0, 540.0), (552.0, 585.0)


def _spectrum(A1, A2, seed, c1=524.0, c2=566.0, w1=9.0, w2=10.0, scale=1.0):
    rng = np.random.default_rng(seed)
    x = np.arange(500.0, 600.0, 0.4)
    y = A1 * pv_unit(x, c1, w1, 0.5) + A2 * pv_unit(x, c2, w2, 0.5)
    y = scale * y + 40.0 + 0.2 * (x - 500.0)
    return x, rng.poisson(y).astype(float)


def test_ratio_recovered_on_clean_synthetic_pair():
    # areas: I1 = 3 * I2 for a clean pair with plenty of counts
    x, y = _spectrum(3000.0, 1000.0, seed=1, scale=30.0)
    rows = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)
    r = rows[0]
    assert not any("poor fit" in f or "not significant" in f for f in r["flags"])
    assert abs(r["ratio"] - 3.0) / 3.0 < 0.03
    assert np.isfinite(r["ratio_err"]) and r["ratio_err"] > 0


def test_error_is_calibrated_over_seeds():
    # z = (ratio - truth) / ratio_err over 60 noise realisations should have spread ~1 (the error is neither too small nor too large)
    z = []
    for seed in range(60):
        x, y = _spectrum(3000.0, 1000.0, seed=100 + seed, scale=8.0)
        r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
        if np.isfinite(r["ratio_err"]):
            z.append((r["ratio"] - 3.0) / r["ratio_err"])
    z = np.asarray(z)
    assert len(z) >= 55
    assert 0.7 < np.std(z) < 1.3
    assert abs(np.mean(z)) < 0.5


def test_only_window_data_are_used():
    # changing the spectrum OUTSIDE both windows must not change the result
    x, y = _spectrum(3000.0, 1000.0, seed=7, scale=30.0)
    r1 = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
    y2 = y.copy()
    outside = (x < W1[0]) | ((x > W1[1]) & (x < W2[0])) | (x > W2[1])
    y2[outside] = 1e6 * np.arange(outside.sum()) % 97 + 5e5
    r2 = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y2}], W1, W2)[0]
    assert r1["ratio"] == pytest.approx(r2["ratio"], rel=1e-12)


def test_series_keeps_order_and_temperature():
    specs = []
    for k, T in enumerate([290.0, 310.0, 330.0]):
        x, y = _spectrum(3000.0 - 500 * k, 1000.0, seed=20 + k, scale=30.0)
        specs.append({"T": T, "x": x, "y": y})
    rows = PW.window_ratio_series(specs, W1, W2)
    assert [r["T"] for r in rows] == [290.0, 310.0, 330.0]
    assert rows[0]["ratio"] > rows[1]["ratio"] > rows[2]["ratio"]


def test_window_cutting_the_band_is_flagged():
    # window 1 ends 4 nm after the band centre: the band tail is cut, the area is an extrapolation
    x, y = _spectrum(3000.0, 1000.0, seed=3, scale=30.0)
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], (512.0, 528.0), W2)[0]
    assert any("I1: tail outside window" in f for f in r["flags"])


def test_too_few_points_is_reported_not_fitted():
    x, y = _spectrum(3000.0, 1000.0, seed=4, scale=30.0)
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], (524.0, 525.0), W2)[0]
    assert any("few points" in f for f in r["flags"])
    assert math.isnan(r["ratio"])
    assert any("ratio undefined" in f for f in r["flags"])


def test_no_signal_in_a_window_gives_no_ratio():
    x, y = _spectrum(3000.0, 0.0, seed=5, scale=30.0)     # band 2 absent
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
    assert math.isnan(r["ratio"])
    assert any("ratio undefined" in f for f in r["flags"])


def test_overlapping_windows_are_refused():
    x, y = _spectrum(3000.0, 1000.0, seed=6, scale=30.0)
    with pytest.raises(ValueError, match="overlap"):
        PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], (520.0, 560.0), (550.0, 580.0))


@pytest.mark.parametrize("bad", [None, (540.0, 512.0, 3.0), (np.nan, 540.0), (530.0, 530.0)])
def test_bad_window_input_is_refused(bad):
    x, y = _spectrum(3000.0, 1000.0, seed=8, scale=30.0)
    with pytest.raises(ValueError):
        PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], bad, W2)


def test_insignificant_band_gives_no_ratio_and_a_flag():
    x, y = _spectrum(3000.0, 0.0, seed=9, scale=30.0)     # band 2 absent: only the tail of band 1 reaches window 2
    a = PW.fit_window(x, y, W2, None, True)
    assert "not significant" in a["flags"] or "no area" in a["flags"]
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
    assert math.isnan(r["ratio"])
