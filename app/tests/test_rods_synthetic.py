"""Synthetic regression test for the rod measurement (mfs.rods.segment_rods, summary).

Synthetic ellipses with known length and width (1 nm per pixel) on a flat background.
Measured biases are pinned per case (within 1 percentage point). There is no general tolerance; the target
of 5 % is an explicit xfail where the current method does not reach it:
  case (length L, width W, angle)   length bias   width bias
  40 x 8 nm,  0 deg                 +6.5 %        -3.0 %   (known: length target of 5 % not met, see xfail)
  40 x 8 nm, 15 deg                 +3.3 %        +8.1 %   (known width error +8.1 %, see xfail below)
  40 x 8 nm, 45 deg                 +3.3 %        +3.4 %
  60 x 10 nm, 30 deg                +0.2 %        +0.4 %
  50 x 12 nm, 75 deg                +0.7 %        +0.3 %
The width error at 40 x 8 nm, 15 deg (+8.1 %) is a known error of the current method; it is pinned at the measured
value and recorded as an xfail for the 5 % target. It is not hidden and no tolerance is set to pass it.
Also checked: summary() counts under key 'n'. (The rename to 'n_rods' on the TEM page is not tested here.)
No real micrograph is used; this file is synthetic only.
"""
import numpy as np
import pytest
from scipy import ndimage as ndi

from mfs import rods as RO


def _ellipse_image(L, W, angle_deg, size=200, centre=100):
    yy, xx = np.mgrid[0:size, 0:size].astype(float)
    img = np.full((size, size), 200.0)
    th = np.deg2rad(angle_deg)
    dx, dy = xx - centre, yy - centre
    u = dx * np.cos(th) + dy * np.sin(th)
    v = -dx * np.sin(th) + dy * np.cos(th)
    img[(u / (L / 2)) ** 2 + (v / (W / 2)) ** 2 <= 1] = 60.0
    return ndi.gaussian_filter(img, 0.7)


def _measure(L, W, angle):
    _, objs = RO.segment_rods(_ellipse_image(L, W, angle), 1.0, min_area_px=25)
    assert len(objs) == 1
    return objs[0]


PINNED = [  # (L, W, angle, length bias %, width bias %)
    (40.0, 8.0, 0.0, 6.5, -3.0),
    (40.0, 8.0, 15.0, 3.3, 8.1),
    (40.0, 8.0, 45.0, 3.3, 3.4),
    (60.0, 10.0, 30.0, 0.2, 0.4),
    (50.0, 12.0, 75.0, 0.7, 0.3),
]


@pytest.mark.parametrize("L,W,angle,len_bias,wid_bias", PINNED)
def test_measured_bias_is_pinned(L, W, angle, len_bias, wid_bias):
    o = _measure(L, W, angle)
    assert abs(100 * (o["length_nm"] - L) / L - len_bias) < 1.0
    assert abs(100 * (o["width_nm"] - W) / W - wid_bias) < 1.0


@pytest.mark.xfail(strict=True, reason="known length overestimate +6.5 % for 40 x 8 nm at 0 deg; the 5 % target is not met")
def test_length_target_5_percent_for_40x8_at_0deg():
    o = _measure(40.0, 8.0, 0.0)
    assert abs(o["length_nm"] - 40.0) / 40.0 < 0.05


@pytest.mark.xfail(strict=True, reason="known width error +8.1 % for 40 x 8 nm at 15 deg; the 5 % target is not met")
def test_width_target_5_percent_for_40x8_at_15deg():
    o = _measure(40.0, 8.0, 15.0)
    assert abs(o["width_nm"] - 8.0) / 8.0 < 0.05


def test_summary_counts_rods_under_key_n():
    o = _measure(60.0, 10.0, 30.0)
    RO.flag_agglomerates([o], 1.0)
    RO.in_circles([o], [])
    s = RO.summary([o])
    assert s["n"] == 1
    assert "length_nm_median" in s and "width_nm_median" in s
