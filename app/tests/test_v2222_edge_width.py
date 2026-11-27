"""Estimate the edge width of the synthetic discs from the radial contrast profile (10 to 90 % rise distance / 2.563).

First measurement on this code (3 noise seeds, true edge width 1, 1.5, 2 px): estimates 1.46, 1.75, 1.98 px. The order is
right and the steps are about 0.25 px apart. At 1 px the pipeline adds about +0.46 px of width (the contrast map is smoothed).
The limits are set above the measured deviation; this is a synthetic check only, not a validation of the instrument's edge."""
import numpy as np
import pytest
from scipy import special, ndimage as ndi

from mfs import particles as PA

H, W = 260, 300
CENT = [(60, 60, 8.0), (60, 180, 10.0), (60, 250, 12.0), (180, 60, 14.0), (180, 180, 9.0), (180, 250, 11.0)]


def _img(seed, s, noise=2.0):
    yy, xx = np.mgrid[:H, :W]
    absorb = np.zeros((H, W))
    for cy, cx, R in CENT:
        d = np.hypot(yy - cy, xx - cx)
        absorb += 0.6 * 0.5 * special.erfc((d - R) / (np.sqrt(2) * s))
    bg = 200 - 20 * xx / W
    return np.clip(bg * np.exp(-absorb) + noise * np.random.default_rng(seed).standard_normal((H, W)), 0, 255).astype(np.uint8)


def edge_width(c, cx, cy, R):
    # radial profile of contrast, averaged over angles, sampled in 0.25 px steps around R
    th = np.linspace(0, 2*np.pi, 180, endpoint=False)
    ds = np.arange(R-5, R+5.01, 0.25)
    yy = cy + np.outer(np.sin(th), ds); xx = cx + np.outer(np.cos(th), ds)
    p = ndi.map_coordinates(c, [yy.ravel(), xx.ravel()], order=1).reshape(len(th), len(ds)).mean(0)
    inner = p[ds < R-4].mean(); outer = p[ds > R+4].mean(); span = inner - outer
    if span <= 0: return np.nan
    q = (p - outer) / span                      # 1 inside, 0 outside
    def cross(level):
        i = np.where(q <= level)[0][0]          # first sample where profile falls to level (moving outward)
        j = i - 1
        return ds[j] + (q[j]-level)/(q[j]-q[i])*(ds[i]-ds[j])
    d10, d90 = cross(0.1), cross(0.9)
    return (d10 - d90) / 2.5631                  # erfc edge: 10-90 % distance = 2.563 s


def _estimates(s):
    ests = []
    for seed in range(3):
        lab, f, hi, info = PA.segment_adaptive(_img(seed, s), sphere="precise", clear_border=True)
        c = info["_contrast"]
        for k, v in info.items():
            if not (isinstance(k, int) and isinstance(v, dict) and v.get("d_px")):
                continue
            cx, cy, r = v["circle"]
            e = edge_width(c, cx, cy, r)
            if np.isfinite(e):
                ests.append(e)
    return np.array(ests)


@pytest.mark.parametrize("s", [1.0, 1.5, 2.0])
def test_estimate_within_limit_of_true_width(s):
    a = _estimates(s)
    assert len(a) >= 10, len(a)
    assert abs(a.mean() - s) < 0.6, (s, a.mean())             # limit: measured deviation 0.46, 0.25, -0.02 px


def test_estimates_are_ordered_with_true_width():
    m = {s: _estimates(s).mean() for s in (1.0, 1.5, 2.0)}
    assert m[1.0] < m[1.5] - 0.15 and m[1.5] < m[2.0] - 0.15, m
