"""Precise sphere mode per edge width, in one table: 1, 1.5 and 2 px (synthetic discs only).

Each row gives the limits for one edge width. They are set above the values measured on this code (3 noise seeds):
edge 1 px: 6 of 6 found, mean error +0.68 px, largest 0.80 px;
edge 1.5 px: 6 of 6 found, mean +1.33 px, largest 1.54 px;
edge 2 px: 3 to 4 of 6 found, mean +1.78 px, largest 2.10 px.
The minimum found is the MEASURED lower bound for each row (not chosen to pass). The bias depends on the edge width;
this is not a validation of the instrument radius. No private measurement data is used or named here."""
import numpy as np
import pytest
from scipy import special

from mfs import particles as PA

H, W = 260, 300
CENT = [(60, 60, 8.0), (60, 180, 10.0), (60, 250, 12.0), (180, 60, 14.0), (180, 180, 9.0), (180, 250, 11.0)]

# edge width (px), minimum discs found, max |error| limit (px), mean error range (px)
TABLE = [
    (1.0, 6, 1.0, (0.4, 1.0)),
    (1.5, 6, 2.0, (1.0, 1.7)),
    (2.0, 3, 2.5, (1.4, 2.2)),
]


def _img(seed, s_blur, noise=2.0):
    yy, xx = np.mgrid[:H, :W]
    absorb = np.zeros((H, W))
    for cy, cx, R in CENT:
        d = np.hypot(yy - cy, xx - cx)
        absorb += 0.6 * 0.5 * special.erfc((d - R) / (np.sqrt(2) * s_blur))
    bg = 200 - 20 * xx / W
    img = bg * np.exp(-absorb) + noise * np.random.default_rng(seed).standard_normal((H, W))
    return np.clip(img, 0, 255).astype(np.uint8)


def _errors(s_blur):
    out, found = [], []
    for seed in range(3):
        lab, f, hi, info = PA.segment_adaptive(_img(seed, s_blur), sphere="precise", clear_border=True)
        n = 0
        for k, v in info.items():
            if not (isinstance(k, int) and isinstance(v, dict) and v.get("d_px")):
                continue
            cx, cy, _ = v["circle"]
            _, _, R = min(CENT, key=lambda c: (c[0] - cy) ** 2 + (c[1] - cx) ** 2)
            out.append(v["d_px"] / 2 - R)
            n += 1
        found.append(n)
    return np.array(out), found


@pytest.mark.parametrize("s_blur, min_found, max_err, mean_range", TABLE)
def test_profile_row_within_its_limits(s_blur, min_found, max_err, mean_range):
    err, found = _errors(s_blur)
    assert min(found) >= min_found, (s_blur, found)
    assert np.abs(err).max() < max_err, (s_blur, np.abs(err).max())
    lo, hi = mean_range
    assert lo < err.mean() < hi, (s_blur, err.mean())


def test_bias_grows_across_the_table():
    means = [_errors(s)[0].mean() for s, *_ in TABLE]
    assert means[0] < means[1] < means[2], means
