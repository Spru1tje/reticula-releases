"""Precise sphere mode across sizes: discs of 8 to 14 px radius in one image, on our own absorption model.

Measured on this code (3 noise seeds, sigma 2): the error (fitted minus true radius) is +0.44 to +0.75 px per size,
with a spread up to 0.27 px at R = 14. The limits below are set above these values. This checks the size dependence
of the precise mode on synthetic discs only; it is not a validation of the instrument radius, and no private
measurement data is used or named here."""
import numpy as np
import pytest
from scipy import special

from mfs import particles as PA

H, W = 260, 300
CENT = [(60, 60, 8.0), (60, 180, 10.0), (60, 250, 12.0), (180, 60, 14.0), (180, 180, 9.0), (180, 250, 11.0)]


def _img(seed, noise=2.0):
    yy, xx = np.mgrid[:H, :W]
    absorb = np.zeros((H, W))
    for cy, cx, R in CENT:
        d = np.hypot(yy - cy, xx - cx)
        absorb += 0.6 * 0.5 * special.erfc((d - R) / (np.sqrt(2) * 1.0))
    bg = 200 - 20 * xx / W
    img = bg * np.exp(-absorb) + noise * np.random.default_rng(seed).standard_normal((H, W))
    return np.clip(img, 0, 255).astype(np.uint8)


@pytest.fixture(scope="module")
def errors_by_size():
    out = {}
    for seed in range(3):
        lab, f, hi, info = PA.segment_adaptive(_img(seed), sphere="precise", clear_border=True)
        for k, v in info.items():
            if not (isinstance(k, int) and isinstance(v, dict) and v.get("d_px")):
                continue
            cx, cy, _ = v["circle"]
            cy0, cx0, R = min(CENT, key=lambda c: (c[0] - cy) ** 2 + (c[1] - cx) ** 2)
            out.setdefault(R, []).append(v["d_px"] / 2 - R)
    return out


def test_every_size_is_found_in_every_seed(errors_by_size):
    assert sorted(errors_by_size) == [8.0, 9.0, 10.0, 11.0, 12.0, 14.0]
    assert all(len(e) == 3 for e in errors_by_size.values()), {k: len(v) for k, v in errors_by_size.items()}


def test_error_within_tolerance_per_size(errors_by_size):
    # tolerance: each fitted radius within 1.0 px of the true radius (measured largest error 0.80 px)
    for R, e in errors_by_size.items():
        assert np.abs(e).max() < 1.0, (R, e)


def test_spread_within_limit_per_size(errors_by_size):
    # limit: standard deviation below 0.5 px (measured up to 0.27 px at R = 14)
    for R, e in errors_by_size.items():
        assert np.std(e) < 0.5, (R, e)
