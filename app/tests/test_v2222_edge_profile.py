"""Precise sphere mode with a second edge profile: the blur of the synthetic disc edge is 2 px instead of 1 px.

Measured on this code (3 noise seeds, sigma 2): with edge width 1 px the mean error is +0.68 px and 6 of 6 discs are found;
with edge width 2 px the mean error is +1.78 px (per size +1.41 to +2.10 px) and only 3 to 4 of 6 discs are found. So the
~0.7 px overestimation of the precise mode is NOT a general property: it depends on the edge profile. The limits below
document the second profile as measured; they are set above the measured values and are not a validation. Own absorption
model, not a validation of the instrument radius. The lower bound of 3 of 6 discs found is the MEASURED minimum at edge
width 2 px (seen: 3, 4, 4), kept as the bound on purpose. No private measurement data is used or named here."""
import numpy as np
import pytest
from scipy import special

from mfs import particles as PA

H, W = 260, 300
CENT = [(60, 60, 8.0), (60, 180, 10.0), (60, 250, 12.0), (180, 60, 14.0), (180, 180, 9.0), (180, 250, 11.0)]


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
    out, found = {}, []
    for seed in range(3):
        lab, f, hi, info = PA.segment_adaptive(_img(seed, s_blur), sphere="precise", clear_border=True)
        n = 0
        for k, v in info.items():
            if not (isinstance(k, int) and isinstance(v, dict) and v.get("d_px")):
                continue
            cx, cy, _ = v["circle"]
            _, _, R = min(CENT, key=lambda c: (c[0] - cy) ** 2 + (c[1] - cx) ** 2)
            out.setdefault(R, []).append(v["d_px"] / 2 - R)
            n += 1
        found.append(n)
    return out, found


@pytest.fixture(scope="module")
def edge2():
    return _errors(2.0)


def test_second_profile_finds_at_least_3_of_6(edge2):
    # limit: 3 of 6 is the measured lower bound for edge width 2 px (seen: 3, 4, 4)
    _, found = edge2
    assert min(found) >= 3, found


def test_second_profile_error_within_documented_limit(edge2):
    # limit: largest error below 2.5 px (measured largest 2.10 px); this documents the profile dependence, it does not
    # claim the precise mode is correct for edge width 2 px
    errors, _ = edge2
    assert max(np.abs(np.asarray(e)).max() for e in errors.values()) < 2.5


def test_bias_is_larger_for_the_second_profile_than_for_the_first(edge2):
    # the fact that matters: the bias is profile-dependent, so the 0.7 px claim cannot be stated as general
    errors2, _ = edge2
    errors1, _ = _errors(1.0)
    mean2 = np.mean(np.concatenate([np.asarray(e) for e in errors2.values()]))
    mean1 = np.mean(np.concatenate([np.asarray(e) for e in errors1.values()]))
    assert mean2 > mean1 + 0.5, (mean1, mean2)
