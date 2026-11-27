"""Precise sphere mode across three edge profiles (blur width 1, 1.5 and 2 px) on synthetic discs.

Measured on this code (3 noise seeds, sigma 2): edge 1 px: 6 of 6 found, mean error +0.68 px, largest 0.80 px.
Edge 1.5 px: 6 of 6 found, mean +1.33 px, largest 1.54 px. Edge 2 px: 3 to 4 of 6 found, mean +1.78 px, largest 2.10 px.
So the bias grows with the edge width: the ~0.7 px overestimation is specific to the 1 px profile. The limits below are set
above the measured values and document this dependence; they are not a validation. The 6 of 6 discs found at edge 1 and
1.5 px is the MEASURED lower bound (every seed found all 6), kept as the bound on purpose. Own absorption model, not a
validation of the instrument radius. No private measurement data is used or named here."""
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


def _flat(errors):
    return np.concatenate([np.asarray(e) for e in errors.values()])


@pytest.mark.parametrize("s_blur, max_err, mean_lo, mean_hi", [
    (1.0, 1.0, 0.4, 1.0),    # measured: largest 0.80, mean +0.68
    (1.5, 2.0, 1.0, 1.7),    # measured: largest 1.54, mean +1.33
])
def test_profile_error_within_measured_limits(s_blur, max_err, mean_lo, mean_hi):
    errors, found = _errors(s_blur)
    assert min(found) >= 6, found                               # measured: all 6 found in every seed at 1 and 1.5 px
    flat = _flat(errors)
    assert np.abs(flat).max() < max_err, (s_blur, np.abs(flat).max())
    assert mean_lo < flat.mean() < mean_hi, (s_blur, flat.mean())


def test_bias_grows_with_edge_width():
    # the fact that matters: the bias is profile-dependent (measured +0.68, +1.33, +1.78 for 1, 1.5, 2 px)
    m1 = _flat(_errors(1.0)[0]).mean()
    m15 = _flat(_errors(1.5)[0]).mean()
    m2 = _flat(_errors(2.0)[0]).mean()
    assert m1 < m15 < m2, (m1, m15, m2)
