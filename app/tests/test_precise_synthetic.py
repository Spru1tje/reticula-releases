"""Synthetic regression test for the precise sphere mode (mfs.precise).

Model: blurred absorbing sphere, A*(1-exp(-k t))/k with t = 2 sqrt(r^2 - d^2), Gaussian blur s, nearest-edge padding.
Known and checked here:
  * the model reproduces a synthetic sphere exactly at the true parameters (residual 0);
  * an isolated sphere is recovered to within 0.01 px from a perturbed start, with free and with fixed s.
Known and NOT solved (open):
  * on grids of touching-spaced spheres the radius bias is within a few percent, but the spread between
    particles is still 4-8 % of the radius. Sizes below ~20 nm must not be published from this mode.
No real micrograph is used; this file is synthetic only.
"""
import numpy as np
import pytest
from scipy import ndimage as ndi
from scipy.optimize import least_squares

from mfs import precise as PR

S, K, A = 1.2, 0.08, 0.6


def _sphere(r, cx, cy, H, W, noise, rng):
    yy, xx = np.mgrid[0:H, 0:W].astype(float)
    t = 2 * np.sqrt(np.clip(r * r - ((yy - cy) ** 2 + (xx - cx) ** 2), 0, None))
    img = ndi.gaussian_filter(A * (-np.expm1(-K * t)) / K, S, mode="nearest")
    return img + rng.normal(0, noise, img.shape), yy, xx


def test_model_reproduces_synthetic_sphere_exactly():
    H = W = 60
    cx, cy, r = W / 2 + 0.3, H / 2 - 0.2, 8.0
    yy, xx = np.mgrid[0:H, 0:W].astype(float)
    t = 2 * np.sqrt(np.clip(r * r - ((yy - cy) ** 2 + (xx - cx) ** 2), 0, None))
    truth = ndi.gaussian_filter(A * (-np.expm1(-K * t)) / K, S, mode="nearest")
    model = PR._model(np.array([cx, cy, r, K, A, 0.0, S]), yy, xx, 0, 0, truth.shape)
    assert np.max(np.abs(model - truth)) < 1e-9


@pytest.mark.parametrize("r", [5.0, 8.0, 12.0])
@pytest.mark.parametrize("free_s", [True, False])
def test_isolated_sphere_recovered_from_perturbed_start(r, free_s):
    rng = np.random.default_rng(0)
    H = W = int(4 * r + 30)
    cx, cy = W / 2 + 0.3, H / 2 - 0.2
    z, yy, xx = _sphere(r, cx, cy, H, W, 0.02, rng)
    p0 = np.array([cx + 0.5, cy - 0.5, r * 0.9, 0.05, 0.5, 0.0, S])

    def res(p):
        q = p if free_s else np.append(p, S)
        return (PR._model(q, yy, xx, 0, 0, z.shape) - z).ravel()

    lo, hi = [-np.inf] * 6 + [0.2], [np.inf] * 6 + [3.0]
    p, lo, hi = (p0, lo, hi) if free_s else (p0[:6], lo[:6], hi[:6])
    sol = least_squares(res, p, bounds=(lo, hi), max_nfev=400)
    assert abs(sol.x[0] - cx) < 0.01 and abs(sol.x[1] - cy) < 0.01
    assert abs(sol.x[2] - r) < 0.01
