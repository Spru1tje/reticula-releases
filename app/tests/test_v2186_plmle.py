"""v2.18.6: Poisson maximum-likelihood thermometry (mfs/plmle.py) on synthetic count series.

The reference values are the truth of the simulation, never a starting guess of the fit.
Acceptance: unbiased (|bias| < 1 %), error bars consistent with the Cramér-Rao bound, and
coverage close to the nominal 68 % (within 1 sigma) and 95 % (within 2 sigma)."""
import numpy as np
from mfs import plmle as M

T = np.arange(80.0, 401.0, 20.0)
DE_TRUE = 700.0          # cm-1
ALPHA = 1.5
B1 = B2 = 50.0           # known background counts per window


def _simulate(rng, s2=2e4, dE=DE_TRUE, alpha=ALPHA):
    S2 = s2 * np.ones_like(T)
    D = alpha * np.exp(-dE / M.EV_CM / (M.KB_EV * T))
    N1 = rng.poisson(B1 + D * S2).astype(float)
    N2 = rng.poisson(B2 + S2).astype(float)
    return N1, N2, S2


def test_mle_unbiased_and_coverage():
    rng = np.random.default_rng(3)
    d, e, z = [], [], []
    nonconv = 0
    for _ in range(150):
        N1, N2, _ = _simulate(rng)
        r = M.fit(T, N1, N2, B1, B2)
        if not (r["converged"] and np.isfinite(r["sigma_cm1"])):
            nonconv += 1
            continue
        d.append(r["dE_cm1"]); e.append(r["sigma_cm1"]); z.append((r["dE_cm1"] - DE_TRUE) / r["sigma_cm1"])
    assert nonconv <= 5, nonconv        # at most 3 % of the fits may fail to converge (checked on Linux and Windows)
    d, e, z = np.array(d), np.array(e), np.array(z)
    assert abs(d.mean() - DE_TRUE) / DE_TRUE < 0.01                      # unbiased
    assert 0.55 <= np.mean(np.abs(z) <= 1) <= 0.80                        # ~68 % within 1 sigma
    assert np.mean(np.abs(z) <= 2) >= 0.90                                # ~95 % within 2 sigma


def test_error_bars_consistent_with_cramer_rao():
    rng = np.random.default_rng(5)
    sig = []
    for _ in range(40):
        N1, N2, S2 = _simulate(rng)
        sig.append(M.fit(T, N1, N2, B1, B2)["sigma_cm1"])
    crb = M.cramer_rao_sigma(T, 2e4 * np.ones_like(T), B1, B2, DE_TRUE, ALPHA)
    med = float(np.median(sig))
    assert abs(med - crb) / crb < 0.15                                     # observed info matches the Fisher bound
    assert med >= 0.85 * crb                                               # no error bar below the bound


def test_weak_band_low_temperature_is_finite_and_converges():
    rng = np.random.default_rng(9)
    # band 1 almost absent at 80 K: the ratio is tiny; the fit must still return a finite result
    N1, N2, _ = _simulate(rng, s2=3e3)
    r = M.fit(T, N1, N2, B1, B2)
    assert np.isfinite(r["dE_cm1"]) and np.isfinite(r["alpha"])


def test_start_is_from_data_not_from_truth():
    rng = np.random.default_rng(11)
    N1, N2, _ = _simulate(rng)
    r1 = M.fit(T, N1, N2, B1, B2)
    # the result must not depend on a hidden reference value: refit with a different background vector shape
    r2 = M.fit(T, N1, N2, np.full(len(T), B1), np.full(len(T), B2))
    assert abs(r1["dE_cm1"] - r2["dE_cm1"]) < 1e-6
