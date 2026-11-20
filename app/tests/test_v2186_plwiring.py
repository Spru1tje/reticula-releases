"""v2.18.6: the Poisson maximum-likelihood window fit is wired to the thermometry results (mfs/pl_thermo._mle_window_block)."""
import numpy as np
from mfs import pl_thermo as TH

KB = 8.617333262e-5
EV_CM = 8065.544
DE_TRUE = 700.0


def _spectra(rng, alpha=1.5, s2=2e4, bg=50.0):
    wl = np.arange(500.0, 800.0, 0.5)
    band1 = np.exp(-0.5 * ((wl - 520) / 4.0) ** 2)
    band2 = np.exp(-0.5 * ((wl - 545) / 4.0) ** 2)
    band1 = band1 / band1.sum()
    band2 = band2 / band2.sum()
    T = np.arange(100.0, 401.0, 20.0)
    out = []
    for t in T:
        D = alpha * np.exp(-DE_TRUE / EV_CM / (KB * t))
        expected = bg + D * s2 * band1 + s2 * band2
        raw = rng.poisson(expected).astype(float)
        out.append({"T": t, "wl": wl, "I": raw, "I_raw": raw, "baseline": np.full(len(wl), bg)})
    return out


def test_mle_window_block_recovers_the_gap():
    rng = np.random.default_rng(21)
    sp = _spectra(rng)
    info = {"windows": [(515.0, 525.0), (540.0, 550.0)]}
    r = TH._mle_window_block(sp, info, swapped=False, counts=True)
    assert r is not None and r["converged"]
    assert abs(r["dE_cm1"] - DE_TRUE) < 4 * r["sigma_cm1"] + 1.0
    assert r["n"] == len(sp)


def test_mle_window_block_swapped_gives_the_same_gap_with_sign_flipped_inputs():
    rng = np.random.default_rng(22)
    sp = _spectra(rng)
    info = {"windows": [(515.0, 525.0), (540.0, 550.0)]}
    r_norm = TH._mle_window_block(sp, info, swapped=False, counts=True)
    # the same spectra, the ratio inverted by the program: the helper must swap the windows back
    r_sw = TH._mle_window_block(sp, {"windows": [(540.0, 550.0), (515.0, 525.0)]}, swapped=True, counts=True)
    assert r_sw is not None
    assert abs(abs(r_sw["dE_cm1"]) - abs(r_norm["dE_cm1"])) < 2.0


def test_mle_window_block_not_applicable_without_counts_or_windows():
    rng = np.random.default_rng(23)
    sp = _spectra(rng)
    assert TH._mle_window_block(sp, {"windows": [(515.0, 525.0), (540.0, 550.0)]}, swapped=False, counts=False) is None
    assert TH._mle_window_block(sp, {"windows": None}, swapped=False, counts=True) is None
