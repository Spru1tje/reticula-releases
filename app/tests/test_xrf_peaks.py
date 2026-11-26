"""Synthetic tests for mfs.xrf_peaks: a known background plus Gaussian peaks at known channels.
No real measurement is used. The 8 sigma default was chosen on this synthetic data, not validated on real spectra."""
import numpy as np

from mfs import xrf_peaks as XP


def _spectrum(seed=0):
    rng = np.random.default_rng(seed)
    ch = np.arange(1024, dtype=float)
    bg = 200 * np.exp(-ch / 150.0) + 30.0
    peaks = [(300.0, 4000.0), (520.0, 2500.0), (800.0, 900.0)]
    y = bg.copy()
    for c, a in peaks:
        y += a * np.exp(-0.5 * ((ch - c) / 3.0) ** 2)
    y = rng.poisson(y).astype(float)
    return y, bg, [c for c, _ in peaks]


def test_background_is_close_to_true_background_away_from_peaks():
    y, bg, _ = _spectrum()
    est = XP.snip_background(y)
    mask = np.ones(len(y), bool)
    for c in (300, 520, 800):
        mask[max(0, c - 15):c + 16] = False
    rel = np.median(np.abs(est[mask] - bg[mask]) / bg[mask])
    assert rel < 0.15


def test_weak_peak_150_counts_found():
    ch = np.arange(1024.0)
    bg = 200 * np.exp(-ch / 150) + 30
    for s in range(10):
        y = np.random.default_rng(s).poisson(bg + 150 * np.exp(-0.5 * ((ch - 520) / 3) ** 2)).astype(float)
        c, _ = XP.find_xrf_peaks(y)
        assert np.min(np.abs(c - 520)) <= 1.5


def test_peaks_found_at_known_channels():
    y, _, truth = _spectrum()
    ch, net = XP.find_xrf_peaks(y)
    for t in truth:
        assert np.min(np.abs(ch - t)) <= 1.5, (t, ch)
    assert len(ch) == len(truth)


def test_flat_spectra_no_false_peaks_at_default_threshold():
    """Regression check on seeds 100-299 (not the tuning seeds). The 8 sigma level was TUNED on seeds 0-19
    (there 7 sigma gave 3 false peaks, 8 sigma gave 0). Synthetic flat Poisson spectra only; not validated on real spectra."""
    total = sum(len(XP.find_xrf_peaks(np.random.default_rng(s).poisson(50.0, 1024).astype(float))[0]) for s in range(100, 300))
    assert total == 0
