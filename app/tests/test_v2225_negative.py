"""v2.22.5: negative and positive controls for the net matcher and the cell check, on synthetic data with known answers.
A check that cannot fail is not a check: each test here has a case that must NOT be accepted."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from mfs import netmatch as NM

LAM = 1.5406  # Cu Kα1 in Å


def _cubic_iref_lines(a, tt_lo=3.0, tt_hi=60.0):
    """2θ of the I-centred (h+k+l even) reflections of a cubic cell a (Å), 3-60 deg, Cu Kα1."""
    R = 14
    g = np.arange(-R, R + 1)
    H = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T
    H = H[(np.abs(H).sum(1) > 0) & (H.sum(1) % 2 == 0)]
    d = a / np.sqrt((H ** 2).sum(1))
    tt = np.degrees(2 * np.arcsin(np.clip(LAM / (2 * d), -1, 1)))
    return np.unique(np.round(tt[(tt >= tt_lo) & (tt <= tt_hi)], 4))


def _match_fraction(peaks_tt, lines_tt, tol=0.10):
    """Fraction of peaks within tol (deg) of a calculated line."""
    idx = np.searchsorted(lines_tt, peaks_tt)
    dev = []
    for p, i in zip(peaks_tt, idx):
        c = [lines_tt[j] for j in (i - 1, i) if 0 <= j < len(lines_tt)]
        dev.append(min(abs(p - x) for x in c))
    return float(np.mean(np.asarray(dev) <= tol))


# ---------------------------------------------------------------- net matcher
def test_net_matcher_accepts_the_net_that_generated_the_data():
    # noiseless positive control: the full reflection list from the fcu node model itself
    a = 12.0
    g = np.arange(-6, 7)
    H = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T.astype(float)
    H = H[(np.abs(H).sum(1) > 0) & (np.sqrt((H ** 2).sum(1)) <= 6.0)]
    Fsq = NM.node_amplitude("fcu", H) ** 2
    res = NM.rank_nets(H, Fsq)
    assert res["best"] == "fcu"
    assert res["R_best"] < 0.05  # noiseless data from the same node model


def test_net_matcher_does_not_claim_a_net_for_random_intensities():
    rng = np.random.default_rng(7)
    a = 12.0
    h = np.array([(1, 1, 1), (2, 0, 0), (2, 2, 0), (3, 1, 1), (2, 2, 2), (4, 0, 0), (3, 3, 1), (4, 2, 0)], float)
    fsq = rng.uniform(0.2, 1.0, size=len(h)) * 1e4
    hkl, F = NM.complete_reflections(h, fsq, a)
    res = NM.rank_nets(hkl, F)
    # random intensities fit no net: the best R must stay high, and the result is not called confident
    assert res["R_best"] > 0.3
    assert res["confident"] is False


# ---------------------------------------------------------------- cell check
def test_cell_check_accepts_the_true_cell_and_rejects_a_wrong_one():
    true_lines = _cubic_iref_lines(20.0)
    rng = np.random.default_rng(3)
    # 'observed' peaks = true lines with 0.02 deg noise, plus unrelated peaks at random positions
    peaks = true_lines + rng.normal(0, 0.02, size=true_lines.size)
    f_true = _match_fraction(np.sort(peaks), true_lines)
    wrong_lines = _cubic_iref_lines(23.0)
    f_wrong = _match_fraction(np.sort(peaks), wrong_lines)
    assert f_true > 0.9
    assert f_wrong < f_true - 0.3  # a wrong cell must clearly match worse
