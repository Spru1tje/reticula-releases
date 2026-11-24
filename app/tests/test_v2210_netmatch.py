"""v2.21: net matcher. Synthetic checks only: metal-pair test, complete reflection list, node-model ranking.
No CIF geometry is used here."""
import numpy as np

from mfs import netmatch as NM

A = 20.0


def _fcu_reflections(n=8):
    """Synthetic fcu intensities: |F|^2 from the 4 fcc nodes of a single-metal net (absent for mixed parity)."""
    g = np.arange(-n, n + 1)
    H = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T.astype(float)
    H = H[np.abs(H).sum(1) > 0]
    fcc = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    F = np.abs(np.exp(2j * np.pi * (H @ fcc.T)).sum(1)) ** 2
    return H, F


def test_single_metal_is_not_paired():
    rng = np.random.default_rng(0)
    frac = rng.random((24, 3))                     # well-separated single metals
    kind, _ = NM.classify_metal_sublattice(frac, np.ones(24), 30.0, K=24)
    assert kind == "single"


def test_metal_pairs_are_detected():
    base = np.array([[0.1, 0.2, 0.3], [0.6, 0.7, 0.2], [0.3, 0.8, 0.9]])
    pair = [base[0], base[0] + np.array([2.6 / A, 0, 0])]  # Cu-Cu 2.6 A along x
    frac = np.vstack([pair, [base[1], base[1] + np.array([2.6 / A, 0, 0])],
                      [base[2], base[2] + np.array([2.6 / A, 0, 0])]])
    kind, meta = NM.classify_metal_sublattice(frac, np.ones(len(frac)), A, K=len(frac))
    assert kind == "paired" and meta["frac_pair"] == 1.0


def test_complete_reflections_fill_absences_with_zero():
    H, F = _fcu_reflections()
    obs_idx = np.where(F > 0)[0]
    Hc, Fc = NM.complete_reflections(H[obs_idx], F[obs_idx], A)
    assert len(Hc) > len(obs_idx)
    assert np.count_nonzero(Fc) == len(obs_idx)


def test_fcu_wins_over_pcu_when_absences_are_in_the_data():
    H, F = _fcu_reflections()
    res = NM.rank_nets(H, F, nets=["fcu", "pcu"])
    assert res["best"] == "fcu"
    assert res["confident"]


def test_pattern_matcher_routes_single_metal_to_single_nets():
    H, F = _fcu_reflections()
    rng = np.random.default_rng(1)
    frac = rng.random((24, 3))
    out = NM.rank_nets_pattern(H, F, frac, np.ones(24), A)
    assert out["metal_type"] == "single"
    assert "tbo" not in out["candidates"]
