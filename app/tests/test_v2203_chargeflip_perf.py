"""v2.20.3: the symmetrise step of the charge-flipping solver is a sparse matrix product; it must equal the orbit-average loop."""
import numpy as np

from mfs import chargeflip as CF


def _naive(P, f):
    pos, ph, valid = P.orb
    acc = np.zeros_like(f)
    cnt = np.zeros(len(f))
    for g in range(len(pos)):
        acc += np.where(valid[g], ph[g] * f[pos[g]], 0.0)
        cnt += valid[g]
    return acc / np.maximum(cnt, 1.0)


def _problem():
    # cubic cell, reflections of a small sphere, cubic point group (the symmetry operations of m-3m, in fractional form)
    cell = (8.0, 8.0, 8.0, 90.0, 90.0, 90.0)
    hkl = np.array([(h, k, l) for h in range(-3, 4) for k in range(-3, 4) for l in range(-3, 4) if (h, k, l) != (0, 0, 0)], int)
    group = np.arange(len(hkl))
    G = np.ones(len(hkl))
    ops = []
    import itertools
    for perm in itertools.permutations(range(3)):
        for sgn in itertools.product([1, -1], repeat=3):
            R = np.zeros((3, 3), int)
            for i, j in enumerate(perm):
                R[i, j] = sgn[i]
            ops.append(R)
    R = np.array(ops)
    t = np.zeros((len(ops), 3))
    return CF.Problem(cell, hkl, group, G, 1.0, sym=(R, t))


def test_symmetrise_matches_the_orbit_loop():
    P = _problem()
    rng = np.random.default_rng(7)
    f = rng.standard_normal(P.M) + 1j * rng.standard_normal(P.M)
    assert np.max(np.abs(P.symmetrise(f) - _naive(P, f))) < 1e-10


def test_symmetrise_without_symmetry_is_identity():
    cell = (8.0, 8.0, 8.0, 90.0, 90.0, 90.0)
    hkl = np.array([(1, 0, 0), (0, 1, 0), (0, 0, 1)], int)
    P = CF.Problem(cell, hkl, np.arange(3), np.ones(3), 1.0)
    f = np.array([1 + 2j, 3j, -1.0])
    assert np.allclose(P.symmetrise(f), f)
