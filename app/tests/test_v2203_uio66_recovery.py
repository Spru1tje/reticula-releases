"""v2.20.3: regression test - charge flipping must still recover the UiO-66 framework (Zr and O) from simulated powder data with
overlapping reflections and 5 % intensity noise. Chance level is about 5-10 %; the threshold is far above it."""
import math
import os

import numpy as np
import pytest

from mfs import chargeflip as CF

CIF = os.path.join(os.path.dirname(__file__), "..", "references", "UiO-66_Zr_COD4512074.cif")


def _simulate(dmin=1.2, noise_pct=5.0, dq=0.01):
    pytest.importorskip("pymatgen")
    from pymatgen.core import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.analysis.diffraction.xrd import ATOMIC_SCATTERING_PARAMS
    s = Structure.from_file(CIF)
    cell = s.lattice.abc
    ops = SpacegroupAnalyzer(s, symprec=0.01).get_symmetry_operations(cartesian=False)
    R = np.array([o.rotation_matrix for o in ops]).round().astype(int)
    t = np.array([o.translation_vector for o in ops])
    Gs = np.linalg.inv(s.lattice.matrix).T
    hm = [int(math.floor(cell[i] / dmin)) for i in range(3)]
    hkl = np.array([(h, k, l) for h in range(-hm[0], hm[0] + 1) for k in range(-hm[1], hm[1] + 1) for l in range(-hm[2], hm[2] + 1)
                    if (h, k, l) != (0, 0, 0) and 1 / np.linalg.norm(Gs @ np.array([h, k, l])) >= dmin], int)
    sites = [(sp.symbol, float(o), site.frac_coords) for site in s for sp, o in site.species.items()]
    dstar = np.linalg.norm(hkl @ Gs.T, axis=1)
    Xs = np.array([x for _, _, x in sites]); occs = np.array([o for _, o, _ in sites]); els = [e for e, _, _ in sites]
    F = np.zeros(len(hkl), complex)
    for c0 in range(0, len(hkl), 1500):
        sl = dstar[c0:c0 + 1500] / 2.0
        ph = np.exp(2j * np.pi * (hkl[c0:c0 + 1500] @ Xs.T))
        fe = {e: np.array([sum(a * np.exp(-b * v ** 2) for a, b in ATOMIC_SCATTERING_PARAMS[e]) for v in sl]) for e in set(els)}
        F[c0:c0 + 1500] = np.sum(np.stack([fe[e] for e in els], axis=1) * occs[None, :] * ph, axis=1)
    Fsq = np.abs(F) ** 2
    keep = Fsq > 1e-6 * Fsq.max()
    hkl, Fsq, dst = hkl[keep], Fsq[keep], dstar[keep]
    q = np.round(dst / dq).astype(int)
    _, group = np.unique(q, return_inverse=True)
    G = np.bincount(group, weights=Fsq)
    G = np.clip(G * (1 + noise_pct / 100 * np.random.default_rng(1).standard_normal(len(G))), 0, None)
    return s, cell, hkl, group, G, R, t


def _matched(s, cell, rho, shift_from_zr=True):
    pk = np.array([p[:3] for p in CF.peaks(rho, tuple(cell) + (90.0, 90.0, 90.0), n_max=400)])
    zr = np.array([site.frac_coords for site in s if site.species_string.startswith("Zr")])
    def dist(a, b):
        d = a - b
        d -= np.round(d)
        return np.linalg.norm(d @ s.lattice.matrix, axis=-1)
    best = (-1.0, None)
    for p in pk[:60]:
        cand = zr + (p - zr[0])
        frac = float(np.mean(np.min(dist(cand[:, None, :], pk[None, :, :]), axis=1) < 0.5))
        if frac > best[0]:
            best = (frac, p - zr[0])
    shift = best[1]
    out = {}
    for el in ("Zr", "O"):
        pos = np.array([site.frac_coords for site in s if site.species_string.startswith(el)]) + shift
        out[el] = float(np.mean(np.min(dist(pos[:, None, :], pk[None, :, :]), axis=1) < 0.7))
    return out


def test_uio66_framework_is_recovered_from_overlapping_noisy_data():
    s, cell, hkl, group, G, R, t = _simulate()
    P = CF.Problem(cell, hkl, group, G, 1.2, sym=(R, t))
    res = CF.run_trials(P, n_trials=4, n_iter=400, seed0=0)
    m = _matched(s, cell, res[0]["rho"])
    assert m["Zr"] >= 0.9, m
    assert m["O"] >= 0.8, m
