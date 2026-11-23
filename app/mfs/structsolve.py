"""v2.20.3: rigid-body structure search for a Zr6-cluster MOF (UiO-66 type) from a powder pattern.

Inputs (from the measurement and the own topology library; no CIF geometry):
  * observed |F|^2 of the reflections (from charge flipping / Le Bail), the lattice and the space-group operations;
  * density peaks: the 24 Zr atoms (grouped into 4 clusters of 6);
  * an ideal BDC linker (mfs/data/bdc_linker.xyz) and the fcu topology (mfs/topology.py).

Search: the linker centre starts on an fcu edge midpoint (origin from the clusters on the fcc nodes). Its orientation is
scanned on a grid (run in parallel), then refined by simulated annealing together with the translation and one radius
per Zr6 face: the mu3-O sites sit on the 8 faces of each octahedron, at a free distance from the face centre.
"""
from __future__ import annotations

import itertools
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.spatial.transform import Rotation as Rot

_DATA = Path(__file__).parent / "data" / "bdc_linker.xyz"


def _periodic(L, a, b):
    d = a - b
    d -= np.round(d)
    return np.linalg.norm(d @ L, axis=-1)


def load_linker(path=_DATA):
    """Ideal BDC (C8H4O4): element list and Cartesian coordinates (A), centred on the carbon centroid."""
    rows = np.loadtxt(path, skiprows=2, dtype=str)
    return list(rows[:, 0]), rows[:, 1:].astype(float)


def zr_clusters(zr_frac, L, cut=4.2):
    """Single-linkage groups of Zr atoms (periodic); returns the cluster label per Zr."""
    n = len(zr_frac)
    D = np.array([[_periodic(L, zr_frac[i], zr_frac[j]) for j in range(n)] for i in range(n)])
    return fcluster(linkage(squareform(D, checks=False), "single"), t=cut, criterion="distance")


def faces_from_clusters(zr_frac, lab, L):
    """The 8 triangular faces of every Zr6 octahedron as (centre, outward normal, circumradius), all in Cartesian A,
    and the fractional cluster centres."""
    faces, centres = [], []
    for k in np.unique(lab):
        q = zr_frac[lab == k]
        ref = q[0]
        qu = np.array([ref + ((p - ref) - np.round(p - ref)) for p in q])
        centres.append(qu.mean(axis=0) % 1)
        Dq = np.array([[np.linalg.norm((qu[a] - qu[b]) @ L) for b in range(len(qu))] for a in range(len(qu))])
        edge_cut = np.sort(Dq[np.triu_indices(len(qu), 1)])[11]           # the 12 shortest pairs are the octahedron edges
        for a, b, c in itertools.combinations(range(len(qu)), 3):
            if max(Dq[a, b], Dq[b, c], Dq[a, c]) <= edge_cut + 1e-6:
                tri = np.array([qu[a], qu[b], qu[c]]) @ L
                cen = tri.mean(axis=0)
                cl = (qu @ L).mean(axis=0)
                nd = cen - cl
                nd = nd / np.linalg.norm(nd)
                faces.append((cen, nd, float(np.linalg.norm(tri - cen, axis=1).mean())))
    return faces, np.array(centres)


def origin_shift(centres, nodes, L):
    """Shift that puts the cluster centres on the fcc nodes of the topology (best cluster/node pair)."""
    best = None
    for c in centres:
        for n in nodes:
            sh = (n - c) % 1
            dev = np.sqrt(np.mean([min(_periodic(L, cc + sh, nn) ** 2 for nn in nodes) for cc in centres]))
            if best is None or dev < best[0]:
                best = (float(dev), sh)
    return best[1], best[0]


def _scattering(elements, sl):
    """Cromer-Mann form factors f(s) for each element at s = sin(theta)/lambda (array sl)."""
    from pymatgen.analysis.diffraction.xrd import ATOMIC_SCATTERING_PARAMS as AS
    return {e: np.array([sum(a * np.exp(-b * v ** 2) for a, b in AS[e]) for v in sl]) for e in set(elements)}


class RigidCost:
    """Scaled R factor between observed and calculated |F|^2: one linker with its symmetry images, plus the Zr and mu3-O atoms."""

    def __init__(self, hkl, Fsq, L, sym_R, sym_t, lk_el, lk_xyz, lk_occ, zr_frac, faces):
        self.hkl = np.asarray(hkl, int)
        self.Fsq = np.asarray(Fsq, float)
        self.L = np.asarray(L, float)
        self.sym_R = np.asarray(sym_R)
        self.sym_t = np.asarray(sym_t)
        self.lk_el = list(lk_el)
        self.lk_xyz = np.asarray(lk_xyz, float)
        self.lk_occ = np.asarray(lk_occ, float)
        self.zr = np.asarray(zr_frac, float)
        self.faces = faces
        self.sl = np.linalg.norm(self.hkl @ np.linalg.inv(self.L).T, axis=1) / 2.0
        self.fe = _scattering(list(self.lk_el) + ["Zr", "O"], self.sl)

    def images(self, frac_lk):
        """All symmetry images of one linker, deduplicated by centroid."""
        out, cents = [], []
        is_c = np.array([e == "C" for e in self.lk_el])
        for Rg, tg in zip(self.sym_R, self.sym_t):
            f2 = frac_lk @ Rg.T + tg
            c = f2[is_c].mean(axis=0) % 1
            if all(np.linalg.norm(((c - co) + 0.5) % 1 - 0.5) > 1e-3 for co in cents):
                cents.append(c)
                out.append(f2)
        return out

    def mu3(self, radii):
        out = []
        for (cen, nd, Rc), r in zip(self.faces, radii):
            h = np.sqrt(max(r ** 2 - Rc ** 2, 0.0))
            out.append(np.linalg.solve(self.L.T, cen + h * nd) % 1)
        return np.array(out)

    def __call__(self, p):
        """p = [centre (3, fractional), rotation vector (3), one radius per face]."""
        centre, rv, radii = p[:3], p[3:6], p[6:6 + len(self.faces)]
        R_ = Rot.from_rotvec(rv).as_matrix()
        cart = self.lk_xyz @ R_.T + centre @ self.L
        ims = self.images(np.linalg.solve(self.L.T, cart.T).T)
        mu = self.mu3(radii)
        frac_all = np.vstack(ims + [self.zr, mu])
        occ = np.concatenate([np.tile(self.lk_occ, len(ims)), np.ones(len(self.zr) + len(mu))])
        F = np.zeros(len(self.hkl), complex)
        n_zr = len(self.zr)
        for c0 in range(0, len(self.hkl), 1500):
            s_ = slice(c0, c0 + 1500)
            ph = np.exp(2j * np.pi * (self.hkl[s_] @ frac_all.T))
            lk = np.tile(np.stack([self.fe[e][s_] for e in self.lk_el], axis=1), (1, len(ims)))
            fa = np.hstack([lk, np.repeat(self.fe["Zr"][s_, None], n_zr, axis=1), np.repeat(self.fe["O"][s_, None], len(mu), axis=1)])
            F[s_] = np.sum(fa * occ[None, :] * ph, axis=1)
        I = np.abs(F) ** 2
        k = np.sum(self.Fsq * I) / max(np.sum(I * I), 1e-30)
        return float(np.sum(np.abs(self.Fsq - k * I)) / np.sum(self.Fsq))


def search(cost, centre0, n_faces, grid=500, steps=2000, seed=0, workers=None, refine_top=3):
    """Orientation grid (parallel), then simulated annealing on the best orientations. Returns (best cost, best parameters)."""
    rng = np.random.default_rng(seed)
    rots = Rot.from_euler("ZYZ", rng.uniform([0, 0, 0], [2 * np.pi, np.pi, 2 * np.pi], size=(grid, 3)))
    base = np.concatenate([np.asarray(centre0, float), np.zeros(3), np.full(n_faces, 2.20)])
    cands = []
    for r in rots:
        q = base.copy()
        q[3:6] = r.as_rotvec()
        cands.append(q)
    workers = workers or min(8, (os.cpu_count() or 2))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        scores = list(ex.map(cost, cands))
    order = np.argsort(scores)[:refine_top]
    best_c, best_p = np.inf, None
    for i in order:
        cur, cc = cands[i].copy(), scores[i]
        bp, bc = cur.copy(), cc
        for it in range(steps):
            T = 0.01 * (0.05) ** (it / steps)
            fr = 1 - it / steps
            cand = cur.copy()
            u = rng.random()
            if u < 0.3:
                cand[:3] += rng.normal(0, 0.02 * fr + 0.002, 3)
            elif u < 0.7:
                cand[3:6] = (Rot.from_rotvec(rng.normal(0, 0.12 * fr + 0.01, 3)) * Rot.from_rotvec(cand[3:6])).as_rotvec()
            else:
                j = 6 + rng.integers(n_faces)
                cand[j] += rng.normal(0, 0.04 * fr + 0.003)
            cn = cost(cand)
            if cn < cc or rng.random() < np.exp(-(cn - cc) / T):
                cur, cc = cand, cn
                if cc < bc:
                    bp, bc = cur.copy(), cc
        if bc < best_c:
            best_c, best_p = bc, bp
    return best_c, best_p
