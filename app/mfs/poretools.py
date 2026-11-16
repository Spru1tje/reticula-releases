"""Geometric pore analysis on a periodic grid (Zeo++-like definitions, grid-based):
Di = largest included sphere (2 × max distance to the vdW surface), Df = largest free sphere (largest sphere that can
percolate through the cell along at least one axis; binary search with periodic union-find), void fraction (helium-
like, probe radius 0 by default), probe-occupiable fraction, and the accessible surface area by Monte Carlo on probe-
inflated atomic spheres (Düren et al. 2007), here with Bondi/Alvarez vdW radii (not UFF)."""
from __future__ import annotations
import math
import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage as ndi
from .cif import VDW, MASS

NA = 6.02214076e23


def cell_atoms(S, part="major"):
    sites, _ = S.grow("cell", part=part)
    F = np.array([s["frac"] for s in sites]); F = F - np.floor(F)
    el = [s["el"] for s in sites]
    occ = np.array([S.atoms[s["i"]]["occ"] for s in sites])
    return F, el, occ


def _tree(S, F, el):
    sh = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    P = (F[None] + sh[:, None]).reshape(-1, 3)
    R = np.tile([VDW.get(e, 1.8) for e in el], len(sh))
    return cKDTree(S.cart(P)), R


def distance_grid(S, spacing=0.4, part="major"):
    F, el, occ = cell_atoms(S, part)
    n = np.maximum(np.ceil(np.array([c[0] for c in S.cell[:3]]) / spacing).astype(int), 8)
    g = np.stack(np.meshgrid(*[(np.arange(k) + 0.5) / k for k in n], indexing="ij"), -1).reshape(-1, 3)
    tree, R = _tree(S, F, el)
    d, idx = tree.query(S.cart(g), k=min(12, len(R)))
    D = np.min(d - R[idx], axis=1)
    return D.reshape(n), n, (F, el, occ)


def _percolates(mask):
    """True if a connected component of `mask` (6-connected, periodic) wraps around the cell along some axis."""
    lab, nl = ndi.label(mask)
    if nl == 0:
        return False
    parent = list(range(nl + 1)); off = [np.zeros(3, int) for _ in range(nl + 1)]

    def find(a):
        if parent[a] == a:
            return a, np.zeros(3, int)
        r, o = find(parent[a]); parent[a] = r; off[a] = off[a] + o; return r, off[a]
    for ax in range(3):
        a = np.take(lab, -1, axis=ax); b = np.take(lab, 0, axis=ax)
        pairs = set(zip(a[(a > 0) & (b > 0)].ravel(), b[(a > 0) & (b > 0)].ravel()))
        sh = np.zeros(3, int); sh[ax] = 1                                    # going from the last layer into the next cell
        for p, q in pairs:
            rp, op = find(p); rq, oq = find(q)
            if rp == rq:
                if np.any(op + sh - oq != 0):
                    return True
            else:
                parent[rq] = rp; off[rq] = op + sh - oq
    return False


def analyse(S, spacing=0.4, probe=1.82, part="major", n_mc=200, seed=0):
    D, n, (F, el, occ) = distance_grid(S, spacing, part)
    Di = 2 * float(D.max())
    lo, hi = 0.0, float(D.max())
    if not _percolates(D > 0):
        Df = 0.0
    else:
        for _ in range(18):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if _percolates(D > mid) else (lo, mid)
        Df = 2 * lo
    vf = float((D > 0).mean())
    # probe-occupiable: centres with D > probe, dilated by the probe radius
    cen = D > probe
    occ_frac = 0.0
    if cen.any():
        idx = np.argwhere(cen); g = (idx + 0.5) / n
        sh = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
        tr = cKDTree(S.cart((g[None] + sh[:, None]).reshape(-1, 3)))
        all_g = np.stack(np.meshgrid(*[(np.arange(k) + 0.5) / k for k in n], indexing="ij"), -1).reshape(-1, 3)
        dd, _ = tr.query(S.cart(all_g), k=1, distance_upper_bound=probe + 1e-6)
        occ_frac = float(np.isfinite(dd).mean())
    # accessible surface area (Monte Carlo)
    rng = np.random.default_rng(seed)
    tree, R = _tree(S, F, el)
    area = 0.0
    u = rng.normal(size=(n_mc, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
    X = S.cart(F)
    for i, e in enumerate(el):
        r = VDW.get(e, 1.8) + probe
        pts = X[i] + u * r
        dd, ii = tree.query(pts, k=min(16, len(R)))
        inside = np.any(dd < (R[ii] + probe) - 1e-6, axis=1)
        area += (1 - inside.mean()) * 4 * math.pi * r * r
    _, el_all, occ_all = cell_atoms(S, "all")                              # mass: every site × its occupancy
    mass = float(sum(MASS.get(e, 50.0) * o for e, o in zip(el_all, occ_all))) / NA
    return {"Di (Å)": Di, "Df (Å)": Df, "void fraction (probe 0)": vf, f"probe-occupiable fraction (r = {probe} Å)": occ_frac,
            "ASA (Å² per cell)": area, "ASA (m² g⁻¹)": area * 1e-20 / mass, "ASA (m² cm⁻³)": area * 1e-20 / (S.V * 1e-24),
            "ρ (g cm⁻³)": mass / (S.V * 1e-24), "grid (Å)": spacing, "atoms in cell": len(el),
            "method": ("grid distance to the vdW surface (Bondi/Alvarez radii), Di = largest included sphere, Df = largest free "
                       "sphere (periodic percolation), void fraction = grid points outside all vdW spheres; ASA = Monte Carlo on "
                       f"vdW+probe spheres (probe r = {probe} Å, N₂ 3.64 Å kinetic diameter), {n_mc} points/atom; "
                       "inaccessible pockets are not excluded. Geometry uses the major disorder part at full radii; the mass uses every site × occupancy.")}
