"""EXPERIMENTAL Hirshfeld surface (point cloud) and 2D fingerprint plot for one molecule in its crystal.
Promolecule densities are single exponentials per atom, ρ_i(r) = Z_i·exp(−r/λ_i) with λ_i ∝ r_vdW (not the Spackman
spherical-atom tables used by CrystalExplorer), so surfaces and d_norm values are approximate; use CrystalExplorer
for publication. The surface is found along rays from each molecule atom (w = ρ_mol/ρ_crystal = 0.5)."""
from __future__ import annotations
import numpy as np
from scipy.spatial import cKDTree
from .cif import VDW

ZN = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "S": 16, "Cl": 17, "Br": 35, "Zr": 40, "Hf": 72, "P": 15}


def _rho(X, El, pts, tree):
    out = np.zeros(len(pts))
    for j, nb in enumerate(tree.query_ball_point(pts, 8.0)):
        if nb:
            r = np.linalg.norm(X[nb] - pts[j], axis=1)
            lam = np.array([VDW.get(El[k], 1.7) for k in nb]) / 5.0
            out[j] = np.sum(np.array([ZN.get(El[k], 6) for k in nb]) * np.exp(-r / lam))
    return out


def surface(S, mol_sites, env_sites, n_dir=120, seed=0):
    """mol_sites: the molecule (list of site dicts); env_sites: surrounding atoms (packing ± 1 cell). Returns dict with
    points, di, de, d_norm, and the element pair of the nearest inside/outside atoms."""
    rng = np.random.default_rng(seed)
    Xm = np.array([s["xyz"] for s in mol_sites]); Em = [s["el"] for s in mol_sites]
    keym = {(s["i"], tuple(np.round(s["frac"], 3))) for s in mol_sites}
    env = [s for s in env_sites if (s["i"], tuple(np.round(s["frac"], 3))) not in keym and
           np.min(np.linalg.norm(Xm - s["xyz"], axis=1)) < 8.0]
    Xe = np.array([s["xyz"] for s in env]) if env else np.zeros((0, 3)); Ee = [s["el"] for s in env]
    tm, te = cKDTree(Xm), (cKDTree(Xe) if len(Xe) else None)
    u = rng.normal(size=(n_dir, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
    pts = []
    for i, x in enumerate(Xm):
        lo = np.full(n_dir, 0.3); hi = np.full(n_dir, 4.5)
        for _ in range(14):
            mid = (lo + hi) / 2; p = x + u * mid[:, None]
            rm = _rho(Xm, Em, p, tm); re_ = _rho(Xe, Ee, p, te) if te is not None else np.zeros(len(p))
            w = rm / np.maximum(rm + re_, 1e-30)
            lo, hi = np.where(w > 0.5, mid, lo), np.where(w > 0.5, hi, mid)
        p = x + u * lo[:, None]
        dmin, jmin = tm.query(p)                    # keep points whose nearest molecule atom is this atom (outer surface)
        pts.append(p[jmin == i])
    P = np.vstack(pts)
    di, ii = tm.query(P)
    if te is not None:
        de, ie = te.query(P)
    else:
        de, ie = np.full(len(P), np.nan), np.zeros(len(P), int)
    ri = np.array([VDW.get(Em[k], 1.7) for k in ii]); re2 = np.array([VDW.get(Ee[k], 1.7) for k in ie]) if te is not None else ri
    dn = (di - ri) / ri + (de - re2) / re2
    pair = [f"{Em[a]}···{Ee[b]}" for a, b in zip(ii, ie)] if te is not None else [""] * len(P)
    return {"points": P, "di": di, "de": de, "d_norm": dn, "pair": np.array(pair)}


def draw_fingerprint(ax, hs, highlight=None, bins=60, rng=(0.6, 2.6)):
    import matplotlib.colors as mc
    H, xe, ye = np.histogram2d(hs["di"], hs["de"], bins=bins, range=[rng, rng])
    ax.imshow(np.ma.masked_equal(H.T, 0), origin="lower", extent=[*rng, *rng], cmap="jet", norm=mc.LogNorm(), aspect="equal")
    if highlight:
        m = np.isin(hs["pair"], highlight)
        frac = m.mean() * 100
        ax.text(0.04, 0.95, f"{', '.join(highlight)}: {frac:.1f} %", transform=ax.transAxes, fontsize=6, va="top")
    ax.set_xlabel(r"$d_\mathrm{i}$ (Å)"); ax.set_ylabel(r"$d_\mathrm{e}$ (Å)")
    ax.text(0.97, 0.03, "EXPERIMENTAL", transform=ax.transAxes, ha="right", fontsize=5, color="#B00")


def contacts(hs):
    import pandas as pd
    p = pd.Series(["···".join(sorted(x.split("···"))) for x in hs["pair"]]).value_counts(normalize=True) * 100
    return p.rename("% of surface").reset_index().rename(columns={"index": "contact"})


def draw_dnorm(ax, hs, V):
    P = hs["points"] @ V.T; o = np.argsort(P[:, 2])
    lim = max(abs(np.nanmin(hs["d_norm"])), 0.1)
    sc = ax.scatter(P[o, 0], P[o, 1], c=hs["d_norm"][o], cmap="bwr_r", vmin=-lim, vmax=lim * 3, s=2, lw=0)
    ax.set_aspect("equal"); ax.set_axis_off()
    ax.text(0.02, 0.02, "d_norm (EXPERIMENTAL point cloud)", transform=ax.transAxes, fontsize=5, color="#B00")
    return sc
