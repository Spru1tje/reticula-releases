"""Geometric (probe-occupiable) pore volume from a CIF - used as the ideal V0 for the porosity-based defect estimate."""
from __future__ import annotations
import numpy as np

VDW = {"H": 1.10, "C": 1.70, "O": 1.52, "N": 1.55, "S": 1.80, "Zr": 2.23, "Hf": 2.23, "Cl": 1.75, "Ce": 2.42, "F": 1.47}  # Bondi / Alvarez


def probe_pore_volume(cif_path, probe=1.82, spacing=0.35, substitute=None):
    """Probe-occupiable pore volume (cm3 g-1) and void fraction: grid points whose distance to every atom centre exceeds
    r_vdW + r_probe are probe-accessible centres; the occupiable volume is obtained by dilating them by r_probe.
    N2 kinetic radius 1.82 Å. Disordered/partial sites are taken at full occupancy (as written in the CIF)."""
    import warnings
    warnings.filterwarnings("ignore")
    from pymatgen.core import Structure
    s = Structure.from_file(cif_path)
    if substitute:
        s.replace_species(substitute)
    lat = s.lattice
    n = np.maximum((np.array(lat.abc) / spacing).astype(int), 8)
    fr = np.stack(np.meshgrid(*[np.arange(k) / k for k in n], indexing="ij"), -1).reshape(-1, 3)
    cart_atoms = s.frac_coords
    radii = np.array([VDW.get(max(site.species.items(), key=lambda kv: kv[1])[0].symbol, 1.8) for site in s])
    free = np.ones(len(fr), bool)
    M = lat.matrix
    for k in range(0, len(cart_atoms)):
        d = fr - cart_atoms[k]
        d -= np.round(d)
        dc = d @ M
        free &= (dc ** 2).sum(1) > (radii[k] + probe) ** 2
    # dilate accessible centres by the probe radius (occupiable volume)
    centers = fr[free]
    occ = np.zeros(len(fr), bool)
    if len(centers):
        # use KD-tree on a 3x3x3 periodic supercell of the centres
        from scipy.spatial import cKDTree
        shifts = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
        sc = (centers[None, :, :] + shifts[:, None, :]).reshape(-1, 3) @ M
        tree = cKDTree(sc)
        dist, _ = tree.query(fr @ M, k=1, distance_upper_bound=probe + 1e-6)
        occ = np.isfinite(dist)
    frac = occ.mean()
    vol_A3 = lat.volume * frac
    mass_g = float(s.composition.weight) / 6.02214076e23
    return {"V_cm3_g": vol_A3 * 1e-24 / mass_g, "void_fraction": frac, "cell_A3": lat.volume, "formula": s.composition.reduced_formula}
