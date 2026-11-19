import numpy as np
from mfs.rietveld import _cell_metric


def compare(true, solved, tol=0.4):
    """Fraction of the atoms of `true` found (same element) in `solved` within tol (A), tolerant to origin shift, inversion and a permutation of the cell axes."""
    Gm = _cell_metric(tuple(true.lattice.abc) + tuple(true.lattice.angles))
    tf = true.frac_coords
    te = np.array([s.symbol for s in true.species])
    sf = solved.frac_coords
    se = np.array([s.symbol for s in solved.species])
    if len(sf) == 0:
        return 0.0
    best = 0.0
    for inv in (1.0, -1.0):
        for j in range(len(tf)):
            for i in range(len(sf)):
                if te[j] != se[i]:
                    continue
                shift = sf[i] - inv * tf[j]
                t2 = (inv * tf + shift) % 1.0
                d = t2[:, None, :] - sf[None, :, :]
                d -= np.round(d)
                dist = np.sqrt(np.einsum("abi,ij,abj->ab", d, Gm, d))
                dist[te[:, None] != se[None, :]] = 99
                best = max(best, float(np.mean(dist.min(axis=1) < tol)))
                if best > 0.999:
                    return best
    return best


def matches(true, solved, ltol=0.2, stol=0.4, angle_tol=5.0):
    """True when the solved structure is the true one: pymatgen's matcher (any cell setting, primitive or conventional, any origin), also for the mirror image
    (powder data cannot tell a structure from its inversion)."""
    from pymatgen.analysis.structure_matcher import StructureMatcher
    from pymatgen.core import Structure
    sm = StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol, primitive_cell=True, scale=True)
    if sm.fit(true, solved):
        return True
    mirror = Structure(solved.lattice, solved.species, -solved.frac_coords, coords_are_cartesian=False)
    return bool(sm.fit(true, mirror))
