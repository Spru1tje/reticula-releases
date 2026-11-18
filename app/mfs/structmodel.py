"""v2.16: symmetry-aware crystal-structure model for refinement of atom positions, atomic displacement parameters and occupancies from powder data.

How it works
  * the space-group operations are found from the structure itself (spglib via pymatgen); every atom of the cell is generated from ONE
    representative of its symmetry class by a stored operation, so any shift of the representative moves all its equivalents consistently;
  * the free coordinates of a representative are restricted to what its SITE SYMMETRY allows (a special position keeps its 0, 1 or 2 degrees
    of freedom: the null space of (R - 1) over the site-symmetry operations), so the symmetry can never be broken by the refinement;
  * reflections are classified into symmetry-equivalent classes (point group + Friedel), so |F|² is computed for ONE reflection per class
    (hundreds, not tens of thousands) - a structure-factor evaluation takes milliseconds;
  * each atom has its own isotropic displacement parameter B (grouped by element by default) and optional occupancy;
  * soft distance restraints from the starting geometry keep molecular fragments (linkers) chemically sensible when the data cannot decide.
Not included: anisotropic displacement, rigid bodies, magnetic scattering, anomalous dispersion."""
from __future__ import annotations

import math

import numpy as np

COVALENT = {}      # filled lazily from pymatgen


def _covalent(el):
    if not COVALENT:
        from pymatgen.core import Element
        for e in Element:
            try:
                COVALENT[e.symbol] = float(e.atomic_radius or 1.5)
            except Exception:
                COVALENT[e.symbol] = 1.5
    return COVALENT.get(el, 1.5)


def _frac_diff(a, b):
    d = a - b
    return d - np.round(d)


class StructureModel:
    def __init__(self, structure, symprec: float = 0.01):
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        self.lattice0 = np.array(structure.lattice.matrix, float)
        sites = list(structure)
        self.n_sites = len(sites)
        self.site_frac = np.array([s.frac_coords for s in sites], float)
        # one "atom" per (site, species): a partially occupied site carries several atoms at the same position
        self.atom_site, self.atom_el, self.atom_occ = [], [], []
        for k, s in enumerate(sites):
            for sp, o in s.species.items():
                self.atom_site.append(k)
                self.atom_el.append(sp.symbol)
                self.atom_occ.append(float(o))
        self.atom_site = np.array(self.atom_site, int)
        self.atom_occ = np.array(self.atom_occ, float)
        self.n_atoms = len(self.atom_site)
        # symmetry from an ordered proxy (dominant species per site)
        proxy = structure.copy()
        try:
            proxy = type(structure)(structure.lattice, [max(s.species.items(), key=lambda kv: kv[1])[0] for s in sites], self.site_frac)
            ds = SpacegroupAnalyzer(proxy, symprec=symprec, angle_tolerance=5.0).get_symmetry_dataset()
            self.rotations = np.array(ds.rotations if hasattr(ds, "rotations") else ds["rotations"], float)
            self.translations = np.array(ds.translations if hasattr(ds, "translations") else ds["translations"], float)
            self.spacegroup = str(getattr(ds, "international", None) or ds["international"])
            self.sg_number = int(getattr(ds, "number", None) or ds["number"])
        except Exception:
            self.rotations = np.eye(3)[None]
            self.translations = np.zeros((1, 3))
            self.spacegroup, self.sg_number = "P1 (symmetry not found)", 1
        self.n_ops = len(self.rotations)
        self._build_orbits(symprec)
        self._build_parameters()

    # ------------------------------------------------------------------------------------------------ orbits and special positions
    def _build_orbits(self, symprec):
        x = self.site_frac
        S = self.n_sites
        self.site_rep = -np.ones(S, int)
        self.site_op = np.zeros((S, 3, 3))
        self.site_t = np.zeros((S, 3))
        self.site_n = np.zeros((S, 3))
        tol = max(3 * symprec, 2e-3)
        for s in range(S):
            if self.site_rep[s] >= 0:
                continue
            # s becomes the representative of every not yet assigned site it maps onto
            for k in range(self.n_ops):
                y = self.rotations[k] @ x[s] + self.translations[k]
                d = np.linalg.norm(_frac_diff(x, y[None, :]) @ self.lattice0, axis=1)
                for j in np.nonzero((d < tol * max(np.linalg.norm(self.lattice0, axis=1)) + 1e-9) & (self.site_rep < 0))[0]:
                    if self._same_species(s, j):
                        self.site_rep[j] = s
                        self.site_op[j] = self.rotations[k]
                        self.site_t[j] = self.translations[k]
                        self.site_n[j] = np.round(x[j] - y)
            if self.site_rep[s] < 0:
                self.site_rep[s] = s
                self.site_op[s] = np.eye(3)
        self.reps = sorted(set(self.site_rep.tolist()))
        # site symmetry of every representative -> allowed displacement directions
        self.basis = {}
        for r in self.reps:
            Rs = []
            for k in range(self.n_ops):
                y = self.rotations[k] @ x[r] + self.translations[k]
                if np.linalg.norm(_frac_diff(x[r][None, :], y[None, :]) @ self.lattice0) < tol * max(np.linalg.norm(self.lattice0, axis=1)) + 1e-9:
                    Rs.append(self.rotations[k])
            M = np.vstack([R - np.eye(3) for R in Rs]) if Rs else np.zeros((1, 3))
            u, sv, vt = np.linalg.svd(M)
            rank = int(np.sum(sv > 1e-6))
            self.basis[r] = vt[rank:].T                      # 3 x d
            self.site_sym_order = getattr(self, "site_sym_order", {})
            self.site_sym_order[r] = len(Rs)

    def _same_species(self, a, b):
        ea = sorted((self.atom_el[i], round(self.atom_occ[i], 4)) for i in np.nonzero(self.atom_site == a)[0])
        eb = sorted((self.atom_el[i], round(self.atom_occ[i], 4)) for i in np.nonzero(self.atom_site == b)[0])
        return ea == eb

    # ------------------------------------------------------------------------------------------------ parameters
    def _build_parameters(self):
        self.pos_params = []              # (rep site, basis column)
        for r in self.reps:
            for c in range(self.basis[r].shape[1]):
                self.pos_params.append((r, c))
        self.n_pos = len(self.pos_params)
        self.elements = sorted(set(self.atom_el))
        self.atom_el_idx = np.array([self.elements.index(e) for e in self.atom_el], int)
        self.atom_rep_site = self.site_rep[self.atom_site]

    def param_names(self):
        return [f"x{r}.{c}" for r, c in self.pos_params]

    def labels(self):
        """Readable label per representative site: element + counter (e.g. Zr1, O3)."""
        cnt, out = {}, {}
        for r in self.reps:
            els = [self.atom_el[i] for i in np.nonzero(self.atom_site == r)[0]]
            e = "/".join(sorted(set(els)))
            cnt[e] = cnt.get(e, 0) + 1
            out[r] = f"{e}{cnt[e]}"
        return out

    # ------------------------------------------------------------------------------------------------ geometry
    def site_positions(self, p: np.ndarray | None = None) -> np.ndarray:
        """Fractional coordinates of all sites for the position parameters p (shifts along the allowed directions of every representative)."""
        xr = self.site_frac.copy()
        if p is not None and len(p):
            for (r, c), v in zip(self.pos_params, p):
                xr[r] = xr[r] + self.basis[r][:, c] * v
        out = np.einsum("sij,sj->si", self.site_op, xr[self.site_rep]) + self.site_t + self.site_n
        return out

    def atom_positions(self, p=None) -> np.ndarray:
        return self.site_positions(p)[self.atom_site]

    # ------------------------------------------------------------------------------------------------ reflection classes
    def classify(self, hkl: np.ndarray):
        """Symmetry-equivalent classes (point group + Friedel) of the reflections hkl (K x 3 ints). Returns (class index per reflection, representative hkl
        per class, number of reflections per class in `hkl`)."""
        rots = np.unique(np.round(self.rotations).astype(int), axis=0)
        h = np.asarray(hkl, int)
        best = None
        for R in rots:
            for sign in (1, -1):
                hp = sign * (h @ R)
                key = ((hp[:, 0] + 128) * 257 + (hp[:, 1] + 128)) * 257 + (hp[:, 2] + 128)
                best = key if best is None else np.minimum(best, key)
        uk, first, inv, cnt = np.unique(best, return_index=True, return_inverse=True, return_counts=True)
        return inv, h[first].astype(float), cnt

    def f2(self, hkl_reps: np.ndarray, s2: np.ndarray, p_pos: np.ndarray | None, B_atoms: np.ndarray, occ: np.ndarray | None = None) -> np.ndarray:
        """|F|² of the class representatives (X-ray, real scattering factors) for positions p_pos, per-atom B (Å²) and occupancies."""
        from pymatgen.analysis.diffraction.xrd import ATOMIC_SCATTERING_PARAMS
        from pymatgen.core import Element
        X = self.atom_positions(p_pos)
        occ = self.atom_occ if occ is None else occ
        fe = np.zeros((len(s2), len(self.elements)))
        for k, e in enumerate(self.elements):
            co = np.asarray(ATOMIC_SCATTERING_PARAMS[e])
            fe[:, k] = Element(e).Z - 41.78214 * s2 * np.sum(co[:, 0] * np.exp(-co[:, 1] * s2[:, None]), axis=1)
        amp = fe[:, self.atom_el_idx] * occ[None, :] * np.exp(-B_atoms[None, :] * s2[:, None])
        ph = np.exp(2j * np.pi * (hkl_reps @ X.T))
        F = np.sum(amp * ph, axis=1)
        return np.abs(F) ** 2

    # ------------------------------------------------------------------------------------------------ restraints
    def bond_list(self, lattice: np.ndarray | None = None, tol: float = 0.45, max_unique: int = 4000):
        """Bonded atom pairs (i, j, integer cell shift n, d0) from the starting geometry, one per symmetry-unique bond (same elements, same length).
        Bonded = distance < sum of covalent radii + tol."""
        L = self.lattice0 if lattice is None else lattice
        X = self.atom_positions()
        pairs, seen = [], set()
        rad = np.array([_covalent(e) for e in self.atom_el])
        for i in range(self.n_atoms):
            d = X - X[i]
            for n in np.array([(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]):
                dd = np.linalg.norm((d + n) @ L, axis=1)
                lim = rad[i] + rad + tol
                for j in np.nonzero((dd < lim) & (dd > 0.4))[0]:
                    if (j, tuple(n)) == (i, (0, 0, 0)):
                        continue
                    a, b = sorted((self.atom_el[i], self.atom_el[j]))
                    key = (a, b, round(float(dd[j]), 3), tuple(sorted((int(self.atom_rep_site[i]), int(self.atom_rep_site[j])))))
                    if key in seen:
                        continue
                    seen.add(key)
                    pairs.append((i, int(j), n.astype(float), float(dd[j])))
                    if len(pairs) >= max_unique:
                        return pairs
        return pairs

    def distances(self, bonds, p_pos, lattice) -> np.ndarray:
        X = self.atom_positions(p_pos)
        out = np.zeros(len(bonds))
        for k, (i, j, n, d0) in enumerate(bonds):
            out[k] = np.linalg.norm(((X[j] + n) - X[i]) @ lattice)
        return out

    def distances_fast(self, bi, bj, bn, p_pos, lattice) -> np.ndarray:
        X = self.atom_positions(p_pos)
        return np.linalg.norm(((X[bj] + bn) - X[bi]) @ lattice, axis=1)

    def distances_metric(self, bi, bj, bn, p_pos, G) -> np.ndarray:
        """Bond lengths for the current cell metric G (3x3, Å²) - the cell is refined together with the atoms."""
        X = self.atom_positions(p_pos)
        d = (X[bj] + bn) - X[bi]
        return np.sqrt(np.maximum(np.einsum("ni,ij,nj->n", d, G, d), 0.0))
