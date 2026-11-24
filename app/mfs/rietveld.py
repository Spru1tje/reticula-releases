"""v2.15: Rietveld and Le Bail refinement of powder patterns - own implementation (numpy / scipy), built to be fast and to need no expert input.

What it does
  * phases from CIF files (any number): structure factors are computed ONCE (X-ray scattering factors as in pymatgen, occupancies, symmetry-
    expanded cell); only an overall Debye-Waller factor is refined, so every refinement step costs milliseconds;
  * peaks: Thompson-Cox-Hastings pseudo-Voigt (GSAS type 3: U, V, W Gaussian; X, Y Lorentzian), optional axial-divergence asymmetry,
    Kα1 / Kα1+Kα2 radiation, zero shift and sample displacement, Lorentz-polarisation (optional monochromator);
  * scale factors and the Chebyshev background are LINEAR parameters: they are solved exactly (non-negative least squares) inside every step
    (variable projection), so only the few non-linear parameters (cell, zero, profile, B) are iterated;
  * AUTOMATIC staged refinement with a start-value search for cell scale and peak width, so that nothing has to be set by the user;
  * Le Bail mode (cell + profile with free intensities, no structure needed);
  * Rwp, Rp, Rexp, GoF, Bragg R, estimated standard deviations (Bérar-corrected), weight fractions (Hill & Howard) and size / strain from X, Y.
Not included (yet, stated in the report): preferred orientation, anomalous-dispersion corrections, anisotropic displacement, atom-position refinement."""
from __future__ import annotations

import io
import math
import os
import tempfile
import time

import numpy as np
from scipy.optimize import least_squares, nnls

RADIATION = {
    "Cu Kα1": [(1.540593, 1.0)],
    "Cu Kα1 + Kα2": [(1.540593, 1.0), (1.544414, 0.5)],
    "Co Kα1": [(1.788965, 1.0)],
    "Co Kα1 + Kα2": [(1.788965, 1.0), (1.792850, 0.5)],
    "Mo Kα1": [(0.709319, 1.0)],
    "Cr Kα1": [(2.289700, 1.0)],
}
LN2 = math.log(2.0)


# ---------------------------------------------------------------------------------------------------------------- crystal structure
def _cell_metric(c):
    a, b, cc, al, be, ga = c
    al, be, ga = (math.radians(v) for v in (al, be, ga))
    return np.array([[a * a, a * b * math.cos(ga), a * cc * math.cos(be)],
                     [a * b * math.cos(ga), b * b, b * cc * math.cos(al)],
                     [a * cc * math.cos(be), b * cc * math.cos(al), cc * cc]])


def cell_volume(c):
    return float(math.sqrt(max(np.linalg.det(_cell_metric(c)), 0.0)))


def detect_system(c, tol_len=2e-3, tol_ang=0.05):
    """Crystal system from the lattice metric alone: 'cubic', 'hexagonal', 'tetragonal', 'orthorhombic', 'monoclinic', 'rhombohedral', 'triclinic'."""
    a, b, cc, al, be, ga = c
    eq = lambda u, v: abs(u - v) <= tol_len * max(u, v)
    r90 = lambda t: abs(t - 90.0) <= tol_ang
    if r90(al) and r90(be) and r90(ga):
        if eq(a, b) and eq(b, cc):
            return "cubic"
        if eq(a, b):
            return "tetragonal"
        return "orthorhombic"
    if r90(al) and r90(be) and abs(ga - 120.0) <= tol_ang and eq(a, b):
        return "hexagonal"
    if eq(a, b) and eq(b, cc) and abs(al - be) <= tol_ang and abs(be - ga) <= tol_ang:
        return "rhombohedral"
    if r90(al) and r90(ga):
        return "monoclinic"
    return "triclinic"


# independent cell parameters per system, and how they map to (a, b, c, alpha, beta, gamma)
CELL_DOF = {"cubic": ["a"], "hexagonal": ["a", "c"], "tetragonal": ["a", "c"], "orthorhombic": ["a", "b", "c"], "rhombohedral": ["a", "alpha"],
            "monoclinic": ["a", "b", "c", "beta"], "triclinic": ["a", "b", "c", "alpha", "beta", "gamma"]}


def cell_from_dof(system, p, base):
    a0, b0, c0, al0, be0, ga0 = base
    p = dict(p)
    if system == "cubic":
        return (p["a"], p["a"], p["a"], 90.0, 90.0, 90.0)
    if system == "hexagonal":
        return (p["a"], p["a"], p["c"], 90.0, 90.0, 120.0)
    if system == "tetragonal":
        return (p["a"], p["a"], p["c"], 90.0, 90.0, 90.0)
    if system == "orthorhombic":
        return (p["a"], p["b"], p["c"], 90.0, 90.0, 90.0)
    if system == "rhombohedral":
        return (p["a"], p["a"], p["a"], p["alpha"], p["alpha"], p["alpha"])
    if system == "monoclinic":
        return (p["a"], p["b"], p["c"], 90.0, p["beta"], 90.0)
    return (p["a"], p["b"], p["c"], p["alpha"], p["beta"], p["gamma"])


def dof_from_cell(system, c):
    names = CELL_DOF[system]
    idx = {"a": 0, "b": 1, "c": 2, "alpha": 3, "beta": 4, "gamma": 5}
    return {n: float(c[idx[n]]) for n in names}


_POOL = None


def _pool_map():
    """v2.20: the finite-difference Jacobian columns are evaluated in parallel threads (same numbers, same steps; the residual only reads its inputs)."""
    global _POOL
    if _POOL is None:
        from concurrent.futures import ThreadPoolExecutor
        import os as _os
        _POOL = ThreadPoolExecutor(max_workers=max(1, min(4, (_os.cpu_count() or 2))))
    return _POOL.map


class Phase:
    """One crystalline phase: cell, symmetry-expanded atoms and the reflection list (grouped by d-spacing) with |F|² at B = 0."""

    def __init__(self, name, cell, system, atoms_frac, occ, elems, mass_cell, formula=""):
        self.name = name
        self.cell0 = tuple(float(v) for v in cell)
        self.system = system
        self.frac = np.asarray(atoms_frac, float)
        self.occ = np.asarray(occ, float)
        self.elems = list(elems)               # element symbol per atom
        self.mass_cell = float(mass_cell)      # g/mol of the cell content (Z·M)
        self.formula = formula
        self.hkl = None
        self.mult = None
        self.F2 = None
        self.po_axis = None            # (h, k, l) of the preferred-orientation axis (March-Dollase), None = random powder
        self.model = None              # StructureModel (symmetry-aware atoms) when the structure is refined
        self.struct_refine = False
        self._structure = None
        self.bond_arr = None

    @classmethod
    def from_cif(cls, src, name=None, substitute=None):
        import warnings
        warnings.filterwarnings("ignore")
        from pymatgen.core import Structure
        if isinstance(src, (bytes, bytearray)):
            with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as t:
                t.write(src)
                path = t.name
        else:
            path = str(src)
        try:
            s = Structure.from_file(path)
        finally:
            if isinstance(src, (bytes, bytearray)):
                try:
                    os.unlink(path)
                except OSError:
                    pass
        if substitute:
            s.replace_species(substitute)
        L = s.lattice
        cell = (L.a, L.b, L.c, L.alpha, L.beta, L.gamma)
        frac, occ, el = [], [], []
        mass = 0.0
        for site in s:
            for sp, o in site.species.items():
                frac.append(site.frac_coords)
                occ.append(float(o))
                el.append(sp.symbol)
                mass += float(o) * float(sp.atomic_mass)
        obj = cls(name or s.composition.reduced_formula, cell, detect_system(cell), frac, occ, el, mass, s.composition.reduced_formula)
        obj._structure = s
        return obj

    @classmethod
    def from_cell(cls, cell, name="unknown phase"):
        """A phase without a structure (e.g. from indexing): every reflection of the cell with intensity 1 - for Le Bail refinement only."""
        obj = cls(name, cell, detect_system(cell), np.zeros((0, 3)), [], [], 0.0, "")
        obj.cell_only = True
        return obj

    def enable_structure(self, symprec: float = 0.01, restraint_tol: float = 0.45):
        """Switches the phase to structure refinement (atom positions, displacement parameters): builds the symmetry-aware model and the bond restraints."""
        from mfs.structmodel import StructureModel
        if self._structure is None:
            raise ValueError("structure refinement needs the phase to be read from a CIF")
        self.model = StructureModel(self._structure, symprec)
        self.struct_refine = True
        bonds = self.model.bond_list(tol=restraint_tol)
        self.bond_arr = (np.array([b[0] for b in bonds], int), np.array([b[1] for b in bonds], int), np.array([b[2] for b in bonds], float),
                         np.array([b[3] for b in bonds], float)) if bonds else None
        if self.hkl is not None:
            self._classify()
        return self

    def _classify(self):
        """Reflection classes for the structure-factor evaluation (one |F|² per symmetry-equivalent class)."""
        inv, reps, cnt = self.model.classify(self.m_hkl)
        self.m_class = inv
        self.cls_hkl = reps
        self.cls_cnt = cnt.astype(float)
        self.cls_group = np.zeros(len(reps), int)
        self.cls_group[inv] = self.m_grp

    def build(self, dmin: float):
        if getattr(self, "cell_only", False):
            return self._build_cell_only(dmin)
        """Reflections down to d = dmin (Å): half-space hkl, structure factors at B = 0, equivalents merged by d-spacing."""
        from pymatgen.analysis.diffraction.xrd import ATOMIC_SCATTERING_PARAMS
        from pymatgen.core import Element
        Gs = np.linalg.inv(_cell_metric(self.cell0))
        a, b, c = self.cell0[:3]
        hm = [int(math.ceil(v / dmin)) for v in (a, b, c)]
        h, k, l = np.meshgrid(*[np.arange(-m, m + 1) for m in hm], indexing="ij")
        hk = np.stack([h.ravel(), k.ravel(), l.ravel()], 1)
        half = (hk[:, 0] > 0) | ((hk[:, 0] == 0) & (hk[:, 1] > 0)) | ((hk[:, 0] == 0) & (hk[:, 1] == 0) & (hk[:, 2] > 0))
        hk = hk[half]
        inv_d2 = np.einsum("ki,ij,kj->k", hk, Gs, hk)
        hk, inv_d2 = hk[inv_d2 <= 1.0 / dmin ** 2], inv_d2[inv_d2 <= 1.0 / dmin ** 2]
        s2 = inv_d2 / 4.0
        F = np.zeros(len(hk), complex)
        for sym in sorted(set(self.elems)):
            sel = np.array([e == sym for e in self.elems])
            co = np.asarray(ATOMIC_SCATTERING_PARAMS[sym])
            Z = Element(sym).Z
            fs = Z - 41.78214 * s2 * np.sum(co[:, 0] * np.exp(-co[:, 1] * s2[:, None]), axis=1)
            fr, oc = self.frac[sel].T, self.occ[sel]
            step = max(1, 4_000_000 // max(1, fr.shape[1]))
            for i0 in range(0, len(hk), step):
                F[i0:i0 + step] += fs[i0:i0 + step] * (np.exp(2j * np.pi * (hk[i0:i0 + step] @ fr)) @ oc)
        F2 = np.abs(F) ** 2
        key = np.round(inv_d2 * 1e6).astype(np.int64)
        uk, first, inv = np.unique(key, return_index=True, return_inverse=True)
        I0 = np.bincount(inv, weights=2.0 * F2)
        keep = I0 > 1e-7 * I0.max()
        newidx = -np.ones(len(I0), int)
        newidx[keep] = np.arange(int(keep.sum()))
        mg = newidx[inv]
        ok = mg >= 0
        self.m_hkl, self.m_F2, self.m_grp = hk[ok].astype(float), 2.0 * F2[ok], mg[ok]      # every single reflection and its group
        self.hkl = hk[first][keep]
        self.I0 = I0[keep]                      # multiplicity · |F|² (Friedel pairs counted)
        self.dmin = dmin
        if self.model is not None:
            self._classify()
        return self

    def _build_cell_only(self, dmin: float):
        Gs = np.linalg.inv(_cell_metric(self.cell0))
        hm = [int(math.ceil(v / dmin)) for v in self.cell0[:3]]
        h, k, l = np.meshgrid(*[np.arange(-m, m + 1) for m in hm], indexing="ij")
        hk = np.stack([h.ravel(), k.ravel(), l.ravel()], 1)
        half = (hk[:, 0] > 0) | ((hk[:, 0] == 0) & (hk[:, 1] > 0)) | ((hk[:, 0] == 0) & (hk[:, 1] == 0) & (hk[:, 2] > 0))
        hk = hk[half]
        inv_d2 = np.einsum("ki,ij,kj->k", hk, Gs, hk)
        sel = inv_d2 <= 1.0 / dmin ** 2
        hk, inv_d2 = hk[sel], inv_d2[sel]
        key = np.round(inv_d2 * 1e6).astype(np.int64)
        uk, first, inv = np.unique(key, return_index=True, return_inverse=True)
        self.hkl = hk[first]
        self.I0 = np.bincount(inv, weights=np.full(len(hk), 2.0))
        self.m_hkl, self.m_F2, self.m_grp = hk.astype(float), np.full(len(hk), 2.0), inv
        self.dmin = dmin
        return self

    @property
    def n_ref(self):
        return 0 if self.hkl is None else len(self.hkl)


# ---------------------------------------------------------------------------------------------------------------- parameters
class _P:
    __slots__ = ("name", "value", "lo", "hi", "free", "scale", "label")

    def __init__(self, name, value, lo, hi, free=False, scale=0.01, label=""):
        self.name, self.value, self.lo, self.hi, self.free, self.scale, self.label = name, float(value), float(lo), float(hi), bool(free), float(scale), label or name


def _tch(Hg2, Hl):
    """Thompson-Cox-Hastings: total FWHM and mixing eta from the Gaussian (squared) and Lorentzian FWHM."""
    Hg = np.sqrt(np.maximum(Hg2, 1e-12))
    H5 = Hg ** 5 + 2.69269 * Hg ** 4 * Hl + 2.42843 * Hg ** 3 * Hl ** 2 + 4.47163 * Hg ** 2 * Hl ** 3 + 0.07842 * Hg * Hl ** 4 + Hl ** 5
    H = H5 ** 0.2
    q = np.clip(Hl / np.maximum(H, 1e-9), 0, 1)
    eta = np.clip(1.36603 * q - 0.47719 * q ** 2 + 0.11116 * q ** 3, 0, 1)
    return H, eta


def _pv(dx, H, eta):
    """Area-normalised pseudo-Voigt (per degree) at offsets dx (…×n) for widths H and mixing eta (…×1).
    v2.21: evaluated in float32 (about 1.9x faster on a full refinement; the refined positions move by < 3e-5 and Rwp by < 3e-6
    against float64). Returns float64."""
    z = (dx / H).astype(np.float32)
    Hf = np.asarray(H, np.float32)
    ef = np.asarray(eta, np.float32)
    G = (np.float32(2.0 * math.sqrt(LN2 / math.pi)) / Hf) * np.exp(np.float32(-4.0 * LN2) * z * z)
    L = (np.float32(2.0 / math.pi) / Hf) / (np.float32(1.0) + np.float32(4.0) * z * z)
    return (ef * L + (np.float32(1.0) - ef) * G).astype(np.float64)


class Refinement:
    """Rietveld / Le Bail refinement of one pattern with one or more phases."""

    def __init__(self, x, y, phases, radiation="Cu Kα1", bg_order=6, mode="rietveld", mono_2theta=0.0, instrument=None, counts=None, progress=None, offset="displacement", glass_bg=True, bg_pattern=None):
        self.x = np.asarray(x, float)
        y = np.asarray(y, float)
        if counts is None:
            counts = bool(np.mean(np.abs(y - np.round(y)) < 1e-6) > 0.9 and y.max() >= 100)
        self.counts = counts
        self.scale_y = 1.0
        if not counts:                           # arbitrary units: rescale to pseudo-counts (weights only matter relatively)
            self.scale_y = 1e4 / max(float(np.max(y)), 1e-12)
        self.y = y * self.scale_y
        self.w = 1.0 / np.maximum(self.y, max(1.0, 1e-3 * float(self.y.max())))
        self.phases = phases
        self.rad = RADIATION[radiation] if isinstance(radiation, str) else list(radiation)
        self.radiation = radiation if isinstance(radiation, str) else "custom"
        self.mode = mode
        self.nb = int(bg_order) + 1
        self.mono = float(mono_2theta)
        # Zero shift and sample displacement are strongly correlated; refining both biases the cell (measured on a real laboratory pattern:
        # -0.06 %). Default = displacement only (reflection geometry); 'zero' for transmission / calibrated instruments; 'both' needs a standard.
        self.offset_mode = offset if offset in ("displacement", "zero", "both") else "displacement"
        self.progress = progress or (lambda frac, text: None)
        self.instrument = instrument
        xr = (self.x - self.x.min()) / max(self.x.max() - self.x.min(), 1e-9) * 2 - 1
        self.nb_poly = self.nb
        T = np.polynomial.chebyshev.chebvander(xr, self.nb - 1)               # (N, nb)
        # v2.18.2: an amorphous hump of a glass sample holder (broad, 15-35 deg 2theta) is modelled by three broad Gaussians next to the
        # polynomial (switch off with glass_bg=False); a measured background pattern (empty holder) can be added as one more column.
        extra = []
        if glass_bg:
            for c_ in (18.0, 24.0, 30.0):
                if self.x.min() - 6 < c_ < self.x.max() + 6 and (self.x.max() - self.x.min()) > 15:
                    extra.append(np.exp(-0.5 * ((self.x - c_) / 4.0) ** 2))
        self.bg_extra = len(extra)
        self.bg_measured = False
        if bg_pattern is not None:
            bx, by = np.asarray(bg_pattern[0], float), np.asarray(bg_pattern[1], float)
            o = np.argsort(bx)
            bgi = np.interp(self.x, bx[o], by[o], left=by[o][0], right=by[o][-1]) * self.scale_y
            if np.ptp(bgi) > 0:
                extra.append(bgi / max(float(np.max(np.abs(bgi))), 1e-12) * float(np.max(self.y)) * 0.1)
                self.bg_measured = True
        if extra:
            T = np.hstack([T, np.stack(extra, 1)])
            self.nb = T.shape[1]
        self.T = T
        self.dx = np.gradient(self.x)
        tt_max = float(self.x.max()) + 2.0
        lam = min(l for l, _ in self.rad)
        dmin = lam / (2 * math.sin(math.radians(min(tt_max, 175.0) / 2)))
        for p in phases:
            if p.hkl is None or getattr(p, "dmin", 99) > dmin * 1.001:
                p.build(dmin * 0.995)
        self._init_params()
        self.I_free = {i: np.array(p.I0 / max(p.I0.max(), 1e-12), float) for i, p in enumerate(phases)}     # Le Bail intensities (start = model)
        self.info: list[str] = []

    # ------------------------------------------------------------------------------------------- parameter table
    def _init_params(self):
        P = {}
        P["zero"] = _P("zero", 0.0, -0.6, 0.6, False, 0.01, "Zero shift (°2θ)")
        P["disp"] = _P("disp", 0.0, -0.6, 0.6, False, 0.01, "Sample displacement (°2θ)")
        P["U"] = _P("U", 0.0, 0.0, 3.0, False, 0.001, "U (°²)")
        P["V"] = _P("V", 0.0, -1.0, 1.0, False, 0.001, "V (°²)")
        P["W"] = _P("W", 0.004, 0.0, 3.0, False, 0.001, "W (°²)")
        P["X"] = _P("X", 0.05, 0.0, 3.0, False, 0.01, "X Lorentzian (°)")
        P["Y"] = _P("Y", 0.0, 0.0, 3.0, False, 0.01, "Y Lorentzian (°)")
        P["asym"] = _P("asym", 0.0, 0.0, 1.0, False, 0.01, "Axial asymmetry")
        for i, ph in enumerate(self.phases):
            base = dof_from_cell(ph.system, ph.cell0)
            for n, v in base.items():
                lo, hi = (v * 0.97, v * 1.03) if n != "alpha" and n != "beta" and n != "gamma" else (v - 3.0, v + 3.0)
                P[f"c{i}.{n}"] = _P(f"c{i}.{n}", v, lo, hi, False, 0.002 * (1 if n in "abc" else 50), f"{ph.name}: {n}")
            P[f"B{i}"] = _P(f"B{i}", 1.0, 0.0, 12.0, False, 0.05, f"{ph.name}: B overall (Å²)")
            if ph.model is not None and ph.struct_refine:
                P[f"B{i}"].value = 0.0                                          # replaced by one B per element (sB)
            P[f"po{i}"] = _P(f"po{i}", 1.0, 0.3, 1.8, False, 0.02, f"{ph.name}: preferred orientation r (March-Dollase)")
            if ph.model is not None and ph.struct_refine:
                lim = 0.35 / max(ph.cell0[:3])                                 # a shift of at most 0.35 Å (fractional units)
                labs = ph.model.labels()
                for k, (r_, c_) in enumerate(ph.model.pos_params):
                    P[f"sx{i}.{k}"] = _P(f"sx{i}.{k}", 0.0, -lim, lim, False, max(lim / 20, 5e-4), f"{ph.name}: {labs[r_]} shift {c_ + 1}")
                for el in ph.model.elements:
                    P[f"sB{i}.{el}"] = _P(f"sB{i}.{el}", 1.0, 0.05, 15.0, False, 0.1, f"{ph.name}: B({el}) (Å²)")
        # v2.18.2: an instrument profile measured on a standard (e.g. Si SRM 640, LaB6 SRM 660) is FIXED: U, V, W, X, Y keep these values
        self.fixed_profile = bool(self.instrument and self.instrument.get("fixed"))
        self.fix_cell = False                                           # v2.18.2: a standard's cell is known (certificate): not refined when True
        if self.fixed_profile:
            for k in ("U", "V", "W", "X", "Y"):
                if k in self.instrument:
                    P[k].value = float(self.instrument[k])
        self.P = P

    def free_names(self):
        return [n for n, p in self.P.items() if p.free]

    def set_free(self, names, on=True):
        for n in names:
            if n in self.P:
                self.P[n].free = on

    # ------------------------------------------------------------------------------------------- model
    def cell(self, i, vals=None):
        ph = self.phases[i]
        g = (lambda n: vals[f"c{i}.{n}"]) if vals else (lambda n: self.P[f"c{i}.{n}"].value)
        return cell_from_dof(ph.system, {n: g(n) for n in CELL_DOF[ph.system]}, ph.cell0)

    def _struct_mode(self, i):
        ph = self.phases[i]
        return ph.model is not None and ph.struct_refine and getattr(ph, "cls_hkl", None) is not None

    def _pos_vec(self, i, v):
        ph = self.phases[i]
        return np.array([v[f"sx{i}.{k}"] for k in range(ph.model.n_pos)])

    def _B_atoms(self, i, v):
        ph = self.phases[i]
        return np.array([v[f"sB{i}.{e}"] for e in ph.model.atom_el])

    def _class_F2(self, i, v):
        """|F|² of the reflection classes. v2.20: cached on the exact cell, positions and B values (profile-only steps of the
        Jacobian reuse it; the result is bit-identical because the key is the exact inputs)."""
        ph = self.phases[i]
        cell = self.cell(i, v)
        pos = self._pos_vec(i, v)
        Bat = self._B_atoms(i, v)
        key = (tuple(np.asarray(cell, float).ravel()), tuple(np.asarray(pos, float).ravel()), tuple(np.asarray(Bat, float).ravel()))
        cache = self.__dict__.setdefault("_f2cache", {})
        hit = cache.get(i)
        if hit is not None and hit[0] == key:
            return hit[1]
        Gs = np.linalg.inv(_cell_metric(cell))
        s2 = np.einsum("ki,ij,kj->k", ph.cls_hkl, Gs, ph.cls_hkl) / 4.0
        F2 = ph.model.f2(ph.cls_hkl, s2, pos, Bat)
        cache[i] = (key, F2)
        return F2

    def _dw(self, i, v, s2):
        """Overall Debye-Waller factor of the intensities; with the structure model every atom carries its own B (already inside |F|²)."""
        if self._struct_mode(i):
            return 1.0
        return np.exp(-2.0 * v[f"B{i}"] * s2)

    def _I0(self, i, v):
        """multiplicity · |F|² of the reflection groups of phase i; with a preferred-orientation axis every single reflection is weighted by the
        March-Dollase function (r² cos²α + sin²α / r)^(-3/2), α = angle between the reflection and the axis (r < 1 plate-like, r > 1 needle-like)."""
        ph = self.phases[i]
        r = v[f"po{i}"]
        struct = self._struct_mode(i)
        if struct:
            F2c = self._class_F2(i, v)
            if ph.po_axis is None:
                return np.bincount(ph.cls_group, weights=2.0 * ph.cls_cnt * F2c, minlength=ph.n_ref)
        if ph.po_axis is None or ph.m_hkl is None:
            return ph.I0
        Gs = np.linalg.inv(_cell_metric(self.cell(i, v)))
        h0 = np.asarray(ph.po_axis, float)
        dot = ph.m_hkl @ Gs @ h0
        n1 = np.einsum("ki,ij,kj->k", ph.m_hkl, Gs, ph.m_hkl)
        n2 = float(h0 @ Gs @ h0)
        cos2 = np.clip(dot ** 2 / np.maximum(n1 * n2, 1e-30), 0.0, 1.0)
        w = (r * r * cos2 + (1.0 - cos2) / r) ** -1.5
        return np.bincount(ph.m_grp, weights=(2.0 * F2c[ph.m_class] if struct else ph.m_F2) * w, minlength=ph.n_ref)

    def _values(self, theta=None):
        v = {n: p.value for n, p in self.P.items()}
        if theta is not None:
            for n, t in zip(self.free_names(), theta):
                v[n] = float(t)
        return v

    def _phase_pattern(self, i, v, want_detail=False):
        """Pattern of phase i with unit scale (all wavelengths) on the data grid; also peak positions and intensities."""
        ph = self.phases[i]
        cell = self.cell(i, v)
        Gs = np.linalg.inv(_cell_metric(cell))
        inv_d2 = np.einsum("ki,ij,kj->k", ph.hkl, Gs, ph.hkl)
        d = 1.0 / np.sqrt(inv_d2)
        s2 = inv_d2 / 4.0
        if self.mode == "lebail":
            I = self.I_free[i].copy()
        else:
            I = self._I0(i, v) * self._dw(i, v, s2)
        out = np.zeros_like(self.x)
        pk_pos, pk_I = [], []
        for lam, wl in self.rad:
            sin_t = lam / (2.0 * d)
            ok = sin_t < 0.9995
            th = np.arcsin(np.clip(sin_t, 0, 0.9995))
            tt = 2.0 * np.degrees(th) + v["zero"] + v["disp"] * np.cos(th)
            c2m = math.cos(math.radians(self.mono)) ** 2
            lp = (1.0 + c2m * np.cos(2 * th) ** 2) / ((1.0 + c2m) * sin_t ** 2 * np.cos(th))
            tan = np.tan(th)
            Hg2 = v["U"] * tan ** 2 + v["V"] * tan + v["W"]
            Hl = v["X"] / np.cos(th) + v["Y"] * tan
            H, eta = _tch(Hg2, Hl)
            Ik = (I * lp * wl if self.mode != "lebail" else I * wl)
            sel = ok & (tt > self.x[0] - 5 * H) & (tt < self.x[-1] + 5 * H)
            if not sel.any():
                continue
            tt_, H_, eta_, Ik_, th_ = tt[sel], H[sel], eta[sel], Ik[sel], th[sel]
            cut = 25.0 * H_
            i0 = np.searchsorted(self.x, tt_ - cut)
            i1 = np.searchsorted(self.x, tt_ + cut)
            nw = int(min(max(int((i1 - i0).max()), 3), 1500))
            idx = i0[:, None] + np.arange(nw)[None, :]
            valid = (idx < len(self.x)) & (idx < i1[:, None])
            idc = np.clip(idx, 0, len(self.x) - 1)
            dxm = self.x[idc] - tt_[:, None]
            if v["asym"] > 0:
                f = 1.0 + v["asym"] / np.maximum(np.tan(th_), 0.05)[:, None] * np.where(dxm < 0, 1.0, -1.0) * 0.5
                dxm = dxm / np.maximum(f, 0.2)
                norm = 1.0
            else:
                norm = 1.0
            prof = _pv(dxm, H_[:, None], eta_[:, None]) * norm
            prof = np.where(valid, prof, 0.0)
            contrib = Ik_[:, None] * prof
            out += np.bincount(idc.ravel(), weights=contrib.ravel(), minlength=len(self.x))
            if want_detail:
                pk_pos.append(tt_)
                pk_I.append(Ik_)
        if want_detail:
            return out, (np.concatenate(pk_pos) if pk_pos else np.array([]), np.concatenate(pk_I) if pk_I else np.array([]))
        return out

    def _design(self, v):
        cols = [self._phase_pattern(i, v) for i in range(len(self.phases))]
        return np.stack(cols, 1) if cols else np.zeros((len(self.x), 0))

    def _solve_linear(self, v, A=None):
        """Scale factors (>= 0) and Chebyshev background (free sign) by weighted non-negative least squares."""
        A = self._design(v) if A is None else A
        sw = np.sqrt(self.w)
        M = np.hstack([A, self.T, -self.T])
        c, _ = nnls(M * sw[:, None], self.y * sw, maxiter=2000)
        n = A.shape[1]
        s = c[:n]
        bg_c = c[n:n + self.nb] - c[n + self.nb:]
        calc = A @ s + self.T @ bg_c
        return s, bg_c, calc, A

    restraint_sigma = 0.0          # Å; > 0 switches the distance restraints on (set by the structure stages)
    restraint_scale = 1.0          # multiplies the restraint widths: < 1 tighter (trust the starting geometry), > 1 looser

    def _restraints(self, v):
        """Soft distance restraints (d - d0) / sigma of the bonded pairs of the starting structure; d0 follows the refined cell (isotropic scale)."""
        out = []
        if self.restraint_sigma <= 0:
            return np.zeros(0)
        for i, ph in enumerate(self.phases):
            if not self._struct_mode(i) or ph.bond_arr is None:
                continue
            bi, bj, bn, d0 = ph.bond_arr
            cell = self.cell(i, v)
            sc = float(np.mean(cell[:3]) / np.mean(ph.cell0[:3]))
            d = ph.model.distances_metric(bi, bj, bn, self._pos_vec(i, v), _cell_metric(cell))
            out.append((d - d0 * sc) / self.restraint_sigma)
        return np.concatenate(out) if out else np.zeros(0)

    def _data_residual(self, v):
        s, bgc, calc, _ = self._solve_linear(v)
        return np.sqrt(self.w) * (self.y - calc)

    def _residual(self, theta):
        v = self._values(theta)
        r = self._data_residual(v)
        rr = self._restraints(v)
        return np.concatenate([r, rr]) if len(rr) else r

    # ------------------------------------------------------------------------------------------- statistics
    def stats(self, calc=None):
        if calc is None:
            v = self._values()
            _, _, calc, _ = self._solve_linear(v)
        r = self.y - calc
        n_par = len(self.free_names()) + len(self.phases) + self.nb
        if self.mode == "lebail":
            n_par += sum(len(a) for a in self.I_free.values())
        N = len(self.x)
        rwp = math.sqrt(float(np.sum(self.w * r ** 2) / np.sum(self.w * self.y ** 2)))
        rp = float(np.sum(np.abs(r)) / np.sum(np.abs(self.y)))
        rexp = math.sqrt(max(N - n_par, 1) / float(np.sum(self.w * self.y ** 2)))
        return {"Rwp": rwp, "Rp": rp, "Rexp": rexp, "chi2": (rwp / rexp) ** 2, "n_par": n_par, "N": N}

    def _fit(self, names, max_nfev=40):
        """Least squares on the named non-linear parameters (others fixed). Keeps the result only if Rwp improves."""
        before = {n: self.P[n].value for n in self.P}
        r0 = float(np.sum(self._residual(np.array([self.P[n].value for n in self.free_names()])) ** 2)) if self.free_names() else None
        self.set_free(list(self.P), False)
        self.set_free(names, True)
        fn = self.free_names()
        if not fn:
            return None
        th0 = np.array([self.P[n].value for n in fn])
        lo = np.array([self.P[n].lo for n in fn])
        hi = np.array([self.P[n].hi for n in fn])
        th0 = np.clip(th0, lo + 1e-12, hi - 1e-12)
        sc = np.array([self.P[n].scale for n in fn])
        f0 = self._residual(th0)
        # v2.18.6: multi-start. The soft bond restraints make the optimum depend on the starting point (small numerical differences between platforms
        # gave different geometries). Two deterministic perturbed starts are tried besides the given one; the lowest cost wins.
        rng_ms = np.random.default_rng(12345)
        starts = [th0] + [np.clip(th0 + rng_ms.normal(0.0, 1.0, th0.size) * sc * 0.5, lo + 1e-12, hi - 1e-12) for _ in range(2)]
        sol = None
        for th_s in starts:
            s_ = least_squares(self._residual, th_s, bounds=(lo, hi), x_scale=sc, diff_step=1e-5, max_nfev=max_nfev, ftol=1e-7, xtol=1e-9, workers=_pool_map())
            if sol is None or float(np.sum(s_.fun ** 2)) < float(np.sum(sol.fun ** 2)):
                sol = s_
        if float(np.sum(sol.fun ** 2)) < float(np.sum(f0 ** 2)):
            for n, t in zip(fn, sol.x):
                self.P[n].value = float(t)
            self._last_jac = (fn, sol.jac, sol.fun)
        else:
            for n in fn:
                self.P[n].value = before[n]
        return sol

    # ------------------------------------------------------------------------------------------- automatic strategy
    def _prealign(self):
        """Global search over an isotropic cell scale and the zero shift with deliberately broad peaks: removes the start-value problem."""
        saved = {n: p.value for n, p in self.P.items()}
        for n, val in (("U", 0.0), ("V", 0.0), ("W", 0.08), ("X", 0.1), ("Y", 0.0)):
            self.P[n].value = val
        for i, ph in enumerate(self.phases):
            base = dof_from_cell(ph.system, ph.cell0)
            best = (1e99, 1.0)
            for f in np.linspace(0.975, 1.025, 51):
                for n, v0 in base.items():
                    if n in ("a", "b", "c"):
                        self.P[f"c{i}.{n}"].value = v0 * f
                v = self._values()
                s, bg, calc, _ = self._solve_linear(v)
                rr = float(np.sum(self.w * (self.y - calc) ** 2))
                if rr < best[0]:
                    best = (rr, f)
            for n, v0 in base.items():
                if n in ("a", "b", "c"):
                    self.P[f"c{i}.{n}"].value = v0 * best[1]
        for n in ("U", "V", "W", "X", "Y"):
            self.P[n].value = saved[n]
        return

    def _width_scan(self):
        """Start value for the peak width: smallest Rwp over a log scale of Gaussian / Lorentzian widths."""
        best = (1e99, None)
        for f in np.geomspace(0.02, 1.0, 14):
            for lor in (0.3, 0.7):
                self.P["W"].value = (f * (1 - lor)) ** 2
                self.P["X"].value = f * lor
                self.P["Y"].value = 0.0
                s, bg, calc, _ = self._solve_linear(self._values())
                rr = float(np.sum(self.w * (self.y - calc) ** 2))
                if rr < best[0]:
                    best = (rr, (self.P["W"].value, self.P["X"].value))
        self.P["W"].value, self.P["X"].value = best[1]

    def _po_scan(self):
        """The March-Dollase objective has two basins (plates r < 1 and needles r > 1): scan r on a grid for every phase with a texture axis
        and start the refinement in the better one."""
        for i, ph in enumerate(self.phases):
            if ph.po_axis is None:
                continue
            best = (1e99, 1.0)
            for r in np.concatenate([np.linspace(0.35, 0.95, 13), np.linspace(1.05, 1.75, 8), [1.0]]):
                self.P[f"po{i}"].value = float(r)
                s_, bg_, calc, _ = self._solve_linear(self._values())
                rr = float(np.sum(self.w * (self.y - calc) ** 2))
                if rr < best[0]:
                    best = (rr, float(r))
            self.P[f"po{i}"].value = best[1]

    def _le_bail_extract(self, cycles=40):
        v = self._values()
        for i in range(len(self.phases)):
            pass
        sw = self.w
        for _ in range(cycles):
            A = self._design(v)
            s, bgc, calc, _ = self._solve_linear(v, A)
            bg = self.T @ bgc
            # partition of the observed intensity over the reflections, in proportion to the current model
            for i, ph in enumerate(self.phases):
                self._lebail_update(i, v, calc, bg, float(s[i]))
        return

    def _lebail_update(self, i, v, calc, bg, scale=1.0):
        """One Le Bail step for phase i: I_k ← Σ_j (I_k φ_kj / y_calc,j) (y_obs,j − bg_j) Δx_j."""
        ph = self.phases[i]
        cell = self.cell(i, v)
        Gs = np.linalg.inv(_cell_metric(cell))
        inv_d2 = np.einsum("ki,ij,kj->k", ph.hkl, Gs, ph.hkl)
        d = 1.0 / np.sqrt(inv_d2)
        ycp = calc - bg
        ratio = np.where(ycp > 1e-9, np.clip(self.y - bg, 0, None) / np.maximum(ycp, 1e-9), 0.0)
        newI = np.zeros_like(self.I_free[i])
        for lam, wl in self.rad:
            sin_t = lam / (2.0 * d)
            th = np.arcsin(np.clip(sin_t, 0, 0.9995))
            tt = 2.0 * np.degrees(th) + v["zero"] + v["disp"] * np.cos(th)
            tan = np.tan(th)
            H, eta = _tch(v["U"] * tan ** 2 + v["V"] * tan + v["W"], v["X"] / np.cos(th) + v["Y"] * tan)
            cut = 25.0 * H
            i0 = np.searchsorted(self.x, tt - cut)
            i1 = np.searchsorted(self.x, tt + cut)
            nw = int(min(max(int((i1 - i0).max()), 3), 1500))
            idx = i0[:, None] + np.arange(nw)[None, :]
            valid = (idx < len(self.x)) & (idx < i1[:, None])
            idc = np.clip(idx, 0, len(self.x) - 1)
            prof = np.where(valid, _pv(self.x[idc] - tt[:, None], H[:, None], eta[:, None]), 0.0)
            share = (prof * (ratio[idc] * self.dx[idc])).sum(1)
            newI += wl * self.I_free[i] * share * scale / sum(w for _, w in self.rad)
        ok = np.isfinite(newI)
        self.I_free[i] = np.where(ok, np.maximum(newI, 0.0), self.I_free[i])

    def run(self, cycles=None):
        """Automatic staged refinement. Returns the result dict."""
        t0 = time.time()
        P = self.P
        prog = self.progress
        n_ph = len(self.phases)
        cell_names = [n for n in P if n.startswith("c")]
        prog(0.02, "Looking for the best start values (cell scale)")
        if not self.fix_cell:                                   # v2.18.2: a standard's cell (certificate) is not rescaled
            self._prealign()
        prog(0.10, "Estimating the peak width")
        if not self.fixed_profile:                              # v2.18.2: a measured instrument profile is not re-estimated
            self._width_scan()
        off = {"displacement": ["disp"], "zero": ["zero"], "both": ["zero", "disp"]}[self.offset_mode]
        stages = [("Scale and background", [], 0.14),
                  ("Instrument offset and cell", off + cell_names, 0.25),
                  ("Peak widths (W, X)", off + ["W", "X"] + cell_names, 0.40),
                  ("All peak-shape parameters", off + ["U", "V", "W", "X", "Y"] + cell_names, 0.55),
                  ]
        if self.mode == "rietveld":
            stages.append(("Atomic displacement (B)", off + ["U", "V", "W", "X", "Y"] + cell_names + [f"B{i}" for i in range(n_ph)], 0.70))
            po_names = [f"po{i}" for i, ph_ in enumerate(self.phases) if ph_.po_axis is not None]
            struct_ph = [i for i in range(n_ph) if self._struct_mode(i)]
            if struct_ph:
                bn_ = [f"sB{i}.{e}" for i in struct_ph for e in self.phases[i].model.elements]
                xn_ = [f"sx{i}.{k}" for i in struct_ph for k in range(self.phases[i].model.n_pos)]
                base_ = off + ["U", "V", "W", "X", "Y"] + cell_names
                stages.append(("Atomic displacement parameters (B per element)", base_ + bn_, 0.72))
                stages.append(("Atom positions (soft bond restraints)", base_ + bn_ + xn_, 0.78))
                stages.append(("Atom positions (tight restraints)", base_ + bn_ + xn_, 0.82))
            if po_names:
                stages.append(("Preferred orientation", off + ["U", "V", "W", "X", "Y"] + cell_names + [f"B{i}" for i in range(n_ph)] + po_names, 0.76))
        low_angle = float(self.x.min()) < 20.0
        stages.append(("Fine tuning" + (" incl. peak asymmetry" if low_angle else ""),
                       off + ["U", "V", "W", "X", "Y"] + (["asym"] if low_angle else []) + cell_names
                       + ([f"B{i}" for i in range(n_ph)] + [f"po{i}" for i, ph_ in enumerate(self.phases) if ph_.po_axis is not None]
                          if self.mode == "rietveld" else []), 0.85))
        if self.fixed_profile or self.fix_cell:                 # v2.18.2: a fixed profile (U, V, W, X, Y) and/or a fixed cell are not free stages
            drop = (("U", "V", "W", "X", "Y") if self.fixed_profile else ()) + (tuple(cell_names) if self.fix_cell else ())
            stages = [(nm_, [n for n in nms if n not in drop], fr_) for nm_, nms, fr_ in stages]
        for nm, names, frac in stages:
            prog(frac, nm)
            self.restraint_sigma = self.restraint_scale * (0.06 if "soft" in nm else (0.025 if ("tight" in nm or "Fine" in nm) else 0.0))
            if "Fine" in nm and any(self._struct_mode(i) for i in range(n_ph)):
                names = names + [f"sB{i}.{e}" for i in range(n_ph) if self._struct_mode(i) for e in self.phases[i].model.elements] \
                    + [f"sx{i}.{k}" for i in range(n_ph) if self._struct_mode(i) for k in range(self.phases[i].model.n_pos)]
            if nm.startswith("Preferred"):
                self._po_scan()
            if self.mode == "lebail":
                self._le_bail_extract(15)
            if names:
                self._fit(names, max_nfev=cycles or 40)
        if self.mode == "lebail":
            for k in range(2):
                self._le_bail_extract(30)
                self._fit(stages[-1][1], max_nfev=cycles or 30)
            self._le_bail_extract(60)
        prog(0.95, "Statistics and uncertainties")
        res = self.result()
        self.restraint_sigma = 0.0
        res["seconds"] = time.time() - t0
        prog(1.0, "Done")
        return res

    # ------------------------------------------------------------------------------------------- results
    def esd(self):
        """Estimated standard deviations of the refined non-linear parameters from the Jacobian, Bérar-corrected for correlated residuals."""
        fn = [n for n, p in self.P.items() if p.free]
        if not fn:
            return {}
        th = np.array([self.P[n].value for n in fn])
        sc = np.array([self.P[n].scale for n in fn])
        f0 = self._residual(th)
        J = np.zeros((len(f0), len(fn)))
        for k in range(len(fn)):
            h = 1e-5 * max(1.0, abs(th[k])) if abs(th[k]) > 1e-8 else 1e-6
            t2 = th.copy()
            t2[k] += h
            J[:, k] = (self._residual(t2) - f0) / h
        try:
            cov = np.linalg.inv(J.T @ J)
        except np.linalg.LinAlgError:
            cov = np.linalg.pinv(J.T @ J)
        n_data = len(self.y)
        f_data = f0[:n_data]
        dof = max(n_data - len(fn) - len(self.phases) - self.nb, 1)
        red = float(f_data @ f_data) / dof
        # Bérar-Lelann: serial correlation of the residuals inflates the errors
        r = f0 - f0.mean()
        ac = float(np.dot(r[:-1], r[1:]) / max(np.dot(r, r), 1e-12))
        infl = max(1.0, math.sqrt(max((1 + max(ac, 0)) / (1 - min(max(ac, 0), 0.95)), 1.0)))
        self._cov = (fn, cov * red * infl ** 2)
        return {n: float(math.sqrt(max(cov[k, k], 0) * red)) * infl for k, n in enumerate(fn)}

    def result(self):
        names_free = [n for n, p in self.P.items() if p.free]
        esd = self.esd()
        v = self._values()
        s, bgc, calc, A = self._solve_linear(v)
        bg = self.T @ bgc
        st = self.stats(calc)
        phases = []
        zmv = []
        for i, ph in enumerate(self.phases):
            cell = self.cell(i)
            vol = cell_volume(cell)
            dens = ph.mass_cell / (0.602214 * vol)
            _, (pk_pos, pk_I) = self._phase_pattern(i, v, True)
            ce = {n: esd.get(f"c{i}.{n}") for n in CELL_DOF[ph.system]}
            veps = None
            if ph.system == "cubic" and ce.get("a"):
                veps = 3 * vol * ce["a"] / cell[0]
            phases.append({"name": ph.name, "system": ph.system, "cell": cell, "cell_esd": ce, "volume": vol, "volume_esd": veps, "density": dens,
                           "B": self.P[f"B{i}"].value, "B_esd": esd.get(f"B{i}"), "po_r": self.P[f"po{i}"].value if ph.po_axis is not None else None,
                           "po_axis": ph.po_axis, "po_esd": esd.get(f"po{i}"), "scale": float(s[i]), "n_reflections": ph.n_ref, "ticks": pk_pos,
                           "formula": ph.formula})
            if self._struct_mode(i):
                phases[-1]["structure"] = self._structure_table(i, v, esd)
            zmv.append(float(s[i]) * ph.mass_cell * vol if self.mode == "rietveld" else 0.0)
        tot = sum(zmv)
        for i, ph in enumerate(phases):
            ph["weight_fraction"] = zmv[i] / tot if tot > 0 else float("nan")
        # Bragg R-factor: observed intensities by partition of the pattern
        if self.mode == "rietveld":
            for i in range(len(self.phases)):
                phases[i]["Rb"] = self._bragg_r(i, v, s, calc, bg)
        glob = {n: {"value": self.P[n].value, "esd": esd.get(n), "label": self.P[n].label, "refined": n in names_free}
                for n in ("zero", "disp", "U", "V", "W", "X", "Y", "asym")}
        lam = self.rad[0][0]
        X_, Y_ = self.P["X"].value, self.P["Y"].value
        size = 0.9 * lam / (math.radians(X_) * 10) if X_ > 1e-6 else float("inf")
        strain = math.radians(Y_) / 4.0
        notes = self._notes(st, phases, size)
        return {"x": self.x, "y": self.y / self.scale_y, "calc": calc / self.scale_y, "bg": bg / self.scale_y, "resid": (self.y - calc) / self.scale_y,
                "stats": st, "phases": phases, "global": glob, "mode": self.mode, "counts": self.counts, "radiation": self.radiation,
                "size_nm_apparent": size, "strain_apparent": strain, "notes": notes, "bg_order": self.nb_poly - 1, "glass_bg": bool(self.bg_extra), "bg_measured": self.bg_measured,
                "free": names_free, "lebail_I": {i: self.I_free[i] for i in self.I_free} if self.mode == "lebail" else None}

    def structure_report(self, i):
        """Bond lengths, angles, coordination, bond-valence sums (with esd from the refinement covariance) and the plausibility verdict of phase i."""
        from mfs import structstats as SST
        ph = self.phases[i]
        mdl = ph.model
        v = self._values()
        names_pos = [f"sx{i}.{k}" for k in range(mdl.n_pos)]
        names_cell = [f"c{i}.{n}" for n in CELL_DOF[ph.system]]
        base = {f"p{k}": v[n] for k, n in enumerate(names_pos)}
        base.update({f"c.{n}": v[f"c{i}.{n}"] for n in CELL_DOF[ph.system]})
        to_base = {**{n: f"p{k}" for k, n in enumerate(names_pos)}, **{f"c{i}.{n}": f"c.{n}" for n in CELL_DOF[ph.system]}}
        cov, names = None, []
        if getattr(self, "_cov", None):
            fn, C = self._cov
            keep = [k for k, n in enumerate(fn) if n in to_base]
            names = [to_base[fn[k]] for k in keep]
            cov = C[np.ix_(keep, keep)] if keep else None
        metric = lambda vals: _cell_metric(cell_from_dof(ph.system, {n: vals[f"c.{n}"] for n in CELL_DOF[ph.system]}, ph.cell0))
        rep = SST.report(mdl, self._pos_vec(i, v), metric, cov, names, base, mdl.labels())
        rep["verdict"] = SST.verdict(rep)
        return rep

    def _structure_table(self, i, v, esd):
        """Refined atoms (label, element, x, y, z ± esd, B ± esd, occupancy, number of equivalent atoms in the cell, site-symmetry order) and the
        restraint statistics."""
        ph = self.phases[i]
        mdl = ph.model
        p = self._pos_vec(i, v)
        xs = mdl.site_frac.copy()
        for (r, c), val in zip(mdl.pos_params, p):
            xs[r] = xs[r] + mdl.basis[r][:, c] * val
        cov = None
        if getattr(self, "_cov", None):
            fn, C = self._cov
            idx = {n: k for k, n in enumerate(fn)}
        labs = mdl.labels()
        atoms = []
        for r in mdl.reps:
            ks = [k for k, (rr, c) in enumerate(mdl.pos_params) if rr == r]
            er = np.zeros(3)
            if ks and getattr(self, "_cov", None):
                names = [f"sx{i}.{k}" for k in ks]
                if all(n in idx for n in names):
                    Cs = C[np.ix_([idx[n] for n in names], [idx[n] for n in names])]
                    Bm = np.array([mdl.basis[r][:, mdl.pos_params[k][1]] for k in ks]).T
                    er = np.sqrt(np.maximum(np.diag(Bm @ Cs @ Bm.T), 0))
            els = [mdl.atom_el[j] for j in np.nonzero(mdl.atom_site == r)[0]]
            e0 = els[0]
            atoms.append({"label": labs[r], "element": "/".join(sorted(set(els))), "x": float(xs[r][0]), "y": float(xs[r][1]), "z": float(xs[r][2]),
                          "x_esd": float(er[0]), "y_esd": float(er[1]), "z_esd": float(er[2]), "B": float(v[f"sB{i}.{e0}"]), "B_esd": esd.get(f"sB{i}.{e0}"),
                          "occupancy": float(mdl.atom_occ[np.nonzero(mdl.atom_site == r)[0][0]]), "n_equivalent": int(np.sum(mdl.site_rep == r)),
                          "site_symmetry_order": int(mdl.site_sym_order[r]), "free_coordinates": len(ks)})
        out = {"spacegroup": mdl.spacegroup, "sg_number": mdl.sg_number, "n_atoms_cell": mdl.n_atoms, "atoms": atoms}
        try:
            out.update(self.structure_report(i))
        except Exception as e:                                       # the statistics must never cost the refinement result
            out["report_error"] = str(e)[:200]
        out["cn"] = {c["atom"]: c["CN"] for c in out.get("coordination", [])}
        out["_state"] = {"model": mdl, "p": p.copy(), "cell": self.cell(i, v), "B": self._B_atoms(i, v), "labels": labs}
        if ph.bond_arr is not None:
            bi, bj, bn, d0 = ph.bond_arr
            cell = self.cell(i, v)
            d = mdl.distances_metric(bi, bj, bn, p, _cell_metric(cell))
            sc = float(np.mean(cell[:3]) / np.mean(ph.cell0[:3]))
            dev = d - d0 * sc
            out["bond_rms_deviation"] = float(np.sqrt(np.mean(dev ** 2)))
            out["bond_max_deviation"] = float(np.max(np.abs(dev)))
            out["n_bonds"] = int(len(d))
        return out

    def _bragg_r(self, i, v, s, calc, bg):
        ph = self.phases[i]
        cell = self.cell(i, v)
        Gs = np.linalg.inv(_cell_metric(cell))
        inv_d2 = np.einsum("ki,ij,kj->k", ph.hkl, Gs, ph.hkl)
        d = 1.0 / np.sqrt(inv_d2)
        s2 = inv_d2 / 4.0
        I = self._I0(i, v) * self._dw(i, v, s2)
        ycp = calc - bg
        ratio = np.where(ycp > 1e-9, np.clip(self.y - bg, 0, None) / np.maximum(ycp, 1e-9), 0.0)
        Ic_tot = 0.0
        Io_tot = 0.0
        Icalc, Iobs = np.zeros(len(d)), np.zeros(len(d))
        for lam, wl in self.rad:
            sin_t = lam / (2.0 * d)
            th = np.arcsin(np.clip(sin_t, 0, 0.9995))
            tt = 2.0 * np.degrees(th) + v["zero"] + v["disp"] * np.cos(th)
            c2m = math.cos(math.radians(self.mono)) ** 2
            lp = (1.0 + c2m * np.cos(2 * th) ** 2) / ((1.0 + c2m) * sin_t ** 2 * np.cos(th))
            tan = np.tan(th)
            H, eta = _tch(v["U"] * tan ** 2 + v["V"] * tan + v["W"], v["X"] / np.cos(th) + v["Y"] * tan)
            Ik = s[i] * I * lp * wl
            cut = 25.0 * H
            i0 = np.searchsorted(self.x, tt - cut)
            i1 = np.searchsorted(self.x, tt + cut)
            nw = int(min(max(int((i1 - i0).max()), 3), 1500))
            idx = i0[:, None] + np.arange(nw)[None, :]
            valid = (idx < len(self.x)) & (idx < i1[:, None])
            idc = np.clip(idx, 0, len(self.x) - 1)
            prof = np.where(valid, _pv(self.x[idc] - tt[:, None], H[:, None], eta[:, None]), 0.0)
            inside = (tt > self.x[0]) & (tt < self.x[-1])
            Iobs += np.where(inside, Ik * (prof * (ratio[idc] * self.dx[idc])).sum(1), 0.0)
            Icalc += np.where(inside, Ik, 0.0)
        tot = float(Iobs.sum())
        return float(np.abs(Iobs - Icalc).sum() / tot) if tot > 0 else float("nan")

    def _notes(self, st, phases, size):
        n = []
        rw = st["Rwp"]
        if self.counts:
            if st["chi2"] > 4:
                n.append(f"GoF² = {st['chi2']:.1f}: the model misses something (a phase, preferred orientation, peak shape or a wrong cell). Look at the difference curve.")
        else:
            n.append("The data are not in counts: Rwp is valid, but Rexp and GoF are not meaningful.")
        if rw > 0.2:
            n.append(f"Rwp = {100 * rw:.0f} %: a poor fit. Check the phase list (a missing phase or wrong CIF is the usual cause).")
        for p in phases:
            if self.mode == "rietveld" and p["weight_fraction"] < 0.01:
                n.append(f"{p['name']}: less than 1 % – hardly present in the pattern.")
            if p.get("Rb") is not None and p["Rb"] > 0.25:
                n.append(f"{p['name']}: Bragg R = {100 * p['Rb']:.0f} % – the intensities do not match the structure (disorder, solvent / guests in the pores, preferred orientation).")
        for k in ("U", "V", "W", "X", "Y"):
            p = self.P[k]
            if p.free and (abs(p.value - p.lo) < 1e-6 or abs(p.value - p.hi) < 1e-6) and k in ("X", "W"):
                n.append(f"{p.label} sits at its limit.")
        if size < 100:
            n.append(f"Apparent size from the Lorentzian width X: {size:.0f} nm (instrument broadening NOT removed → a lower limit; measure a standard to separate it).")
        for p_ in phases:
            if p_.get("po_r") is not None and (p_["po_r"] < 0.45 or p_["po_r"] > 1.6):
                n.append(f"{p_['name']}: preferred-orientation r = {p_['po_r']:.2f} is far from 1 – a strong texture; check the sample preparation (back-loading, spinning, capillary).")
        n.append("Not refined: atom positions, anisotropic broadening, dispersion corrections – results are for phase identification, "
                 "lattice parameters and weight fractions, not for structure determination.")
        return n


def quick_refine(x, y, cifs, radiation="Cu Kα1", mode="rietveld", tt_range=None, bg_order=6, progress=None, mono_2theta=0.0, names=None, offset="displacement", po_axes=None, refine_structure=False, restraint_scale=1.0, instrument=None, glass_bg=True, bg_pattern=None):
    """One call: cifs = list of CIF paths / bytes (or Phase objects). Returns the result dict of Refinement.run().
    instrument: {"U", "V", "W", "X", "Y", "fixed": True} keeps a measured instrument profile fixed (see instrument_profile_from_standard)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if tt_range:
        m = (x >= tt_range[0]) & (x <= tt_range[1])
        x, y = x[m], y[m]
    phases = [c if isinstance(c, Phase) else Phase.from_cif(c, (names or [None] * len(cifs))[k]) for k, c in enumerate(cifs)]
    if refine_structure:
        for ph_ in phases:
            if ph_.model is None:
                ph_.enable_structure()
    for k, ax_ in enumerate(po_axes or []):
        if ax_ is not None and k < len(phases):
            phases[k].po_axis = tuple(int(v) for v in ax_)
    R = Refinement(x, y, phases, radiation, bg_order, mode, mono_2theta, progress=progress, offset=offset, instrument=instrument, glass_bg=glass_bg, bg_pattern=bg_pattern)
    R.restraint_scale = float(restraint_scale)
    return R.run()


# Silicon, NIST SRM 640 (lattice parameter a = 5.431195 Å at 22.5 °C; the value used for the cell start)
SI_CELL = (5.4312, 5.4312, 5.4312, 90.0, 90.0, 90.0)


def instrument_profile_from_standard(x, y, radiation="Cu Kα1", tt_range=None, bg_order=6, standard="Si", offset="displacement"):
    """v2.18.2: the instrument line profile (Caglioti-free TCH: U, V, W Gaussian and X, Y Lorentzian, in °2θ) measured on a line-profile
    standard: a Le Bail fit of the standard's own cell (Si SRM 640 by default). The returned values are then passed as
    instrument={..., 'fixed': True} to a sample refinement, which keeps them fixed. Returns {"U","V","W","X","Y","zero","rwp","cell_a"}."""
    if standard != "Si":
        raise ValueError("Only Si (SRM 640) is built in; for LaB6 give its cell with Phase.from_cell and call quick_refine(mode='lebail').")
    ph = Phase.from_cell(SI_CELL, "Si standard")
    if tt_range:
        m = (np.asarray(x, float) >= tt_range[0]) & (np.asarray(x, float) <= tt_range[1])
        x, y = np.asarray(x, float)[m], np.asarray(y, float)[m]
    R = Refinement(np.asarray(x, float), np.asarray(y, float), [ph], radiation, bg_order, "lebail", 0.0, offset=offset)
    R.fix_cell = True                                            # the standard's lattice parameter is known; only the line profile is fitted
    res = R.run()
    g = res["global"]
    out = {n: float(g[n]["value"]) for n in ("U", "V", "W", "X", "Y")}
    out["zero"] = float(g["zero"]["value"])                      # the instrument zero shift (°2θ) when offset="zero" (a separate standard sample: its own displacement is not the sample's)
    out["disp"] = float(g["disp"]["value"])
    out["rwp"] = float(res["stats"]["Rwp"])
    out["cell_a"] = float(res["phases"][0]["cell"][0])
    return out


PHASE_COLOURS = ["#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00"]


def draw(ax, res, rax=None, phase_colours=None, show_bg=True):
    """The standard Rietveld figure: observed (dots), calculated (line), background, difference (offset below, or in `rax`) and a row of
    reflection ticks per phase. Everything on ONE axis when `rax` is None (so it can go into the combined figure)."""
    x, y, c, bg, r = res["x"], res["y"], res["calc"], res["bg"], res["resid"]
    top = float(max(y.max(), c.max()))
    cols = phase_colours or PHASE_COLOURS
    ax.plot(x, y, "o", ms=1.6, mfc="none", mec="#222222", mew=0.4, label="observed")
    ax.plot(x, c, "-", color="#D62728", lw=0.9, label="calculated")
    if show_bg:
        ax.plot(x, bg, "--", color="#7F7F7F", lw=0.6, label="background")
    n_ph = len(res["phases"])
    tick_h = 0.035 * top
    if rax is None:
        off = -0.06 * top - tick_h * (n_ph + 1)
        span = float(np.max(np.abs(r))) if len(r) else 1.0
        scale = min(1.0, 0.18 * top / max(span, 1e-12))
        ax.plot(x, off + r * scale, "-", color="#0072B2", lw=0.6, label="difference")
        ax.axhline(off, color="#0072B2", lw=0.3)
        y0 = -0.04 * top
        ax.set_ylim(off - span * scale - 0.03 * top, top * 1.06)
    else:
        rax.plot(x, r, "-", color="#0072B2", lw=0.6)
        rax.axhline(0, color="k", lw=0.3)
        rax.set_ylabel("Δ")
        y0 = -0.04 * top
        ax.set_ylim(y0 - tick_h * (n_ph + 1.5), top * 1.06)
    for i, ph in enumerate(res["phases"]):
        t = np.asarray(ph["ticks"])
        t = t[(t >= x.min()) & (t <= x.max())]
        yy = (-0.045 * top - tick_h * i) if rax is None else (y0 - tick_h * (i + 1))
        ax.vlines(t, yy - tick_h * 0.4, yy + tick_h * 0.4, color=cols[i % len(cols)], lw=0.6)
        ax.plot([], [], "|", color=cols[i % len(cols)], label=ph["name"])
    ax.set_xlim(float(x.min()), float(x.max()))
    ax.set_ylabel("Intensity (counts)" if res.get("counts") else "Intensity (a.u.)")
    if rax is None:
        ax.set_xlabel(r"2$\theta$ (°)")
    else:
        rax.set_xlabel(r"2$\theta$ (°)")
    st_ = res["stats"]
    txt = f"R$_{{wp}}$ = {100 * st_['Rwp']:.1f} %" + (f",  GoF = {math.sqrt(st_['chi2']):.2f}" if res.get("counts") else "")
    ax.text(0.99, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=6)
    ax.legend(frameon=False, fontsize=6, loc="upper right", bbox_to_anchor=(1.0, 0.9))


def report_xlsx(res, sample="") -> bytes:
    """Excel report: parameters, phases, pattern (obs / calc / background / difference)."""
    import pandas as pd
    st_ = res["stats"]
    rows = [{"quantity": "mode", "value": res["mode"]}, {"quantity": "radiation", "value": res["radiation"]}, {"quantity": "sample", "value": sample},
            {"quantity": "Rwp (%)", "value": 100 * st_["Rwp"]}, {"quantity": "Rp (%)", "value": 100 * st_["Rp"]},
            {"quantity": "Rexp (%)", "value": 100 * st_["Rexp"] if res["counts"] else None}, {"quantity": "GoF", "value": math.sqrt(st_["chi2"]) if res["counts"] else None},
            {"quantity": "free parameters", "value": st_["n_par"]}, {"quantity": "points", "value": st_["N"]}]
    ph = []
    for p in res["phases"]:
        row = {"phase": p["name"], "system": p["system"], "weight fraction (%)": 100 * p["weight_fraction"] if res["mode"] == "rietveld" else None,
               "Bragg R (%)": 100 * p["Rb"] if p.get("Rb") is not None else None, "volume (Å³)": p["volume"], "density (g/cm³)": p["density"],
               "B overall (Å²)": p["B"] if res["mode"] == "rietveld" else None}
        for n, v in zip(("a", "b", "c", "alpha", "beta", "gamma"), p["cell"]):
            row[n] = v
        for n, e in (p.get("cell_esd") or {}).items():
            row[f"esd {n}"] = e
        ph.append(row)
    gl = pd.DataFrame([{"parameter": v["label"], "value": v["value"], "esd": v["esd"], "refined": v["refined"]} for v in res["global"].values()])
    pat = pd.DataFrame({"two_theta": res["x"], "observed": res["y"], "calculated": res["calc"], "background": res["bg"], "difference": res["resid"]})
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as xw:
        pd.DataFrame(rows).to_excel(xw, sheet_name="Summary", index=False)
        pd.DataFrame(ph).to_excel(xw, sheet_name="Phases", index=False)
        gl.to_excel(xw, sheet_name="Profile", index=False)
        pat.to_excel(xw, sheet_name="Pattern", index=False)
        pd.DataFrame({"note": res["notes"]}).to_excel(xw, sheet_name="Notes", index=False)
    return buf.getvalue()
