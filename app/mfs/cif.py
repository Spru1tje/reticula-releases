"""Small, dependency-free CIF reader + crystal-structure geometry for the Structure (CIF) page: cell, symmetry, atoms
with ADPs, bonds by covalent radii, growing (asymmetric unit, complete molecules, packing, supercell), H-bonds with
symmetry codes, distances/angles/torsions with approximate ESDs."""
from __future__ import annotations
import math, re
from fractions import Fraction
import numpy as np

# ---------------------------------------------------------------- element data (covalent radii: Cordero 2008; vdW: Bondi/Alvarez)
COV = {"H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11,
       "P": 1.07, "S": 1.05, "Cl": 1.02, "K": 2.03, "Ca": 1.76, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32,
       "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Br": 1.20, "Zr": 1.75, "Mo": 1.54, "Ag": 1.45, "Cd": 1.44, "I": 1.39,
       "Ce": 2.04, "Hf": 1.75, "La": 2.07, "Y": 1.90, "Al ": 1.21}
VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98, "P": 1.80,
       "Zr": 2.23, "Hf": 2.23, "Ti": 2.11, "Zn": 2.01, "Cu": 1.96, "Ce": 2.42, "B": 1.92, "Si": 2.10, "Na": 2.27, "Al": 1.84}
COLOR = {"H": "#FFFFFF", "C": "#555555", "N": "#3050F8", "O": "#E8221F", "F": "#90E050", "S": "#E7C21C", "P": "#FF8000",
         "Cl": "#1FF01F", "Br": "#A62929", "I": "#940094", "Zr": "#3FA9C9", "Hf": "#4DC2FF", "Ti": "#BFC2C7", "Zn": "#7D80B0",
         "Cu": "#C88033", "Ce": "#E5C07B", "Na": "#AB5CF2", "Al": "#BFA6A6", "Si": "#F0C8A0", "B": "#FFB5B5"}
METALS = {"Zr", "Hf", "Ti", "Ce", "Zn", "Cu", "Al", "Fe", "Co", "Ni", "Mn", "Cr", "V", "Mg", "Ca", "Y", "La", "Cd", "Ag", "Mo"}
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "S": 32.06, "P": 30.974, "Cl": 35.45, "Br": 79.904,
        "Zr": 91.224, "Hf": 178.49, "Ti": 47.867, "Zn": 65.38, "Cu": 63.546, "Ce": 140.116, "Na": 22.990, "Al": 26.982}


def num(v):
    """'1.234(5)' -> (1.234, 0.005); '?' / '.' -> (nan, nan)."""
    s = str(v).strip()
    m = re.match(r"^([-+]?\d*\.?\d*(?:[eE][-+]?\d+)?)(?:\((\d+)\))?$", s)
    if not m or m.group(1) in ("", "-", "+", "."):
        return float("nan"), float("nan")
    x = float(m.group(1))
    if m.group(2) is None:
        return x, float("nan")
    dec = len(m.group(1).split(".")[1]) if "." in m.group(1) and "e" not in m.group(1).lower() else 0
    return x, int(m.group(2)) * 10 ** (-dec)


def fmt_esd(x, e, maxdec=6):
    """1.53712, 0.0047 -> '1.537(5)'; nan esd -> plain."""
    if not np.isfinite(x):
        return "–"
    if not np.isfinite(e) or e <= 0:
        return f"{x:.{3 if abs(x) < 100 else 1}f}"
    d = max(0, min(maxdec, -int(math.floor(math.log10(e))) + (1 if int(e / 10 ** math.floor(math.log10(e))) == 1 else 0)))
    ei = int(round(e * 10 ** d))
    if ei >= 20 and d > 0:
        d -= 1; ei = int(round(e * 10 ** d))
    return f"{x:.{d}f}({ei})"


# ---------------------------------------------------------------- parser
def _tokens(text):
    lines = text.splitlines(); i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith(";"):
            buf = [ln[1:]]; i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                buf.append(lines[i]); i += 1
            yield ("\n".join(buf)).strip(), True
            i += 1; continue
        for m in re.finditer(r"""'((?:[^']|'(?=\S))*)'(?=\s|$)|"((?:[^"]|"(?=\S))*)"(?=\s|$)|(#.*$)|(\S+)""", ln):
            if m.group(3):
                break
            if m.group(1) is not None or m.group(2) is not None:
                yield (m.group(1) if m.group(1) is not None else m.group(2)), True
            else:
                yield m.group(4), False
        i += 1


def parse(text):
    """Returns {block name: {'items': {tag: value}, 'loops': [ {tag: [values]} ]}} (tags lower-case)."""
    blocks, cur = {}, None
    toks = list(_tokens(text)); i = 0; n = len(toks)
    while i < n:
        t, q = toks[i]
        tl = t.lower() if not q else ""
        if not q and tl.startswith("data_"):
            cur = {"items": {}, "loops": []}; blocks[t[5:]] = cur; i += 1; continue
        if cur is None:
            i += 1; continue
        if not q and tl == "loop_":
            i += 1; tags = []
            while i < n and not toks[i][1] and toks[i][0].startswith("_"):
                tags.append(toks[i][0].lower()); i += 1
            vals = []
            while i < n and not (not toks[i][1] and (toks[i][0].startswith("_") or toks[i][0].lower() in ("loop_",) or toks[i][0].lower().startswith("data_"))):
                vals.append(toks[i][0]); i += 1
            k = len(tags)
            if k:
                cur["loops"].append({tg: vals[j::k][: len(vals) // k] for j, tg in enumerate(tags)})
            continue
        if not q and t.startswith("_"):
            if i + 1 < n:
                cur["items"][tl] = toks[i + 1][0]; i += 2
            else:
                i += 1
            continue
        i += 1
    return blocks


def _loop(b, tag):
    for lp in b["loops"]:
        if tag in lp:
            return lp
    return None


# ---------------------------------------------------------------- symmetry
def parse_symop(s):
    s = s.replace(" ", "").lower(); R = np.zeros((3, 3)); t = np.zeros(3)
    for i, part in enumerate(s.split(",")):
        for sign, term in re.findall(r"([+-]?)([^+-]+)", part):
            sg = -1 if sign == "-" else 1
            m = re.match(r"^(\d+(?:\.\d+)?(?:/\d+)?)?\*?([xyz])$", term)
            if m:
                c = float(Fraction(m.group(1))) if m.group(1) else 1.0
                R[i, "xyz".index(m.group(2))] += sg * c
            else:
                m2 = re.match(r"^([xyz])/(\d+)$", term)
                if m2:
                    R[i, "xyz".index(m2.group(1))] += sg / float(m2.group(2))
                else:
                    t[i] += sg * float(Fraction(term))
    return R, t


def symop_str(R, t):
    """(R, t) -> '1-x, 1/2+y, 1-z' (CIF/Mercury style)."""
    out = []
    for i in range(3):
        terms = ""
        tt = Fraction(t[i]).limit_denominator(12)
        if tt != 0:
            terms += str(tt)
        for j, v in enumerate("xyz"):
            c = R[i, j]
            if abs(c) < 1e-9:
                continue
            terms += ("+" if c > 0 and terms else "") + ("-" if c < 0 else "") + v
        out.append(terms or "0")
    return ", ".join(out)


# ---------------------------------------------------------------- structure
class Structure:
    def __init__(self, block: dict, name=""):
        self.name = name; self.items = block["items"]; self.block = block
        g = lambda k: num(self.items.get(k, "nan"))
        self.cell = [g(f"_cell_length_{x}") for x in "abc"] + [g(f"_cell_angle_{x}") for x in ("alpha", "beta", "gamma")]
        a, b, c, al, be, ga = [v for v, _ in self.cell]
        al, be, ga = map(math.radians, (al, be, ga))
        V = a * b * c * math.sqrt(1 - math.cos(al) ** 2 - math.cos(be) ** 2 - math.cos(ga) ** 2 + 2 * math.cos(al) * math.cos(be) * math.cos(ga))
        self.V = V
        self.M = np.array([[a, b * math.cos(ga), c * math.cos(be)],
                           [0, b * math.sin(ga), c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)],
                           [0, 0, V / (a * b * math.sin(ga))]])
        self.Minv = np.linalg.inv(self.M)
        self.recip = np.array([b * c * math.sin(al), a * c * math.sin(be), a * b * math.sin(ga)]) / V   # a*, b*, c*
        lp = _loop(block, "_space_group_symop_operation_xyz") or _loop(block, "_symmetry_equiv_pos_as_xyz")
        ops = (lp.get("_space_group_symop_operation_xyz") or lp.get("_symmetry_equiv_pos_as_xyz")) if lp else ["x,y,z"]
        self.ops = [parse_symop(o) for o in ops]
        self.sg = self.items.get("_space_group_name_h-m_alt") or self.items.get("_symmetry_space_group_name_h-m") or "?"
        al_ = _loop(block, "_atom_site_label") or _loop(block, "_atom_site_fract_x")
        if al_ is None:
            raise ValueError("no _atom_site loop")
        nA = len(al_["_atom_site_fract_x"])
        labels = al_.get("_atom_site_label") or [f"X{i + 1}" for i in range(nA)]
        types = al_.get("_atom_site_type_symbol") or [re.match(r"[A-Za-z]{1,2}", l_).group(0) for l_ in labels]
        el = [re.match(r"([A-Z][a-z]?)", t_).group(1) if re.match(r"([A-Z][a-z]?)", t_) else t_ for t_ in types]
        el = [e if e in COV else e[0] for e in el]
        self.atoms = []
        for i in range(nA):
            f = [num(al_[f"_atom_site_fract_{x}"][i]) for x in "xyz"]
            gv = lambda k, d=".": (al_.get(k) or [d] * nA)[i]
            self.atoms.append({"label": labels[i], "el": el[i], "frac": np.array([v for v, _ in f]), "sfrac": np.array([e for _, e in f]),
                               "uiso": num(gv("_atom_site_u_iso_or_equiv", "nan"))[0] if gv("_atom_site_u_iso_or_equiv", "nan") not in (".", "?") else
                               (num(gv("_atom_site_b_iso_or_equiv", "nan"))[0] / (8 * math.pi ** 2) if gv("_atom_site_b_iso_or_equiv", "nan") not in (".", "?") else 0.02),
                               "occ": num(gv("_atom_site_occupancy", "1"))[0] if gv("_atom_site_occupancy", "1") not in (".", "?") else 1.0,
                               "calc": gv("_atom_site_calc_flag"), "dis_assembly": gv("_atom_site_disorder_assembly"),
                               "dis_group": gv("_atom_site_disorder_group"), "U": None})
        for a_ in self.atoms:
            if not np.isfinite(a_["uiso"]):
                a_["uiso"] = 0.02
        an = _loop(block, "_atom_site_aniso_label")
        if an:
            idx = {a_["label"]: a_ for a_ in self.atoms}
            N = np.diag(self.recip)
            for j, lab in enumerate(an["_atom_site_aniso_label"]):
                if lab not in idx:
                    continue
                u = {k: num(an[f"_atom_site_aniso_u_{k}"][j])[0] for k in ("11", "22", "33", "12", "13", "23")}
                Uc = np.array([[u["11"], u["12"], u["13"]], [u["12"], u["22"], u["23"]], [u["13"], u["23"], u["33"]]])
                idx[lab]["U"] = self.M @ N @ Uc @ N @ self.M.T          # Cartesian U (Å²)
        self.idx = {a_["label"]: i for i, a_ in enumerate(self.atoms)}

    # ---- basics
    def cart(self, f):
        return np.asarray(f) @ self.M.T

    def frac(self, x):
        return np.asarray(x) @ self.Minv.T

    @property
    def Z(self):
        return num(self.items.get("_cell_formula_units_z", "nan"))[0]

    def formula(self):
        f = self.items.get("_chemical_formula_sum")
        if f and f not in ("?", "."):
            return f
        from collections import Counter
        c = Counter()
        for a_ in self.atoms:
            c[a_["el"]] += a_["occ"]
        return " ".join(f"{e}{'' if abs(n - 1) < 1e-6 else f'{n:g}'}" for e, n in sorted(c.items(), key=lambda kv: (kv[0] != "C", kv[0] != "H", kv[0])))

    def density(self):
        mw = num(self.items.get("_chemical_formula_weight", "nan"))[0]
        z = self.Z
        if np.isfinite(mw) and np.isfinite(z):
            return mw * z / (self.V * 0.602214)
        return float("nan")

    # ---- atom images
    def image(self, i, op, T=(0, 0, 0)):
        R, t = self.ops[op]; return R @ self.atoms[i]["frac"] + t + np.asarray(T, float)

    def _all_images(self, ncell=1):
        """All symmetry images of all ASU atoms in cells -ncell..ncell around the origin cell (fractional coords + keys)."""
        F, K = [], []
        rng = range(-ncell, ncell + 1)
        for oi, (R, t) in enumerate(self.ops):
            base = np.array([R @ a_["frac"] + t for a_ in self.atoms])
            sh = np.floor(base + 1e-6)                                      # bring into [0,1)
            for T in [(i, j, k) for i in rng for j in rng for k in rng]:
                F.append(base - sh + T)
                K += [(ai, oi, tuple((np.array(T) - sh[ai]).astype(int))) for ai in range(len(self.atoms))]
        return np.vstack(F), K

    def keep(self, part="all"):
        """Atom indices to use: all, or disorder part filter ('major' = group blank/1, 'minor' = blank/2)."""
        out = []
        for i, a_ in enumerate(self.atoms):
            g = str(a_["dis_group"]).strip(".?-")
            if part == "all" or g == "" or (part == "major" and g in ("1",)) or (part == "minor" and g not in ("1",)):
                out.append(i)
        return out

    def bond_ok(self, e1, e2, d, tol=0.45):
        if e1 == "H" and e2 == "H":
            return False
        return 0.4 < d < COV.get(e1, 1.5) + COV.get(e2, 1.5) + tol

    def grow(self, mode="molecules", part="all", supercell=(1, 1, 1), max_atoms=6000, include_h=True):
        """Returns (sites, bonds): sites = list of dict(i, op, T, xyz, label, el, U, uiso, dis_group); bonds = [(s1, s2)].
        modes: 'asu', 'molecules' (complete molecules from the ASU), 'packing' (complete molecules with centroid in the
        cell, every symmetry copy), 'cell' (all atoms in [0,1)^3, bonds across cell edges not completed), 'supercell'."""
        use = set(self.keep(part))
        if not include_h:
            use = {i for i in use if self.atoms[i]["el"] != "H"}
        F, K = self._all_images(1)
        sel = np.array([k[0] in use for k in K])
        F = F[sel]; K = [k for k, s in zip(K, sel) if s]
        X = self.cart(F); EL = [self.atoms[k[0]]["el"] for k in K]
        key_index = {}
        for j, k in enumerate(K):
            key_index.setdefault((k[0], k[1], k[2]), j)
        from scipy.spatial import cKDTree
        tree = cKDTree(X)

        def neighbors(j):
            out = []
            for m in tree.query_ball_point(X[j], 3.2):
                if m != j and self.bond_ok(EL[j], EL[m], float(np.linalg.norm(X[j] - X[m]))):
                    out.append(m)
            return out
        if mode in ("cell", "supercell"):
            sc = supercell if mode == "supercell" else (1, 1, 1)
            sites = []
            for oi in range(len(self.ops)):
                for i in use:
                    f0 = self.image(i, oi); f0 = f0 - np.floor(f0 + 1e-6)
                    for T in np.ndindex(*sc):
                        f = f0 + T; sites.append(self._site(i, oi, tuple(np.array(T) - np.floor(self.image(i, oi) + 1e-6).astype(int)), f))
            sites = _dedup(sites)
            return sites, self._bonds(sites)
        # asu / molecules
        start = [key_index[(i, 0, (0, 0, 0))] for i in sorted(use) if (i, 0, (0, 0, 0)) in key_index]
        if mode == "asu":
            sites = [self._site(K[j][0], K[j][1], K[j][2], F[j]) for j in start]
            return sites, self._bonds(sites)
        seen = set(start); queue = list(start); polymeric = False
        while queue:
            j = queue.pop()
            for m in neighbors(j):
                if m not in seen:
                    seen.add(m); queue.append(m)
                    if len(seen) > max_atoms:
                        polymeric = True; queue = []; break
            # atoms at the edge of the ±1 cell block: the fragment is (probably) polymeric
            if np.any(np.abs(F[j]) > 1.6):
                polymeric = True
        mol = [self._site(K[j][0], K[j][1], K[j][2], F[j]) for j in sorted(seen)]
        mol = _dedup(mol)
        if mode == "molecules":
            s_, b_ = mol, self._bonds(mol)
            self.polymeric = polymeric
            return s_, b_
        # packing: each symmetry copy of each complete molecule, translated so its centroid lies in the cell
        comps = _components(mol, self._bonds(mol))
        sites = []
        for comp in comps:
            fr = np.array([mol[c]["frac"] for c in comp])
            for oi, (R, t) in enumerate(self.ops):
                g = fr @ R.T + t; shift = -np.floor(g.mean(0) + 1e-9)
                for c, f in zip(comp, g + shift):
                    s0 = mol[c]; sites.append(self._site(s0["i"], oi, tuple(np.round(f - (R @ self.atoms[s0["i"]]["frac"] + t)).astype(int)), f))
        sites = _dedup(sites)
        if supercell != (1, 1, 1):
            base = list(sites); sites = []
            for T in np.ndindex(*supercell):
                for s0 in base:
                    sites.append(self._site(s0["i"], s0["op"], tuple(np.array(s0["T"]) + T), s0["frac"] + T))
        self.polymeric = polymeric
        return sites, self._bonds(sites)

    def _site(self, i, op, T, f):
        a_ = self.atoms[i]; R = self.ops[op][0]
        Rc = self.M @ R @ self.Minv                                       # Cartesian rotation part of the op
        U = None if a_["U"] is None else Rc @ a_["U"] @ Rc.T
        return {"i": i, "op": op, "T": tuple(int(v) for v in T), "frac": np.asarray(f, float), "xyz": self.cart(f),
                "label": a_["label"], "el": a_["el"], "U": U, "uiso": a_["uiso"], "dis_group": a_["dis_group"]}

    def _bonds(self, sites, tol=0.45):
        if not sites:
            return []
        from scipy.spatial import cKDTree
        X = np.array([s["xyz"] for s in sites]); tree = cKDTree(X); out = []
        for a, b in tree.query_pairs(3.2):
            sa, sb = sites[a], sites[b]
            ga, gb = str(sa["dis_group"]).strip(".?"), str(sb["dis_group"]).strip(".?")
            if ga and gb and ga != gb:                                     # different disorder parts never bond
                continue
            if self.bond_ok(sa["el"], sb["el"], float(np.linalg.norm(X[a] - X[b])), tol):
                out.append((a, b))
        return out

    def symcode(self, op, T):
        R, t = self.ops[op]; return symop_str(R, t + np.asarray(T, float))

    # ---- geometry with approximate ESDs (coordinate su only, correlations and cell su neglected)
    def _sx(self, s):
        sf = self.atoms[s["i"]]["sfrac"]; sf = np.where(np.isfinite(sf), sf, 0.0)
        return np.abs(self.M) @ np.abs(self.ops[s["op"]][0]) @ sf            # crude Cartesian su (Å)

    def dist(self, s1, s2):
        d = s2["xyz"] - s1["xyz"]; L = float(np.linalg.norm(d)); u = d / (L or 1)
        e = math.sqrt(float((u * self._sx(s1)) @ (u * self._sx(s1)) + (u * self._sx(s2)) @ (u * self._sx(s2))))
        return L, e

    @staticmethod
    def angle_of(p1, p2, p3):
        a, b = p1 - p2, p3 - p2
        return math.degrees(math.acos(max(-1, min(1, a @ b / (np.linalg.norm(a) * np.linalg.norm(b))))))

    @staticmethod
    def torsion_of(p1, p2, p3, p4):
        b0, b1, b2 = p1 - p2, p3 - p2, p4 - p3
        b1n = b1 / np.linalg.norm(b1); v = b0 - (b0 @ b1n) * b1n; w = b2 - (b2 @ b1n) * b1n
        return math.degrees(math.atan2(np.cross(b1n, v) @ w, v @ w))


def _dedup(sites, tol=1e-3):
    out, seen = [], set()
    for s in sites:
        k = (s["i"],) + tuple(np.round(s["frac"] / tol).astype(int))
        if k not in seen:
            seen.add(k); out.append(s)
    return out


def _components(sites, bonds):
    adj = {i: [] for i in range(len(sites))}
    for a, b in bonds:
        adj[a].append(b); adj[b].append(a)
    seen, comps = set(), []
    for i in range(len(sites)):
        if i in seen:
            continue
        st_, comp = [i], []
        seen.add(i)
        while st_:
            j = st_.pop(); comp.append(j)
            for m in adj[j]:
                if m not in seen:
                    seen.add(m); st_.append(m)
        comps.append(comp)
    return comps


def read(src, name=""):
    """Path, bytes or file-like -> list of Structure (one per data block with atoms)."""
    if hasattr(src, "read"):
        raw = src.read()
    elif isinstance(src, (bytes, bytearray)):
        raw = src
    else:
        raw = open(src, "rb").read()
    text = raw.decode("utf-8", "replace") if isinstance(raw, (bytes, bytearray)) else raw
    out = []
    for bn, b in parse(text).items():
        if _loop(b, "_atom_site_fract_x"):
            out.append(Structure(b, bn))
    if not out:
        raise ValueError("no data block with atom sites")
    return out


# ---------------------------------------------------------------- hydrogen bonds
DONORS, ACCEPTORS = {"O", "N"}, {"O", "N", "F", "Cl"}


def _hb_ok(S, d, h, a, dHA, ang, cut):
    return dHA < cut(a) and ang >= 120 and d is not a


def hbonds(S, part="all", delta=0.0, min_angle=120.0):
    """D–H···A from the ASU H atoms: H···A < r_vdW(H) + r_vdW(A) − delta and D–H···A ≥ min_angle (Mercury-like).
    Returns rows with symmetry codes of the acceptor (and donor, always the ASU)."""
    use = S.keep(part)
    F, K = S._all_images(1)
    X = S.cart(F)
    from scipy.spatial import cKDTree
    tree = cKDTree(X)
    rows = []
    ident = {i: j for j, k in enumerate(K) if k[1] == 0 and k[2] == (0, 0, 0) for i in [k[0]]}
    for hi in use:
        if S.atoms[hi]["el"] != "H" or hi not in ident:
            continue
        jh = ident[hi]; xh = X[jh]
        # donor = nearest bonded N/O
        cand = [m for m in tree.query_ball_point(xh, 1.3) if S.atoms[K[m][0]]["el"] in DONORS]
        if not cand:
            continue
        jd = min(cand, key=lambda m: np.linalg.norm(X[m] - xh))
        for ja in tree.query_ball_point(xh, 3.2):
            ea = S.atoms[K[ja][0]]["el"]
            if ea not in ACCEPTORS or ja == jd or K[ja][0] not in use:
                continue
            dHA = float(np.linalg.norm(X[ja] - xh))
            if dHA >= VDW["H"] + VDW.get(ea, 1.55) - delta:
                continue
            ang = Structure.angle_of(X[jd], xh, X[ja])
            if ang < min_angle:
                continue
            sD, sH, sA = (S._site(K[j][0], K[j][1], K[j][2], F[j]) for j in (jd, jh, ja))
            dDH = S.dist(sD, sH)[0]; dDA, eDA = S.dist(sD, sA)
            rows.append({"D": sD["label"], "H": sH["label"], "A": sA["label"], "D–H (Å)": dDH, "H···A (Å)": dHA,
                         "D···A (Å)": dDA, "σ D···A": eDA, "D–H···A (°)": ang, "symmetry (A)": S.symcode(K[ja][1], K[ja][2]),
                         "_A": (K[ja][0], K[ja][1], K[ja][2])})
    rows.sort(key=lambda r: (r["D"], r["H···A (Å)"]))
    return rows


def hbond_pairs(S, sites, delta=0.0, min_angle=120.0):
    """(i_H, i_A) site-index pairs inside a drawn site list (packing / molecules)."""
    if not sites:
        return []
    from scipy.spatial import cKDTree
    X = np.array([s["xyz"] for s in sites]); tree = cKDTree(X); out = []
    for i, s in enumerate(sites):
        if s["el"] != "H":
            continue
        cand = [m for m in tree.query_ball_point(X[i], 1.3) if sites[m]["el"] in DONORS]
        if not cand:
            continue
        d = min(cand, key=lambda m: np.linalg.norm(X[m] - X[i]))
        for a in tree.query_ball_point(X[i], 3.2):
            if a == d or sites[a]["el"] not in ACCEPTORS:
                continue
            if np.linalg.norm(X[a] - X[i]) < VDW["H"] + VDW.get(sites[a]["el"], 1.55) - delta and \
                    Structure.angle_of(X[d], X[i], X[a]) >= min_angle:
                out.append((i, a))
    return out


# ---------------------------------------------------------------- tables of bonds / angles / torsions
def geom_tables(S, include_h=False, publ_only=False):
    """From the CIF _geom loops (values with the refinement ESDs). If a loop is missing, values are computed from the
    coordinates (approximate ESDs, correlations neglected)."""
    import pandas as pd
    out = {}
    specs = {"bonds": ("_geom_bond_atom_site_label_1", ["_geom_bond_atom_site_label_1", "_geom_bond_atom_site_label_2"], "_geom_bond_distance", "_geom_bond_site_symmetry_2", "Distance (Å)"),
             "angles": ("_geom_angle_atom_site_label_1", ["_geom_angle_atom_site_label_1", "_geom_angle_atom_site_label_2", "_geom_angle_atom_site_label_3"], "_geom_angle", None, "Angle (°)"),
             "torsions": ("_geom_torsion_atom_site_label_1", [f"_geom_torsion_atom_site_label_{i}" for i in range(1, 5)], "_geom_torsion", None, "Torsion (°)")}
    for key, (tag, labs, val, sym, col) in specs.items():
        lp = _loop(S.block, tag)
        if lp:
            n = len(lp[tag]); rows = []
            for j in range(n):
                names = [lp[l_][j] for l_ in labs]
                if not include_h and any(S.atoms[S.idx[nm_]]["el"] == "H" for nm_ in names if nm_ in S.idx):
                    continue
                pf = (lp.get(val.replace("_distance", "") + "_publ_flag") or lp.get(val + "_publ_flag") or ["?"] * n)[j]
                if publ_only and pf.lower() not in ("y", "yes"):
                    continue
                v, e = num(lp[val][j])
                syms = []
                for c in [k_ for k_ in lp if "site_symmetry" in k_]:
                    if lp[c][j] not in (".", "?"):
                        syms.append(lp[c][j])
                rows.append({"atoms": "–".join(names), col: fmt_esd(v, e), "value": v, "esd": e, "symmetry": " ".join(syms) or "",
                             "source": "CIF"})
            out[key] = pd.DataFrame(rows)
        else:
            out[key] = _computed(S, key, include_h, col)
    return out


def _computed(S, key, include_h, col):
    import pandas as pd
    sites, bonds = S.grow("asu", include_h=include_h)
    adj = {i: set() for i in range(len(sites))}
    for a, b in bonds:
        adj[a].add(b); adj[b].add(a)
    rows = []
    X = [s["xyz"] for s in sites]
    if key == "bonds":
        for a, b in bonds:
            d, e = S.dist(sites[a], sites[b]); rows.append({"atoms": f"{sites[a]['label']}–{sites[b]['label']}", col: fmt_esd(d, e) + "*", "value": d, "esd": e})
    elif key == "angles":
        for j in adj:
            nb = sorted(adj[j])
            for x in range(len(nb)):
                for y in range(x + 1, len(nb)):
                    v = Structure.angle_of(X[nb[x]], X[j], X[nb[y]])
                    rows.append({"atoms": f"{sites[nb[x]]['label']}–{sites[j]['label']}–{sites[nb[y]]['label']}", col: f"{v:.1f}*", "value": v, "esd": float("nan")})
    else:
        for b, c in bonds:
            for a in adj[b] - {c}:
                for d in adj[c] - {b, a}:
                    v = Structure.torsion_of(X[a], X[b], X[c], X[d])
                    rows.append({"atoms": "–".join(sites[k]["label"] for k in (a, b, c, d)), col: f"{v:.1f}*", "value": v, "esd": float("nan")})
    df = pd.DataFrame(rows)
    if len(df):
        df["source"] = "computed (* approx. ESD, correlations neglected)"
    return df
