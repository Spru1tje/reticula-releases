"""
Icon generators: make house-style icons from data instead of drawing them.

    icon_from_cif("MOF-801.cif", mode="net")        # SBU-node net (MOF topology), or mode="atoms"
    icon_from_xy(x, y, kind="line")                 # the user's own spectrum / isotherm / TGA as an icon
    icon_from_molecule(atoms, bonds)                # 2D molecule (e.g. from RDKit coordinates)
    icon_from_smiles("OC(=O)c1ccc(cc1)C(O)=O")      # needs RDKit (optional dependency)
    animate(svg, kind="auto")                       # SMIL animation: bubbles / heat / pulse / stir

Pure Python (numpy optional for speed is not required).
"""
from __future__ import annotations

import math
import re
import shlex

from .core import Icon, Op
from .library import axes

METALS = set("Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr Nb Mo Ru Rh Pd Ag Cd In Sn Cs Ba La Ce "
             "Pr Nd Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Pb Bi Th U".split())
COV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "S": 1.05, "P": 1.07, "Cl": 1.02, "Br": 1.20,
       "B": 0.84, "Si": 1.11}


# ----------------------------------------------------------------- CIF
def _num(s: str) -> float:
    return float(re.sub(r"\(.*?\)", "", s))


def parse_cif(text: str) -> dict:
    """Minimal CIF reader: cell, symmetry operators, atom sites (label, element, fract xyz)."""
    lines = [l.rstrip() for l in text.splitlines()]
    cell = {}
    for k in ("a", "b", "c", "alpha", "beta", "gamma"):
        m = re.search(rf"^_cell_(?:length|angle)_{k}\s+(\S+)", text, re.M)
        cell[k] = _num(m.group(1)) if m else (90.0 if k in ("alpha", "beta", "gamma") else 10.0)
    loops, i = [], 0
    while i < len(lines):
        if lines[i].strip().lower() == "loop_":
            heads, rows = [], []
            i += 1
            while i < len(lines) and lines[i].strip().startswith("_"):
                heads.append(lines[i].strip().split()[0].lower())
                i += 1
            buf = []
            while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith(("_", "loop_", "#")) \
                    and not lines[i].startswith(";"):
                try:
                    buf += shlex.split(lines[i], posix=True)
                except ValueError:
                    buf += lines[i].split()
                while len(buf) >= len(heads):
                    rows.append(buf[:len(heads)])
                    buf = buf[len(heads):]
                i += 1
            loops.append((heads, rows))
            continue
        i += 1
    ops, atoms = [], []
    for heads, rows in loops:
        key = next((h for h in heads if h in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz")), None)
        if key:
            ops = [r[heads.index(key)] for r in rows]
        if "_atom_site_fract_x" in heads:
            ix = heads.index("_atom_site_fract_x")
            it = heads.index("_atom_site_type_symbol") if "_atom_site_type_symbol" in heads else heads.index("_atom_site_label")
            for r in rows:
                el = re.match(r"[A-Z][a-z]?", r[it]).group(0)
                atoms.append((el, _num(r[ix]), _num(r[ix + 1]), _num(r[ix + 2])))
    return {"cell": cell, "ops": ops or ["x,y,z"], "atoms": atoms}


def _symop(op: str):
    parts = op.lower().replace(" ", "").split(",")

    def f(x, y, z):
        env = {"x": x, "y": y, "z": z}
        return tuple(eval(re.sub(r"(\d)\/(\d)", r"(\1/\2)", p), {"__builtins__": {}}, env) for p in parts)  # noqa: S307
    return f


def expand(cif: dict, tol=1e-3):
    """Apply symmetry, wrap into [0,1), remove duplicates. Returns [(el, fx, fy, fz)]."""
    fns = [_symop(o) for o in cif["ops"]]
    seen, out = set(), []
    for el, x, y, z in cif["atoms"]:
        for f in fns:
            p = tuple(v % 1.0 for v in f(x, y, z))
            k2 = (el,) + tuple(round(v, 2) % 1.0 for v in p)
            if k2 in seen:
                continue
            seen.add(k2)
            out.append((el,) + p)
    return out


def _cart(cell):
    a, b, c = cell["a"], cell["b"], cell["c"]
    al, be, ga = (math.radians(cell[k]) for k in ("alpha", "beta", "gamma"))
    v = math.sqrt(1 - math.cos(al) ** 2 - math.cos(be) ** 2 - math.cos(ga) ** 2 + 2 * math.cos(al) * math.cos(be) * math.cos(ga))
    M = ((a, b * math.cos(ga), c * math.cos(be)),
         (0, b * math.sin(ga), c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)),
         (0, 0, c * v / math.sin(ga)))
    return lambda f: tuple(sum(M[i][j] * f[j] for j in range(3)) for i in range(3))


def _view(yaw=35, pitch=25):
    cy, sy, cp, sp = math.cos(math.radians(yaw)), math.sin(math.radians(yaw)), math.cos(math.radians(pitch)), math.sin(math.radians(pitch))

    def r(p):
        x, y, z = p
        x, z = cy * x + sy * z, -sy * x + cy * z
        y, z = cp * y - sp * z, sp * y + cp * z
        return x, -y, z
    return r


def _fit(points, margin=7):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    s = (64 - 2 * margin) / max(max(xs) - min(xs), max(ys) - min(ys), 1e-6)
    ox = margin - min(xs) * s + ((64 - 2 * margin) - (max(xs) - min(xs)) * s) / 2
    oy = margin - min(ys) * s + ((64 - 2 * margin) - (max(ys) - min(ys)) * s) / 2
    return lambda p: (p[0] * s + ox, p[1] * s + oy), s


def icon_from_cif(path_or_text: str, name: str = "cif_structure", label: str = "Crystal Structure",
                  mode: str = "net", yaw: float = 35, pitch: float = 25, cell_box: bool = True) -> Icon:
    """mode="net": metal clusters (SBUs) as nodes, linked to their nearest neighbours (MOF topology view).
    mode="atoms": ball-and-stick of the unit-cell content (metals accent, O/N secondary, C body; H hidden)."""
    text = open(path_or_text, encoding="utf-8", errors="ignore").read() if "\n" not in path_or_text else path_or_text
    cif = parse_cif(text)
    frac = expand(cif)
    to_cart = _cart(cif["cell"])
    rot = _view(yaw, pitch)
    ic = Icon(name, label, "materials", ["CIF", "crystal structure", "generated"])
    corners = [to_cart((i, j, k)) for i in (0, 1) for j in (0, 1) for k in (0, 1)]

    if mode == "net":
        metals = [to_cart(p[1:]) for p in frac if p[0] in METALS]
        # periodic images so clusters on the cell edge stay whole
        imgs = [to_cart((p[1] + i, p[2] + j, p[3] + k)) for p in frac if p[0] in METALS
                for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
        clusters, used = [], [False] * len(imgs)
        for i, p in enumerate(imgs):
            if used[i]:
                continue
            grp, stack = [], [i]
            used[i] = True
            while stack:
                q = stack.pop()
                grp.append(imgs[q])
                for j, r in enumerate(imgs):
                    if not used[j] and math.dist(imgs[q], r) < 4.2:
                        used[j] = True
                        stack.append(j)
            c = tuple(sum(v[k] for v in grp) / len(grp) for k in range(3))
            clusters.append(c)
        # keep cluster centres inside (or on) the cell
        inv_ok = []
        lo = [min(cc[k] for cc in corners) - 0.5 for k in range(3)]
        hi = [max(cc[k] for cc in corners) + 0.5 for k in range(3)]
        for c in clusters:
            if all(lo[k] <= c[k] <= hi[k] for k in range(3)):
                inv_ok.append(c)
        nodes = inv_ok or clusters[:8]
        ds = sorted(math.dist(a, b) for i, a in enumerate(nodes) for b in nodes[i + 1:])
        dmin = ds[0] if ds else 1
        edges = [(a, b) for i, a in enumerate(nodes) for b in nodes[i + 1:] if math.dist(a, b) < dmin * 1.12]
        P = [rot(p) for p in nodes] + [rot(p) for p in corners]
        fit, s = _fit(P)
        if cell_box:
            cc = [fit(rot(p)) for p in corners]
            for i in range(8):
                for j in range(i + 1, 8):
                    if bin(i ^ j).count("1") == 1:
                        ic.line(*cc[i], *cc[j], 1.2)
        for a, b in sorted(edges, key=lambda e: rot(e[0])[2] + rot(e[1])[2]):
            ic.sec().line(*fit(rot(a)), *fit(rot(b)), 2.4)
        for p in sorted(nodes, key=lambda p: rot(p)[2]):
            x, y = fit(rot(p))
            ic.acc().circle(x, y, 3.6)
        return ic

    # atoms mode
    cart = [(el, to_cart(p)) for el, *p in frac if el != "H"]
    if len(cart) > 900:
        cart = cart[:900]
    P = [rot(c) for _, c in cart] + [rot(p) for p in corners]
    fit, s = _fit(P)
    if cell_box:
        cc = [fit(rot(p)) for p in corners]
        for i in range(8):
            for j in range(i + 1, 8):
                if bin(i ^ j).count("1") == 1:
                    ic.line(*cc[i], *cc[j], 1)
    bonds = []
    for i in range(len(cart)):
        ei, pi = cart[i]
        for j in range(i + 1, len(cart)):
            ej, pj = cart[j]
            lim = (1.75 if ei in METALS else COV.get(ei, 0.8)) + (1.75 if ej in METALS else COV.get(ej, 0.8)) - 0.25 \
                if (ei in METALS or ej in METALS) else COV.get(ei, 0.8) + COV.get(ej, 0.8) + 0.35
            if math.dist(pi, pj) < lim:
                bonds.append((i, j))
    for i, j in bonds:
        ic.line(*fit(rot(cart[i][1])), *fit(rot(cart[j][1])), 0.9)
    for el, p in sorted(cart, key=lambda t: rot(t[1])[2]):
        x, y = fit(rot(p))
        if el in METALS:
            ic.acc().circle(x, y, 1.9)
        elif el in ("O", "N", "S"):
            ic.sec().circle(x, y, 1.0)
        else:
            ic.circle(x, y, 0.8)
    return ic


# ----------------------------------------------------------------- data -> icon
def icon_from_xy(x, y, kind: str = "line", name: str = "my_data", label: str = "My Data", invert: bool = False) -> Icon:
    """kind: line | sticks | dots. invert=True for transmission-type data (dips)."""
    n = len(x)
    if n > 400:
        step = n // 400 + 1
        x, y = x[::step], y[::step]
    x0, x1 = min(x), max(x)
    y0, y1 = min(y), max(y)
    sx = 44 / ((x1 - x0) or 1)
    sy = 40 / ((y1 - y0) or 1)
    pts = [(12 + (xi - x0) * sx, (12 + (yi - y0) * sy) if invert else (52 - (yi - y0) * sy)) for xi, yi in zip(x, y)]
    ic = Icon(name, label, "data", ["generated", "my data"])
    axes(ic)
    if kind == "sticks":
        for px, py in pts:
            ic.acc().line(px, 54, px, py, 1.8)
    elif kind == "dots":
        for px, py in pts[:: max(1, len(pts) // 25)]:
            ic.acc().circle(px, py, 1.8)
    else:
        ic.acc().spoly(pts, 2.4)
    return ic


# ----------------------------------------------------------------- molecules
EL_ROLE = {"C": "m", "H": "m", "O": "s", "N": "a", "S": "a", "P": "a"}
# CPK-inspired atom colours for colour/illustration styles (solid stays monochrome)
CPK = {"H": "#C9D1DB", "O": "#E5484D", "N": "#3E63DD", "S": "#F5C518", "P": "#F08A24", "Cl": "#30A46C",
       "F": "#8BD45C", "Br": "#A0522D", "I": "#7B4FB5", "B": "#F4A6A6", "Si": "#D6B48C", "C": "#3B4552"}
CPK_TEXT = {"O": "#D13438", "N": "#2F55C9", "S": "#B58900", "Cl": "#1F8A55", "H": "#64748B"}


def icon_from_molecule(atoms, bonds, name: str = "molecule", label: str = "Molecule", show_h: bool = False,
                       mode: str = "skeletal") -> Icon:
    """atoms: [(element, x, y)], bonds: [(i, j, order)] in any 2D units.
    mode: "skeletal" (carbons implicit), "ball" (ball-and-stick, all atoms), "space" (space-filling),
    "label" (skeletal with element letters on heteroatoms)."""
    if mode in ("ball", "space"):
        show_h = True
    keep = [i for i, a in enumerate(atoms) if show_h or a[0] != "H"]
    if mode == "label":
        keep = [i for i, a in enumerate(atoms) if a[0] != "H" or len(atoms) <= 5]
    idx = {i: k for k, i in enumerate(keep)}
    pts = [(atoms[i][1], -atoms[i][2]) for i in keep]
    fit, s = _fit(pts, margin=12 if mode == "space" else 8)
    P = [fit(p) for p in pts]
    ic = Icon(name, label, "materials", ["molecule", "generated"])
    ic.ink = True
    if mode == "space":
        rad = {"H": 0.55, "C": 0.85, "N": 0.8, "O": 0.78, "S": 1.0, "Cl": 0.98}
        order = sorted(range(len(keep)), key=lambda k: (P[k][1], -P[k][0]))     # back (top) to front
        for k in order:
            el = atoms[keep[k]][0]
            r = max(3.2, min(13.0, s * rad.get(el, 0.9)))
            ic.circle(*P[k], r, role=EL_ROLE.get(el, "a"))
            ic.fixed(CPK.get(el, "#4B5563"))
            ic.acc().circle(P[k][0] - r * 0.32, P[k][1] - r * 0.34, r * 0.28)
            ic.fixed("#FFFFFF")
            ic.ops[-1].attrs["opacity"] = 0.5
        return ic
    for i, j, order in bonds:
        if i not in idx or j not in idx:
            continue
        (xa, ya), (xb, yb) = P[idx[i]], P[idx[j]]
        if order >= 2:
            dx, dy = yb - ya, xa - xb
            L = math.hypot(dx, dy) or 1
            o = 2.0
            ic.line(xa + dx / L * o, ya + dy / L * o, xb + dx / L * o, yb + dy / L * o, 2.2)
            ic.line(xa - dx / L * o, ya - dy / L * o, xb - dx / L * o, yb - dy / L * o, 2.2)
        else:
            ic.line(xa, ya, xb, yb, 3.0)
    for k, i in enumerate(keep):
        el = atoms[i][0]
        if mode == "ball":
            r = 3.6 if el == "H" else 5.6
            ic.circle(*P[k], r, role=EL_ROLE.get(el, "a"))
            ic.fixed(CPK.get(el, "#4B5563"))
            ic.acc().circle(P[k][0] - r * 0.3, P[k][1] - r * 0.32, r * 0.3)
            ic.fixed("#FFFFFF")
            ic.ops[-1].attrs["opacity"] = 0.55
            continue
        if el == "C":
            continue
        if mode == "label":
            ic.cut().circle(*P[k], 5.6)
            ic.text(P[k][0], P[k][1] + 3.4, el, 9.5 if len(el) == 1 else 8, "bold", role="m", fill=CPK_TEXT.get(el, "#1F2328"))
            continue
        ic.circle(*P[k], 3.4 if el == "H" else 4.4, role=EL_ROLE.get(el, "a"))
        if el in CPK:
            ic.fixed(CPK[el])
    return ic


def icon_from_smiles(smiles: str, name: str = "molecule", label: str | None = None) -> Icon:
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
    except ImportError as e:  # pragma: no cover
        raise RuntimeError("icon_from_smiles needs RDKit: pip install rdkit") from e
    m = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(m)
    conf = m.GetConformer()
    atoms = [(a.GetSymbol(), conf.GetAtomPosition(a.GetIdx()).x, conf.GetAtomPosition(a.GetIdx()).y) for a in m.GetAtoms()]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), 2 if b.GetBondTypeAsDouble() >= 1.5 else 1) for b in m.GetBonds()]
    return icon_from_molecule(atoms, bonds, name, label or smiles)


# ----------------------------------------------------------------- animation (SMIL)
_SELF = re.compile(r'<(circle|path|rect|line|polyline|polygon|ellipse)\b([^>]*class="(rt-a[^"]*)"[^>]*?)/>')


def animate(svg: str, kind: str = "auto", dur: float = 1.6) -> str:
    """Add SMIL animation to accent shapes. kind: bubbles | heat | pulse | stir | auto.
    Works in browsers and most SVG viewers; Office/PDF show the still frame."""
    k = [0]

    def sub(m):
        tag, attrs, cls = m.groups()
        k[0] += 1
        delay = f"{(k[0] % 5) * dur / 5:.2f}s"
        kk = kind
        if kk == "auto":
            kk = "bubbles" if tag == "circle" else ("heat" if "-line" in cls and tag == "path" else "pulse")
        if kk == "bubbles":
            anim = (f'<animateTransform attributeName="transform" type="translate" values="0 0;0 -6" dur="{dur}s" '
                    f'begin="{delay}" repeatCount="indefinite"/><animate attributeName="opacity" values="1;0" '
                    f'dur="{dur}s" begin="{delay}" repeatCount="indefinite"/>')
        elif kk == "heat":
            anim = (f'<animateTransform attributeName="transform" type="translate" values="0 2;0 -3" dur="{dur}s" '
                    f'begin="{delay}" repeatCount="indefinite"/><animate attributeName="opacity" values="0.2;1;0.2" '
                    f'dur="{dur}s" begin="{delay}" repeatCount="indefinite"/>')
        elif kk == "stir":
            anim = (f'<animateTransform attributeName="transform" type="scale" values="1 1;0.3 1;1 1" dur="{dur / 2}s" '
                    f'repeatCount="indefinite" additive="sum"/>')
        else:
            anim = f'<animate attributeName="opacity" values="1;0.35;1" dur="{dur}s" begin="{delay}" repeatCount="indefinite"/>'
        return f'<{tag}{attrs}>{anim}</{tag}>'
    return _SELF.sub(sub, svg)
