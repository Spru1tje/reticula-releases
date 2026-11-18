"""v2.16: data for the 3D structure viewer (mfs/components/struct_view/viewer.html) and a static publication figure of a structure.

The viewer is a self-contained page (no external libraries, works offline): rotate, zoom, click an atom for its bonds, thermal spheres, labels, cell,
views along a / b / c, PNG export. This module builds the atoms (a supercell block with the atoms just outside the cell so that no bond is cut) and the bonds."""
from __future__ import annotations

import json
import math
import os

import numpy as np

from mfs.structmodel import is_bond, _covalent

HERE = os.path.dirname(os.path.abspath(__file__))
# Jmol-like colours
COLOURS = {"H": "#E8E8E8", "C": "#707070", "N": "#3050F8", "O": "#E8222A", "F": "#90E050", "P": "#FF8000", "S": "#FFD123", "Cl": "#1FF01F", "Br": "#A62929",
           "Zr": "#2BC4C4", "Hf": "#4DC2FF", "Ti": "#BFC2C7", "Zn": "#7D80B0", "Cu": "#C88033", "Fe": "#E06633", "Co": "#F090A0", "Ni": "#50D050", "Mn": "#9C7AC7",
           "Cr": "#8A99C7", "V": "#A6A6AB", "Al": "#BFA6A6", "Si": "#F0C8A0", "Ca": "#3DFF00", "Mg": "#8AFF00", "Na": "#AB5CF2", "K": "#8F40D4", "Ce": "#FFFFC7",
           "La": "#70D4FF", "Y": "#94FFFF", "B": "#FFB5B5", "Li": "#CC80FF", "Ba": "#00C900", "Sr": "#00FF00", "Cd": "#FFD98F", "Ag": "#C0C0C0", "Au": "#FFD123",
           "Pt": "#D0D0E0", "Pd": "#006985", "Ru": "#248F8F", "Rh": "#0A7D8C", "Mo": "#54B5B5", "W": "#2194D6", "Nb": "#73C2C9", "Ga": "#C28F8F", "In": "#A67573"}


def lattice_matrix(cell) -> np.ndarray:
    """Rows = a, b, c vectors (Å) for the cell (a, b, c, alpha, beta, gamma): a along x, b in the xy plane."""
    a, b, c, al, be, ga = cell
    al, be, ga = (math.radians(v) for v in (al, be, ga))
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return np.array([[a, 0, 0], [b * math.cos(ga), b * math.sin(ga), 0], [cx, cy, cz]], float)


def payload(model, p, cell, B=None, labels=None, repeat=(1, 1, 1), margin: float = 0.06, max_atoms: int = 4000, cn=None) -> dict:
    """Atoms, bonds and cell of the structure for the viewer. `model` = StructureModel, `p` = position parameters, `cell` = refined cell, `B` = B per atom of the
    model (or None). The block is repeat x the cell, plus every atom within `margin` (fractional) outside it, so bonds across the faces are complete."""
    L = lattice_matrix(cell)
    X = model.atom_positions(p) % 1.0
    X[np.abs(X - 1.0) < 1e-6] = 0.0
    labs = labels or model.labels()
    rep = tuple(int(v) for v in repeat)
    atoms, frac = [], []
    for i in range(rep[0] + 2):
        for j in range(rep[1] + 2):
            for k in range(rep[2] + 2):
                sh = np.array([i - 1, j - 1, k - 1], float)
                Y = X + sh
                keep = np.all((Y >= -margin) & (Y <= np.array(rep) + margin), axis=1)
                for a in np.nonzero(keep)[0]:
                    atoms.append((int(a), Y[a]))
    if len(atoms) > max_atoms:
        atoms = atoms[:max_atoms]
    cart = np.array([y for _, y in atoms]) @ L
    out_atoms = []
    for (a, y), c in zip(atoms, cart):
        el = model.atom_el[a]
        occ = float(model.atom_occ[a])
        u = None if B is None else float(B[a]) / (8 * math.pi ** 2)
        out_atoms.append({"x": float(c[0]), "y": float(c[1]), "z": float(c[2]), "fx": float(y[0]), "fy": float(y[1]), "fz": float(y[2]), "e": el,
                          "l": labs[int(model.atom_rep_site[a])], "r": 0.36 * _covalent(el) if el != "H" else 0.2,
                          "rt": (1.538 * math.sqrt(u) if u is not None and u > 0 else None), "B": None if B is None else float(B[a]), "occ": occ})
    # bonds between the displayed atoms (the chemical bond rule of the refinement)
    bonds = []
    n = len(out_atoms)
    pos = cart
    els = [a["e"] for a in out_atoms]
    for i in range(n):
        d = np.linalg.norm(pos[i + 1:] - pos[i], axis=1)
        for jj in np.nonzero(d < 3.2)[0]:
            j = i + 1 + int(jj)
            if is_bond(els[i], els[j], float(d[jj])):
                bonds.append([i, j, float(d[jj])])
    return {"atoms": out_atoms, "bonds": bonds, "cell": L.tolist(), "origin": [0.0, 0.0, 0.0] if rep == (1, 1, 1) else [0.0, 0.0, 0.0],
            "colors": {e: COLOURS.get(e, "#B0B0B0") for e in set(els)}, "repeat": list(rep), "cn": cn or {}}


def html(data: dict) -> str:
    """The viewer page with the data embedded."""
    with open(os.path.join(HERE, "components", "struct_view", "viewer.html"), encoding="utf-8") as fh:
        tpl = fh.read()
    return tpl.replace("__DATA__", json.dumps(data, separators=(",", ":")))


def show(st, data: dict, height: int = 560):
    import streamlit.components.v1 as components
    components.html(html(data), height=height, scrolling=False)


# --------------------------------------------------------------------------------------------------- static figure
def draw_structure(ax, data: dict, azim: float = 35.0, elev: float = 25.0, along: str | None = None, show_cell: bool = True, size: float = 1.0,
                   labels: bool = False, hide_h: bool = True):
    """Orthographic ball-and-stick picture (publication figure): depth-sorted spheres with shading, bonds, cell edges. `along` = 'a' | 'b' | 'c' views down that axis."""
    A = data["atoms"]
    L = np.array(data["cell"])
    pts = np.array([[a["x"], a["y"], a["z"]] for a in A])
    keep = np.array([not (hide_h and a["e"] == "H") for a in A])
    if along:
        z = L["abc".index(along)]
        z = z / np.linalg.norm(z)
        u = np.array([0, 0, 1.0]) if abs(z[2]) < 0.9 else np.array([0, 1.0, 0])
        x = np.cross(u, z)
        x /= np.linalg.norm(x)
        y = np.cross(z, x)
        Rm = np.array([x, y, z])
    else:
        az, el = math.radians(azim), math.radians(elev)
        Rz = np.array([[math.cos(az), -math.sin(az), 0], [math.sin(az), math.cos(az), 0], [0, 0, 1]])
        Rx = np.array([[1, 0, 0], [0, math.cos(el), -math.sin(el)], [0, math.sin(el), math.cos(el)]])
        Rm = Rx @ Rz
    cen = pts.mean(axis=0)
    Q = (pts - cen) @ Rm.T
    ext = float(np.max(np.linalg.norm(Q, axis=1))) or 1.0
    order = np.argsort(Q[:, 2])
    from matplotlib.patches import Circle
    cols = data["colors"]
    pt_per_ang = 0.0
    fig = ax.figure
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    ax.set_xlim(-ext * 1.05, ext * 1.05)
    ax.set_ylim(-ext * 1.05, ext * 1.05)
    ax.set_aspect("equal")
    ax.axis("off")
    if show_cell:
        o = np.array(data["origin"])
        corners = {(i, j, k): o + i * L[0] + j * L[1] + k * L[2] for i in (0, 1) for j in (0, 1) for k in (0, 1)}
        for (a, b) in [((0, 0, 0), (1, 0, 0)), ((0, 0, 0), (0, 1, 0)), ((0, 0, 0), (0, 0, 1)), ((1, 0, 0), (1, 1, 0)), ((1, 0, 0), (1, 0, 1)), ((0, 1, 0), (1, 1, 0)),
                       ((0, 1, 0), (0, 1, 1)), ((0, 0, 1), (1, 0, 1)), ((0, 0, 1), (0, 1, 1)), ((1, 1, 0), (1, 1, 1)), ((1, 0, 1), (1, 1, 1)), ((0, 1, 1), (1, 1, 1))]:
            p1, p2 = (corners[a] - cen) @ Rm.T, (corners[b] - cen) @ Rm.T
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], "-", color="#888888", lw=0.5, zorder=-1e6 + (p1[2] + p2[2]) / 2)
    bond_by_atom = {}
    for i, j, d in data["bonds"]:
        if keep[i] and keep[j]:
            ax.plot([Q[i, 0], Q[j, 0]], [Q[i, 1], Q[j, 1]], "-", color="#505050", lw=1.2 * size, solid_capstyle="round", zorder=(Q[i, 2] + Q[j, 2]) / 2 - 0.01)
    for i in order:
        if not keep[i]:
            continue
        a = A[i]
        r = a["r"] * size
        c = cols.get(a["e"], "#B0B0B0")
        ax.add_patch(Circle((Q[i, 0], Q[i, 1]), r, fc=c, ec="#202020", lw=0.4, zorder=Q[i, 2]))
        ax.add_patch(Circle((Q[i, 0] - 0.3 * r, Q[i, 1] + 0.3 * r), 0.35 * r, fc="white", ec="none", alpha=0.35, zorder=Q[i, 2] + 1e-3))
        if labels and a["l"]:
            ax.text(Q[i, 0] + r * 0.7, Q[i, 1] + r * 0.7, a["l"], fontsize=5, zorder=1e6)
    return {"extent": ext}
