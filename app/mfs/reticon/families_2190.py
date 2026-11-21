"""Reticon 2.19.0: drawing families for the material-science catalogue. Each family draws one schematic motif in five shape
variants (base, dense, open, stacked, tilted). A term picks its family by keywords (family_for) and gets a small, term-specific
offset so that no two icons are identical. Roles: m outline, s soft fill, a accent (the one element), x cut.
General practice: these are schematic figure icons, not crystallographic drawings; check the source before publication."""
from __future__ import annotations

import math
import zlib

from .core import Icon

VARIANTS = ["base", "dense", "open", "stacked", "tilted"]


def _offset(name: str, span: float = 4.0) -> tuple[float, float]:
    h = zlib.crc32(name.encode("utf-8"))
    return ((h % 1000) / 1000 - 0.5) * span, (((h // 1000) % 1000) / 1000 - 0.5) * span


def _hex(cx, cy, r, rot=30):
    return [(cx + r * math.cos(math.radians(60 * k + rot)), cy + r * math.sin(math.radians(60 * k + rot))) for k in range(6)]


def _cube(ic, x, y, s, dx=14, dy=-10, role="m", corners=True):
    f = [(x, y), (x + s, y), (x + s, y + s), (x, y + s)]
    b = [(p[0] + dx, p[1] + dy) for p in f]
    for i in range(4):
        ic.line(*f[i], *f[(i + 1) % 4], w=2.2, role=role)
        ic.line(*b[i], *b[(i + 1) % 4], w=2.2, role=role)
        ic.line(*f[i], *b[i], w=2.2, role=role)
    if corners:
        for px, py in f + b:
            ic.circle(px, py, 3.0, role="a")
    return f, b


def fam_cube(ic, v, n):
    # crystal cells and cubic structures: the cube with atoms, variants change density and arrangement
    dx, dy = _offset(n, 6)
    if v == "base":
        _cube(ic, 14 + dx, 22 + dy, 26)
    elif v == "dense":
        _cube(ic, 14 + dx, 22 + dy, 26)
        for px, py in [(27, 35), (27, 22), (14, 35), (40, 35), (27, 48)]:
            ic.circle(px + dx, py + dy, 2.2, role="s")
    elif v == "open":
        _cube(ic, 14 + dx, 22 + dy, 26, corners=False)
        for px, py in [(14, 22), (40, 48), (54, 12)]:
            ic.circle(px + dx, py + dy, 3.0, role="a")
    elif v == "stacked":
        _cube(ic, 8 + dx, 26 + dy, 22, dx=10, dy=-8)
        _cube(ic, 22 + dx, 18 + dy, 22, dx=10, dy=-8, role="s")
    else:
        _cube(ic, 12 + dx, 24 + dy, 24, dx=-12, dy=-12)


def fam_spheres(ic, v, n):
    # clusters and single atoms: a group of spheres
    dx, dy = _offset(n, 6)
    pts = [(20, 22), (34, 18), (44, 30), (28, 38), (40, 46), (18, 44)]
    k = {"base": 6, "dense": 6, "open": 4, "stacked": 6, "tilted": 5}[v]
    for i, (x, y) in enumerate(pts[:k]):
        ic.circle(x + dx, y + dy, 6 if i == 0 else 4.5, role="a" if i == 0 else "s")
    if v == "stacked":
        ic.circle(32 + dx, 32 + dy, 8, role="m")


def fam_rods(ic, v, n):
    dx, dy = _offset(n, 6)
    cnt = {"base": 3, "dense": 5, "open": 1, "stacked": 4, "tilted": 2}[v]
    for i in range(cnt):
        x = 14 + i * (36 / max(cnt, 1)) + dx
        if v == "tilted":
            ic.line(x, 50 + dy, x + 14, 14 + dy, w=6, role="m")
        else:
            ic.srect(x, 14 + dy, 6, 36, r=3, sw=2.2, role="m" if i % 2 == 0 else "s")
    ic.circle(32 + dx, 56 + dy, 2.5, role="a")


def fam_tube(ic, v, n):
    dx, dy = _offset(n, 6)
    ic.circle(32 + dx, 32 + dy, 20, role="m")
    if v == "open":
        ic.cut().circle(32 + dx, 32 + dy, 9, role="x")
    ic.circle(32 + dx, 32 + dy, 14, role="m")
    if v != "open":
        ic.circle(32 + dx, 32 + dy, 8, role="m")
    if v == "dense":
        ic.circle(32 + dx, 32 + dy, 3, role="a")
    if v == "stacked":
        ic.circle(32 + dx, 32 + dy, 20, role="s")
    if v == "tilted":
        ic.ellipse(32 + dx, 32 + dy, 20, 11, role="m")


def fam_sheet(ic, v, n):
    dx, dy = _offset(n, 6)
    if v == "base":
        ic.poly(_hex(32 + dx, 32 + dy, 20), role="m")
    elif v == "dense":
        ic.poly(_hex(32 + dx, 32 + dy, 20), role="m"); ic.poly(_hex(32 + dx, 32 + dy, 9), role="a")
    elif v == "open":
        for cx, cy in [(22, 26), (42, 26), (32, 40)]:
            ic.poly(_hex(cx + dx, cy + dy, 9), role="m")
    elif v == "stacked":
        for k in range(3):
            ic.poly(_hex(30 + k * 3 + dx, 30 + k * 5 + dy, 17), role="s" if k < 2 else "m")
    else:
        ic.ellipse(32 + dx, 32 + dy, 20, 9, role="m"); ic.ellipse(32 + dx, 32 + dy, 20, 9, role="s")


def fam_pores(ic, v, n):
    dx, dy = _offset(n, 4)
    ic.circle(32 + dx, 32 + dy, 24, role="s"); ic.circle(32 + dx, 32 + dy, 24, role="m")
    r = {"base": 3.4, "dense": 4.0, "open": 2.6, "stacked": 3.4, "tilted": 3.4}[v]
    grid = [(x, y) for x in (18, 28, 38, 48) for y in (18, 28, 38, 48)]
    for x, y in (grid if v != "open" else grid[::2]):
        ic.cut().circle(x + dx, y + dy, r, role="x")


def fam_channels(ic, v, n):
    dx, dy = _offset(n, 4)
    gap = {"base": 16, "dense": 10, "open": 22, "stacked": 14, "tilted": 16}[v]
    top = 32 - gap / 2 + dy
    ic.srect(10 + dx, 10, 44, top - 10, r=1, sw=2, role="m")
    ic.srect(10 + dx, top + gap, 44, top - 10, r=1, sw=2, role="m")
    for x in range(16, 52, 9):
        ic.circle(x + dx, 32 + dy, 2.2, role="a")


def fam_foam(ic, v, n):
    dx, dy = _offset(n, 4)
    cells = [(20, 20), (36, 20), (28, 34), (44, 34), (20, 46), (36, 46)]
    if v == "dense":
        cells = [(16 + 12 * i, 18 + 12 * j) for i in range(4) for j in range(4)]
    for cx, cy in cells[: (16 if v == "dense" else 6)]:
        ic.poly(_hex(cx + dx, cy + dy, 7), role="s" if v == "stacked" else "m")
    if v == "open":
        ic.line(8, 56, 56, 56, w=1.8, role="m")


def fam_shell(ic, v, n):
    dx, dy = _offset(n, 4)
    ic.circle(32 + dx, 32 + dy, 22, role="s"); ic.circle(32 + dx, 32 + dy, 22, role="m")
    rr = {"base": 12, "dense": 9, "open": 14, "stacked": 12, "tilted": 12}[v]
    ic.circle(32 + dx, 32 + dy, rr, role="a")
    if v == "stacked":
        ic.circle(32 + dx, 32 + dy, 16, role="m")
    if v == "open":
        ic.cut().circle(32 + dx, 32 + dy, 6, role="x")


def fam_polyhedron(ic, v, n):
    dx, dy = _offset(n, 4)
    pts = [(32, 8), (54, 32), (32, 56), (10, 32)]
    ic.poly([(x + dx, y + dy) for x, y in pts], role="m")
    ic.line(32 + dx, 8 + dy, 32 + dx, 56 + dy, w=1.8, role="m")
    ic.line(10 + dx, 32 + dy, 54 + dx, 32 + dy, w=1.8, role="m")
    if v in ("dense", "stacked"):
        ic.poly([(x + dx, y + dy) for x, y in [(32, 16), (46, 32), (32, 48), (18, 32)]], role="a")
    if v == "tilted":
        ic.poly([(x + dx + 6, y + dy - 4) for x, y in pts], role="s")


def fam_film(ic, v, n):
    dx, dy = _offset(n, 4)
    ic.poly([(8 + dx, 44), (48 + dx, 44), (56 + dx, 36), (16 + dx, 36)], role="m")
    thick = {"base": 6, "dense": 9, "open": 3, "stacked": 12, "tilted": 6}[v]
    ic.srect(16 + dx, 36 - thick + dy, 40, thick, r=1, sw=2, role="s")
    for x in range(20, 52, 8):
        ic.circle(x + dx, 34 - thick / 2 + dy, 2.4, role="a")


def fam_curve(ic, v, n):
    dx, dy = _offset(n, 4)
    ic.line(8, 54, 56, 54, w=1.8, role="m"); ic.line(8, 54, 8, 8, w=1.8, role="m")
    peaks = {"base": 1, "dense": 3, "open": 1, "stacked": 2, "tilted": 1}[v]
    pts = []
    for i in range(0, 49):
        x = 8 + i
        y = 50
        for p in range(peaks):
            cx = 16 + p * (34 / max(peaks, 1)) + dx
            y -= 30 / (1 + ((x - cx) / 3.5) ** 2) / peaks ** 0.5
        pts.append((x, y + dy))
    ic.spoly(pts, w=2.4, role="a")


def fam_chain(ic, v, n):
    dx, dy = _offset(n, 4)
    pts = [(10, 40), (20, 22), (30, 36), (40, 20), (50, 34)]
    k = {"base": 5, "dense": 5, "open": 3, "stacked": 5, "tilted": 4}[v]
    ic.spoly([(x + dx, y + dy) for x, y in pts[:k]], w=3.2, role="s")
    for x, y in pts[:k]:
        ic.circle(x + dx, y + dy, 3.6, role="a" if (x, y) == pts[0] else "m")


def fam_layers(ic, v, n):
    dx, dy = _offset(n, 4)
    count = {"base": 3, "dense": 5, "open": 2, "stacked": 4, "tilted": 3}[v]
    for i in range(count):
        y = 16 + i * (36 / count) + dy
        ic.poly([(10 + dx, y), (50 + dx, y), (56 + dx, y + 6), (16 + dx, y + 6)], role="a" if i == count // 2 else "m")


def fam_wave(ic, v, n):
    dx, dy = _offset(n, 4)
    f = {"base": 2, "dense": 4, "open": 1, "stacked": 3, "tilted": 2}[v]
    pts = [(8 + i, 32 + 10 * math.sin(f * (8 + i) / 9 + dx) + dy) for i in range(49)]
    ic.spoly(pts, w=2.2, role="a")
    ic.line(8, 52, 56, 52, w=1.5, role="m")


def fam_vessel(ic, v, n):
    dx, dy = _offset(n, 4)
    ic.srect(14 + dx, 12, 36, 42, r=6, sw=2.2, role="m")
    fill = {"base": 18, "dense": 28, "open": 10, "stacked": 22, "tilted": 18}[v]
    ic.srect(16 + dx, 54 - fill, 32, fill, r=2, sw=0, role="s")
    ic.line(20 + dx, 12, 44 + dx, 12, w=2.2, role="m")


def fam_grains(ic, v, n):
    dx, dy = _offset(n, 4)
    k = {"base": 7, "dense": 11, "open": 4, "stacked": 9, "tilted": 6}[v]
    for i in range(k):
        ang = 2 * math.pi * i / k
        cx = 32 + 16 * math.cos(ang) * (0.5 + (i % 3) / 4) + dx
        cy = 32 + 16 * math.sin(ang) * (0.5 + (i % 2) / 3) + dy
        ic.poly(_hex(cx, cy, 6 + (i % 2) * 2, rot=15 * i), role="a" if i == 0 else ("s" if i % 2 else "m"))


def fam_defect(ic, v, n):
    dx, dy = _offset(n, 4)
    for x in range(14, 52, 12):
        for y in range(14, 52, 12):
            if v == "open" and (x, y) == (26, 26):
                continue
            ic.circle(x + dx, y + dy, 3.0, role="m")
    ic.circle(26 + dx, 26 + dy, 4.2, role="a")


def fam_device(ic, v, n):
    dx, dy = _offset(n, 4)
    ic.srect(10 + dx, 36, 44, 14, r=1.5, sw=2.2, role="m")
    ic.srect(20 + dx, 22, 24, 14, r=1.5, sw=2.2, role="s")
    ic.line(32 + dx, 22 + dy, 32 + dx, 10 + dy, w=2.2, role="a")
    if v in ("dense", "stacked"):
        ic.srect(14 + dx, 26, 36, 6, r=1, sw=1.8, role="a")


def fam_droplet(ic, v, n):
    dx, dy = _offset(n, 4)
    ic.circle(32 + dx, 36 + dy, 16 if v != "open" else 10, role="s")
    ic.circle(32 + dx, 36 + dy, 16 if v != "open" else 10, role="m")
    k = {"base": 1, "dense": 4, "open": 0, "stacked": 2, "tilted": 3}[v]
    for i in range(k):
        ic.circle(20 + i * 8 + dx, 14 + (i % 2) * 6 + dy, 2.6, role="a")


FAMILIES = {"cube": fam_cube, "spheres": fam_spheres, "rods": fam_rods, "tube": fam_tube, "sheet": fam_sheet,
            "pores": fam_pores, "channels": fam_channels, "foam": fam_foam, "shell": fam_shell,
            "polyhedron": fam_polyhedron, "film": fam_film, "curve": fam_curve, "chain": fam_chain,
            "layers": fam_layers, "wave": fam_wave, "vessel": fam_vessel, "grains": fam_grains,
            "defect": fam_defect, "device": fam_device, "droplet": fam_droplet}

# keyword rules: first match wins
RULES = [
    (("unit_cell", "lattice", "cubic", "crystal_structure", "fluorite", "rock_salt", "spinel", "perovskite", "wurtzite", "zinc_blende", "rutile",
      "diamond", "gyroid", "kagome", "honeycomb_lattice", "lieb", "chiral", "structure"), "cube"),
    (("nanowire", "rod", "fiber", "fibre", "wire", "ribbon", "pillar", "column"), "rods"),
    (("foam", "aerogel", "xerogel", "cryogel", "sponge", "scaffold", "honeycomb", "mat", "textile", "woven"), "foam"),
    (("nanotube", "tube", "nanofiber"), "tube"),
    (("sheet", "flake", "mxene", "graphene", "boron_nitride", "dichalcogenide", "phosphorus", "nanoplate", "nanodisc", "layer_double", "clay", "mica", "kaolin"), "sheet"),
    (("pore", "porous", "mesopor", "micropor", "macropor", "adsorb", "isotherm", "sorption", "bottleneck", "ink_bottle", "channel"), "pores"),
    (("slit", "cylindrical", "diffusion", "permeation", "transport"), "channels"),
    (("core", "shell", "yolk", "hollow", "capsule", "micelle", "vesicle", "liposome", "bead", "particle", "colloid", "microsphere", "dot", "nanoflower", "pickering"), "shell"),
    (("octahedr", "dodecahedr", "polyhedr", "habit", "cubic_nano", "pyramid", "truncated", "crystal_growth", "nucleation"), "polyhedron"),
    (("film", "coating", "deposition", "sputter", "layer_by", "langmuir", "dip", "electroplat", "anodis", "passive", "oxide_scale", "pvd", "spray", "monolayer"), "film"),
    (("spectrum", "spectra", "peak", "curve", "plot", "tauc", "kubelka", "voltammogram", "impedance", "tafel", "isotherm_type", "kinetic",
      "emission", "absorption", "edge", "hysteresis", "conversion", "loop", "exafs", "xanes", "epr", "mossbauer", "nmr", "diffraction", "scatter", "magnetization", "hysteresis_magnetic",
      "heat_capacity", "thermogravimetric", "differential", "dilatometry", "stress", "density_of_states", "band_structure", "dispersion",
      "refractive", "birefringence", "reflectance", "fermi", "bragg", "scherrer", "pair_distribution", "photoelectron"), "curve"),
    (("chain", "polymer", "brush", "block", "crosslink", "hydrogel", "linker", "ligand", "cyclodextrin", "crown", "supramolecular", "switch",
      "molecule", "complex", "fullerene", "c60", "nitride", "organic", "cage", "coordination"), "chain"),
    (("layer", "stack", "van_der", "heterostructure", "junction", "bilayer", "lamella", "block"), "layers"),
    (("wave", "phonon", "thermal", "conduct", "photon", "plasmon", "resonance", "raman", "infrared", "electron", "neutron", "x_ray"), "wave"),
    (("reactor", "autoclave", "furnace", "synthesis", "sinter", "ball", "milling", "calcin", "press", "vapour", "vapor", "process", "sol_gel",
      "hydrothermal", "solvothermal", "microwave", "drying", "electrospin", "melt", "freeze"), "vessel"),
    (("grain", "twin", "boundary", "domain", "phase", "eutectic", "solid_solution", "texture", "mosaic", "strain", "fracture", "crack", "creep",
      "fatigue", "hardness", "tensile", "yield", "work_hardening", "corrosion", "elastic", "young", "poisson", "shape_memory"), "grains"),
    (("defect", "vacancy", "interstitial", "dopant", "antisite", "frenkel", "schottky_defect", "dislocation", "stacking_fault"), "defect"),
    (("battery", "electrode", "cell", "anode", "cathode", "electrolyte", "separator", "collector", "supercapacitor", "capacitor", "sensor",
      "device", "memristor", "resistive", "junction", "contact", "transistor", "dielectric", "spin_valve", "tunnel", "solar", "thermoelectric",
      "piezo", "ferro", "magnet", "skyrmion", "domain_wall", "vortex", "water_splitting", "evolution", "reduction", "photocatal", "exciton"), "device"),
    (("droplet", "emulsion", "aerosol", "ferrofluid", "dust", "droplet", "suspension", "solution"), "droplet"),
]


def family_for(name: str) -> str:
    for keys, fam in RULES:
        if any(k in name for k in keys):
            return fam
    return "spheres"


def draw(name: str, v: str, fam: str | None = None) -> Icon:
    fam = fam or family_for(name)
    label = name.replace("_", " ").capitalize()
    ic = Icon(name, label, "materials", name.split("_"))
    FAMILIES[fam](ic, v, name)
    return ic
