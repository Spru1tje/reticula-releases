"""Reticon 2.19.0: ten material-science icons, five shape variants each (docs/icons/STYLE.md).

Each variant is a schematic drawing on the 64-unit grid, built with the same Icon builder as the rest of the library.
Roles: m = main outline, s = secondary (soft fill), a = accent (the one element the icon is about), x = cut (knocks out what was drawn before).
Known-answer counts (checked in tests/test_v2190_reticon.py): e.g. simple cubic has 8 corner atoms; face-centred adds 6 face atoms.
General practice: schematic icons for figures; they are not crystallographic drawings. Check against the source before you publish."""
from __future__ import annotations

import math

from .core import Icon

VARIANTS = ["base", "dense", "open", "stacked", "tilted"]


def _cube(ic, x, y, s, dx=14, dy=-10, edges=True, corners=True, role="m"):
    """Oblique cube: front face square at (x, y) with side s, depth vector (dx, dy)."""
    f = [(x, y), (x + s, y), (x + s, y + s), (x, y + s)]
    b = [(p[0] + dx, p[1] + dy) for p in f]
    if edges:
        for i in range(4):
            ic.line(*f[i], *f[(i + 1) % 4], w=2.2, role=role)
            ic.line(*b[i], *b[(i + 1) % 4], w=2.2, role=role)
            ic.line(*f[i], *b[i], w=2.2, role=role)
    pts = f + b
    if corners:
        for px, py in pts:
            ic.circle(px, py, 3.2, role="a")
    return f, b


def simple_cubic(v):
    ic = Icon("unit_cell_simple_cubic", "Simple cubic unit cell", "materials", ["unit cell", "simple cubic", "lattice"])
    if v == "base":
        _cube(ic, 14, 22, 26)
    elif v == "dense":
        f, b = _cube(ic, 14, 22, 26)
        for px, py in [(27, 35), (27, 22), (14, 35)]:
            ic.circle(px, py, 2.0, role="s")
    elif v == "open":
        _cube(ic, 14, 22, 26, corners=False)
        for px, py in [(14, 22), (40, 22), (40, 48), (14, 48), (54, 12), (27, 9)]:
            ic.circle(px, py, 3.2, role="a")
    elif v == "stacked":
        _cube(ic, 8, 26, 22, dx=10, dy=-8)
        _cube(ic, 22, 18, 22, dx=10, dy=-8, role="s")
    else:  # tilted: the cube turned, drawn with a larger depth
        _cube(ic, 12, 24, 24, dx=-12, dy=-12)
    return ic


def face_centred(v):
    ic = Icon("unit_cell_face_centred", "Face-centred cubic unit cell", "materials", ["unit cell", "fcc", "face centred", "lattice"])
    f, b = _cube(ic, 14, 22, 26, corners=True)
    if v in ("base", "dense", "open"):
        centres = [(27, 35), (27, 22), (14, 35), (40, 35), (27, 48), (27, 35)]
        for px, py in centres[:5] if v == "base" else centres:
            ic.circle(px, py, 2.4, role="s")
    elif v == "stacked":
        _cube(ic, 22, 16, 26, dx=12, dy=-9, role="s")
    else:
        for px, py in [(27, 35), (27, 22), (14, 35), (40, 35), (27, 48)]:
            ic.circle(px, py, 2.4, role="s")
    return ic


def nanowire(v):
    ic = Icon("nanowire", "Nanowire", "materials", ["nanowire", "nanorod", "one-dimensional"])
    if v == "base":
        for i, x in enumerate([24, 34, 44]):
            ic.srect(x, 14, 5, 36, r=2.5, sw=2.2, role="m")
    elif v == "dense":
        for i in range(5):
            ic.srect(12 + i * 8, 14, 5, 36, r=2.5, sw=2.2, role="m")
    elif v == "open":
        ic.srect(28, 12, 6, 40, r=3, sw=2.2, role="m")
        ic.circle(31, 12, 3, role="a")
    elif v == "stacked":
        for i in range(4):
            ic.srect(16 + i * 8, 14 + i * 4, 5, 34 - i * 6, r=2.5, sw=2.2, role="m")
    else:
        ic.line(12, 48, 52, 20, w=6, role="m")
        ic.circle(52, 20, 4, role="a")
    return ic


def multiwall_nanotube(v):
    ic = Icon("nanotube_multiwall", "Multi-walled nanotube", "materials", ["nanotube", "carbon nanotube", "multi-walled", "CNT"])
    ic.circle(32, 32, 20, role="m")
    if v == "base":
        ic.circle(32, 32, 14, role="m"); ic.circle(32, 32, 8, role="m")
    elif v == "dense":
        ic.circle(32, 32, 14, role="m"); ic.circle(32, 32, 8, role="m"); ic.circle(32, 32, 3, role="a")
    elif v == "open":
        ic.cut().circle(32, 32, 9, role="x")
        ic.circle(32, 32, 14, role="m")
    elif v == "stacked":
        ic.circle(32, 32, 14, role="m"); ic.circle(32, 32, 8, role="m")
        ic.circle(32, 32, 20, role="s")
    else:
        ic.ellipse(32, 32, 20, 12, role="m"); ic.ellipse(32, 32, 12, 7, role="m")
    return ic


def two_d_flake(v):
    ic = Icon("two_dimensional_flake", "Two-dimensional flake", "materials", ["2D", "nanosheet", "flake", "graphene", "MXene"])
    def hexagon(cx, cy, r):
        return [(cx + r * math.cos(math.radians(60 * k + 30)), cy + r * math.sin(math.radians(60 * k + 30))) for k in range(6)]
    if v == "base":
        ic.poly(hexagon(32, 32, 20), role="m")
    elif v == "dense":
        ic.poly(hexagon(32, 32, 20), role="m"); ic.poly(hexagon(32, 32, 9), role="a")
    elif v == "open":
        for cx, cy in [(22, 26), (42, 26), (32, 40)]:
            ic.poly(hexagon(cx, cy, 9), role="m")
    elif v == "stacked":
        for k in range(3):
            ic.poly(hexagon(30 + k * 3, 30 + k * 5, 17), role="s" if k < 2 else "m")
    else:
        ic.ellipse(32, 32, 20, 9, role="m")
        ic.ellipse(32, 32, 20, 9, role="s")
    return ic


def mesoporous(v):
    ic = Icon("mesoporous_material", "Mesoporous material", "materials", ["mesoporous", "pores", "2-50 nm", "silica"])
    ic.circle(32, 32, 24, role="s")
    ic.circle(32, 32, 24, role="m")
    grid = [(x, y) for x in (18, 28, 38, 48) for y in (18, 28, 38, 48)]
    r = {"base": 3.4, "dense": 4.0, "open": 2.6, "stacked": 3.4, "tilted": 3.4}[v]
    pts = grid if v != "open" else grid[::2]
    for x, y in pts:
        ic.cut().circle(x, y, r, role="x")
    return ic


def open_foam(v):
    ic = Icon("open_cell_foam", "Open-cell foam", "materials", ["foam", "open cell", "porous", "scaffold"])
    cells = [(20, 20), (36, 20), (28, 34), (44, 34), (20, 46), (36, 46)] if v != "dense" else [(16 + 12 * i, 18 + 12 * j) for i in range(4) for j in range(4)]
    for cx, cy in cells[:6 if v != "dense" else 16]:
        ic.poly([(cx + 7 * math.cos(math.radians(60 * k + 30)), cy + 7 * math.sin(math.radians(60 * k + 30))) for k in range(6)],
                role="m" if v != "stacked" else "s")
    if v == "open":
        ic.line(8, 56, 56, 56, w=1.8, role="m")
    return ic


def slit_pore(v):
    ic = Icon("slit_pore", "Slit pore", "materials", ["slit pore", "micropore", "2D pore", "adsorption"])
    gap = {"base": 16, "dense": 10, "open": 22, "stacked": 14, "tilted": 16}[v]
    top = 32 - gap / 2
    ic.srect(10, 10, 44, top - 10, r=1, sw=2, role="m")
    ic.srect(10, top + gap, 44, top - 10, r=1, sw=2, role="m")
    for x in range(16, 52, 9):
        ic.circle(x, 32, 2.2, role="a")
    return ic


def cylindrical_pore(v):
    ic = Icon("cylindrical_pore", "Cylindrical pore", "materials", ["cylindrical pore", "pore", "channel", "mesopore"])
    ic.srect(8, 16, 48, 32, r=3, sw=2.2, role="s")
    rr = {"base": 6, "dense": 4, "open": 8, "stacked": 6, "tilted": 6}[v]
    for x in ([22, 42] if v != "dense" else [16, 26, 36, 46]):
        ic.cut().circle(x, 32, rr, role="x")
        ic.circle(x, 32, rr, role="m")
    return ic


def core_shell(v):
    ic = Icon("core_shell_particle", "Core-shell particle", "materials", ["core shell", "core-shell", "nanoparticle", "coating"])
    ic.circle(32, 32, 22, role="s")
    ic.circle(32, 32, 22, role="m")
    rr = {"base": 12, "dense": 9, "open": 14, "stacked": 12, "tilted": 12}[v]
    ic.circle(32, 32, rr, role="a")
    if v == "stacked":
        ic.circle(32, 32, 16, role="m")
    return ic


BUILDERS = {
    "unit_cell_simple_cubic": simple_cubic,
    "unit_cell_face_centred": face_centred,
    "nanowire": nanowire,
    "nanotube_multiwall": multiwall_nanotube,
    "two_dimensional_flake": two_d_flake,
    "mesoporous_material": mesoporous,
    "open_cell_foam": open_foam,
    "slit_pore": slit_pore,
    "cylindrical_pore": cylindrical_pore,
    "core_shell_particle": core_shell,
}


def build_all() -> dict:
    """{icon name: {variant: spec dict}} with spec = Icon.to_spec()."""
    out = {}
    for name, fn in BUILDERS.items():
        out[name] = {v: fn(v).to_spec() for v in VARIANTS}
    return out
