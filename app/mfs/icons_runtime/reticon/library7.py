"""
v5 batch 3: nanomaterials & photophysics — quantum dots, PL, nanoparticle morphologies
(core-shell, yolk-shell, Janus, hollow, nanorods/cubes/stars/plates/cages), MOF & COF motifs,
composites and catalysis. All original geometry (CC0).
"""
from __future__ import annotations

import math
import random

from .core import Op, polar, smooth_path, lorentz_trace, _f
from .library import icon, axes
from .library2 import rot_rect, gauss_trace
from .library4 import arrowhead, hexagon

EMIT = ["#3B5BDB", "#15AABF", "#40C057", "#FAB005", "#FD7E14", "#E03131"]   # blue → red
GOLD, SILVER, NIR = "#E0A526", "#B9C2CC", "#9C1C1C"


def ligands(ic, cx, cy, r0, r1, n=12, w=1.8, role="a", dot=1.6, start=15):
    for k in range(n):
        a = start + k * 360 / n
        x1, y1 = polar(cx, cy, r0, a)
        x2, y2 = polar(cx, cy, r1, a)
        ic.line(x1, y1, x2, y2, w, role=role)
        if dot:
            ic.circle(x2, y2, dot, role=role)


def glow(ic, cx, cy, r, colour, rings=3):
    for k in range(rings, 0, -1):
        ic.ops.append(Op("a", "circle", {"cx": cx, "cy": cy, "r": r + k * 4, "opacity": 0.12 + 0.06 * (rings - k),
                                         "_fill": colour}))


# ==========================================================================
# QUANTUM DOTS & PHOTOPHYSICS
# ==========================================================================

@icon("qd_core_shell", "materials", label="Core–Shell Quantum Dot", tags=["quantum dot", "CdSe/ZnS", "core-shell", "QD", "ligands"])
def _(ic):
    ligands(ic, 32, 32, 19, 27)
    ic.circle(32, 32, 19)
    ic.cut().circle(32, 32, 14)
    ic.sec().circle(32, 32, 12)
    ic.cut().circle(28, 28, 2.6)


@icon("qd_size_series", "materials", label="QD Size → Emission", tags=["quantum confinement", "size-dependent emission", "quantum dots", "band gap"])
def _(ic):
    for k, (x, r) in enumerate(((7, 3.2), (17, 4.2), (29, 5.4), (43, 6.6), (57, 7.8))):
        glow(ic, x, 32, r, EMIT[[0, 1, 2, 4, 5][k]], rings=1)
        ic.circle(x, 32, r).fixed(EMIT[[0, 1, 2, 4, 5][k]])
    ic.line(6, 50, 56, 50, 1.8)
    arrowhead(ic, 59, 50, 0, 3.5, "m")


@icon("qd_emission", "materials", label="Emitting Quantum Dot", tags=["photoluminescence", "emission", "quantum dot", "glow"])
def _(ic):
    glow(ic, 36, 34, 12, EMIT[2])
    ic.circle(36, 34, 12).fixed(EMIT[2])
    ic.cut().circle(32, 30, 2.6)
    ic.acc().spath("M4 10 q3 -4 6 0 t6 0 t6 0", 2)
    arrowhead(ic, 24, 12, 30, 3.5)


@icon("carbon_dots", "materials", label="Carbon Dots", tags=["carbon dots", "CDs", "graphene quantum dots", "fluorescent carbon"])
def _(ic):
    ic.circle(32, 32, 20)
    for cx, cy in ((26, 28), (36, 28), (31, 37), (41, 37), (21, 37)):
        hexagon(ic, cx, cy, 5, 1.4, "s", rot=0)
    for a in range(0, 360, 40):
        x2, y2 = polar(32, 32, 26, a)
        x1, y1 = polar(32, 32, 20, a)
        ic.acc().line(x1, y1, x2, y2, 1.6)
        ic.acc().circle(x2, y2, 2.2) if a % 80 == 0 else ic.sec().circle(x2, y2, 2)


@icon("perovskite_nanocrystal", "materials", label="Perovskite Nanocrystal", tags=["CsPbBr3", "perovskite", "nanocube", "nanocrystal"])
def _(ic):
    for k in range(12):
        x, y = (10 + k * 4, 12) if k < 6 else (52, 16 + (k - 6) * 6)
        ic.acc().line(x, y, x - 3, y - 6, 1.4) if k < 6 else ic.acc().line(x, y, x + 6, y - 2, 1.4)
    ic.sec().poly([(10, 24), (40, 24), (40, 54), (10, 54)])
    ic.poly([(10, 24), (22, 14), (52, 14), (40, 24)])
    ic.ops[-1].role = "s"
    ic.sec().poly([(40, 24), (52, 14), (52, 44), (40, 54)])
    ic.cut().spoly([(10, 24), (40, 24), (52, 14)], 1)
    ic.cut().line(40, 24, 40, 54, 1)


@icon("perovskite_structure", "materials", label="Perovskite ABX₃", tags=["perovskite", "ABX3", "octahedra", "crystal structure"])
def _(ic):
    for (cx, cy) in ((16, 16), (48, 16), (16, 48), (48, 48)):
        ic.poly([(cx, cy - 11), (cx + 11, cy), (cx, cy + 11), (cx - 11, cy)], role="s")
        ic.circle(cx, cy, 2.6)
    ic.acc().circle(32, 32, 6.5)


@icon("upconversion_np", "materials", label="Upconversion NP", tags=["upconversion", "UCNP", "NaYF4", "NIR to visible", "lanthanide"])
def _(ic):
    ic.circle(32, 34, 14)
    ic.cut().circle(32, 34, 10)
    ic.sec().circle(32, 34, 8)
    for k, y in enumerate((22, 30)):
        ic.spath(f"M2 {y} q2.5 -4 5 0 t5 0 t5 0", 2).fixed(NIR)
    for k, (a, col) in enumerate(((-40, EMIT[2]), (0, EMIT[0]), (40, EMIT[5]))):
        x1, y1 = polar(32, 34, 16, a)
        x2, y2 = polar(32, 34, 26, a)
        ic.line(x1, y1, x2, y2, 2.2).fixed(col)


@icon("jablonski_diagram", "data", label="Jablonski Diagram", tags=["Jablonski", "fluorescence", "phosphorescence", "ISC", "excited states"])
def _(ic):
    for y, x0, x1 in ((54, 6, 34), (48, 6, 34), (18, 6, 34), (13, 6, 34), (30, 40, 60), (26, 40, 60)):
        ic.line(x0, y, x1, y, 2 if y in (54, 18, 30) else 1)
    ic.line(12, 53, 12, 20, 2.2).fixed(EMIT[0])
    arrowhead(ic, 12, 18, -90, 3.5)
    ic.ops[-1].attrs["_fill"] = EMIT[0]
    ic.line(24, 19, 24, 52, 2.2).fixed(EMIT[2])
    arrowhead(ic, 24, 54, 90, 3.5)
    ic.ops[-1].attrs["_fill"] = EMIT[2]
    ic.spath("M30 19 C36 19 40 24 44 29", 1.6)
    ic.ops[-1].attrs["stroke-dasharray"] = "2.5 2"
    ic.line(52, 31, 52, 52, 2.2).fixed(EMIT[5])
    arrowhead(ic, 52, 54, 90, 3.5)
    ic.ops[-1].attrs["_fill"] = EMIT[5]


@icon("stokes_shift", "data", label="Stokes Shift", tags=["Stokes shift", "absorption", "emission", "spectra"])
def _(ic):
    axes(ic)
    ic.spoly(gauss_trace(12, 56, 54, [(26, 34, 5)], n=120), 2.4).fixed(EMIT[0])
    ic.spoly(gauss_trace(12, 56, 54, [(40, 30, 5)], n=120), 2.4).fixed(EMIT[4])
    ic.line(26, 14, 40, 14, 1.4)
    arrowhead(ic, 42, 14, 0, 3, "m")


@icon("pl_decay", "data", label="PL Lifetime Decay", tags=["TCSPC", "lifetime", "decay", "time-resolved PL"])
def _(ic):
    axes(ic)
    rnd = random.Random(2)
    for i in range(0, 44, 2):
        y = 12 + 40 * (1 - math.exp(-i / 7))
        ic.sec().circle(12 + i, y + rnd.uniform(-1.2, 1.2), 1.1)
    ic.acc().spoly([(12 + i, 12 + 40 * (1 - math.exp(-i / 7))) for i in range(0, 44)], 2.2)
    ic.spath("M12 12 L13 10 L14 14", 1.2)


@icon("fret_pair", "materials", label="FRET / Energy Transfer", tags=["FRET", "energy transfer", "donor acceptor", "sensing"])
def _(ic):
    ic.circle(16, 36, 10).fixed(EMIT[1])
    ic.circle(48, 36, 10).fixed(EMIT[4])
    ic.spath("M26 32 C32 26 34 26 38 32", 2).fixed("#868E96")
    arrowhead(ic, 39, 34, 55, 3.5, "m")
    ic.spath("M6 10 q2 -4 4 0 t4 0", 1.8).fixed(EMIT[0])
    arrowhead(ic, 15, 16, 60, 3)
    ic.ops[-1].attrs["_fill"] = EMIT[0]
    ic.spath("M52 22 q3 -4 6 -1", 1.8).fixed(EMIT[5])


@icon("luminescence_thermometry", "materials", label="Luminescence Thermometry", tags=["thermometry", "ratiometric", "Er3+", "temperature sensing", "PL ratio"])
def _(ic):
    ic.spath("M8 6 V42", 2)
    ic.circle(8, 48, 5).fixed(EMIT[5])
    ic.rect(6.5, 20, 3, 26, 1.5).fixed(EMIT[5])
    ic.line(16, 54, 60, 54, 2)
    ic.spoly(gauss_trace(18, 58, 54, [(30, 30, 2.4)], n=80), 2.2).fixed(EMIT[2])
    ic.spoly(gauss_trace(18, 58, 54, [(44, 18, 2.8)], n=80), 2.2).fixed(EMIT[4])


@icon("led_phosphor", "materials", label="Phosphor-Converted LED", tags=["LED", "phosphor", "white light", "wLED"])
def _(ic):
    ic.rect(14, 46, 36, 8, 2)
    ic.rect(20, 54, 4, 8)
    ic.rect(40, 54, 4, 8)
    ic.ops.append(Op("s", "path", {"d": "M16 46 C16 28 48 28 48 46 Z", "_fill": "#FFE066"}))
    ic.rect(28, 40, 8, 6, 1).fixed(EMIT[0])
    for a in (-60, -30, 0, 30, 60):
        x1, y1 = polar(32, 40, 16, a - 90)
        x2, y2 = polar(32, 40, 24, a - 90)
        ic.acc().line(x1, y1, x2, y2, 2)


@icon("lanthanide_mof", "materials", label="Luminescent Ln-MOF", tags=["lanthanide", "Eu3+", "Tb3+", "LnMOF", "luminescent MOF", "antenna effect"])
def _(ic):
    f = [(13, 25), (39, 25), (39, 51), (13, 51)]
    b = [(x + 12, y - 12) for x, y in f]
    for i in range(4):
        for p, q in ((f[i], f[(i + 1) % 4]), (b[i], b[(i + 1) % 4]), (f[i], b[i])):
            ic.sec().line(*p, *q, 2.4)
    for p in f + b:
        glow(ic, *p, 3.5, EMIT[5], rings=1)
        ic.circle(*p, 3.8).fixed(EMIT[5])


# ==========================================================================
# NANOPARTICLE MORPHOLOGIES
# ==========================================================================

@icon("gold_nanoparticle", "materials", label="Gold Nanoparticle", tags=["AuNP", "gold", "plasmonic", "LSPR", "nanoparticle"])
def _(ic):
    glow(ic, 32, 32, 16, GOLD, rings=2)
    ic.circle(32, 32, 16).fixed(GOLD)
    ic.cut().spath("M22 28 A11 11 0 0 1 30 20", 2.4)


@icon("gold_nanorod", "materials", label="Gold Nanorods", tags=["nanorod", "AuNR", "plasmonic", "anisotropic"])
def _(ic):
    for cx, cy, a, l in ((22, 24, -30, 34), (42, 40, -30, 34), (24, 46, 40, 22)):
        ic.ops.append(Op("m", "rect", {"x": cx - l / 2, "y": cy - 5.5, "width": l, "height": 11, "rx": 5.5,
                                       "transform": f"rotate({a} {cx} {cy})", "_fill": GOLD}))
        ic.ops.append(Op("x", "rect", {"x": cx - l / 2 + 4, "y": cy - 3.5, "width": l * 0.4, "height": 2, "rx": 1,
                                       "transform": f"rotate({a} {cx} {cy})"}))


@icon("nanocube", "materials", label="Nanocube", tags=["nanocube", "Ag nanocube", "shape-controlled", "facets"])
def _(ic):
    ic.poly([(12, 24), (38, 24), (38, 50), (12, 50)]).fixed(SILVER)
    ic.poly([(12, 24), (24, 12), (50, 12), (38, 24)]).fixed("#DDE3E9")
    ic.poly([(38, 24), (50, 12), (50, 38), (38, 50)]).fixed("#97A1AC")


@icon("nanostar", "materials", label="Nanostar", tags=["nanostar", "gold nanostar", "branched", "SERS", "hot spots"])
def _(ic):
    pts = []
    for k in range(18):
        r = 25 if k % 2 == 0 else 10
        pts.append(polar(32, 32, r * (1 - 0.12 * ((k // 2) % 3 == 1)), k * 20))
    ic.poly(pts).fixed(GOLD)
    ic.cut().circle(28, 28, 2.6)


@icon("nanoplate", "materials", label="Nanoplates / Prisms", tags=["nanoplate", "triangular prism", "nanodisk", "2D nanocrystal"])
def _(ic):
    ic.poly([(6, 44), (28, 8), (38, 30)]).fixed(SILVER)
    ic.poly([(6, 44), (38, 30), (38, 36), (6, 50)]).fixed("#97A1AC")
    ic.sec().ellipse(46, 46, 14, 5)
    ic.sec().rect(32, 42, 28, 4)
    ic.ops[-1].attrs["opacity"] = 0.8
    ic.sec().ellipse(46, 42, 14, 5)


@icon("nanocage", "materials", label="Nanocage / Nanoframe", tags=["nanocage", "nanoframe", "hollow", "galvanic replacement"])
def _(ic):
    f = [(12, 24), (38, 24), (38, 50), (12, 50)]
    b = [(x + 12, y - 12) for x, y in f]
    for i in range(4):
        for p, q in ((f[i], f[(i + 1) % 4]), (b[i], b[(i + 1) % 4]), (f[i], b[i])):
            ic.line(*p, *q, 4.5).fixed(GOLD)
    ic.cut().circle(25, 37, 5)


@icon("hollow_sphere", "materials", label="Hollow Sphere", tags=["hollow", "nanoshell", "capsule", "Kirkendall"])
def _(ic):
    ic.sec().circle(32, 32, 24)
    ic.cut().circle(32, 32, 17)
    ic.cut().spath("M15 28 A18 18 0 0 1 26 13", 2)


@icon("yolk_shell", "materials", label="Yolk–Shell", tags=["yolk-shell", "rattle", "nanorattle", "void", "core in shell"])
def _(ic):
    ic.sec().circle(32, 32, 25)
    ic.cut().circle(32, 32, 19)
    ic.circle(36, 37, 10)
    ic.cut().circle(33, 34, 2.4)


@icon("multi_yolk_shell", "materials", label="Multi-Yolk Shell", tags=["multiyolk", "encapsulated nanoparticles", "NP@MOF", "yolk-shell"])
def _(ic):
    ic.sec().circle(32, 32, 25)
    ic.cut().circle(32, 32, 19.5)
    for x, y, r in ((26, 26, 5), (38, 28, 4.5), (30, 39, 5), (41, 39, 3.5)):
        ic.acc().circle(x, y, r)


@icon("janus_particle", "materials", label="Janus Particle", tags=["Janus", "asymmetric", "two-faced", "anisotropic particle"])
def _(ic):
    ic.sec().path("M32 8 A24 24 0 0 0 32 56 Z")
    ic.acc().path("M32 8 A24 24 0 0 1 32 56 Z")
    ic.cut().spath("M15 26 A18 18 0 0 1 24 14", 2)


@icon("mesoporous_silica", "materials", label="Mesoporous Silica (MSN)", tags=["MSN", "MCM-41", "SBA-15", "mesoporous", "silica"])
def _(ic):
    ic.circle(32, 32, 25)
    for j in range(-3, 4):
        for i in range(-3, 4):
            x = 32 + i * 7.2 + (3.6 if j % 2 else 0)
            y = 32 + j * 6.2
            if (x - 32) ** 2 + (y - 32) ** 2 < 20 ** 2:
                ic.cut().poly([polar(x, y, 2.6, 30 + 60 * m) for m in range(6)])


@icon("magnetic_np", "materials", label="Magnetic Nanoparticle", tags=["Fe3O4", "magnetite", "superparamagnetic", "magnetic separation"])
def _(ic):
    for r in (20, 26):
        ic.acc().spath(f"M{32 - r * .3} {32 - 14} C{32 - r} {32 - 26} {32 - r} {32 + 26} {32 - r * .3} {32 + 14}", 1.4)
        ic.acc().spath(f"M{32 + r * .3} {32 - 14} C{32 + r} {32 - 26} {32 + r} {32 + 26} {32 + r * .3} {32 + 14}", 1.4)
    ic.circle(32, 32, 13).fixed("#4A4F57")
    ic.rect(29, 18, 6, 10, 1).fixed("#E03131")
    ic.rect(29, 36, 6, 10, 1).fixed("#3B5BDB")


@icon("carbon_nanotube", "materials", label="Carbon Nanotube", tags=["CNT", "nanotube", "SWCNT", "MWCNT"])
def _(ic):
    ic.ops.append(Op("s", "rect", {"x": 6, "y": 22, "width": 52, "height": 20, "rx": 10, "transform": "rotate(-25 32 32)"}))
    for k in range(7):
        x = 12 + k * 7
        ic.ops.append(Op("x", "polygon", {"points": " ".join(f"{_f(px)},{_f(py)}" for px, py in
                                                           [polar(x, 32, 3.4, 30 + 60 * m) for m in range(6)]),
                                          "transform": "rotate(-25 32 32)"}, 1.1))
    ic.ops.append(Op("m", "ellipse", {"cx": 58, "cy": 32, "rx": 4, "ry": 10, "transform": "rotate(-25 32 32)"}))


@icon("graphene_oxide", "materials", label="Graphene Oxide", tags=["GO", "rGO", "graphene oxide", "2D carbon"])
def _(ic):
    ic.path("M4 40 L20 24 H60 L44 40 Z")
    for k in range(6):
        for j in range(2):
            cx, cy = 18 + k * 6.5 - j * 6, 28.5 + j * 6
            ic.cut().poly([polar(cx, cy, 2.4, 30 + 60 * m) for m in range(6)])
    for x, y in ((14, 22), (30, 18), (46, 20), (24, 46), (40, 46)):
        ic.acc().line(x, y + (4 if y < 30 else -4), x, y, 1.4)
        ic.acc().circle(x, y, 2.4)


@icon("nanowire_array", "materials", label="Nanowire Array", tags=["nanowires", "vertical array", "ZnO nanorods", "forest"])
def _(ic):
    ic.rect(4, 50, 56, 8, 1.5)
    for k, x in enumerate(range(9, 58, 7)):
        h = 30 + (k % 3) * 6
        ic.sec().rect(x - 2.2, 50 - h, 4.4, h, 2.2)


@icon("np_at_mof", "materials", label="NP@MOF", tags=["NP@MOF", "encapsulation", "nanoparticles in MOF", "composite"])
def _(ic):
    top, bot, l, r, f = (32, 4), (32, 60), (6, 30), (58, 30), (27, 38)
    ic.sec().poly([top, l, f])
    ic.poly([top, f, r])
    ic.poly([bot, l, f])
    ic.sec().poly([bot, f, r])
    for x, y in ((24, 26), (36, 22), (40, 34), (30, 44), (20, 38)):
        ic.circle(x, y, 2.8).fixed(GOLD)


@icon("mof_shell_np", "materials", label="Core–Shell NP@MOF", tags=["core-shell", "MOF shell", "NP@MOF", "coating"])
def _(ic):
    for k in range(10):
        a = k * 36
        ic.sec().poly([polar(32, 32, 24, a - 9), polar(32, 32, 27, a), polar(32, 32, 24, a + 9)])
    ic.sec().circle(32, 32, 24)
    for k in range(18):
        ic.cut().circle(*polar(32, 32, 20.5, k * 20), 1.3)
    ic.circle(32, 32, 14).fixed(GOLD)
    ic.cut().circle(28, 28, 2.4)


@icon("polymer_coated_np", "materials", label="Polymer-Coated NP", tags=["polymer brush", "PEGylated", "stealth coating", "nanoparticle"])
def _(ic):
    for k in range(14):
        a = k * 360 / 14
        p0 = polar(32, 32, 13, a)
        p1 = polar(32, 32, 20, a + 10)
        p2 = polar(32, 32, 27, a - 6)
        ic.sec().spath(f"M{_f(p0[0])} {_f(p0[1])} Q{_f(p1[0])} {_f(p1[1])} {_f(p2[0])} {_f(p2[1])}", 1.8)
    ic.circle(32, 32, 13)
    ic.cut().circle(28.5, 28.5, 2.2)


@icon("single_atom_catalyst", "materials", label="Single-Atom Catalyst", tags=["SAC", "single atoms", "atomically dispersed", "M-N-C"])
def _(ic):
    ic.path("M4 46 L18 32 H60 L46 46 Z")
    ic.path("M4 46 H46 V52 H4 Z")
    ic.path("M46 46 L60 32 V38 L46 52 Z")
    for x, y in ((18, 40), (30, 37), (42, 40), (34, 43), (24, 44), (48, 35)):
        ic.acc().circle(x, y - 4, 2.6)
        ic.sec().line(x, y - 1.4, x, y, 1.2)


@icon("photocatalysis", "materials", label="Photocatalysis (e⁻/h⁺)", tags=["photocatalysis", "electron-hole", "semiconductor", "solar fuels"])
def _(ic):
    ic.sec().circle(32, 36, 18)
    ic.cut().rect(16, 26, 32, 3, 1.5)
    ic.cut().rect(16, 44, 32, 3, 1.5)
    ic.acc().circle(26, 27.5, 2.4)
    ic.circle(38, 45.5, 2.4)
    ic.acc().line(30, 43, 30, 31, 1.6)
    arrowhead(ic, 30, 29, -90, 3)
    ic.acc().spath("M4 6 q3 4 6 0 t6 0 t6 0", 2)
    arrowhead(ic, 23, 9, 30, 3.5)


@icon("heterojunction", "data", label="Heterojunction (Type II)", tags=["heterojunction", "band alignment", "type II", "charge separation"])
def _(ic):
    ic.rect(6, 12, 22, 7, 1)
    ic.sec().rect(6, 44, 22, 7, 1)
    ic.rect(36, 20, 22, 7, 1)
    ic.ops[-1].role = "a"
    ic.sec().rect(36, 50, 22, 7, 1)
    ic.ops[-1].role = "a"
    ic.line(26, 15, 40, 23, 1.6)
    arrowhead(ic, 42, 24, 30, 3.5, "m")
    ic.line(40, 53, 26, 48, 1.6)
    arrowhead(ic, 24, 47.5, 200, 3.5, "m")


# ==========================================================================
# MOF & COF MOTIFS
# ==========================================================================

@icon("cof_2d_sheet", "materials", label="2D COF (hcb)", tags=["COF", "2D COF", "hcb", "honeycomb", "covalent organic framework"])
def _(ic):
    r = 11
    centres = [(32, 32), (32 - 1.732 * r, 32), (32 + 1.732 * r, 32), (32 - .866 * r, 32 - 1.5 * r),
               (32 + .866 * r, 32 - 1.5 * r), (32 - .866 * r, 32 + 1.5 * r), (32 + .866 * r, 32 + 1.5 * r)]
    nodes = set()
    for cx, cy in centres:
        hp = [polar(cx, cy, r, 30 + 60 * k) for k in range(6)]
        ic.sec().spoly(hp + [hp[0]], 2.6)
        nodes.update((round(x, 1), round(y, 1)) for x, y in hp)
    for x, y in nodes:
        ic.poly([polar(x, y, 2.6, 90 + 120 * m) for m in range(3)])


@icon("cof_stacking_aa", "materials", label="COF AA Stacking", tags=["AA stacking", "eclipsed", "COF layers", "1D channels"])
def _(ic):
    for k, y in enumerate((14, 28, 42)):
        pts = [(cx + dx, y + dy) for cx, dx, dy in ((22, 0, 0),)]
        ic.path(f"M6 {y + 8} L20 {y} H58 L44 {y + 8} Z", "s" if k % 2 else "m")
        for cx in (22, 38):
            ic.cut().ellipse(cx + 3, y + 4, 4.5, 2)


@icon("cof_stacking_ab", "materials", label="COF AB Stacking", tags=["AB stacking", "staggered", "COF layers", "offset"])
def _(ic):
    for k, y in enumerate((14, 28, 42)):
        off = 6 if k % 2 else 0
        ic.path(f"M{6 + off} {y + 8} L{20 + off} {y} H{52 + off} L{38 + off} {y + 8} Z", "s" if k % 2 else "m")
        for cx in (20, 34):
            ic.cut().ellipse(cx + 3 + off, y + 4, 4, 2)


@icon("imine_linkage", "materials", label="Imine Linkage (C=N)", tags=["imine", "Schiff base", "C=N", "COF linkage"])
def _(ic):
    hexagon(ic, 13, 32, 8, 2.4)
    hexagon(ic, 51, 32, 8, 2.4)
    ic.line(21, 32, 27, 32, 2.4)
    ic.line(27, 30, 37, 30, 2.2)
    ic.line(27, 34, 37, 34, 2.2)
    ic.circle(27, 32, 3)
    ic.acc().circle(37, 32, 3.6)
    ic.line(37, 32, 43, 32, 2.4)


@icon("boronate_ester", "materials", label="Boronate Ester Linkage", tags=["boronate ester", "boroxine", "B-O", "COF-5", "COF linkage"])
def _(ic):
    hexagon(ic, 12, 32, 8, 2.4)
    ic.line(20, 32, 28, 32, 2.4)
    ic.acc().circle(29, 32, 3.6)
    ic.spoly([(29, 32), (37, 24), (45, 28), (45, 36), (37, 40), (29, 32)], 2.2)
    ic.sec().circle(37, 24, 2.8)
    ic.sec().circle(37, 40, 2.8)
    hexagon(ic, 52, 32, 7, 2.2)


@icon("triazine_node", "materials", label="Triazine Node", tags=["triazine", "CTF", "covalent triazine framework", "C3 node"])
def _(ic):
    pts = hexagon(ic, 32, 32, 9, 2.6)
    for k, p in enumerate(pts):
        if k % 2 == 0:
            ic.acc().circle(*p, 2.6)
    for a in (-90, 30, 150):
        x1, y1 = polar(32, 32, 9, a)
        x2, y2 = polar(32, 32, 20, a)
        ic.line(x1, y1, x2, y2, 2.4)
        hexagon(ic, *polar(32, 32, 26, a), 5.5, 2, "s")


@icon("cof_3d_dia", "materials", label="3D COF / dia Net", tags=["3D COF", "dia", "diamondoid", "tetrahedral node"])
def _(ic):
    nodes = [(32, 10), (14, 22), (50, 22), (32, 32), (14, 44), (50, 44), (32, 54)]
    for a, b in ((0, 1), (0, 2), (1, 3), (2, 3), (3, 4), (3, 5), (4, 6), (5, 6), (1, 4), (2, 5)):
        ic.sec().line(*nodes[a], *nodes[b], 2.4)
    for p in nodes:
        ic.circle(*p, 4)
    ic.acc().circle(32, 32, 4.4)


@icon("zif_cage", "materials", label="ZIF Sodalite Cage", tags=["ZIF-8", "sodalite", "sod", "imidazolate", "cage"])
def _(ic):
    sq = [(18, 14), (46, 14), (58, 32), (46, 50), (18, 50), (6, 32)]
    for i in range(6):
        ic.sec().line(*sq[i], *sq[(i + 1) % 6], 2.4)
    inner = [(24, 24), (40, 24), (46, 32), (40, 40), (24, 40), (18, 32)]
    for i in range(6):
        ic.sec().line(*inner[i], *inner[(i + 1) % 6], 2)
        ic.sec().line(*sq[i], *inner[i], 1.8)
    for p in sq + inner:
        ic.circle(*p, 3)
    ic.acc().circle(32, 32, 5)


@icon("mof_channel", "materials", label="1D Channel with Guests", tags=["channel", "1D pore", "MIL-53", "confinement", "guest"])
def _(ic):
    ic.rect(4, 12, 56, 10, 2)
    ic.rect(4, 42, 56, 10, 2)
    for x in range(10, 58, 8):
        ic.acc().circle(x, 17, 2.2)
        ic.acc().circle(x, 47, 2.2)
    for x, y in ((14, 30), (26, 34), (38, 29), (50, 33)):
        ic.sec().circle(x, y, 3.2)
        ic.sec().circle(x + 4, y + 1, 2)


@icon("open_metal_site", "materials", label="Open Metal Site", tags=["open metal site", "coordinatively unsaturated", "Lewis acid", "binding"])
def _(ic):
    for a in (90, 180, 270, 0):
        x2, y2 = polar(32, 36, 16, a)
        ic.sec().line(32, 36, x2, y2, 2.6)
        ic.sec().circle(x2, y2, 3)
    ic.acc().circle(32, 36, 6)
    ic.cut().line(32, 30, 32, 21, 1.4)
    ic.ops[-1].attrs["stroke-dasharray"] = "2 2"
    ic.circle(32, 12, 4.4)
    ic.acc().circle(26, 12, 2.6)
    ic.acc().circle(38, 12, 2.6)


@icon("mixed_linker", "materials", label="Mixed-Linker MOF", tags=["mixed linker", "multivariate", "MTV-MOF", "functional groups"])
def _(ic):
    nodes = [(12, 12), (32, 12), (52, 12), (12, 32), (32, 32), (52, 32), (12, 52), (32, 52), (52, 52)]
    edges = [(0, 1), (1, 2), (3, 4), (4, 5), (6, 7), (7, 8), (0, 3), (3, 6), (1, 4), (4, 7), (2, 5), (5, 8)]
    for k, (i, j) in enumerate(edges):
        ic.line(*nodes[i], *nodes[j], 3, role="s" if k % 3 else "a")
    for p in nodes:
        ic.circle(*p, 4.6)


@icon("hierarchical_pores", "materials", label="Hierarchical Porosity", tags=["hierarchical", "micro-meso", "mesopores", "macropores"])
def _(ic):
    ic.circle(32, 32, 26)
    ic.cut().circle(24, 26, 8)
    ic.cut().circle(42, 40, 7)
    rnd = random.Random(4)
    for _ in range(26):
        x, y = 10 + rnd.random() * 44, 10 + rnd.random() * 44
        if (x - 32) ** 2 + (y - 32) ** 2 < 22 ** 2 and (x - 24) ** 2 + (y - 26) ** 2 > 100 and (x - 42) ** 2 + (y - 40) ** 2 > 81:
            ic.cut().circle(x, y, 1.1)


@icon("mof_glass", "materials", label="MOF Glass / Amorphous", tags=["MOF glass", "amorphous", "melt-quenched", "disordered network"])
def _(ic):
    rnd = random.Random(9)
    pts = [(8 + rnd.random() * 48, 8 + rnd.random() * 48) for _ in range(16)]
    for i, p in enumerate(pts):
        ds = sorted(range(len(pts)), key=lambda j: math.dist(p, pts[j]))[1:4]
        for j in ds:
            if j > i:
                ic.sec().line(*p, *pts[j], 2)
    for p in pts:
        ic.circle(*p, 3.2)


# ==========================================================================
# COMPOSITES & SOFT MATTER
# ==========================================================================

@icon("fiber_composite", "materials", label="Fibre-Reinforced Composite", tags=["fibre composite", "reinforcement", "matrix", "laminate"])
def _(ic):
    ic.sec().rect(6, 14, 52, 36, 4)
    ic.ops[-1].attrs["opacity"] = 0.6
    for y in (20, 28, 36, 44):
        ic.rect(4, y - 2, 56, 4, 2)


@icon("nanocomposite", "materials", label="Nanocomposite", tags=["nanocomposite", "dispersed nanoparticles", "polymer matrix", "filler"])
def _(ic):
    ic.sec().rect(6, 10, 52, 44, 6)
    ic.ops[-1].attrs["opacity"] = 0.6
    rnd = random.Random(6)
    for _ in range(14):
        ic.circle(12 + rnd.random() * 40, 16 + rnd.random() * 32, 2 + rnd.random() * 1.6)


@icon("mixed_matrix_membrane", "materials", label="Mixed-Matrix Membrane", tags=["MMM", "MOF-polymer", "membrane", "composite"])
def _(ic):
    ic.sec().rect(4, 22, 56, 20, 3)
    ic.ops[-1].attrs["opacity"] = 0.65
    for x in (12, 26, 40, 52):
        ic.poly([(x, 26), (x + 5, 32), (x, 38), (x - 5, 32)])
    for x in (16, 34, 50):
        ic.acc().circle(x, 12, 2.4)
        ic.acc().line(x, 15, x, 19, 1.4)
    ic.acc().circle(26, 52, 2.4)


@icon("laminate", "materials", label="Layered Composite", tags=["laminate", "multilayer", "heterostructure", "stack"])
def _(ic):
    for k, (y, role) in enumerate(((12, "m"), (20, "s"), (28, "a"), (36, "s"), (44, "m"))):
        ic.path(f"M6 {y + 8} L18 {y} H58 L46 {y + 8} Z", role)


@icon("hydrogel_network", "materials", label="Hydrogel Network", tags=["hydrogel", "crosslinked", "polymer network", "swelling"])
def _(ic):
    rnd = random.Random(8)
    chains = []
    for k in range(5):
        y0 = 10 + k * 11
        pts = [(4 + i * 7, y0 + rnd.uniform(-4, 4)) for i in range(9)]
        chains.append(pts)
        ic.sec().spath(smooth_path(pts), 2)
    for k in range(4):
        i = 2 + (k * 3) % 6
        ic.acc().line(*chains[k][i], *chains[k + 1][i], 1.6)
        ic.acc().circle(*chains[k][i], 2)
    for x, y in ((20, 30), (40, 18), (34, 44)):
        ic.circle(x, y, 1.6)


@icon("core_shell_fiber", "materials", label="Core–Shell Fibre", tags=["core-shell fibre", "coaxial electrospinning", "fibre", "coating"])
def _(ic):
    ic.ops.append(Op("s", "rect", {"x": 2, "y": 22, "width": 56, "height": 20, "rx": 10, "transform": "rotate(-20 32 32)"}))
    ic.ops.append(Op("m", "rect", {"x": 2, "y": 27, "width": 50, "height": 10, "rx": 5, "transform": "rotate(-20 32 32)"}))
    ic.ops.append(Op("s", "ellipse", {"cx": 56, "cy": 32, "rx": 4, "ry": 10, "transform": "rotate(-20 32 32)"}))
    ic.ops.append(Op("m", "ellipse", {"cx": 56, "cy": 32, "rx": 2, "ry": 5, "transform": "rotate(-20 32 32)"}))
