"""
v4 content: MOF building blocks, instruments, sample holders, catalysis / proteomics schemes,
arrows & annotations, gloved hands, extra plot types, and badges that tell look-alike plots apart.
All original geometry (CC0).
"""
from __future__ import annotations

import math

from .core import Op, polar, smooth_path, lorentz_trace, _f
from .library import icon, axes, REGISTRY
from .library2 import rot_rect, gauss_trace


def hexagon(ic, cx, cy, r, w=2.8, role="m", rot=90):
    pts = [polar(cx, cy, r, rot + 60 * k) for k in range(6)]
    ic.spoly(pts + [pts[0]], w, role=role)
    return pts


def arrowhead(ic, x, y, deg, s=5, role="a"):
    p1 = polar(x, y, s, deg + 150)
    p2 = polar(x, y, s, deg - 150)
    ic.poly([(x, y), p1, p2], role=role)


def carboxylate(ic, x, y, deg, role="s"):
    """C with two O at the end of a linker arm pointing in direction deg."""
    for d in (-35, 35):
        ox, oy = polar(x, y, 5.5, deg + d)
        ic.line(x, y, ox, oy, 2.4)
        ic.circle(ox, oy, 2.6, role=role)


# ==========================================================================
# MOF BUILDING BLOCKS & MATERIALS
# ==========================================================================

@icon("paddlewheel_sbu", "materials", label="Paddlewheel SBU", tags=["Cu2", "paddlewheel", "HKUST-1", "SBU", "MOF"])
def _(ic):
    for a in (0, 90, 180, 270):
        x1, y1 = polar(32, 32, 7, a)
        x2, y2 = polar(32, 32, 22, a)
        ic.sec().line(x1, y1, x2, y2, 3)
        ic.sec().circle(x2, y2, 3.4)
    for (x, y) in ((25, 25), (39, 25), (25, 39), (39, 39)):
        ic.circle(x, y, 3)
    ic.acc().circle(29, 32, 5)
    ic.acc().circle(35, 32, 5)
    ic.line(29, 32, 35, 32, 2.4)


@icon("zn4o_sbu", "materials", label="Zn4O SBU", tags=["Zn4O", "MOF-5", "tetrahedral", "SBU"])
def _(ic):
    for a in range(0, 360, 60):
        x1, y1 = polar(32, 32, 10, a)
        x2, y2 = polar(32, 32, 24, a)
        ic.sec().line(x1, y1, x2, y2, 3)
        ic.sec().circle(x2, y2, 3.2)
    tet = [(32, 21), (22.5, 38), (41.5, 38), (32, 34)]
    for i in range(4):
        for j in range(i + 1, 4):
            ic.line(*tet[i], *tet[j], 1.6)
    for p in tet[:3]:
        ic.acc().circle(*p, 4.6)
    ic.acc().circle(*tet[3], 3.6)
    ic.circle(32, 31, 2.6)


@icon("rod_sbu", "materials", label="Rod / Chain SBU", tags=["rod", "chain", "Al-MOF", "MIL-53", "infinite SBU"])
def _(ic):
    for k in range(4):
        cx = 11 + k * 14
        ic.sec().line(cx, 32, cx, 10, 2.6)
        ic.sec().line(cx, 32, cx, 54, 2.6)
        ic.sec().circle(cx, 9, 2.8)
        ic.sec().circle(cx, 55, 2.8)
    for k in range(4):
        cx = 11 + k * 14
        ic.poly([(cx, 22), (cx + 9, 32), (cx, 42), (cx - 9, 32)], role="m")
        ic.cut().line(cx, 22, cx, 42, 1.2)
        ic.acc().circle(cx, 32, 2.6)


@icon("linker_btc", "materials", label="BTC Linker", tags=["BTC", "trimesate", "tritopic", "linker", "HKUST-1"])
def _(ic):
    hexagon(ic, 32, 32, 9)
    ic.scircle(32, 32, 4.2, 2)
    for a in (-90, 30, 150):
        x1, y1 = polar(32, 32, 9, a)
        x2, y2 = polar(32, 32, 18, a)
        ic.line(x1, y1, x2, y2, 2.8)
        carboxylate(ic, x2, y2, a)


@icon("linker_fumarate", "materials", label="Fumarate Linker", tags=["fumarate", "fumaric acid", "MOF-801", "linker"])
def _(ic):
    pts = [(16, 38), (26, 28), (38, 36), (48, 26)]
    ic.spoly(pts, 3)
    ic.spoly([(26.5, 32.5), (36.5, 40)], 2)
    carboxylate(ic, 16, 38, 200)
    carboxylate(ic, 48, 26, 20)


@icon("linker_nh2_bdc", "materials", label="NH2-BDC Linker", tags=["amino", "NH2-BDC", "UiO-66-NH2", "linker"])
def _(ic):
    hexagon(ic, 32, 32, 10)
    ic.scircle(32, 32, 4.8, 2)
    for a in (-90, 90):
        x1, y1 = polar(32, 32, 10, a)
        x2, y2 = polar(32, 32, 18, a)
        ic.line(x1, y1, x2, y2, 2.8)
        carboxylate(ic, x2, y2, a)
    x1, y1 = polar(32, 32, 10, -30)
    x2, y2 = polar(32, 32, 18, -30)
    ic.line(x1, y1, x2, y2, 2.6)
    ic.acc().circle(x2, y2, 4)


@icon("linker_bipyridine", "materials", label="Bipyridine Linker", tags=["bipyridine", "bpy", "pillar", "N-donor", "linker"])
def _(ic):
    for cx in (20, 44):
        pts = hexagon(ic, cx, 32, 9, 2.6, rot=0)
    ic.line(29, 32, 35, 32, 2.6)
    ic.acc().circle(11, 32, 3.2)
    ic.acc().circle(53, 32, 3.2)


@icon("net_pcu", "materials", label="pcu Net", tags=["topology", "pcu", "primitive cubic", "net"])
def _(ic):
    xs = (12, 32, 52)
    for x in xs:
        ic.sec().line(x, 10, x, 54, 2.4)
        ic.sec().line(10, x, 54, x, 2.4)
    for x in xs:
        for y in xs:
            ic.circle(x, y, 4.2)


@icon("net_fcu", "materials", label="fcu Net", tags=["topology", "fcu", "12-connected", "UiO-66", "MOF-801"])
def _(ic):
    outer = [polar(32, 32, 24, 15 + 30 * k) for k in range(12)]
    for p in outer:
        ic.sec().line(32, 32, *p, 2)
    for i in range(12):
        ic.sec().line(*outer[i], *outer[(i + 1) % 12], 1.6)
    for p in outer:
        ic.circle(*p, 3)
    ic.acc().circle(32, 32, 5)


@icon("crystal_octahedron", "materials", label="Octahedral Crystal", tags=["morphology", "octahedral", "crystal", "UiO-66"])
def _(ic):
    top, bot, l, r, f, b = (32, 4), (32, 60), (8, 30), (56, 30), (27, 38), (39, 24)
    ic.sec().poly([top, l, f])
    ic.poly([top, f, r])
    ic.poly([bot, l, f])
    ic.sec().poly([bot, f, r])


@icon("crystal_cube", "materials", label="Cubic Crystal", tags=["morphology", "cube", "crystal", "ZIF-8"])
def _(ic):
    ic.poly([(12, 22), (36, 22), (36, 56), (12, 56)])
    ic.sec().poly([(12, 22), (28, 8), (52, 8), (36, 22)])
    ic.sec().poly([(36, 22), (52, 8), (52, 42), (36, 56)])
    ic.cut().spoly([(12, 22), (36, 22), (52, 8)], 1.2)
    ic.cut().line(36, 22, 36, 56, 1.2)


@icon("crystal_rods", "materials", label="Rod Crystals", tags=["morphology", "needles", "hexagonal prism", "rods"])
def _(ic):
    for cx, cy, ang, ln in ((22, 30, -60, 40), (38, 34, -100, 44), (46, 22, -30, 30)):
        rot_rect(ic, cx, cy, ln, 9, 1.5, ang)
        rot_rect(ic, cx, cy, ln, 3, 0.5, ang, "s")


@icon("thin_film", "materials", label="Thin Film", tags=["film", "coating", "substrate", "layer"])
def _(ic):
    ic.path("M4 44 L18 34 H60 L46 44 Z")
    ic.rect(4, 44, 42, 8, 1)
    ic.path("M46 44 L60 34 V42 L46 52 Z")
    ic.sec().path("M8 41 L19 33 H55 L44 41 Z")
    for x in (18, 26, 34, 42):
        ic.acc().poly([(x, 36), (x + 3, 32), (x + 6, 36)])


@icon("membrane", "materials", label="Membrane", tags=["membrane", "separation", "permeation", "MMM"])
def _(ic):
    ic.rect(28, 4, 8, 56, 2)
    for y in range(9, 58, 8):
        ic.cut().rect(30, y, 4, 4, 1.5)
    for y in (14, 30, 46):
        ic.sec().circle(14, y, 3.6)
        ic.acc().line(19, y, 25, y, 2.2)
    ic.sec().circle(50, 22, 3.6)
    ic.acc().line(39, 22, 45, 22, 2.2)


@icon("sc_co2_dryer", "synthesis", label="Supercritical CO2 Dryer", tags=["scCO2", "supercritical", "activation", "critical point dryer"])
def _(ic):
    ic.rect(6, 20, 46, 28, 10)
    ic.rect(52, 26, 8, 16, 2)
    ic.cut().circle(28, 34, 8)
    ic.sec().circle(28, 34, 5)
    ic.acc().circle(44, 12, 6)
    ic.cut().line(44, 12, 47, 9, 1.4)
    ic.line(44, 18, 44, 20, 2.4)
    ic.rect(12, 48, 6, 8, 1.5)
    ic.rect(40, 48, 6, 8, 1.5)


@icon("solvent_exchange", "synthesis", label="Solvent Exchange", tags=["solvent exchange", "washing", "activation", "soaking"])
def _(ic):
    ic.srect(20, 18, 24, 38, 3, 3, glass=True)
    ic.sec().rect(23.5, 30, 17, 22.5, 2)
    ic.acc().spath("M14 26 C6 18 14 6 26 8", 2.4)
    arrowhead(ic, 26, 8, -10, 4.5)
    ic.acc().spath("M50 48 C58 56 52 62 42 60", 2.4)
    arrowhead(ic, 42, 60, 170, 4.5)
    for x, y in ((28, 47), (34, 44), (31, 50)):
        ic.circle(x, y, 2)


@icon("post_synthetic_modification", "synthesis", label="Post-Synthetic Modification", tags=["PSM", "functionalisation", "grafting"])
def _(ic):
    ic.circle(14, 32, 7)
    ic.sec().line(21, 32, 26, 32, 2.6)
    ic.acc().line(30, 32, 40, 32, 2.6)
    arrowhead(ic, 42, 32, 0, 5)
    ic.circle(52, 26, 6)
    ic.sec().line(52, 32, 52, 40, 2.6)
    ic.acc().circle(52, 44, 4.5)


# ==========================================================================
# INSTRUMENTS (generic, brand-free)
# ==========================================================================

@icon("tga_instrument", "instruments", label="TGA Instrument", tags=["TGA", "thermogravimetric analyser", "furnace", "microbalance"])
def _(ic):
    ic.rect(6, 40, 52, 18, 4)
    ic.cut().rect(12, 46, 16, 7, 1.5)
    ic.rect(38, 8, 14, 32, 4)
    for y in (14, 20, 26, 32):
        ic.acc().line(40, y, 50, y, 1.6)
    ic.line(22, 22, 45, 22, 2)
    ic.line(22, 22, 22, 40, 2)
    ic.sec().path("M41 33 H49 L48 37 H42 Z")


@icon("dsc_instrument", "instruments", label="DSC Instrument", tags=["DSC", "calorimeter", "pans", "heat flow"])
def _(ic):
    ic.rect(6, 26, 52, 32, 4)
    ic.ellipse(32, 26, 22, 7)
    ic.cut().ellipse(32, 26, 18, 5)
    ic.sec().circle(24, 26, 3.2)
    ic.circle(40, 26, 3.2)
    ic.cut().rect(12, 38, 18, 10, 1.5)
    ic.cut().circle(46, 43, 4)
    ic.path("M12 20 L18 6 H46 L52 20 Z")


@icon("sorption_analyzer", "instruments", label="Gas Sorption Analyser", tags=["BET", "N2 sorption", "porosimeter", "3Flex", "dewar"])
def _(ic):
    ic.rect(4, 4, 56, 26, 3)
    ic.cut().rect(40, 9, 15, 9, 1.5)
    for x in (14, 24, 34):
        ic.spath(f"M{x - 2} 30 V44 A2 2 0 0 0 {x + 2} 44 V30", 2, glass=True)
        ic.sec().circle(x, 43.5, 1.6)
    ic.rect(8, 40, 32, 20, 5)
    ic.cut().rect(12, 44, 24, 3, 1.5)
    ic.acc().circle(48, 24, 2.4)


@icon("icp_oes", "instruments", label="ICP-OES", tags=["ICP", "OES", "elemental analysis", "plasma", "metal content"])
def _(ic):
    ic.rect(4, 14, 56, 42, 4)
    ic.cut().rect(9, 20, 22, 28, 2)
    ic.acc().path("M20 22 C26 28 27 34 25 38 C23.5 41 16.5 41 15 38 C13 34 14 28 20 22 Z")
    ic.cut().path("M20 30 C22.5 33 22.6 35.5 21.6 37 C20.8 38 19.2 38 18.4 37 C17.4 35.5 17.5 33 20 30 Z")
    ic.rect(17, 40, 6, 7, 1)
    ic.cut().rect(37, 20, 18, 10, 1.5)
    for x in (39, 44, 49, 54):
        ic.sec().line(x, 36, x, 46 - (x % 3) * 3, 2)


@icon("xrf_spectrometer", "instruments", label="XRF Spectrometer", tags=["XRF", "X-ray fluorescence", "elemental"])
def _(ic):
    ic.rect(6, 28, 52, 30, 4)
    ic.cut().rect(12, 36, 18, 10, 1.5)
    ic.rect(22, 18, 20, 10, 2)
    ic.sec().rect(26, 12, 12, 7, 2)
    ic.acc().spath("M8 6 l4 4 l-3 3 l4 4 l-3 3 l4 3", 2.2)
    ic.acc().spath("M50 6 l-3 4 l3 3 l-3 4 l3 3", 2)
    ic.cut().circle(46, 41, 4)


@icon("dls_zeta", "instruments", label="DLS / Zeta Analyser", tags=["DLS", "zeta potential", "particle size", "light scattering"])
def _(ic):
    ic.rect(4, 18, 56, 38, 5)
    ic.cut().rect(30, 10, 12, 20, 2)
    ic.srect(32.5, 4, 7, 22, 1.5, 2, glass=True)
    ic.sec().rect(34, 14, 4, 11, 0.8)
    ic.cut().rect(10, 26, 16, 10, 1.5)
    ic.acc().line(10, 46, 30, 46, 2.4)
    ic.cut().circle(48, 44, 5)


@icon("raman_microscope", "instruments", label="Raman Microscope", tags=["Raman", "confocal", "laser", "microscope"])
def _(ic):
    ic.rect(8, 54, 44, 6, 3)
    ic.rect(40, 14, 10, 40, 2)
    ic.rect(18, 8, 32, 10, 2.5)
    ic.rect(22, 18, 8, 14, 1.5)
    ic.rect(14, 38, 30, 4, 1.5)
    ic.sec().rect(19, 35.5, 14, 2.5, 0.8)
    ic.acc().line(26, 32, 26, 35, 2.4)
    ic.acc().line(50, 13, 60, 4, 2.2)
    ic.acc().circle(60, 4, 2.2)


@icon("ftir_atr", "instruments", label="FTIR-ATR", tags=["FTIR", "ATR", "diamond", "infrared"])
def _(ic):
    ic.rect(4, 30, 56, 28, 4)
    ic.cut().rect(10, 38, 16, 10, 1.5)
    ic.acc().poly([(34, 30), (44, 30), (41, 25), (37, 25)])
    ic.rect(36, 6, 6, 13, 1.5)
    ic.rect(32, 4, 14, 4, 1.5)
    ic.rect(34, 18, 10, 5, 1.5)
    ic.sec().rect(35, 23.5, 8, 1.6, 0.6)
    ic.line(50, 6, 50, 30, 3)
    ic.line(42, 6, 50, 6, 3)


@icon("integrating_sphere", "instruments", label="Integrating Sphere", tags=["DRS", "diffuse reflectance", "UV-Vis", "band gap"])
def _(ic):
    ic.circle(34, 32, 22)
    ic.cut().circle(34, 32, 16)
    ic.cut().rect(8, 28, 8, 8)
    ic.acc().line(2, 32, 30, 32, 2.6)
    ic.rect(30, 52, 8, 8, 1.5)
    ic.sec().rect(30.5, 50.5, 7, 2.5, 0.6)
    ic.rect(52, 26, 8, 12, 1.5)
    for a in (40, 70, 110, 150):
        x, y = polar(34, 32, 13, -a + 200)
        ic.acc().line(34, 49, x, y, 1.2)


@icon("fluorimeter", "instruments", label="Fluorimeter", tags=["fluorescence", "PL", "spectrofluorometer", "emission"])
def _(ic):
    ic.rect(4, 4, 56, 56, 4)
    ic.cut().rect(10, 10, 44, 44, 3)
    ic.srect(26, 26, 12, 12, 1.2, 2, glass=True)
    ic.sec().rect(28, 29, 8, 7, 0.6)
    ic.acc().line(12, 32, 25, 32, 2.6)
    ic.sec().line(32, 39, 32, 52, 2.6)
    arrowhead(ic, 32, 53, 90, 4, "s")


@icon("gas_chromatograph", "instruments", label="Gas Chromatograph", tags=["GC", "chromatography", "column oven", "GC-MS"])
def _(ic):
    ic.rect(6, 14, 52, 44, 4)
    ic.cut().rect(11, 20, 30, 32, 2)
    ic.scircle(26, 36, 10, 2.2, role="s")
    ic.scircle(26, 36, 6.5, 2.2, role="s")
    ic.rect(18, 6, 6, 8, 1.5)
    ic.acc().line(21, 2, 21, 6, 1.8)
    ic.cut().rect(45, 20, 9, 7, 1.2)
    ic.cut().circle(49.5, 38, 3)


@icon("orbitrap", "instruments", label="Orbitrap Analyser", tags=["Orbitrap", "high resolution", "mass spectrometer", "proteomics"])
def _(ic):
    ic.path("M6 32 C14 18 50 18 58 32 C50 46 14 46 6 32 Z")
    ic.cut().path("M11 32 C18 22 46 22 53 32 C46 42 18 42 11 32 Z")
    ic.path("M14 32 C20 29 44 29 50 32 C44 35 20 35 14 32 Z")
    ic.ops.append(Op("a", "ellipse", {"cx": 32, "cy": 32, "rx": 17, "ry": 6.5}, 1.6))
    ic.acc().circle(49, 32, 2.2)


@icon("potentiostat", "instruments", label="Potentiostat", tags=["electrochemistry", "CV", "potentiostat", "3-electrode"])
def _(ic):
    ic.rect(4, 8, 30, 22, 3)
    ic.cut().rect(8, 12, 14, 8, 1.2)
    ic.acc().spath("M9 17 q2 -4 4 0 t4 0", 1.4)
    for x, role in ((26, "a"), (30, "s")):
        ic.circle(x, 24, 1.8, role=role)
    ic.spath("M36 42 V56 Q36 60 40 60 H56 Q60 60 60 56 V42", 2.6, glass=True)
    ic.sec().rect(38.5, 48, 19, 9.5, 1.5)
    for x, role in ((42, "a"), (48, "s"), (54, "m")):
        ic.rect(x - 1.2, 36, 2.4, 18, 1, role)
    ic.acc().spath("M26 26 C26 34 42 30 42 36", 1.6)
    ic.sec().spath("M30 26 C30 32 48 30 48 36", 1.6)
    ic.spath("M22 30 C22 36 54 30 54 36", 1.6)


@icon("photoreactor", "instruments", label="Photoreactor", tags=["photocatalysis", "LED", "irradiation", "photochemistry"])
def _(ic):
    ic.rect(6, 10, 52, 48, 6)
    ic.cut().circle(32, 34, 17)
    ic.srect(27, 18, 10, 30, 2, 2.4, glass=True)
    ic.sec().rect(29.5, 30, 5, 15.5, 1)
    for a in range(0, 360, 45):
        x1, y1 = polar(32, 34, 13, a)
        x2, y2 = polar(32, 34, 16, a)
        ic.acc().line(x1, y1, x2, y2, 1.8)


@icon("flow_reactor", "instruments", label="Flow Reactor", tags=["continuous flow", "microreactor", "coil", "flow chemistry"])
def _(ic):
    ic.rect(18, 14, 28, 36, 4)
    ic.cut().rect(22, 18, 20, 28, 2)
    for k in range(5):
        ic.sec().ellipse(32, 22 + k * 5.2, 8, 2.2)
    ic.acc().line(2, 22, 18, 22, 2.4)
    ic.acc().line(2, 34, 18, 22, 2.4)
    ic.sec().line(46, 44, 62, 44, 2.4)
    arrowhead(ic, 62, 44, 0, 4, "s")


@icon("pressure_reactor", "instruments", label="High-Pressure Reactor", tags=["Parr", "autoclave", "hydrogenation", "pressure vessel"])
def _(ic):
    ic.rect(14, 24, 36, 34, 5)
    ic.rect(10, 16, 44, 9, 2.5)
    for x in (16, 26, 38, 48):
        ic.cut().circle(x, 20.5, 1.8)
    ic.rect(30, 8, 4, 8)
    ic.acc().circle(44, 8, 6)
    ic.cut().line(44, 8, 47, 5, 1.4)
    ic.line(34, 10, 38, 10, 2)
    ic.cut().rect(22, 32, 20, 18, 2)
    ic.sec().rect(25, 39, 14, 9, 1.5)


@icon("uv_lamp", "instruments", label="UV Lamp", tags=["UV", "254 nm", "365 nm", "TLC lamp", "irradiation"])
def _(ic):
    ic.rect(6, 10, 52, 14, 4)
    ic.acc().rect(10, 20, 44, 5, 2)
    for x in range(12, 56, 8):
        ic.acc().spath(f"M{x} 30 q2 4 0 8 t0 8 t0 8", 1.8)
    ic.rect(28, 4, 8, 6, 1.5)


# ==========================================================================
# SAMPLE HOLDERS & PREP
# ==========================================================================

@icon("pxrd_holder", "glassware", label="PXRD Sample Holder", tags=["flat plate", "sample holder", "PXRD", "well"])
def _(ic):
    ic.rect(6, 20, 52, 24, 3)
    ic.cut().circle(32, 32, 9)
    ic.sec().circle(32, 32, 7.5)
    ic.cut().line(10, 24, 18, 24, 1.4)


@icon("zero_background_plate", "glassware", label="Zero-Background Plate", tags=["Si", "zero background", "low background", "PXRD"])
def _(ic):
    ic.circle(32, 32, 24)
    ic.cut().circle(32, 32, 20)
    ic.circle(32, 32, 18)
    ic.sec().ellipse(32, 32, 8, 6)
    ic.cut().spath("M20 22 A14 14 0 0 1 28 17", 2)


@icon("capillary", "glassware", label="Capillary", tags=["capillary", "transmission", "Debye-Scherrer", "synchrotron"])
def _(ic):
    rot_rect(ic, 32, 32, 6, 56, 3, 45)
    ic.ops[-1].glass = True
    ic.ops[-1].stroke = 2
    ic.ops.append(Op("s", "rect", {"x": 29.8, "y": 32, "width": 4.4, "height": 22, "rx": 2,
                                   "transform": "rotate(45 32 32)"}))
    ic.ops.append(Op("m", "rect", {"x": 28.5, "y": 2, "width": 7, "height": 7, "rx": 1.5,
                                   "transform": "rotate(45 32 32)"}))


@icon("nmr_tube", "glassware", label="NMR Tube", tags=["NMR tube", "5 mm", "solution NMR"])
def _(ic):
    ic.sec().rect(28, 2, 8, 7, 2)
    ic.spath("M29.5 9 V56 A2.5 2.5 0 0 0 34.5 56 V9", 2.4, glass=True)
    ic.sec().rect(30.8, 42, 2.4, 15, 1.2)


@icon("ssnmr_rotor", "glassware", label="MAS Rotor", tags=["solid-state NMR", "MAS", "rotor", "ZrO2 rotor"])
def _(ic):
    ic.rect(24, 14, 16, 44, 3)
    ic.cut().rect(27, 20, 10, 32, 2)
    ic.sec().rect(28.5, 26, 7, 24, 1)
    ic.acc().rect(22, 6, 20, 9, 2)
    for x in (26, 30, 34, 38):
        ic.cut().line(x, 7.5, x - 2, 13.5, 1.4)
    ic.acc().spath("M10 34 A22 8 0 0 0 54 34", 1.8)


@icon("tem_grid", "glassware", label="TEM Grid", tags=["TEM grid", "copper grid", "holey carbon"])
def _(ic):
    ic.circle(32, 32, 26)
    ic.cut().circle(32, 32, 21)
    ic.circle(32, 32, 20)
    for k in range(-3, 4):
        for j in range(-3, 4):
            x, y = 32 + k * 5.6, 32 + j * 5.6
            if (x - 32) ** 2 + (y - 32) ** 2 < 17.5 ** 2:
                ic.cut().rect(x - 2, y - 2, 4, 4, 0.6)
    ic.rect(54, 28, 8, 8, 1.5)


@icon("sem_stub", "glassware", label="SEM Stub", tags=["SEM stub", "pin stub", "carbon tape", "sample mount"])
def _(ic):
    ic.ellipse(32, 24, 24, 8)
    ic.rect(8, 24, 48, 6)
    ic.ellipse(32, 30, 24, 8)
    ic.rect(29, 34, 6, 26, 1)
    ic.ops.append(Op("m", "ellipse", {"cx": 32, "cy": 24, "rx": 24, "ry": 8}))
    ic.ops.append(Op("s", "ellipse", {"cx": 32, "cy": 24, "rx": 15, "ry": 5}))
    for x, y in ((28, 23), (34, 25), (36, 22)):
        ic.acc().circle(x, y, 1.6)


@icon("sputter_coater", "instruments", label="Sputter Coater", tags=["sputter", "gold coating", "SEM prep", "plasma"])
def _(ic):
    ic.rect(6, 44, 52, 14, 3)
    ic.cut().rect(12, 48, 14, 6, 1.2)
    ic.spath("M14 44 V20 A18 12 0 0 1 50 20 V44", 2.6, glass=True)
    ic.rect(22, 12, 20, 4, 1.5)
    ic.acc().path("M24 18 C28 26 36 26 40 18 Z")
    ic.rect(24, 38, 16, 4, 1.5)
    for x in (27, 32, 37):
        ic.acc().circle(x, 31, 1.3)


# ==========================================================================
# CATALYSIS / PROTEOMICS SCHEMES
# ==========================================================================

@icon("dipeptide_hydrolysis", "proteomics", label="Dipeptide Hydrolysis", tags=["GlyGly", "dipeptide", "hydrolysis", "peptide bond"])
def _(ic):
    ic.circle(8, 22, 5)
    ic.sec().circle(20, 22, 5)
    ic.line(13, 22, 15, 22, 2.6)
    ic.acc().line(28, 22, 38, 22, 2.6)
    arrowhead(ic, 40, 22, 0, 5)
    ic.circle(48, 16, 5)
    ic.sec().circle(56, 28, 5)
    ic.sec().circle(32, 40, 4)
    ic.acc().circle(27, 37, 2.4)
    ic.acc().circle(37, 37, 2.4)


@icon("substrate_to_product", "proteomics", label="Substrate → Product", tags=["conversion", "reaction", "catalysis", "S to P"])
def _(ic):
    ic.sec().circle(12, 32, 8)
    ic.acc().line(23, 32, 40, 32, 3)
    arrowhead(ic, 43, 32, 0, 6)
    ic.circle(32, 22, 4.5)
    hexagon(ic, 53, 32, 8, 3, "s")


@icon("enzyme_vs_nanozyme", "proteomics", label="Enzyme vs Nanozyme", tags=["nanozyme", "artificial enzyme", "mimic", "comparison"])
def _(ic):
    ic.sec().path(smooth_path([(6, 34), (8, 22), (18, 16), (26, 22), (28, 34), (20, 44), (10, 44), (6, 34)]))
    ic.cut().spath("M14 26 q4 4 8 0", 1.6)
    ic.acc().line(31, 28, 37, 28, 2.2)
    ic.acc().line(31, 36, 37, 36, 2.2)
    top, bot, l, r, f = (50, 16), (50, 50), (38, 32), (62, 32), (48, 35)
    ic.poly([top, l, f])
    ic.sec().poly([top, f, r])
    ic.poly([bot, l, f])
    ic.sec().poly([bot, f, r])


@icon("protein_globule", "proteomics", label="Globular Protein", tags=["BSA", "myoglobin", "lysozyme", "protein", "substrate"])
def _(ic):
    ic.sec().path(smooth_path([(12, 34), (14, 18), (28, 10), (44, 14), (54, 28), (50, 46), (34, 54), (18, 50), (12, 34)]))
    ic.cut().spath(smooth_path([(20, 30), (26, 22), (32, 30), (38, 22), (44, 30)]), 2.2)
    ic.cut().spath(smooth_path([(22, 40), (30, 44), (38, 40)]), 2.2)


@icon("heavy_water", "proteomics", label="D2O", tags=["D2O", "deuterium", "pD", "NMR solvent", "heavy water"])
def _(ic):
    ic.line(32, 30, 16, 44, 3)
    ic.line(32, 30, 48, 44, 3)
    ic.circle(32, 28, 12)
    ic.acc().circle(16, 46, 8)
    ic.acc().circle(48, 46, 8)
    ic.cut().circle(16, 46, 2.6)
    ic.cut().circle(48, 46, 2.6)


@icon("selective_cleavage", "proteomics", label="Selective Cleavage", tags=["selectivity", "specific site", "protease", "cleavage"])
def _(ic):
    xs = [6, 15, 24, 40, 49, 58]
    for i in range(len(xs) - 1):
        if i != 2:
            ic.sec().line(xs[i], 40, xs[i + 1], 40, 2.4)
    for i, x in enumerate(xs):
        ic.circle(x, 40, 3.8) if i not in (2, 3) else ic.acc().circle(x, 40, 3.8)
    ic.acc().path("M32 8 C38 8 41 12 41 17 C41 24 32 32 32 32 C32 32 23 24 23 17 C23 12 26 8 32 8 Z")
    ic.cut().circle(32, 17, 3.4)


@icon("metal_leaching", "materials", label="Metal Leaching", tags=["leaching", "stability", "ICP", "hot filtration"])
def _(ic):
    ic.circle(26, 36, 16)
    ic.cut().circle(21, 31, 3)
    for (x, y, a) in ((48, 16, -40), (54, 34, 0), (46, 52, 40)):
        ic.acc().circle(x, y, 3.2)
        sx, sy = polar(26, 36, 18.5, a)
        ex, ey = polar(x, y, 5, a + 180)
        ic.acc().line(sx, sy, ex, ey, 1.6)


@icon("recycling_runs", "materials", label="Recyclability", tags=["reuse", "recycling", "cycles", "stability"])
def _(ic):
    for a0 in (0, 120, 240):
        p1 = polar(32, 32, 22, a0 + 10)
        p2 = polar(32, 32, 22, a0 + 95)
        ic.acc().spath(f"M{_f(p1[0])} {_f(p1[1])} A22 22 0 0 1 {_f(p2[0])} {_f(p2[1])}", 3)
        tip = polar(32, 32, 22, a0 + 108)
        ic.acc().poly([tip, polar(32, 32, 17.5, a0 + 96), polar(32, 32, 26.5, a0 + 96)])
    ic.circle(32, 32, 10)
    ic.cut().circle(28.5, 28.5, 2.4)


# ==========================================================================
# ARROWS & ANNOTATIONS (for schemes)
# ==========================================================================

@icon("arrow_forward", "annotations", label="Reaction Arrow", tags=["arrow", "reaction", "yields", "forward"])
def _(ic):
    ic.line(6, 32, 50, 32, 3.6)
    arrowhead(ic, 58, 32, 0, 9, "m")


@icon("arrow_equilibrium", "annotations", label="Equilibrium Arrows", tags=["equilibrium", "harpoons", "reversible reaction"])
def _(ic):
    ic.line(8, 27, 52, 27, 3)
    ic.poly([(56, 27), (46, 20), (46, 27)])
    ic.line(12, 37, 56, 37, 3)
    ic.poly([(8, 37), (18, 44), (18, 37)])


@icon("arrow_double", "annotations", label="Double Arrow", tags=["reversible", "both directions", "exchange"])
def _(ic):
    ic.line(14, 32, 50, 32, 3.6)
    arrowhead(ic, 58, 32, 0, 9, "m")
    arrowhead(ic, 6, 32, 180, 9, "m")


@icon("arrow_curved", "annotations", label="Curved Arrow", tags=["curved arrow", "mechanism", "process step"])
def _(ic):
    ic.spath("M10 46 C14 16 42 12 52 30", 5)
    arrowhead(ic, 54, 36, 70, 8, "m")


@icon("arrow_branch", "annotations", label="Branching Arrow", tags=["split", "branch", "two products", "workflow"])
def _(ic):
    ic.line(4, 32, 24, 32, 3.4)
    ic.spath("M24 32 C32 32 34 16 48 16", 3.4)
    ic.spath("M24 32 C32 32 34 48 48 48", 3.4)
    arrowhead(ic, 56, 16, 0, 8, "m")
    arrowhead(ic, 56, 48, 0, 8, "m")


@icon("heat_symbol", "annotations", label="Heat (Δ)", tags=["delta", "heat", "heating", "thermal", "Δ"])
def _(ic):
    ic.acc().spoly([(32, 10), (54, 52), (10, 52), (32, 10)], 4.5)


@icon("light_symbol", "annotations", label="Light (hν)", tags=["hv", "photon", "light", "irradiation"])
def _(ic):
    ic.acc().spath("M4 40 q4 -10 8 0 t8 0 t8 0 t8 0 t8 0", 3)
    ic.acc().line(44, 40, 50, 40, 3)
    arrowhead(ic, 58, 40, 0, 8)
    ic.acc().circle(14, 16, 6)
    for a in range(0, 360, 45):
        x1, y1 = polar(14, 16, 8.5, a)
        x2, y2 = polar(14, 16, 11.5, a)
        ic.acc().line(x1, y1, x2, y2, 1.8)


@icon("ultrasound_symbol", "annotations", label="Ultrasound", tags=["ultrasound", "sonication", ")))", "US"])
def _(ic):
    for r in (8, 16, 24):
        p1 = polar(12, 32, r, -45)
        p2 = polar(12, 32, r, 45)
        ic.acc().spath(f"M{_f(p1[0])} {_f(p1[1])} A{r} {r} 0 0 1 {_f(p2[0])} {_f(p2[1])}", 3.4)


@icon("pressure_symbol", "annotations", label="Pressure", tags=["pressure", "bar", "gauge", "vacuum"])
def _(ic):
    ic.circle(32, 34, 22)
    ic.cut().circle(32, 34, 17)
    for k in range(7):
        x1, y1 = polar(32, 34, 13, 150 + k * 40)
        x2, y2 = polar(32, 34, 15.5, 150 + k * 40)
        ic.line(x1, y1, x2, y2, 1.6)
    ic.acc().line(32, 34, 42, 24, 2.6)
    ic.circle(32, 34, 2.4)
    ic.rect(28, 6, 8, 6, 1.5)


@icon("time_symbol", "annotations", label="Time", tags=["time", "clock", "duration", "hours"])
def _(ic):
    ic.circle(32, 32, 24)
    ic.cut().circle(32, 32, 19.5)
    for k in range(12):
        x1, y1 = polar(32, 32, 16, k * 30)
        x2, y2 = polar(32, 32, 18, k * 30)
        ic.line(x1, y1, x2, y2, 1.4)
    ic.acc().line(32, 32, 32, 19, 2.6)
    ic.acc().line(32, 32, 42, 36, 2.6)


@icon("bracket", "annotations", label="Bracket", tags=["bracket", "group", "brace"])
def _(ic):
    ic.spath("M22 6 H14 V58 H22", 3.4)
    ic.spath("M42 6 H50 V58 H42", 3.4)


@icon("zoom_inset", "annotations", label="Zoom Inset", tags=["zoom", "magnify", "inset", "detail", "close-up"])
def _(ic):
    ic.scircle(16, 46, 7, 2.6)
    ic.scircle(42, 22, 17, 3)
    ic.line(21, 41, 30, 34, 1.6)
    ic.line(22.5, 47, 47, 38, 1.6)
    ic.acc().circle(42, 22, 5)


@icon("callout", "annotations", label="Callout", tags=["label", "note", "speech bubble", "annotation"])
def _(ic):
    ic.path("M6 10 H58 V40 H30 L18 54 L20 40 H6 Z")
    ic.cut().rect(12, 17, 40, 3.5, 1.75)
    ic.cut().rect(12, 25, 30, 3.5, 1.75)


@icon("scale_bar", "annotations", label="Scale Bar", tags=["scale bar", "length", "micrometre", "SEM"])
def _(ic):
    ic.rect(8, 34, 48, 7, 1)
    ic.line(8, 26, 8, 46, 2)
    ic.line(56, 26, 56, 46, 2)
    ic.cut().rect(30, 34, 2, 7)


# ==========================================================================
# GLOVED HANDS
# ==========================================================================

def _fist(ic, x, y):
    ic.sec().rect(x, y, 22, 20, 8)
    for k in range(3):
        ic.cut().line(x + 4, y + 6 + k * 5, x + 14, y + 6 + k * 5, 1.4)
    ic.ops.append(Op("s", "rect", {"x": x - 5, "y": y + 2, "width": 9, "height": 16, "rx": 4.5,
                                   "transform": f"rotate(-25 {x} {y + 10})"}))
    ic.sec().rect(x + 3, y + 18, 16, 12, 3)
    ic.cut().line(x + 3, y + 22, x + 19, y + 22, 1.2)


@icon("hand_pipette", "safety", label="Hand with Pipette", tags=["hand", "gloved hand", "pipetting", "pipette"])
def _(ic):
    ic.rect(29, 2, 8, 6, 2)
    ic.path("M28 8 H38 V40 L35 46 H31 L28 40 Z")
    ic.sec().path("M31.4 48 H34.6 L33.4 60 H32.6 Z")
    _fist(ic, 26, 14)


@icon("hand_vial", "safety", label="Hand with Vial", tags=["hand", "gloved hand", "holding vial", "sample"])
def _(ic):
    ic.rect(25, 2, 16, 7, 2)
    ic.srect(26.5, 9, 13, 24, 2, 2.4, glass=True)
    ic.sec().rect(29, 18, 8, 12, 1)
    _fist(ic, 22, 26)


@icon("hand_spatula", "safety", label="Hand with Spatula", tags=["hand", "gloved hand", "spatula", "weighing"])
def _(ic):
    ic.line(30, 30, 56, 6, 2.8)
    rot_rect(ic, 57, 5, 10, 6, 3, -43)
    ic.acc().circle(57, 5, 2)
    _fist(ic, 12, 26)


# ==========================================================================
# PLOTS
# ==========================================================================

@icon("pxrd_stacked", "data", label="Stacked PXRD", tags=["stacked", "PXRD", "comparison", "offset patterns"])
def _(ic):
    axes(ic)
    for k, role in enumerate(("a", "s", "m")):
        base = 50 - k * 13
        pts = lorentz_trace(12, 56, base, [(18, 10, .5), (24, 5, .5), (33, 7, .5), (42, 4, .5), (50, 5, .5)], n=160)
        ic.spoly(pts, 2, role=role) if role != "m" else ic.spoly(pts, 2)


@icon("rietveld_fit", "data", label="Rietveld Fit", tags=["Rietveld", "refinement", "difference curve", "Yobs Ycalc"])
def _(ic):
    axes(ic)
    peaks = [(18, 26, .6), (26, 12, .6), (34, 18, .6), (44, 8, .6), (51, 10, .6)]
    for x, y in lorentz_trace(12, 56, 40, peaks, n=44):
        ic.sec().circle(x, y, 0.9)
    ic.acc().spoly(lorentz_trace(12, 56, 40, peaks, n=200), 1.8)
    for x, _, _ in peaks:
        ic.line(x, 44, x, 47, 1.4)
    ic.spoly([(12 + i, 51 + 1.2 * math.sin(i * 1.7)) for i in range(0, 45, 2)], 1.4)


def _iso(ic, d, d2=None):
    axes(ic)
    ic.acc().spath(d, 3)
    if d2:
        ic.sec().spath(d2, 2.4)


@icon("isotherm_type_i", "data", label="Isotherm (Type I)", tags=["microporous", "Langmuir", "IUPAC type I", "BET"])
def _(ic):
    _iso(ic, "M12 52 C13 22 18 16 56 14")


@icon("isotherm_type_ii", "data", label="Isotherm (Type II)", tags=["nonporous", "macroporous", "IUPAC type II"])
def _(ic):
    _iso(ic, "M12 52 C14 40 18 36 30 34 C42 32 50 26 56 8")


@icon("isotherm_type_iii", "data", label="Isotherm (Type III)", tags=["weak interaction", "IUPAC type III"])
def _(ic):
    _iso(ic, "M12 52 C30 50 46 42 56 10")


@icon("isotherm_type_v", "data", label="Isotherm (Type V)", tags=["water adsorption", "IUPAC type V", "hysteresis"])
def _(ic):
    _iso(ic, "M12 52 C24 51 30 48 34 34 C38 20 44 16 56 14", "M56 14 C48 15 44 18 40 30 C38 40 34 46 28 48")


@icon("isotherm_type_vi", "data", label="Isotherm (Type VI)", tags=["stepwise", "layer by layer", "IUPAC type VI"])
def _(ic):
    _iso(ic, "M12 52 C16 52 16 44 20 44 H24 C28 44 28 36 32 36 H36 C40 36 40 28 44 28 H48 C52 28 52 20 56 20")


@icon("tga_dtg", "data", label="TGA + DTG", tags=["TGA", "DTG", "derivative", "mass loss"])
def _(ic):
    axes(ic)
    ic.acc().spath("M12 12 H22 C28 12 26 22 32 22 H38 C44 22 42 32 56 32", 3)
    ic.sec().spoly(gauss_trace(12, 56, 46, [(27, 9, 2.2), (42, 7, 2.6)], n=140, sign=1), 2.2)


@icon("band_diagram", "data", label="Band Diagram", tags=["band gap", "VB", "CB", "semiconductor", "Eg"])
def _(ic):
    ic.rect(8, 6, 48, 14, 2)
    ic.sec().rect(8, 44, 48, 14, 2)
    ic.acc().line(32, 24, 32, 40, 2.4)
    arrowhead(ic, 32, 22, -90, 4.5)
    arrowhead(ic, 32, 42, 90, 4.5)
    ic.acc().circle(22, 13, 2.6)
    ic.cut().circle(42, 51, 2.6)


# ==========================================================================
# FIXES: clearer ESI and nanosheet; badges for look-alike plots
# ==========================================================================

@icon("esi_source", "proteomics", label="ESI Source", tags=["electrospray", "ionisation", "mass spectrometry", "Taylor cone"])
def _(ic):
    ic.rect(2, 29, 18, 6, 2)
    ic.poly([(20, 29), (27, 31.5), (27, 32.5), (20, 35)])
    ic.acc().poly([(27, 32), (34, 27), (34, 37)])
    for (x, y), r in zip(((37, 31), (39, 35.5), (41, 28), (43, 33), (45, 37), (46, 29.5), (48, 33.5)),
                         (1.7, 1.4, 1.3, 1.1, 1, 0.9, 0.8)):
        ic.sec().circle(x, y, r)
    ic.path("M52 18 H62 V46 H52 V36 H55 V28 H52 Z")


@icon("nanosheet", "materials", label="Nanosheet", tags=["2D", "nanosheet", "monolayer", "exfoliated", "MOF nanosheet"])
def _(ic):
    ic.path("M4 38 L22 22 H60 L42 38 Z")
    ic.path("M4 38 H42 V42 H4 Z")
    ic.path("M42 38 L60 22 V26 L42 42 Z")
    for k in range(5):
        for j in range(3):
            cx, cy = 18 + k * 7 - j * 4.6, 26 + j * 4
            ic.cut().poly([polar(cx, cy, 2, 30 + 60 * m) for m in range(6)])
    for x, y in ((14, 50), (30, 52), (46, 48)):
        ic.acc().circle(x, y, 2.4)


BADGES = {
    "ftir": "ir", "raman": "laser", "hplc": "column", "edx_spectrum": "electron",
    "xrd_pattern": "lattice", "uv_vis_spectrum": "sun", "pl_spectrum": "glow", "mass_spec": "ion",
}


def _badge(ic, kind):
    x, y = 53, 10
    if kind == "ir":
        ic.acc().spath("M44 50 q2 -4 4 0 t4 0 t4 0", 1.8)
    elif kind == "laser":
        ic.acc().line(46, 10, 58, 10, 2)
        ic.acc().circle(58, 10, 2.2)
    elif kind == "column":
        ic.rect(47, 4, 10, 4, 1)
        ic.sec().rect(48.5, 8, 7, 8, 1)
    elif kind == "electron":
        ic.acc().circle(x, y, 3)
        ic.ops.append(Op("s", "ellipse", {"cx": x, "cy": y, "rx": 7, "ry": 2.5, "transform": f"rotate(-30 {x} {y})"}, 1.2))
    elif kind == "lattice":
        for dx in (-4, 0, 4):
            for dy in (-3, 3):
                ic.acc().circle(x + dx, y + dy, 1.4)
    elif kind == "sun":
        ic.acc().circle(x, y, 3)
        for a in range(0, 360, 60):
            x1, y1 = polar(x, y, 4.5, a)
            x2, y2 = polar(x, y, 6.5, a)
            ic.acc().line(x1, y1, x2, y2, 1.2)
    elif kind == "glow":
        ic.sec().circle(x, y, 4.5)
        ic.acc().circle(x, y, 2.4)
    elif kind == "ion":
        ic.acc().circle(x - 3, y, 2.6)
        ic.acc().line(x + 1, y, x + 6, y, 1.4)
        ic.acc().line(x + 3.5, y - 2.5, x + 3.5, y + 2.5, 1.4)


def apply_badges():
    for name, kind in BADGES.items():
        if name not in REGISTRY or getattr(REGISTRY[name], "_badged", False):
            continue
        orig = REGISTRY[name]

        def build(orig=orig, kind=kind):
            ic = orig()
            _badge(ic, kind)
            return ic
        build._badged = True
        REGISTRY[name] = build


apply_badges()
