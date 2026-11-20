"""
v5 batch 2: sample states, lab infrastructure + original GHS pictograms, cells & biology
(viability, uptake, kinetics plots), computation & data. All original geometry (CC0).
"""
from __future__ import annotations

import math
import random

from .core import Op, polar, smooth_path
from .library import icon, axes
from .library4 import arrowhead

RED, BLACK, GREEN, PURPLE = "#E3000F", "#1F2328", "#2EAD5B", "#7B3FA6"


# ==========================================================================
# 5. SAMPLE STATES
# ==========================================================================

@icon("powder_vial", "glassware", label="Powder in Vial", tags=["powder", "sample", "MOF powder", "solid"])
def _(ic):
    ic.rect(19, 6, 26, 10, 2.5)
    ic.srect(20.75, 18.75, 22.5, 38.5, 4, 3, glass=True)
    ic.sec().path("M24.5 44 Q29 40 32 42 Q36 39 39.5 43 V52.5 Q39.5 54 38 54 H26 Q24.5 54 24.5 52.5 Z")
    for x, y in ((28, 46), (33, 47), (36, 45), (30, 50)):
        ic.acc().circle(x, y, 0.9)


@icon("pressed_pellet", "materials", label="Pressed Pellet", tags=["pellet", "KBr", "pressed disc", "sample"])
def _(ic):
    ic.sec().ellipse(32, 38, 24, 9)
    ic.sec().rect(8, 30, 48, 8)
    ic.ellipse(32, 30, 24, 9)
    ic.ops[-1].role = "s"
    ic.cut().spath("M16 28 A16 5 0 0 1 30 24", 1.4)


@icon("monolith", "materials", label="Monolith", tags=["monolith", "shaped body", "pellet", "extrudate"])
def _(ic):
    ic.sec().poly([(14, 22), (40, 22), (40, 56), (14, 56)])
    ic.poly([(14, 22), (26, 10), (52, 10), (40, 22)])
    ic.ops[-1].role = "s"
    ic.sec().poly([(40, 22), (52, 10), (52, 44), (40, 56)])
    for x in (20, 27, 34):
        for y in (30, 38, 46):
            ic.cut().circle(x, y, 2)


@icon("gel_vial", "materials", label="Gel (Inversion Test)", tags=["gel", "inverted vial", "gelation", "hydrogel"])
def _(ic):
    ic.rect(19, 48, 26, 10, 2.5)
    ic.srect(20.75, 8.75, 22.5, 38.5, 4, 3, glass=True)
    ic.sec().rect(24.5, 30, 15, 15, 2)
    arrowhead(ic, 54, 20, 60, 4)
    ic.acc().spath("M12 26 A22 22 0 0 1 52 14", 2)


@icon("aerogel", "materials", label="Aerogel", tags=["aerogel", "ultralight", "porous solid", "supercritical drying"])
def _(ic):
    ic.sec().rect(12, 22, 40, 28, 6)
    ic.ops[-1].attrs["opacity"] = 0.55
    for k in range(18):
        x = 15 + (k * 37) % 34
        y = 25 + (k * 23) % 22
        ic.cut().circle(x, y, 1.4 + (k % 3) * 0.5)
    ic.acc().poly([(28, 14), (36, 14), (34, 22), (30, 22)])
    ic.acc().ellipse(32, 12, 6, 2.4)


@icon("suspension_vial", "glassware", label="Suspension", tags=["suspension", "dispersion", "colloid", "cloudy"])
def _(ic):
    ic.rect(19, 6, 26, 10, 2.5)
    ic.srect(20.75, 18.75, 22.5, 38.5, 4, 3, glass=True)
    ic.sec().rect(24.5, 26, 15, 27.5, 2)
    ic.ops[-1].attrs["opacity"] = 0.75
    rnd = random.Random(7)
    for _ in range(16):
        ic.acc().circle(26 + rnd.random() * 11, 28 + rnd.random() * 24, 0.7 + rnd.random() * 0.6)


@icon("single_crystal", "materials", label="Single Crystal", tags=["single crystal", "faceted", "crystal habit"])
def _(ic):
    ic.sec().poly([(20, 16), (44, 16), (54, 32), (44, 50), (20, 50), (10, 32)])
    ic.poly([(20, 16), (44, 16), (38, 28), (26, 28)])
    ic.ops[-1].role = "s"
    ic.cut().spoly([(26, 28), (38, 28), (44, 50)], 1.2)
    ic.cut().spoly([(26, 28), (20, 50)], 1.2)
    ic.cut().spoly([(10, 32), (26, 28)], 1.2)
    ic.cut().spoly([(54, 32), (38, 28)], 1.2)


@icon("film_on_glass", "materials", label="Film on Substrate", tags=["film", "coating", "substrate", "glass slide", "FTO"])
def _(ic):
    ic.srect(6, 30, 52, 12, 1.5, 2.2, glass=True)
    ic.sec().rect(16, 26, 36, 4, 1)
    ic.acc().line(20, 18, 20, 25, 1.4)
    ic.acc().line(48, 18, 48, 25, 1.4)


@icon("stability_series", "materials", label="Stability: Fresh → Degraded", tags=["stability", "degradation", "activation", "as-synthesised"])
def _(ic):
    def oct(cx, s, broken=False, pores=False):
        top, bot, l, r = (cx, 32 - s), (cx, 32 + s), (cx - s * .8, 32), (cx + s * .8, 32)
        if broken:
            ic.sec().poly([top, (cx - 2, 30), l])
            ic.sec().poly([(cx + 2, 31), r, bot])
            ic.sec().poly([(cx - 3, 36), l, (cx - 1, 44)])
        else:
            ic.sec().poly([top, r, bot, l])
            if pores:
                for dy in (-5, 0, 5):
                    ic.cut().circle(cx, 32 + dy, 1.4)
    oct(9, 9)
    ic.line(19, 32, 22, 32, 1.8)
    arrowhead(ic, 24, 32, 0, 3, "m")
    oct(33, 9, pores=True)
    ic.line(43, 32, 46, 32, 1.8)
    arrowhead(ic, 48, 32, 0, 3, "m")
    oct(56, 8, broken=True)


@icon("powder_heap", "materials", label="Powder", tags=["powder", "heap", "solid sample", "grains"])
def _(ic):
    ic.sec().path("M6 52 C14 34 24 22 32 22 C40 22 50 34 58 52 Z")
    rnd = random.Random(3)
    for _ in range(14):
        x = 14 + rnd.random() * 36
        y = 34 + rnd.random() * 16
        ic.cut().circle(x, y, 0.9)
    ic.rect(4, 52, 56, 4, 2)


# ==========================================================================
# 6. LAB INFRASTRUCTURE & GHS
# ==========================================================================

@icon("lab_sink", "infrastructure", label="Sink", tags=["sink", "tap", "washing", "lab bench"])
def _(ic):
    ic.rect(4, 30, 56, 8, 2)
    ic.rect(8, 38, 48, 20, 3)
    ic.cut().rect(14, 30, 36, 5, 2)
    ic.spath("M42 30 V16 Q42 10 36 10 H30 V16", 3)
    ic.sec().path("M30 18 C31.5 21 31.5 23 30 23 C28.5 23 28.5 21 30 18 Z")
    ic.rect(46, 24, 6, 6, 1.5)


@icon("lab_shelves", "infrastructure", label="Shelves", tags=["shelves", "storage", "reagents", "lab"])
def _(ic):
    ic.rect(6, 4, 4, 56, 1)
    ic.rect(54, 4, 4, 56, 1)
    for y in (22, 40, 58):
        ic.rect(6, y - 3, 52, 4, 1)
    for x, w, h, role in ((14, 7, 12, "s"), (24, 6, 9, "a"), (34, 8, 14, "s"), (46, 6, 10, "m")):
        ic.rect(x, 19 - h, w, h, 1.5, role)
    for x, w, h, role in ((14, 10, 10, "m"), (28, 6, 14, "s"), (40, 9, 12, "a")):
        ic.rect(x, 37 - h, w, h, 1.5, role)


@icon("chemical_cabinet", "infrastructure", label="Safety Cabinet", tags=["flammables cabinet", "storage", "safety cabinet", "solvents"])
def _(ic):
    ic.rect(10, 4, 44, 56, 3)
    ic.ops[-1].attrs["_fill"] = "#F4C542"
    ic.cut().line(32, 6, 32, 58, 1.6)
    ic.rect(27, 26, 3, 10, 1)
    ic.rect(34, 26, 3, 10, 1)
    ic.poly([(32, 42), (40, 50), (32, 58), (24, 50)]).fixed(RED)
    ic.cut().poly([(32, 45), (37, 50), (32, 55), (27, 50)])
    ic.path("M32 46.5 C34 49 34 51 32 52.5 C30 51 30 49 32 46.5 Z").fixed(BLACK)


@icon("waste_container", "infrastructure", label="Solvent Waste", tags=["waste", "solvent waste", "jerrycan", "disposal"])
def _(ic):
    ic.rect(36, 4, 10, 8, 2)
    ic.path("M12 14 H52 Q56 14 56 18 V56 Q56 60 52 60 H12 Q8 60 8 56 V18 Q8 14 12 14 Z")
    ic.rect(14, 8, 14, 8, 3)
    ic.cut().rect(16, 26, 32, 20, 2)
    ic.sec().rect(19, 29, 26, 14, 1)


@icon("waste_halogenated", "infrastructure", label="Halogenated Waste", tags=["halogenated waste", "chlorinated solvents", "DCM", "disposal"])
def _(ic):
    ic.rect(36, 4, 10, 8, 2)
    ic.path("M12 14 H52 Q56 14 56 18 V56 Q56 60 52 60 H12 Q8 60 8 56 V18 Q8 14 12 14 Z")
    ic.rect(14, 8, 14, 8, 3)
    ic.cut().rect(16, 26, 32, 20, 2)
    ic.acc().rect(19, 29, 26, 14, 1)
    ic.cut().line(26, 32, 38, 40, 2.4)
    ic.cut().line(38, 32, 26, 40, 2.4)


@icon("safety_shower", "infrastructure", label="Safety Shower", tags=["safety shower", "emergency", "decontamination"])
def _(ic):
    ic.rect(10, 4, 5, 56, 1.5)
    ic.rect(10, 4, 30, 5, 1.5)
    ic.ellipse(40, 12, 12, 4)
    for x in (32, 37, 42, 47):
        ic.sec().line(x, 18, x - 2, 30, 1.6)
    ic.rect(15, 32, 10, 3, 1)
    ic.acc().rect(23, 30, 3, 18, 1)
    ic.rect(6, 58, 40, 3, 1)


@icon("eye_wash", "infrastructure", label="Eye Wash", tags=["eye wash", "emergency", "safety station"])
def _(ic):
    ic.rect(29, 34, 6, 26, 1.5)
    ic.rect(20, 56, 24, 4, 1.5)
    ic.path("M12 26 H52 L46 36 H18 Z")
    for x in (24, 40):
        ic.rect(x - 3, 20, 6, 6, 1.5)
        ic.sec().spath(f"M{x - 2} 18 Q{x} 10 {x + 2} 18", 1.8)
    ic.acc().rect(42, 40, 10, 4, 2)


@icon("ln2_dewar", "infrastructure", label="Liquid-N₂ Dewar", tags=["liquid nitrogen", "dewar", "cryogenic", "LN2"])
def _(ic):
    ic.rect(14, 16, 36, 44, 8)
    ic.rect(22, 8, 20, 9, 3)
    ic.rect(18, 4, 28, 5, 2)
    ic.cut().rect(20, 28, 24, 18, 2)
    ic.acc().line(32, 31, 32, 43, 2)
    ic.acc().line(26.8, 34, 37.2, 40, 2)
    ic.acc().line(26.8, 40, 37.2, 34, 2)
    ic.sec().spath("M30 4 q-2 -2 0 -4 M36 4 q2 -2 0 -4", 1.4)


@icon("cryostat", "infrastructure", label="Cryostat", tags=["cryostat", "low temperature", "helium", "variable temperature"])
def _(ic):
    ic.rect(20, 4, 24, 10, 3)
    ic.rect(22, 14, 20, 30, 2)
    ic.rect(26, 44, 12, 12, 2)
    ic.acc().rect(28, 50, 8, 4, 1)
    ic.line(44, 9, 56, 9, 2.4)
    ic.line(56, 9, 56, 30, 2.4)
    ic.rect(50, 30, 12, 10, 2)
    for y in (20, 26, 32, 38):
        ic.cut().line(24, y, 40, y, 1.4)


@icon("manometer", "infrastructure", label="Pressure Gauge", tags=["manometer", "pressure gauge", "vacuum gauge", "bar"])
def _(ic):
    ic.rect(28, 50, 8, 12, 1.5)
    ic.circle(32, 30, 22)
    ic.cut().circle(32, 30, 18)
    for k in range(9):
        x1, y1 = polar(32, 30, 14, 135 + k * 33.75)
        x2, y2 = polar(32, 30, 16.5, 135 + k * 33.75)
        ic.line(x1, y1, x2, y2, 1.4)
    ic.acc().path("M32 30 L45 20 L33.5 32 Z")
    ic.circle(32, 30, 2.6)
    ic.acc().spath("M20 42 A15 15 0 0 0 44 42", 1.6)


def ghs(ic):
    ic.poly([(32, 2.5), (61.5, 32), (32, 61.5), (2.5, 32)]).fixed(RED)
    ic.cut().poly([(32, 8.5), (55.5, 32), (32, 55.5), (8.5, 32)])


def k(ic, op_builder):
    op_builder.fixed(BLACK)


@icon("ghs_explosive", "safety", label="GHS01 Explosive", tags=["GHS01", "explosive", "hazard", "pictogram"])
def _(ic):
    ghs(ic)
    ic.circle(30, 36, 8).fixed(BLACK)
    ic.line(35, 30, 39, 25, 2).fixed(BLACK)
    for a in (-80, -40, 0, 40):
        x1, y1 = polar(40, 24, 3, a)
        x2, y2 = polar(40, 24, 8, a)
        ic.line(x1, y1, x2, y2, 1.8).fixed(BLACK)
    for x, y in ((20, 24), (44, 42), (20, 46)):
        ic.poly([(x, y - 2.5), (x + 2.5, y), (x, y + 2.5), (x - 2.5, y)]).fixed(BLACK)


@icon("ghs_flammable", "safety", label="GHS02 Flammable", tags=["GHS02", "flammable", "fire", "hazard", "pictogram"])
def _(ic):
    ghs(ic)
    ic.path("M32 16 C35 24 42 28 42 37 A10 10 0 0 1 22 37 C22 31 26 29 27 24 C29 28 30 29 30 31 C33 27 33 22 32 16 Z").fixed(BLACK)
    ic.rect(20, 46, 24, 3, 1).fixed(BLACK)


@icon("ghs_oxidizing", "safety", label="GHS03 Oxidising", tags=["GHS03", "oxidiser", "oxidizing", "hazard", "pictogram"])
def _(ic):
    ghs(ic)
    ic.path("M32 14 C35 20 40 23 40 30 A8 8 0 0 1 24 30 C24 25 27 24 28 20 C30 23 31 24 31 26 C33 23 33 19 32 14 Z").fixed(BLACK)
    ic.scircle(32, 40, 7, 3).fixed(BLACK)
    ic.rect(20, 49, 24, 3, 1).fixed(BLACK)


@icon("ghs_gas", "safety", label="GHS04 Gas Under Pressure", tags=["GHS04", "compressed gas", "gas cylinder", "hazard", "pictogram"])
def _(ic):
    ghs(ic)
    ic.rect(18, 30, 26, 9, 4.5).fixed(BLACK)
    ic.rect(44, 32, 3, 5, 1).fixed(BLACK)
    ic.rect(47, 31, 3, 7, 1).fixed(BLACK)


@icon("ghs_corrosive", "safety", label="GHS05 Corrosive", tags=["GHS05", "corrosive", "acid", "hazard", "pictogram"])
def _(ic):
    ghs(ic)
    for x0, rot in ((20, -30), (36, -30)):
        ic.ops.append(Op("m", "rect", {"x": x0, "y": 16, "width": 5, "height": 12, "rx": 1,
                                       "transform": f"rotate({rot} {x0 + 2.5} 22)", "_fill": BLACK}))
    ic.path("M25 30 C26 32 26 33 25 33.5 C24 33 24 32 25 30 Z").fixed(BLACK)
    ic.path("M41 30 C42 32 42 33 41 33.5 C40 33 40 32 41 30 Z").fixed(BLACK)
    ic.rect(16, 40, 13, 4, 1).fixed(BLACK)
    ic.path("M34 38 H46 V44 H38 L34 41 Z").fixed(BLACK)
    ic.rect(14, 46, 36, 3, 1).fixed(BLACK)


@icon("ghs_toxic", "safety", label="GHS06 Acute Toxicity", tags=["GHS06", "toxic", "skull", "poison", "hazard", "pictogram"])
def _(ic):
    ghs(ic)
    ic.line(19, 44, 45, 52, 3).fixed(BLACK)
    ic.line(45, 44, 19, 52, 3).fixed(BLACK)
    ic.path("M32 15 C40 15 44 20 44 27 C44 31 41 33 40 34 V38 H24 V34 C23 33 20 31 20 27 C20 20 24 15 32 15 Z").fixed(BLACK)
    ic.cut().circle(27, 27, 3)
    ic.cut().circle(37, 27, 3)
    ic.cut().poly([(32, 31), (34, 34), (30, 34)])


@icon("ghs_harmful", "safety", label="GHS07 Harmful / Irritant", tags=["GHS07", "exclamation mark", "irritant", "harmful", "hazard", "pictogram"])
def _(ic):
    ghs(ic)
    ic.path("M28.5 15 H35.5 L34 37 H30 Z").fixed(BLACK)
    ic.circle(32, 44, 3.6).fixed(BLACK)


@icon("ghs_health", "safety", label="GHS08 Health Hazard", tags=["GHS08", "health hazard", "carcinogen", "CMR", "hazard", "pictogram"])
def _(ic):
    ghs(ic)
    ic.circle(32, 18, 4).fixed(BLACK)
    ic.path("M23 50 V30 Q23 24 29 24 H35 Q41 24 41 30 V50 Z").fixed(BLACK)
    pts = []
    for kk in range(16):
        r = 7 if kk % 2 == 0 else 3.2
        pts.append(polar(32, 35, r, kk * 22.5))
    ic.cut().poly(pts)


@icon("ghs_environment", "safety", label="GHS09 Environment", tags=["GHS09", "aquatic toxicity", "environment", "hazard", "pictogram"])
def _(ic):
    ghs(ic)
    ic.line(22, 46, 22, 22, 2.4).fixed(BLACK)
    ic.line(22, 30, 16, 24, 1.8).fixed(BLACK)
    ic.line(22, 34, 28, 28, 1.8).fixed(BLACK)
    ic.rect(14, 46, 36, 2.5, 1).fixed(BLACK)
    ic.path("M32 42 C36 37 44 37 47 42 C44 46 36 46 32 42 Z").fixed(BLACK)
    ic.poly([(47, 42), (51, 38.5), (51, 45.5)]).fixed(BLACK)
    ic.cut().line(36, 40, 40, 44, 1)
    ic.cut().line(40, 40, 36, 44, 1)


# ==========================================================================
# 7. CELLS & BIOLOGY (proteomics, nanozymes, nanomedicine)
# ==========================================================================

def cell(ic, cx, cy, r, role="s", nucleus=True, fill=None):
    pts = [polar(cx, cy, r * (1 + 0.08 * math.sin(k * 1.7)), k * 30) for k in range(12)]
    ic.path(smooth_path(pts + [pts[0], pts[1]]) + " Z", role=role)
    if fill:
        ic.fixed(fill)
    if nucleus:
        ic.cut().circle(cx + r * 0.15, cy - r * 0.05, r * 0.42)
        ic.circle(cx + r * 0.15, cy - r * 0.05, r * 0.3, role="a")


@icon("animal_cell", "biology", label="Cell", tags=["cell", "eukaryotic", "mammalian cell", "cell culture"])
def _(ic):
    cell(ic, 32, 32, 25)
    for x, y in ((20, 24), (22, 42), (44, 46)):
        ic.ellipse(x, y, 3.5, 2)


@icon("bacterium", "biology", label="Bacterium", tags=["bacteria", "E. coli", "antibacterial", "microbe"])
def _(ic):
    for a in range(0, 360, 30):
        x1, y1 = polar(32, 32, 16, a)
        x2, y2 = polar(32, 32, 22, a)
        ic.line(x1, y1, x2, y2, 1.2)
    ic.ops.append(Op("s", "rect", {"x": 12, "y": 22, "width": 40, "height": 20, "rx": 10, "transform": "rotate(-25 32 32)"}))
    ic.spath("M52 22 C58 18 56 12 62 8", 1.6)
    ic.ops.append(Op("a", "path", {"d": smooth_path([(22, 34), (28, 28), (34, 34), (40, 28)]), "transform": "rotate(-25 32 32)"}, 1.6))


@icon("lipid_bilayer", "biology", label="Lipid Bilayer", tags=["membrane", "bilayer", "phospholipid", "liposome"])
def _(ic):
    for x in range(8, 60, 7):
        for y, d in ((14, 1), (50, -1)):
            ic.line(x - 1, y + d * 5, x - 1, y + d * 15, 1.4)
            ic.line(x + 1, y + d * 5, x + 1, y + d * 15, 1.4)
            ic.sec().circle(x, y, 3.2)
    ic.acc().rect(28, 10, 8, 44, 3)


@icon("antibody", "biology", label="Antibody", tags=["antibody", "IgG", "immunoassay", "ELISA"])
def _(ic):
    ic.sec().rect(29, 32, 6, 26, 3)
    ic.ops.append(Op("s", "rect", {"x": 29, "y": 6, "width": 6, "height": 28, "rx": 3, "transform": "rotate(-38 32 32)"}))
    ic.ops.append(Op("s", "rect", {"x": 29, "y": 6, "width": 6, "height": 28, "rx": 3, "transform": "rotate(38 32 32)"}))
    ic.ops.append(Op("a", "rect", {"x": 23, "y": 4, "width": 5, "height": 16, "rx": 2.5, "transform": "rotate(-38 32 32)"}))
    ic.ops.append(Op("a", "rect", {"x": 36, "y": 4, "width": 5, "height": 16, "rx": 2.5, "transform": "rotate(38 32 32)"}))
    ic.circle(32, 32, 3.5)


@icon("enzyme_substrate", "biology", label="Enzyme–Substrate Complex", tags=["enzyme", "substrate", "lock and key", "active site", "binding"])
def _(ic):
    ic.sec().path("M6 50 C4 34 10 18 24 14 C30 13 30 24 36 24 C42 24 40 12 48 14 C58 18 60 36 56 50 C50 58 14 58 6 50 Z")
    ic.acc().path("M30 28 C30 16 42 16 42 28 C42 30 30 30 30 28 Z")
    ic.acc().circle(36, 18, 1.8)


@icon("protein_digestion", "biology", label="Protein → Peptides", tags=["digestion", "proteolysis", "bottom-up proteomics", "peptides"])
def _(ic):
    ic.sec().path(smooth_path([(4, 32), (6, 20), (16, 14), (24, 20), (26, 32), (18, 42), (8, 42), (4, 32)]))
    ic.line(30, 30, 36, 30, 2)
    arrowhead(ic, 39, 30, 0, 3.5, "m")
    for k, (x, y) in enumerate(((44, 14), (54, 22), (44, 32), (55, 40), (46, 50))):
        ic.sec().line(x - 3, y, x + 3, y, 2)
        ic.circle(x - 3, y, 2.2)
        ic.acc().circle(x + 3, y, 2.2)


@icon("charged_peptides", "biology", label="Peptide Ions", tags=["peptide ions", "charge states", "ESI", "MS"])
def _(ic):
    for k, (y, z) in enumerate(((14, 1), (32, 2), (50, 3))):
        for j in range(4 + k):
            ic.sec().circle(8 + j * 7, y, 3) if j % 2 else ic.circle(8 + j * 7, y, 3)
        x = 8 + (4 + k) * 7
        for c in range(z):
            ic.acc().circle(x + c * 6, y - 4, 2.2)
            ic.cut().line(x + c * 6 - 1.2, y - 4, x + c * 6 + 1.2, y - 4, 0.8)
            ic.cut().line(x + c * 6, y - 5.2, x + c * 6, y - 2.8, 0.8)


@icon("mof_drug_release", "biology", label="Drug Release", tags=["drug delivery", "release", "cargo", "nanocarrier", "MOF"])
def _(ic):
    ic.circle(26, 34, 18)
    for x, y in ((18, 28), (26, 24), (34, 30), (22, 40), (30, 42)):
        ic.cut().circle(x, y, 3)
        ic.acc().circle(x, y, 1.8)
    for (x, y), a in (((50, 18), -30), ((54, 36), 0), ((48, 52), 30)):
        ic.acc().circle(x, y, 2.4)
        sx, sy = polar(26, 34, 21, a)
        ic.sec().line(sx, sy, x - 3, y + (34 - y) * 0.05, 1.4)


@icon("cell_uptake", "biology", label="Cellular Uptake", tags=["endocytosis", "internalisation", "nanoparticle uptake", "cell"])
def _(ic):
    ic.sec().path("M4 60 V30 C14 26 18 38 26 40 C34 42 38 30 46 28 C54 26 58 34 60 34 V60 Z")
    ic.cut().circle(26, 46, 5.5)
    ic.acc().circle(26, 46, 3.2)
    ic.cut().circle(44, 50, 4)
    ic.acc().circle(44, 50, 2.4)
    ic.acc().circle(26, 18, 3.2)
    ic.line(26, 23, 26, 30, 1.6)
    arrowhead(ic, 26, 33, 90, 3, "m")


@icon("cell_viability", "biology", label="Cell Viability", tags=["viability", "cytotoxicity", "MTT", "biocompatibility"])
def _(ic):
    axes(ic)
    for x, h in ((16, 40), (26, 38), (36, 34), (46, 20)):
        ic.sec().rect(x - 3.5, 55 - h, 7, h - 1.75, 1)
    ic.line(12, 14, 56, 14, 1.2)
    ic.ops[-1].attrs["stroke-dasharray"] = "3 2"


@icon("live_dead_staining", "biology", label="Live / Dead Staining", tags=["live/dead", "calcein", "propidium iodide", "fluorescence"])
def _(ic):
    ic.rect(4, 4, 56, 56, 6)
    ic.ops[-1].attrs["_fill"] = "#14161A"
    for k, (x, y, r) in enumerate(((16, 16, 6), (34, 14, 5), (48, 24, 6), (20, 36, 5.5), (40, 40, 6), (26, 52, 4.5), (50, 50, 5))):
        cell(ic, x, y, r, nucleus=False, fill=GREEN if k not in (2, 5) else RED)


@icon("mtt_assay", "biology", label="MTT / Colorimetric Assay", tags=["MTT", "colorimetric assay", "formazan", "well plate"])
def _(ic):
    ic.rect(4, 13, 56, 38, 4)
    shades = ["#F5EFFA", "#E2CFF0", "#C9A8E3", "#A97BD0", "#8A52BA", "#6E2FA0"]
    for r in range(4):
        for c in range(6):
            x, y = 12 + c * 8, 21 + r * 8
            ic.cut().circle(x, y, 3.2)
            ic.circle(x, y, 2.6).fixed(shades[c])


@icon("ros_generation", "biology", label="ROS Generation", tags=["ROS", "reactive oxygen species", "radicals", "oxidative stress"])
def _(ic):
    ic.circle(22, 34, 12)
    ic.cut().circle(19, 31, 2.4)
    for (x, y) in ((44, 14), (52, 30), (46, 48)):
        ic.acc().circle(x, y, 4.5)
        ic.sec().circle(x + 5, y - 3, 3.2)
        ic.circle(x + 7.5, y - 7.5, 1.3)
    for a in (-40, 0, 40):
        x1, y1 = polar(22, 34, 15, a)
        x2, y2 = polar(22, 34, 20, a)
        ic.acc().line(x1, y1, x2, y2, 1.6)


@icon("fluorescence_imaging", "biology", label="Fluorescence Imaging", tags=["confocal", "fluorescence microscopy", "cell imaging", "staining"])
def _(ic):
    ic.rect(4, 4, 56, 56, 6)
    ic.ops[-1].attrs["_fill"] = "#14161A"
    cell(ic, 30, 32, 20, nucleus=False, fill="#2B7BE4")
    ic.circle(32, 30, 8).fixed("#33D17A")
    for x, y in ((18, 22), (40, 44), (22, 44)):
        ic.circle(x, y, 2).fixed("#FF4FA3")


@icon("lineweaver_burk", "data", label="Lineweaver–Burk Plot", tags=["Lineweaver-Burk", "double reciprocal", "Km", "Vmax", "enzyme kinetics"])
def _(ic):
    ic.line(8, 46, 58, 46, 2)
    ic.line(24, 6, 24, 58, 2)
    ic.acc().line(10, 50, 56, 18, 2.6)
    for x in (30, 38, 46, 54):
        ic.circle(x, 50 - (x - 10) * 32 / 46, 2)


@icon("dose_response", "data", label="Dose–Response / IC₅₀", tags=["dose response", "IC50", "EC50", "inhibition", "sigmoid"])
def _(ic):
    axes(ic)
    pts = [(12 + i, 14 + 36 / (1 + math.exp(-(i - 22) / 4))) for i in range(0, 45)]
    ic.acc().spoly(pts, 2.6)
    ic.sec().line(34, 32, 34, 55, 1.2)
    ic.ops[-1].attrs["stroke-dasharray"] = "2.5 2"
    for i in (4, 12, 20, 26, 34, 42):
        ic.circle(12 + i, 14 + 36 / (1 + math.exp(-(i - 22) / 4)), 1.8)


@icon("liposome", "biology", label="Liposome / Vesicle", tags=["liposome", "vesicle", "nanocarrier", "bilayer"])
def _(ic):
    for k in range(20):
        a = k * 18
        ic.sec().circle(*polar(32, 32, 25, a), 2.6)
        ic.sec().circle(*polar(32, 32, 15, a + 9), 2.2)
        x1, y1 = polar(32, 32, 22, a)
        x2, y2 = polar(32, 32, 18, a + 9)
        ic.line(x1, y1, x2, y2, 1)
    ic.acc().circle(30, 30, 3)
    ic.acc().circle(36, 36, 2.2)


# ==========================================================================
# 8. COMPUTATION & DATA
# ==========================================================================

@icon("dft_unit_cell", "computation", label="DFT / Electron Density", tags=["DFT", "electron density", "simulation", "quantum chemistry"])
def _(ic):
    f = [(12, 24), (40, 24), (40, 52), (12, 52)]
    b = [(x + 12, y - 12) for x, y in f]
    for i in range(4):
        for p, q in ((f[i], f[(i + 1) % 4]), (b[i], b[(i + 1) % 4]), (f[i], b[i])):
            ic.line(*p, *q, 1.4)
    for r, op in ((15, 0.25), (10, 0.45), (5.5, 0.8)):
        ic.ops.append(Op("s", "ellipse", {"cx": 32, "cy": 32, "rx": r, "ry": r * 0.8, "opacity": op}))
    ic.acc().circle(32, 32, 3)


@icon("supercell", "computation", label="Supercell", tags=["supercell", "periodic boundary", "2x2x2", "simulation cell"])
def _(ic):
    for i in range(3):
        ic.line(8 + i * 24, 8, 8 + i * 24, 56, 1.6)
        ic.line(8, 8 + i * 24, 56, 8 + i * 24, 1.6)
    for i in range(2):
        for j in range(2):
            cx, cy = 20 + i * 24, 20 + j * 24
            ic.acc().circle(cx, cy, 4.5) if (i, j) != (0, 0) else ic.sec().circle(cx, cy, 4.5)
            ic.circle(cx - 7, cy + 6, 2.4)


@icon("computing_cluster", "computation", label="Computing Cluster", tags=["HPC", "server", "cluster", "calculation"])
def _(ic):
    for k, x in enumerate((6, 24, 42)):
        ic.rect(x, 6, 16, 52, 2)
        for y in range(10, 54, 8):
            ic.cut().rect(x + 2, y, 12, 5, 1)
            ic.acc().circle(x + 4.5, y + 2.5, 1) if (y // 8 + k) % 2 else ic.sec().circle(x + 4.5, y + 2.5, 1)


@icon("neural_network", "computation", label="Neural Network / ML", tags=["machine learning", "neural network", "AI", "model"])
def _(ic):
    layers = [(10, (16, 32, 48)), (32, (10, 24, 40, 54)), (54, (24, 40))]
    for (x1, ys1), (x2, ys2) in zip(layers, layers[1:]):
        for y1 in ys1:
            for y2 in ys2:
                ic.line(x1, y1, x2, y2, 1)
    for k, (x, ys) in enumerate(layers):
        for y in ys:
            ic.circle(x, y, 4.5, role=("s", "m", "a")[k])


@icon("database", "computation", label="Database", tags=["database", "data", "repository", "storage"])
def _(ic):
    for y in (44, 30, 16):
        ic.ellipse(32, y + 8, 22, 7)
        ic.rect(10, y, 44, 8)
        ic.ellipse(32, y, 22, 7)
        ic.ops[-1].role = "s"
    ic.acc().circle(46, 46, 1.8)


@icon("spreadsheet", "computation", label="Spreadsheet / Table", tags=["spreadsheet", "table", "data", "Excel"])
def _(ic):
    ic.rect(6, 8, 52, 48, 4)
    ic.sec().rect(6, 8, 52, 10, 4)
    ic.sec().rect(6, 14, 52, 4)
    for x in (20, 34, 46):
        ic.cut().line(x, 18, x, 56, 1.2)
    for y in (26, 34, 42, 50):
        ic.cut().line(6, y, 58, y, 1.2)


@icon("sample_label", "computation", label="Sample Label (QR)", tags=["sample ID", "QR code", "label", "tracking"])
def _(ic):
    ic.rect(4, 12, 56, 40, 4)
    ic.cut().rect(8, 16, 32, 32, 2)
    rnd = random.Random(11)
    for r in range(8):
        for c in range(8):
            corner = (r < 3 and c < 3) or (r < 3 and c > 4) or (r > 4 and c < 3)
            if corner or rnd.random() < 0.45:
                ic.rect(10 + c * 3.5, 18 + r * 3.5, 3.2, 3.2)
    for y in (20, 28, 36, 44):
        ic.cut().rect(44, y, 12, 2.6, 1.3)


@icon("ml_workflow", "computation", label="Data → Model → Prediction", tags=["workflow", "machine learning", "prediction", "pipeline"])
def _(ic):
    ic.rect(4, 22, 14, 20, 2)
    for y in (26, 31, 36):
        ic.cut().line(6, y, 16, y, 1.2)
    ic.line(19, 32, 23, 32, 1.8)
    arrowhead(ic, 25, 32, 0, 3, "m")
    for (x, y) in ((30, 24), (30, 40), (38, 32)):
        ic.sec().circle(x, y, 3.6)
    ic.line(30, 24, 38, 32, 1)
    ic.line(30, 40, 38, 32, 1)
    ic.line(43, 32, 47, 32, 1.8)
    arrowhead(ic, 49, 32, 0, 3, "m")
    ic.acc().poly([polar(56, 32, 7, -90 + 72 * k + (36 if j else 0)) if not j else polar(56, 32, 3, -90 + 72 * k + 36)
                   for k in range(5) for j in (0, 1)])


@icon("molecular_dynamics", "computation", label="Molecular Dynamics", tags=["MD simulation", "trajectory", "force field", "dynamics"])
def _(ic):
    ic.srect(4, 4, 56, 56, 4, 1.6)
    rnd = random.Random(5)
    for k in range(7):
        x, y = 10 + rnd.random() * 44, 10 + rnd.random() * 44
        pts = [(x, y)]
        for _ in range(5):
            x += rnd.uniform(-5, 5)
            y += rnd.uniform(-5, 5)
            pts.append((max(6, min(58, x)), max(6, min(58, y))))
        ic.acc().spoly(pts, 1.3) if k % 2 else ic.sec().spoly(pts, 1.3)
        ic.circle(*pts[-1], 2.8, role="a" if k % 2 else "s")
