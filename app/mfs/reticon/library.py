"""
Built-in icon set: labware, synthesis equipment, characterisation techniques,
materials and proteomics. Every icon is a plain function returning an Icon, so
adding one is ~10 lines. Register with @icon.
"""
from __future__ import annotations

import math

from .core import Icon, polar, smooth_path, lorentz_trace, _f

REGISTRY: dict[str, callable] = {}


def icon(name, category, label="", tags=()):
    def deco(fn):
        def build():
            ic = Icon(name, label, category, tags)
            fn(ic)
            return ic
        REGISTRY[name] = build
        return fn
    return deco


def axes(ic, w=3.5):
    ic.spath("M9 7 V55 H58", w)


# ==========================================================================
# GLASSWARE
# ==========================================================================

@icon("round_bottom_flask", "glassware", tags=["flask", "synthesis", "reflux"])
def _(ic):
    ic.spath("M27 6 V22.7 A18 18 0 1 0 37 22.7 V6", 3.5, glass=True)
    ic.line(23.5, 6, 40.5, 6, 3.5)
    ic.sec().path("M18.65 42 A13.5 13.5 0 0 0 45.35 42 Z")
    ic.acc().circle(28, 35.5, 2)
    ic.acc().circle(35, 31, 1.5)


@icon("erlenmeyer_flask", "glassware", tags=["flask", "conical"])
def _(ic):
    ic.spath("M27 6 V22 L12.5 52.5 Q10.5 58 16.5 58 H47.5 Q53.5 58 51.5 52.5 L37 22 V6", 3.5, glass=True)
    ic.line(23.5, 6, 40.5, 6, 3.5)
    ic.sec().path("M20.6 40 H43.4 L48 50.6 Q49 53.5 46 53.5 H18 Q15 53.5 16 50.6 Z")
    ic.acc().circle(29, 33, 1.8)
    ic.acc().circle(34, 28, 1.3)


@icon("beaker", "glassware", tags=["beaker", "solution"])
def _(ic):
    ic.spath("M11 8 Q15.5 8 16.5 12 V52 Q16.5 58 22.5 58 H41.5 Q47.5 58 47.5 52 V8", 3.5, glass=True)
    ic.sec().rect(21, 30, 22, 23.5, 2.5)
    for y in (16, 22):
        ic.line(22, y, 28, y, 2.5)


@icon("test_tube", "glassware", tags=["tube", "sample"])
def _(ic):
    ic.spath("M24 7 V47 A8 8 0 0 0 40 47 V7", 3.5, glass=True)
    ic.line(20.5, 7, 43.5, 7, 3.5)
    ic.sec().path("M28.5 30 V47 A3.5 3.5 0 0 0 35.5 47 V30 Z")


@icon("vial", "glassware", tags=["vial", "sample", "storage"])
def _(ic):
    ic.rect(17, 5, 30, 13, 3)
    ic.cut().line(23, 8.5, 23, 14.5, 2)
    ic.cut().line(29, 8.5, 29, 14.5, 2)
    ic.cut().line(35, 8.5, 35, 14.5, 2)
    ic.cut().line(41, 8.5, 41, 14.5, 2)
    ic.srect(19.75, 21.75, 24.5, 34.5, 4, 3.5, glass=True)
    ic.sec().rect(24, 34, 16, 18, 2)


@icon("micropipette", "glassware", tags=["pipette", "dosing", "liquid handling"])
def _(ic):
    ic.rect(27.5, 3, 9, 7, 2.5)
    ic.rect(30.5, 9, 3, 4)
    ic.path("M23 13 H41 V36 L35.5 46 H28.5 L23 36 Z")
    ic.rect(39, 15, 8, 5, 2.5)
    ic.cut().rect(27, 18, 10, 9, 1.5)
    ic.sec().path("M29.2 48.5 H34.8 L32.8 61 H31.2 Z")


# ==========================================================================
# SYNTHESIS & SAMPLE PREP
# ==========================================================================

@icon("autoclave", "synthesis", tags=["solvothermal", "hydrothermal", "teflon", "MOF synthesis"])
def _(ic):
    ic.rect(9, 5, 46, 11, 2.5)
    for x in (17, 32, 47):
        ic.cut().circle(x, 10.5, 2.2)
    ic.rect(13, 18, 38, 41, 4)
    ic.cut().rect(20, 23, 24, 31, 3)
    ic.sec().rect(23.5, 26.5, 17, 24, 2)


@icon("oven", "synthesis", tags=["drying", "activation", "heating"])
def _(ic):
    ic.rect(5, 9, 54, 43, 5)
    ic.cut().rect(11, 15, 31, 31, 3)
    ic.sec().rect(14.5, 34, 24, 3.5, 1)
    ic.sec().rect(17, 27, 8, 7, 1.5)
    ic.sec().rect(28, 29, 8, 5, 1.5)
    ic.cut().circle(50.5, 21, 3.5)
    ic.cut().circle(50.5, 33, 3.5)
    ic.rect(11, 52, 7, 5, 1.5)
    ic.rect(46, 52, 7, 5, 1.5)


@icon("tube_furnace", "synthesis", tags=["calcination", "annealing", "pyrolysis"])
def _(ic):
    ic.rect(11, 20, 42, 24, 4)
    for x in (19, 25, 31, 37, 43):
        ic.cut().line(x, 37, x, 41, 2)
    ic.cut().rect(1, 25, 62, 11, 5.5)
    ic.sec().rect(2.5, 27.5, 59, 6, 3)
    ic.rect(15, 44, 5, 12, 1.5)
    ic.rect(44, 44, 5, 12, 1.5)
    for x in (24, 32, 40):
        ic.acc().spath(f"M{x} 16 q-3 -3 0 -6 q3 -3 0 -6", 2.5)


@icon("centrifuge", "synthesis", tags=["separation", "washing", "pellet"])
def _(ic):
    ic.path("M9 24 Q32 2 55 24 Z")
    ic.cut().ellipse(32, 17.5, 10, 3.6)
    ic.rect(7, 27, 50, 29, 5)
    ic.cut().rect(14, 35, 16, 9, 2)
    ic.cut().circle(44, 40, 6)
    ic.circle(44, 40, 2.4)
    ic.sec().rect(16.5, 37.5, 11, 4, 1)


@icon("glovebox", "synthesis", tags=["inert", "air-free", "argon", "nitrogen"])
def _(ic):
    ic.rect(4, 10, 56, 39, 4)
    ic.cut().rect(9, 15, 46, 15, 2.5)
    ic.cut().circle(21, 40, 5.5)
    ic.cut().circle(43, 40, 5.5)
    ic.sec().circle(21, 40, 3)
    ic.sec().circle(43, 40, 3)
    ic.sec().rect(13, 22, 9, 6, 1.5)
    ic.rect(8, 49, 5, 9, 1.5)
    ic.rect(51, 49, 5, 9, 1.5)


@icon("analytical_balance", "synthesis", tags=["weighing", "mass"])
def _(ic):
    ic.path("M7 58 L11.5 44 H52.5 L57 58 Z")
    ic.cut().rect(23, 49, 18, 5, 1.2)
    ic.rect(30, 38, 4, 7)
    ic.rect(13, 34.5, 38, 4.5, 2.25)
    ic.sec().path("M22 33 L24.5 27 H39.5 L42 33 Z")


@icon("hotplate_stirrer", "synthesis", tags=["stirring", "heating", "magnetic"])
def _(ic):
    ic.spath("M20 6 V27.5 Q20 31 23.5 31 H40.5 Q44 31 44 27.5 V6", 3.5, glass=True)
    ic.sec().rect(23.6, 15, 16.8, 12.6, 1.6)
    ic.cut().rect(26.5, 23, 11, 2.6, 1.3)
    ic.sec().rect(9, 34, 46, 4, 1.5)
    ic.rect(7, 40, 50, 18, 4)
    ic.cut().circle(20, 49, 4)
    ic.cut().circle(44, 49, 4)


@icon("buchner_filtration", "synthesis", tags=["filtration", "vacuum", "washing"])
def _(ic):
    ic.path("M9 6 H55 V18 Q55 23 50 23 H38 V29 H26 V23 H14 Q9 23 9 18 Z")
    for x in range(17, 50, 6):
        ic.cut().circle(x, 13, 1.7)
    ic.rect(29, 28, 6, 18)
    ic.spath("M26 36 V43 L12.5 58 H51.5 L38 43 V36", 3.5, glass=True)
    ic.spath("M38 39.5 H47 V45", 3.5)


@icon("sonicator", "synthesis", tags=["ultrasound", "dispersion", "exfoliation"])
def _(ic):
    ic.spath("M7 26 V50 Q7 57 14 57 H50 Q57 57 57 50 V26", 3.5, glass=True)
    ic.sec().rect(11, 33, 42, 20, 3)
    for y in (39, 46):
        ic.cut().spath(f"M13 {y} q3 -2.5 6 0 t6 0 t6 0 t6 0 t6 0 t6 0 t6 0", 1.6)
    ic.cut().rect(23.5, 7.5, 17, 36, 4)
    ic.srect(27, 11, 10, 29, 2.5, 3.5, glass=True)
    ic.acc().spath("M45 8 q4 5 0 10", 2.5)
    ic.acc().spath("M50 5 q6 8 0 16", 2.5)


@icon("mortar_pestle", "synthesis", tags=["grinding", "mechanochemistry", "milling"])
def _(ic):
    ic.path("M7 30 H57 Q55.5 50 41 54 H23 Q8.5 50 7 30 Z")
    ic.rect(22, 55.5, 20, 4.5, 2)
    ic.cut().line(33, 33, 52, 8, 13)
    ic.sec().line(33, 33, 52, 8, 7.5)


@icon("gas_cylinder", "synthesis", tags=["gas", "N2", "Ar", "H2"])
def _(ic):
    ic.rect(24, 3, 16, 3.5, 1.75)
    ic.rect(29, 6, 6, 7)
    ic.rect(21, 14, 22, 46, 9)
    ic.cut().rect(21, 28, 22, 3)
    ic.cut().rect(21, 34, 22, 3)
    ic.acc().circle(32, 8, 1.6)


# ==========================================================================
# CHARACTERISATION
# ==========================================================================

@icon("xrd", "characterisation", label="XRD", tags=["diffraction", "PXRD", "crystallinity"])
def _(ic):
    cx, cy = 32, 46
    # goniometer arc
    p1 = polar(cx, cy, 25, 200)
    p2 = polar(cx, cy, 25, 340)
    ic.sec().spath(f"M{_f(p1[0])} {_f(p1[1])} A25 25 0 0 1 {_f(p2[0])} {_f(p2[1])}", 2.2)
    # tube + detector, rotated blocks
    for ang, role in ((210, "m"), (330, "m")):
        bx, by = polar(cx, cy, 25, ang)
        rot = ang - 180 if ang < 270 else ang
        ic.ops.append(_rot_rect(bx, by, 14, 10, 2.5, rot, role))
    # beams
    s = polar(cx, cy, 17, 210)
    d = polar(cx, cy, 17, 330)
    ic.acc().line(s[0], s[1], cx - 1, cy - 1.5, 3)
    ic.acc().line(cx + 1, cy - 1.5, d[0], d[1], 3)
    # sample stage
    ic.rect(22, 47, 20, 4.5, 1.5)
    ic.rect(29, 51.5, 6, 6)
    ic.rect(20, 57.5, 24, 3.5, 1.5)


def _rot_rect(cx, cy, w, h, r, deg, role):
    from .core import Op
    return Op(role, "rect", {"x": cx - w / 2, "y": cy - h / 2, "width": w, "height": h, "rx": r,
                             "transform": f"rotate({_f(deg)} {_f(cx)} {_f(cy)})"})


@icon("xrd_pattern", "characterisation", label="XRD Pattern", tags=["diffractogram", "PXRD"])
def _(ic):
    axes(ic)
    pts = lorentz_trace(12, 56, 50, [(17, 34, .45), (23, 16, .45), (29, 9, .45),
                                     (36, 22, .45), (44, 7, .45), (50, 11, .45)], n=220)
    ic.acc().spoly(pts, 2.4)


@icon("electron_microscope", "characterisation", label="SEM / TEM",
      tags=["SEM", "TEM", "STEM", "EDX", "morphology"])
def _(ic):
    ic.rect(23, 3, 18, 9, 3)
    ic.rect(26, 11, 12, 27)
    ic.rect(19, 16, 26, 6, 2.5)
    ic.rect(19, 27, 26, 6, 2.5)
    ic.cut().line(32, 7, 32, 36, 2.4)
    ic.rect(8, 38, 48, 20, 3.5)
    ic.cut().rect(17, 43, 30, 10, 2)
    ic.acc().poly([(32, 38), (28, 48), (36, 48)])
    ic.sec().rect(22, 48.5, 20, 3, 1)


@icon("gas_sorption", "characterisation", label="Gas Sorption (BET)",
      tags=["BET", "N2 isotherm", "porosity", "surface area"])
def _(ic):
    axes(ic)
    ic.acc().spath("M12 51 C13.5 22 20 17 56 14", 3)
    ic.sec().spath("M30 21 C40 21.5 48 22 56 22", 2.4)
    for x, y in ((14, 33), (19, 22), (28, 17.5), (40, 15.5), (51, 14.5)):
        ic.circle(x, y, 2.3)


@icon("tga", "characterisation", label="TGA", tags=["thermogravimetry", "stability", "mass loss"])
def _(ic):
    axes(ic)
    ic.acc().spath("M12 13 H22 C29 13 27 27 34 27 H41 C47 27 46 43 56 43", 3)
    ic.sec().spath("M12 49 L56 44", 2.2)


@icon("ftir", "characterisation", label="FTIR", tags=["IR", "infrared", "ATR", "functional groups"])
def _(ic):
    axes(ic)
    pts = lorentz_trace(12, 56, 13, [(18, 9, 1.3), (27, 22, .9), (31, 18, .9), (40, 12, 1.4),
                                     (48, 26, .7), (53, 8, 1.0)], n=180, sign=+1)
    ic.acc().spoly(pts, 2.6)


@icon("raman", "characterisation", label="Raman", tags=["raman", "vibrational", "defects"])
def _(ic):
    axes(ic)
    pts = lorentz_trace(12, 56, 50, [(24, 22, 1.6), (40, 34, 1.4), (50, 8, 2.2)], n=180)
    ic.acc().spoly(pts, 2.8)


@icon("uv_vis", "characterisation", label="UV-Vis", tags=["absorbance", "spectrophotometry", "band gap"])
def _(ic):
    ic.acc().line(3, 32, 17, 32, 5)
    ic.acc().poly([(17, 26), (23, 32), (17, 38)])
    ic.acc().line(45, 32, 56, 32, 2.6)
    ic.acc().poly([(56, 28), (61, 32), (56, 36)])
    ic.cut().rect(23, 7, 18, 51, 3)
    ic.srect(25.5, 9.5, 13, 46, 2, 3, glass=True)
    ic.sec().rect(29, 24, 6, 28.5, 1)


@icon("nmr", "characterisation", label="NMR", tags=["NMR", "solid-state", "digestion", "magnet"])
def _(ic):
    ic.path("M13 16 A19 6.5 0 0 1 51 16 V52 Q51 56 47 56 H17 Q13 56 13 52 Z")
    ic.cut().ellipse(32, 16, 19, 6.5)
    ic.path("M15.5 16.5 A16.5 4.6 0 0 1 48.5 16.5 A16.5 4.6 0 0 1 15.5 16.5 Z")
    ic.cut().ellipse(32, 16.3, 5.5, 2.1)
    ic.cut().line(13, 36, 51, 36, 2)
    ic.sec().line(32, 2, 32, 16, 3)
    ic.rect(17, 56, 6, 4, 1)
    ic.rect(41, 56, 6, 4, 1)


@icon("xps", "characterisation", label="XPS", tags=["photoelectron", "surface", "oxidation state"])
def _(ic):
    cx, cy = 37, 35
    ic.path(f"M{cx-21} {cy} A21 21 0 0 1 {cx+21} {cy} H{cx+12} A12 12 0 0 0 {cx-12} {cy} Z")
    ic.rect(cx + 11, cy + 1, 11, 6, 1.5)
    ic.rect(4, 51, 26, 5, 1.5)
    ic.acc().spath("M4 15 l3 4 l-4 3 l3 4 l-4 3 l3 4 l-4 3 l3 4 l1.5 2.5", 2.4)
    ic.acc().spath("M19 50 L24 37", 2.4)
    ic.acc().circle(24.4, 36, 2.1)


@icon("mass_spec", "characterisation", label="Mass Spectrum", tags=["MS", "LC-MS", "MALDI", "proteomics", "peptides"])
def _(ic):
    axes(ic)
    for x, h in ((15, 12), (20, 30), (24, 8), (29, 40), (33, 16), (38, 22), (43, 6), (47, 27), (52, 10)):
        ic.acc().line(x, 55, x, 55 - h, 2.8)


@icon("hplc", "characterisation", label="HPLC / LC", tags=["chromatography", "LC", "separation"])
def _(ic):
    axes(ic)
    pts = lorentz_trace(12, 56, 50, [(19, 18, 1.2), (28, 36, 1.4), (34, 12, 1.1), (45, 26, 1.8)], n=200)
    ic.acc().spoly(pts, 2.6)


# ==========================================================================
# MATERIALS & STRUCTURES
# ==========================================================================

@icon("mof_framework", "materials", label="MOF Framework", tags=["MOF", "framework", "SBU", "lattice", "UiO", "MOF-801"])
def _(ic):
    f = [(13, 25), (39, 25), (39, 51), (13, 51)]
    b = [(x + 12, y - 12) for x, y in f]
    edges = [(f[i], f[(i + 1) % 4]) for i in range(4)] + \
            [(b[i], b[(i + 1) % 4]) for i in range(4)] + [(f[i], b[i]) for i in range(4)]
    for (x1, y1), (x2, y2) in edges:
        ic.sec().line(x1, y1, x2, y2, 3)
    for x, y in f + b:
        ic.circle(x, y, 5)
        ic.cut().circle(x, y, 1.8)


@icon("octahedral_sbu", "materials", label="Zr6/Hf6 Cluster", tags=["SBU", "Zr6", "Hf6", "cluster", "octahedron"])
def _(ic):
    c = (32, 32)
    top, bot = (32, 8), (32, 56)
    eq = [(10, 34), (26, 40), (54, 30), (38, 24)]
    # back faces lighter, front faces solid
    ic.sec().poly([top, eq[3], eq[0]])
    ic.sec().poly([top, eq[2], eq[3]])
    ic.poly([top, eq[0], eq[1]])
    ic.poly([top, eq[1], eq[2]])
    ic.poly([bot, eq[0], eq[1]])
    ic.poly([bot, eq[1], eq[2]])
    ic.cut().spoly([eq[0], eq[1], eq[2]], 1.6)
    ic.cut().line(top[0], top[1], eq[1][0], eq[1][1], 1.6)
    ic.cut().line(bot[0], bot[1], eq[1][0], eq[1][1], 1.6)
    for p in (top, bot, *eq):
        ic.acc().circle(*p, 3.2)


@icon("nanoparticle", "materials", tags=["nanoparticle", "core-shell", "ligand", "colloid"])
def _(ic):
    for k in range(12):
        a = k * 30 + 15
        x1, y1 = polar(32, 32, 21, a)
        x2, y2 = polar(32, 32, 27.5, a)
        ic.acc().line(x1, y1, x2, y2, 2.2)
        ic.acc().circle(x2, y2, 2.2)
    ic.sec().circle(32, 32, 20)
    ic.cut().circle(32, 32, 15.5)
    ic.circle(32, 32, 13)
    ic.cut().spath("M24.5 30 A8 8 0 0 1 30 24.5", 2.6)


@icon("linker_bdc", "materials", label="Organic Linker", tags=["linker", "ligand", "BDC", "carboxylate", "fumarate"])
def _(ic):
    cx, cy, r = 32, 32, 11
    hexpts = [polar(cx, cy, r, 90 + 60 * k) for k in range(6)]
    ic.spoly(hexpts + [hexpts[0]], 3.2)
    ic.scircle(cx, cy, 5.2, 2.4)
    for sgn in (-1, 1):
        y0 = cy + sgn * r
        y1 = cy + sgn * 19
        ic.line(cx, y0, cx, y1, 3.2)
        for dx in (-7, 7):
            ox, oy = cx + dx, y1 + sgn * 6
            ic.line(cx, y1, ox, oy, 3.2)
            ic.sec().circle(ox, oy, 3.6)


@icon("porous_crystal", "materials", tags=["crystal", "pores", "particle", "octahedral crystal"])
def _(ic):
    top, bot = (32, 5), (32, 59)
    l, r, f, b = (7, 30), (57, 30), (28, 37), (38, 25)
    ic.sec().poly([top, l, f])
    ic.poly([top, f, r])
    ic.sec().poly([bot, l, f])
    ic.poly([bot, f, r])
    for (x, y) in ((40, 25), (35, 31), (45, 30), (40, 37), (37, 44), (45, 42)):
        ic.cut().circle(x, y, 2)
    ic.cut().line(*f, *top, 1.4)
    ic.cut().line(*f, *bot, 1.4)


# ==========================================================================
# PROTEOMICS / BIOCATALYSIS
# ==========================================================================

@icon("peptide_cleavage", "proteomics", label="Peptide Cleavage",
      tags=["protease", "hydrolysis", "peptide bond", "nanozyme", "artificial protease"])
def _(ic):
    def y(x):
        return 34 + 8 * math.sin((x - 6) / 8)
    xs_l = [7, 14, 21]
    xs_r = [43, 50, 57]
    ic.spath(smooth_path([(x, y(x)) for x in range(5, 25, 2)]), 3)
    ic.sec().spath(smooth_path([(x, y(x)) for x in range(41, 60, 2)]), 3)
    for x in xs_l:
        ic.circle(x, y(x), 4.6)
    for x in xs_r:
        ic.sec().circle(x, y(x), 4.6)
    ic.acc().poly([(34, 14), (27, 32), (32.5, 32), (29, 50), (38, 29), (32.5, 29), (36, 14)])


@icon("protein", "proteomics", tags=["protein", "enzyme", "substrate", "BSA"])
def _(ic):
    pts = [(14, 46), (10, 34), (16, 24), (26, 28), (24, 40), (34, 44), (40, 34), (34, 24),
           (40, 14), (50, 18), (54, 30), (48, 40), (52, 50)]
    ic.spath(smooth_path(pts), 4.2)
    ic.sec().circle(14, 46, 3.6)
    ic.acc().circle(52, 50, 3.6)


@icon("digestion_tube", "proteomics", label="Digestion", tags=["digestion", "incubation", "eppendorf", "37 C"])
def _(ic):
    ic.rect(21, 4, 22, 7, 2)
    ic.rect(41, 6, 9, 4, 2)
    ic.spath("M23.5 14 V28 L29.5 56 Q32 59 34.5 56 L40.5 28 V14", 3.5, glass=True)
    ic.sec().path("M27.5 32 H36.5 L32.6 50.5 Q32 52 31.4 50.5 Z")
    ic.acc().spath("M51 22 v20", 2.5)
    ic.acc().circle(51, 46, 4)
