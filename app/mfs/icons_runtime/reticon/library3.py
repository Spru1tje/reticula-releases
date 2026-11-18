"""
Original redraws ("paraphrases") of common lab subjects, so the library does
not depend on any icon set that needs attribution. Geometry is drawn from
scratch on the Reticon grid; nothing is traced from another set.
"""
from __future__ import annotations

import math

from .core import Op, polar, smooth_path, _f
from .library import icon, axes
from .library2 import rot_rect, gauss_trace


def ell(ic, cx, cy, rx, ry, deg, w, role="s"):
    ic.ops.append(Op(role, "ellipse", {"cx": cx, "cy": cy, "rx": rx, "ry": ry,
                                       "transform": f"rotate({_f(deg)} {_f(cx)} {_f(cy)})"}, w))


# ==========================================================================
# GLASSWARE & CONSUMABLES
# ==========================================================================

@icon("measuring_cylinder", "glassware", tags=["graduated cylinder", "volume"])
def _(ic):
    ic.spath("M24 6 Q24 4 22 4 M24 6 V54 H40 V6", 3, glass=True)
    ic.sec().rect(26.5, 26, 11, 25.5, 1)
    for y in range(12, 52, 6):
        ic.line(34, y, 38, y, 1.4)
    ic.rect(16, 54, 32, 5, 2)


@icon("volumetric_flask", "glassware", tags=["volumetric", "standard solution", "dilution"])
def _(ic):
    ic.rect(27, 1, 10, 5, 1.5)
    ic.spath("M29 6 V29.3 A14 14 0 1 0 35 29.3 V6", 3, glass=True)
    ic.line(27.5, 15, 36.5, 15, 1.6)
    ic.sec().path("M21.55 44 A10.5 10.5 0 0 0 42.45 44 Z")


@icon("separatory_funnel", "glassware", tags=["extraction", "phase separation", "work-up"])
def _(ic):
    ic.rect(28, 1, 8, 5, 1.5)
    ic.spath("M28 6 V10 C14 18 14 34 29 45 V50 H35 V45 C50 34 50 18 36 10 V6", 3, glass=True)
    ic.acc().path("M19.6 24 H44.4 C44.7 26 44.6 28 44.1 30 H19.9 C19.4 28 19.3 26 19.6 24 Z")
    ic.sec().path("M20 30.6 H44 C42.8 35.5 38.5 40 32 44 C25.5 40 21.2 35.5 20 30.6 Z")
    ic.rect(24, 50, 16, 4, 2)
    ic.line(32, 54, 32, 61, 2.4)


@icon("liebig_condenser", "glassware", label="Condenser", tags=["condenser", "distillation", "reflux", "cooling"])
def _(ic):
    ic.srect(23, 8, 18, 46, 5, 3, glass=True)
    ic.sec().rect(26, 11, 12, 40, 3)
    ic.srect(29.5, 2, 5, 60, 2.5, 2.4, glass=True)
    ic.line(41, 14, 50, 14, 3)
    ic.line(23, 48, 14, 48, 3)


@icon("burette", "glassware", tags=["titration", "burette", "stand"])
def _(ic):
    ic.rect(8, 56, 34, 5, 2)
    ic.rect(12, 4, 4, 52, 1.5)
    ic.rect(16, 18, 14, 3.5, 1.5)
    ic.srect(28, 5, 7, 40, 2.5, 2.6, glass=True)
    ic.sec().rect(30, 14, 3, 30, 1)
    ic.rect(25, 44, 13, 3.5, 1.75)
    ic.spath("M30 47.5 L31 56 H32 L33 47.5", 2, glass=True)


@icon("cuvette", "glassware", tags=["cuvette", "UV-Vis", "quartz", "spectroscopy"])
def _(ic):
    ic.srect(23.5, 5.5, 17, 53, 2, 3, glass=True)
    ic.sec().rect(27, 22, 10, 33, 1)
    ic.rect(23.5, 5.5, 4, 53, 1)


@icon("falcon_tube", "glassware", label="Centrifuge Tube", tags=["falcon", "50 mL", "15 mL", "conical tube"])
def _(ic):
    ic.rect(22, 3, 20, 9, 2.5)
    ic.spath("M24.5 12 V49 L30 59 H34 L39.5 49 V12", 3, glass=True)
    ic.sec().path("M28 32 H36 V48.5 L32.8 54.5 H31.2 L28 48.5 Z")
    for y in (18, 24, 30, 36, 42):
        ic.line(34, y, 37, y, 1.4)


@icon("microtube", "glassware", label="Microtube", tags=["eppendorf", "1.5 mL", "sample tube"])
def _(ic):
    ic.spath("M23 12 V28 L29.5 57 Q32 60 34.5 57 L41 28 V12", 3, glass=True)
    ic.rect(21, 8, 22, 5, 2)
    ic.spath("M42 10 C50 10 52 4 46 3 H36", 2.6)
    ic.sec().path("M26.5 33 H37.5 L33 53 Q32 55 31 53 Z")


@icon("reagent_bottle", "glassware", tags=["bottle", "amber", "reagent", "chemical"])
def _(ic):
    ic.rect(25, 3, 14, 9, 2)
    ic.sec().path("M26 12 H38 V17 Q50 20 50 31 V54 Q50 59 45 59 H19 Q14 59 14 54 V31 Q14 20 26 17 Z")
    ic.cut().rect(20, 31, 24, 15, 2)


@icon("wash_bottle", "glassware", tags=["wash bottle", "squeeze bottle", "solvent"])
def _(ic):
    ic.spath("M19 27 Q19 21 25 21 H39 Q45 21 45 27 V54 Q45 59 40 59 H24 Q19 59 19 54 Z", 3, glass=True)
    ic.sec().rect(22.5, 36, 19, 19.5, 2)
    ic.rect(26, 14, 12, 7, 1.5)
    ic.spath("M32 14 V9 Q32 5 37 5 H52", 2.4)


@icon("pasteur_pipette", "glassware", tags=["pipette", "dropper", "transfer"])
def _(ic):
    ic.acc().ellipse(32, 10, 6, 8)
    ic.spath("M29 17 V40 L31.3 58 H32.7 L35 40 V17", 2.6, glass=True)
    ic.sec().path("M30.3 34 H33.7 V40 L32.4 50 H31.6 L30.3 40 Z")
    ic.sec().path("M32 59 C34 62 34 63.5 32 63.5 C30 63.5 30 62 32 59 Z")


@icon("petri_dish", "glassware", tags=["petri dish", "agar", "culture", "crystallisation dish"])
def _(ic):
    ic.scircle(32, 32, 25, 3, glass=True)
    ic.sec().circle(32, 32, 20.5)
    for x, y, r in ((24, 26, 3), (38, 24, 2.2), (36, 38, 3.4), (25, 39, 2), (44, 32, 1.8)):
        ic.acc().circle(x, y, r)


@icon("glass_slide", "glassware", label="Microscope Slide", tags=["slide", "coverslip", "sample"])
def _(ic):
    ic.srect(5, 21, 54, 22, 2, 2.4, glass=True)
    ic.rect(7.5, 23.5, 12, 17, 1)
    ic.sec().circle(41, 32, 5.5)
    ic.srect(32, 24, 18, 16, 1, 1.4, glass=True)


@icon("tube_rack", "glassware", tags=["rack", "test tubes", "samples"])
def _(ic):
    hts = (30, 24, 36, 28)
    for x, h in zip((17, 27, 37, 47), hts):
        ic.spath(f"M{x - 3.5} 8 V45 A3.5 3.5 0 0 0 {x + 3.5} 45 V8", 2.4, glass=True)
        ic.sec().rect(x - 2, h, 4, 47 - h, 1.6)
    ic.rect(6, 30, 52, 6, 2)
    ic.rect(6, 50, 52, 6, 2)
    ic.rect(8, 30, 4, 26, 1)
    ic.rect(52, 30, 4, 26, 1)


@icon("well_plate_96", "glassware", label="96-Well Plate", tags=["microplate", "96-well", "assay", "screening"])
def _(ic):
    ic.rect(4, 13, 56, 38, 4)
    for r in range(5):
        for c in range(8):
            x, y = 11 + c * 6.3, 20 + r * 6
            ic.cut().circle(x, y, 2.3)
    for r, c in ((0, 0), (0, 1), (1, 0), (1, 1), (2, 3), (2, 4), (3, 6), (4, 6), (4, 7), (3, 7)):
        ic.sec().circle(11 + c * 6.3, 20 + r * 6, 1.8)


@icon("well_plate_6", "glassware", label="6-Well Plate", tags=["well plate", "cell culture", "6-well"])
def _(ic):
    ic.rect(5, 14, 54, 36, 5)
    for x in (18, 32, 46):
        for y in (25, 39):
            ic.cut().circle(x, y, 6)
            ic.sec().circle(x, y, 4.5)


@icon("pipette_tip_box", "glassware", label="Pipette Tips", tags=["tips", "consumables", "box"])
def _(ic):
    for x in range(12, 56, 7):
        ic.sec().poly([(x - 2.6, 24), (x + 2.6, 24), (x + 1.8, 15), (x - 1.8, 15)])
    ic.rect(6, 23, 52, 6, 2)
    ic.rect(9, 29, 46, 28, 3)
    ic.cut().rect(17, 37, 30, 10, 2)


@icon("tweezers", "glassware", tags=["forceps", "tweezers", "handling"])
def _(ic):
    ic.spath("M24 6 L31 56", 3.4)
    ic.spath("M40 6 L33 56", 3.4)
    ic.rect(22, 3, 20, 5, 2.5)


@icon("spatula", "glassware", tags=["spatula", "weighing", "powder"])
def _(ic):
    ic.line(14, 50, 46, 18, 3.2)
    rot_rect(ic, 50, 14, 14, 8, 4, -45)
    ic.ops.append(Op("m", "ellipse", {"cx": 11, "cy": 53, "rx": 6, "ry": 4, "transform": "rotate(-45 11 53)"}))
    ic.sec().circle(10, 54, 2)


@icon("multichannel_pipette", "glassware", tags=["multichannel", "pipette", "plates"])
def _(ic):
    ic.rect(27, 2, 10, 6, 2)
    ic.rect(30.5, 7, 3, 4)
    ic.rect(24, 11, 16, 20, 3)
    ic.cut().rect(28, 15, 8, 7, 1.2)
    ic.rect(8, 30, 48, 8, 2.5)
    for x in range(13, 54, 6):
        ic.sec().poly([(x - 1.8, 38), (x + 1.8, 38), (x + 0.5, 52), (x - 0.5, 52)])


# ==========================================================================
# INSTRUMENTS & BENCH EQUIPMENT
# ==========================================================================

@icon("ph_meter", "instruments", label="pH Meter", tags=["pH", "electrode", "acidity"])
def _(ic):
    ic.spath("M5 34 V54 Q5 58 9 58 H24 Q28 58 28 54 V34", 2.8, glass=True)
    ic.sec().rect(7.5, 42, 18, 14, 1.5)
    ic.spath("M16.5 8 C16.5 1 30 2 33 22", 2)
    ic.rect(14, 8, 5, 42, 2.5)
    ic.rect(31, 20, 28, 36, 4)
    ic.cut().rect(35, 25, 20, 10, 1.5)
    for x in (38, 45, 52):
        ic.cut().circle(x, 45, 2.4)


@icon("stopwatch", "instruments", tags=["timer", "time", "reaction time"])
def _(ic):
    ic.rect(28, 4, 8, 6, 1.5)
    rot_rect(ic, 48, 16, 6, 5, 1.5, 45)
    ic.circle(32, 36, 22)
    ic.cut().circle(32, 36, 17)
    for k in range(12):
        x1, y1 = polar(32, 36, 14.5, k * 30)
        x2, y2 = polar(32, 36, 16, k * 30)
        ic.line(x1, y1, x2, y2, 1.4)
    ic.acc().line(32, 36, 32, 23, 2.6)
    ic.acc().line(32, 36, 40, 40, 2.6)
    ic.circle(32, 36, 2.4)


@icon("hourglass", "instruments", label="Hourglass Timer", tags=["time", "duration", "aging"])
def _(ic):
    ic.spath("M20 11 H44 C44 22 34.5 28 34.5 32 C34.5 36 44 42 44 53 H20 C20 42 29.5 36 29.5 32 "
             "C29.5 28 20 22 20 11 Z", 2.8, glass=True)
    ic.sec().path("M24.5 17 H39.5 C38.5 22 34.5 25.5 32 28 C29.5 25.5 25.5 22 24.5 17 Z")
    ic.sec().path("M23 50.5 C23 46 28 43 32 42 C36 43 41 46 41 50.5 Z")
    ic.sec().line(32, 28, 32, 42, 1.2)
    ic.rect(15, 5, 34, 6, 2.5)
    ic.rect(15, 53, 34, 6, 2.5)


@icon("thermometer", "instruments", tags=["temperature", "thermometer", "heating"])
def _(ic):
    ic.spath("M27 9 A5 5 0 0 1 37 9 V43.5 A9 9 0 1 1 27 43.5 Z", 3, glass=True)
    ic.sec().rect(30.5, 18, 3, 30, 1.5)
    ic.sec().circle(32, 51, 6)
    for y in (14, 20, 26, 32, 38):
        ic.line(40, y, 44, y, 1.6)


@icon("bar_magnet", "instruments", label="Magnet", tags=["magnet", "magnetic separation", "Fe3O4"])
def _(ic):
    ic.path("M12 8 H25 V34 A7 7 0 0 0 39 34 V8 H52 V34 A20 20 0 0 1 12 34 Z")
    ic.acc().rect(12, 8, 13, 9, 1.5)
    ic.sec().rect(39, 8, 13, 9, 1.5)


@icon("water_bath", "instruments", tags=["water bath", "incubation", "heating"])
def _(ic):
    for x in (20, 32, 44):
        ic.spath(f"M{x - 3.5} 6 V26 A3.5 3.5 0 0 0 {x + 3.5} 26 V6", 2.4, glass=True)
    ic.sec().rect(8, 20, 48, 11, 2)
    ic.rect(6, 30, 52, 26, 4)
    ic.cut().rect(12, 38, 16, 9, 1.5)
    ic.cut().circle(46, 42.5, 4)


@icon("incubator_co2", "instruments", label="Incubator", tags=["CO2 incubator", "cell culture", "37 C"])
def _(ic):
    ic.rect(9, 5, 46, 52, 4)
    ic.cut().rect(13, 9, 31, 44, 2)
    for y in (22, 34, 46):
        ic.line(15, y, 42, y, 2)
        ic.sec().rect(18, y - 6, 8, 5, 1)
        ic.sec().rect(30, y - 5, 8, 4, 1)
    ic.cut().rect(47, 11, 5, 8, 1)
    ic.rect(47, 26, 4, 14, 2)
    ic.rect(13, 57, 6, 4, 1)
    ic.rect(45, 57, 6, 4, 1)


@icon("freezer", "instruments", label="-80 °C Freezer", tags=["freezer", "-80", "storage", "cold"])
def _(ic):
    ic.rect(12, 3, 40, 56, 4)
    ic.cut().rect(17, 8, 18, 8, 1.5)
    ic.rect(44, 20, 4, 18, 2)
    for k in range(3):
        a = 90 + k * 60
        x1, y1 = polar(29, 40, 10, a)
        x2, y2 = polar(29, 40, 10, a + 180)
        ic.acc().line(x1, y1, x2, y2, 2.6)
    ic.rect(15, 59, 6, 3, 1)
    ic.rect(43, 59, 6, 3, 1)


@icon("fume_hood", "instruments", tags=["fume hood", "ventilation", "safety"])
def _(ic):
    ic.rect(5, 3, 54, 58, 3)
    ic.cut().rect(10, 10, 44, 32, 2)
    ic.spath("M28 26 V29 L22 39 H38 L32 29 V26", 2, glass=True)
    ic.sec().path("M25.2 34 H34.8 L36.6 37.2 H23.4 Z")
    ic.srect(10, 10, 44, 13, 1, 2, glass=True)
    ic.cut().rect(10, 46, 20, 11, 1.5)
    ic.cut().rect(34, 46, 20, 11, 1.5)
    ic.acc().spath("M20 6 h24", 1.6)


@icon("pcr_thermocycler", "instruments", label="Thermocycler", tags=["PCR", "thermal cycler", "heating block"])
def _(ic):
    ic.path("M8 27 L14 9 H50 L56 27 Z")
    ic.cut().path("M18 22 L21 14 H43 L46 22 Z")
    ic.rect(6, 28, 52, 28, 4)
    ic.cut().rect(12, 35, 18, 11, 1.5)
    ic.sec().rect(14.5, 37.5, 13, 6, 1)
    for x in (40, 48):
        ic.cut().circle(x, 40.5, 3)


@icon("plate_reader", "instruments", tags=["plate reader", "absorbance", "fluorescence", "ELISA"])
def _(ic):
    ic.rect(6, 14, 52, 40, 4)
    ic.cut().rect(40, 20, 13, 9, 1.5)
    for y in (36, 41):
        ic.cut().line(42, y, 52, y, 1.6)
    ic.cut().rect(10, 42, 28, 8, 1.5)
    ic.sec().rect(12, 36, 24, 5, 1)
    ic.rect(10, 41, 28, 4, 1)


@icon("spectrophotometer", "instruments", tags=["spectrophotometer", "UV-Vis", "absorbance"])
def _(ic):
    ic.srect(20, 6, 7, 14, 1, 2, glass=True)
    ic.sec().rect(21.5, 12, 4, 7, 0.8)
    ic.path("M6 30 L10 18 H54 L58 30 V54 Q58 58 54 58 H10 Q6 58 6 54 Z")
    ic.cut().rect(14, 21, 20, 5, 1.5)
    ic.cut().rect(36, 35, 16, 10, 1.5)
    for x in (13, 20, 27):
        ic.cut().rect(x - 2.5, 38, 5, 4, 1)
        ic.cut().rect(x - 2.5, 46, 5, 4, 1)


@icon("gel_box", "instruments", label="Gel Electrophoresis", tags=["agarose", "electrophoresis", "DNA", "gel tank"])
def _(ic):
    ic.srect(6, 22, 52, 28, 3, 2.6, glass=True)
    ic.sec().rect(9, 30, 46, 17, 2)
    ic.rect(15, 31.5, 34, 13, 1)
    for x in (19, 26, 33, 40):
        ic.cut().rect(x, 33, 4.5, 2, 0.6)
    ic.acc().spath("M10 22 V12 H18", 2.6)
    ic.sec().spath("M54 22 V12 H46", 2.6)


@icon("blot_membrane", "instruments", label="Western Blot", tags=["western blot", "membrane", "bands", "immunoblot"])
def _(ic):
    ic.rect(12, 6, 40, 52, 2.5)
    for k, y in enumerate((14, 20, 27, 35, 44, 52)):
        ic.sec().rect(15, y, 7, 2, 1)
    for x, ys in ((27, (22, 36)), (36, (22, 36, 46)), (45, (22, 30))):
        for y in ys:
            ic.acc().rect(x - 3.5, y, 8, 3, 1.5)


@icon("light_microscope", "instruments", label="Light Microscope", tags=["microscope", "optical", "brightfield"])
def _(ic):
    ic.rect(12, 54, 40, 6, 3)
    ic.path("M40 54 V40 C40 30 35 24 27 20 L30 14 C41 19 48 28 48 40 V54 Z")
    rot_rect(ic, 25, 19, 9, 24, 2.5, -28)
    rot_rect(ic, 19.5, 7.5, 7, 7, 1.5, -28)
    ic.rect(28, 29, 5, 8, 1.2)
    ic.rect(18, 38, 26, 4, 1.5)
    ic.sec().rect(24, 35.5, 12, 2.5, 0.8)
    ic.acc().circle(30.5, 47, 3)


@icon("lc_ms", "instruments", label="LC-MS System", tags=["LC-MS", "proteomics", "Orbitrap", "mass spectrometer"])
def _(ic):
    ic.rect(4, 8, 22, 50, 3)
    for y in (12, 24, 36):
        ic.cut().rect(7, y, 16, 9, 1.5)
    ic.sec().rect(9, 47, 12, 6, 1)
    ic.acc().spath("M26 44 H32", 2.6)
    ic.rect(32, 22, 28, 36, 3)
    ic.cut().rect(36, 27, 13, 9, 1.5)
    ic.cut().line(36, 44, 56, 44, 1.6)
    ic.cut().line(36, 50, 56, 50, 1.6)


@icon("bunsen_burner", "instruments", tags=["burner", "flame", "sterilisation"])
def _(ic):
    ic.acc().path("M32 3 C41 12 41 21 37 26 H27 C23 21 23 12 32 3 Z")
    ic.cut().path("M32 13 C36 18 36 22 34.5 25 H29.5 C28 22 28 18 32 13 Z")
    ic.rect(28, 26, 8, 27, 2)
    ic.rect(26, 40, 12, 4, 1.5)
    ic.line(36, 47, 52, 47, 3)
    ic.rect(16, 52, 32, 7, 3)


@icon("magnifier", "instruments", label="Magnifier", tags=["inspect", "zoom", "magnifying glass", "detail"])
def _(ic):
    rot_rect(ic, 46, 46, 9, 24, 4, -45)
    ic.scircle(27, 27, 17, 4, glass=True)
    ic.cut().spath("M17 24 A11 11 0 0 1 24 16", 2)


@icon("ruler", "instruments", tags=["ruler", "scale bar", "length", "measure"])
def _(ic):
    rot_rect(ic, 32, 32, 58, 16, 2, -35)
    ux, uy = math.cos(math.radians(-35)), math.sin(math.radians(-35))
    nx, ny = -uy, ux
    for k in range(-5, 6):
        cx, cy = 32 + ux * k * 4.6, 32 + uy * k * 4.6
        ln = 6 if k % 5 == 0 else 3.5
        x1, y1 = cx - nx * 8, cy - ny * 8
        x2, y2 = cx - nx * (8 - ln), cy - ny * (8 - ln)
        ic.ops.append(Op("x", "line", {"x1": x1, "y1": y1, "x2": x2, "y2": y2}, 1.4))


@icon("lab_coat", "safety", tags=["PPE", "lab coat", "safety"])
def _(ic):
    ic.path("M22 5 L32 14 L42 5 L54 12 L59 44 L51 46 L48 30 V59 H16 V30 L13 46 L5 44 L10 12 Z")
    ic.cut().line(32, 16, 32, 57, 1.4)
    ic.cut().spath("M22 5 L29 20 L32 14 L35 20 L42 5", 1.4)
    ic.cut().rect(36, 38, 8, 7, 1)


@icon("safety_goggles", "safety", tags=["goggles", "PPE", "eye protection"])
def _(ic):
    ic.acc().line(2, 30, 8, 30, 3)
    ic.acc().line(56, 30, 62, 30, 3)
    ic.path("M6 26 Q6 19 13 19 H51 Q58 19 58 26 V35 Q58 44 50 44 H41 Q37 44 35 40 Q32 36 29 40 "
            "Q27 44 23 44 H14 Q6 44 6 35 Z")
    ic.spath("M11 24 H28.5 V33 Q28.5 39 23 39 H15 Q11 39 11 34 Z", 2, glass=True)
    ic.spath("M53 24 H35.5 V33 Q35.5 39 41 39 H49 Q53 39 53 34 Z", 2, glass=True)


@icon("nitrile_gloves", "safety", label="Gloves", tags=["gloves", "PPE", "nitrile"])
def _(ic):
    for x, y in ((18, 12), (25.5, 7), (33, 9), (40.5, 14)):
        ic.sec().rect(x, y, 6.5, 26, 3.25)
    rot_rect(ic, 15, 37, 7, 20, 3.5, -40, "s")
    ic.sec().rect(17, 26, 30.5, 22, 9)
    ic.sec().rect(20, 44, 25, 15, 2)
    ic.cut().line(20, 47.5, 45, 47.5, 1.2)


@icon("lab_notebook", "safety", label="Lab Notebook", tags=["notebook", "ELN", "protocol", "records"])
def _(ic):
    ic.sec().rect(14, 5, 38, 54, 3)
    ic.cut().rect(21, 12, 24, 11, 1.5)
    for y in range(10, 58, 8):
        ic.rect(9, y, 9, 3, 1.5)
    rot_rect(ic, 52, 42, 6, 30, 2, 28, "a")


# ==========================================================================
# SYMBOLS
# ==========================================================================

@icon("atom", "symbols", tags=["atom", "element", "physics", "nucleus"])
def _(ic):
    for a in (0, 60, 120):
        ell(ic, 32, 32, 26, 9.5, a, 2.6, "s")
    ic.acc().circle(32, 32, 5.5)
    for a in (0, 60, 120):
        ic.circle(*polar(32, 32, 26, a + 180 * (a == 60)), 3)


@icon("dna_helix", "symbols", label="DNA", tags=["DNA", "genetics", "helix", "nucleic acid"])
def _(ic):
    a, b = [], []
    for i in range(0, 53):
        y = 6 + i
        t = (y - 6) / 52 * 2 * math.pi * 1.25
        a.append((32 + 14 * math.sin(t), y))
        b.append((32 - 14 * math.sin(t), y))
    for i in range(3, 52, 5):
        ic.line(a[i][0], a[i][1], b[i][0], b[i][1], 2)
    ic.sec().spath(smooth_path(a[::2]), 3.4)
    ic.acc().spath(smooth_path(b[::2]), 3.4)


@icon("molecule", "symbols", tags=["molecule", "ball and stick", "chemistry", "compound"])
def _(ic):
    pts = [(32, 30), (14, 18), (50, 16), (46, 50), (18, 48)]
    for p in pts[1:]:
        ic.line(*pts[0], *p, 3.4)
    ic.circle(*pts[0], 8)
    ic.sec().circle(*pts[1], 6)
    ic.acc().circle(*pts[2], 5)
    ic.sec().circle(*pts[3], 6)
    ic.acc().circle(*pts[4], 5)


@icon("droplet", "symbols", tags=["liquid", "drop", "solvent", "water"])
def _(ic):
    ic.sec().path("M32 5 C42 19 49 28 49 39 A17 17 0 0 1 15 39 C15 28 22 19 32 5 Z")
    ic.cut().spath("M22 40 A10 10 0 0 0 28 50", 2.6)


@icon("flame", "symbols", tags=["fire", "heat", "combustion", "flammable"])
def _(ic):
    ic.acc().path("M32 3 C36 14 48 20 48 38 A16 16 0 0 1 16 38 C16 28 22 24 24 16 C28 22 30 24 30 28 "
                  "C34 22 34 12 32 3 Z")
    ic.cut().path("M32 30 C36 36 40 40 40 46 A8 8 0 0 1 24 46 C24 40 28 36 32 30 Z")


@icon("snowflake", "symbols", tags=["cold", "cryo", "freezing", "low temperature"])
def _(ic):
    for k in range(6):
        a = k * 60 - 90
        x, y = polar(32, 32, 25, a)
        ic.sec().line(32, 32, x, y, 3)
        for d in (-1, 1):
            bx, by = polar(32, 32, 15, a)
            ex, ey = polar(bx, by, 7, a + d * 45)
            ic.sec().line(bx, by, ex, ey, 2.6)
    ic.circle(32, 32, 3.5)


@icon("lightning", "symbols", label="Electricity", tags=["electricity", "voltage", "energy", "electrochemistry"])
def _(ic):
    ic.acc().poly([(36, 3), (14, 36), (29, 36), (24, 61), (50, 24), (35, 24), (42, 3)])


@icon("laser", "symbols", tags=["laser", "light source", "photon", "irradiation"])
def _(ic):
    ic.rect(3, 25, 23, 14, 3)
    ic.rect(26, 28, 4, 8, 1)
    ic.cut().circle(10, 32, 2.5)
    ic.acc().line(31, 32, 54, 32, 3.4)
    for a in (0, 45, 90, 135, 180, 225, 270, 315):
        x1, y1 = polar(56, 32, 4.5, a)
        x2, y2 = polar(56, 32, 7.5, a)
        ic.acc().line(x1, y1, x2, y2, 1.8)


@icon("periodic_table", "symbols", tags=["periodic table", "elements", "chemistry"])
def _(ic):
    s = 5.2
    def cell(c, r, role="m"):
        ic.rect(4 + c * s + .5, 10 + r * s + .5, s - 1, s - 1, 0.8, role)
    for r in range(7):
        cell(0, r, "s")
    for r in range(1, 7):
        cell(1, r, "s")
    for c in range(6, 11):
        for r in range(1, 7):
            cell(c, r, "a" if c in (9, 10) and r < 4 else "m")
    for c in range(2, 6):
        for r in range(3, 7):
            cell(c, r)
    cell(10, 0, "a")
    for c in range(2, 10):
        cell(c, 8, "s")


# ==========================================================================
# GENERIC CHARTS
# ==========================================================================

@icon("line_chart", "data", tags=["line chart", "trend", "time series"])
def _(ic):
    axes(ic)
    pts = [(14, 46), (22, 36), (30, 40), (38, 24), (46, 28), (54, 14)]
    ic.acc().spoly(pts, 3)
    for p in pts:
        ic.circle(*p, 2.6)


@icon("bar_chart", "data", tags=["bar chart", "comparison", "columns"])
def _(ic):
    axes(ic)
    for x, h, role in ((15, 18, "s"), (24, 32, "a"), (33, 24, "s"), (42, 40, "a"), (51, 30, "s")):
        ic.rect(x - 3, 55 - h, 6, h - 1.75, 1, role)


@icon("scatter_plot", "data", tags=["scatter", "correlation", "points"])
def _(ic):
    axes(ic)
    for x, y in ((15, 46), (19, 40), (24, 42), (27, 34), (31, 36), (36, 28), (39, 31), (44, 22), (49, 24), (53, 15)):
        ic.sec().circle(x, y, 2.6)
    for x, y in ((18, 30), (30, 20), (46, 40)):
        ic.acc().circle(x, y, 2.6)


@icon("histogram", "data", tags=["histogram", "distribution", "frequency"])
def _(ic):
    axes(ic)
    for i, h in enumerate((6, 14, 26, 38, 30, 18, 9, 4)):
        ic.sec().rect(12 + i * 5.6, 55 - h, 5.6, h - 1.75, 0)


@icon("heatmap", "data", tags=["heatmap", "matrix", "correlation"])
def _(ic):
    roles = "sasmxasxmsxmsaxsmaxsamxsa"
    for i, role in enumerate(roles):
        r, c = divmod(i, 5)
        ic.rect(7 + c * 10.4, 7 + r * 10.4, 9.6, 9.6, 1.2, role)


@icon("box_plot", "data", tags=["box plot", "statistics", "quartiles"])
def _(ic):
    axes(ic)
    for x, lo, q1, med, q3, hi, role in ((20, 44, 36, 30, 24, 14, "s"), (40, 50, 44, 40, 34, 26, "a")):
        ic.line(x, hi, x, lo, 2)
        ic.line(x - 4, hi, x + 4, hi, 2)
        ic.line(x - 4, lo, x + 4, lo, 2)
        ic.rect(x - 7, q3, 14, q1 - q3, 1.5, role)
        ic.line(x - 7, med, x + 7, med, 2)


@icon("violin_plot", "data", tags=["violin", "distribution", "statistics"])
def _(ic):
    axes(ic)
    for cx, role, w in ((22, "s", 1.0), (42, "a", 0.8)):
        right = [(cx + w * 10 * math.exp(-((y - 30) ** 2) / 90) + w * 3 * math.exp(-((y - 44) ** 2) / 20), y)
                 for y in range(12, 52, 2)]
        left = [(2 * cx - x, y) for x, y in right[::-1]]
        d = "M" + " L".join(f"{_f(x)} {_f(y)}" for x, y in right + left) + " Z"
        ic.path(d, role)
        ic.line(cx, 26, cx, 34, 2)


@icon("gaussian", "data", label="Normal Distribution", tags=["gaussian", "bell curve", "statistics", "fit"])
def _(ic):
    axes(ic)
    ic.acc().spoly(gauss_trace(12, 56, 52, [(34, 40, 7)], n=120), 3)
    ic.sec().line(34, 12, 34, 52, 1.6)


@icon("exponential_decay", "data", tags=["decay", "kinetics", "first-order", "half-life"])
def _(ic):
    axes(ic)
    pts = [(12 + i, 12 + 40 * (1 - math.exp(-i / 9))) for i in range(0, 45)]
    ic.acc().spoly(pts, 3)
    for i in (2, 8, 15, 24, 36):
        ic.circle(12 + i, 12 + 40 * (1 - math.exp(-i / 9)), 2.4)
