"""Reticula 2.20 icon set: professional scientific vector icons (CC0), parametric families, prefix r22_.

Style (the same for every family): one outline width (2.2 on the 64 grid), slate outlines, one steel-blue accent,
a light fill only for liquids and samples, no faces, no toys, no decoration without meaning. Every icon is pure
vector (no masks, images or text inside the artwork), so the app adds labels and every export stays vector.
"""
import math
import random

from reticon.library import icon, REGISTRY, axes

_TAKEN = set(REGISTRY)          # built-in names, registered before this module loads
_NEW = []
SW = 2.2                        # outline width
AW = 1.4                        # accent detail width


def _reg(name, cat, label, tags, fn):
    if name in _TAKEN or name in _NEW:
        raise ValueError("icon name already in use: " + name)
    _NEW.append(name)
    icon(name, cat, label=label, tags=tags)(fn)


def _pd(pts, close=True):
    d = "M" + " L".join("%.2f %.2f" % p for p in pts)
    return d + (" Z" if close else "")


def _scale(pts, s, c=32.0):
    return [(c + (x - c) * s, c + (y - c) * s) for x, y in pts]


# ================================================================ 1. lab vessels (36 x 5 contents x 4 details = 720)
VESSEL_PROFILES = [(nw, bw, nb) for nw in (3.0, 5.0, 8.0) for bw in (6.0, 9.0, 12.0, 16.0) for nb in (14.0, 24.0, 34.0)]
CONTENTS = ["empty", "half", "full", "crystals", "gel"]
VDETAILS = ["plain", "scale", "label", "foot"]


def _vessel(nw, bw, nb, content, detail):
    cx, top, bb = 32.0, 7.0, 56.0
    r = min(7.0, bw - 1.0)

    def fn(ic):
        ic.line(cx - nw - 2.2, top, cx + nw + 2.2, top, SW + 0.6)
        wall = ("M%.2f %.2f V%.2f L%.2f %.2f V%.2f Q%.2f %.2f %.2f %.2f H%.2f Q%.2f %.2f %.2f %.2f V%.2f L%.2f %.2f V%.2f"
                % (cx - nw, top, nb, cx - bw, nb + 4, bb - r, cx - bw, bb, cx - bw + r, bb, cx + bw - r,
                   cx + bw, bb, cx + bw, bb - r, nb + 4, cx + nw, nb, top))
        ic.spath(wall, SW, glass=True)
        if content != "empty":
            level = {"half": bb - 16, "full": nb + 7, "crystals": bb - 12, "gel": bb - 16}[content]
            ic.sec().path(_pd([(cx - bw + 2.4, level), (cx + bw - 2.4, level), (cx + bw - 2.4, bb - 3), (cx - bw + 2.4, bb - 3)]))
        if content == "crystals":
            for dx, dy in ((-3, 2), (2, 0), (-1, 6), (4, 5), (-4, 8)):
                ic.acc().rect(cx + dx - 1.2, bb - 14 + dy, 2.4, 2.4, 0.3)
        if content == "gel":
            ic.acc().line(cx - bw + 3, bb - 16, cx + bw - 3, bb - 16, AW)
        if detail == "scale":
            for i in range(4):
                y = nb + 9 + i * (bb - nb - 16) / 4.0
                ic.acc().line(cx - bw + 1.6, y, cx - bw + 4.4, y, AW)
        if detail == "label":
            ic.acc().rect(cx - bw + 1.8, 36.0, 2 * bw - 3.6, 5.0, 0.8)
        if detail == "foot":
            ic.line(cx - bw - 3, bb + 1.4, cx + bw + 3, bb + 1.4, SW)
    return fn


for _i, (_nw, _bw, _nb) in enumerate(VESSEL_PROFILES):
    for _c in CONTENTS:
        for _d in VDETAILS:
            _reg("r22_vessel_%02d_%s_%s" % (_i, _c, _d), "glassware",
                 "Lab vessel %d (%s, %s)" % (_i + 1, _c, _d), ["vessel", _c, _d, "glass"],
                 _vessel(_nw, _bw, _nb, _c, _d))


# ================================================================ 2. instruments (24 bodies x 8 features x 3 details = 576)
def _instrument(w, h, feat, det, prefix_tags):
    x, y = 32.0 - w / 2, 32.0 - h / 2

    def fn(ic):
        ic.srect(x, y, w, h, 2.5, SW)
        if det == "feet":
            ic.line(x + 3, y + h + 1.4, x + 3, y + h + 4, SW)
            ic.line(x + w - 3, y + h + 1.4, x + w - 3, y + h + 4, SW)
        if det == "hood":
            ic.line(x - 1.5, y - 1.5, x + w + 1.5, y - 1.5, SW)
        if feat == "screen":
            ic.cut().rect(x + w * 0.18, y + h * 0.2, w * 0.64, h * 0.3, 1.2)
            ic.acc().line(x + w * 0.24, y + h * 0.36, x + w * 0.72, y + h * 0.36, AW)
        elif feat == "knob":
            ic.acc().circle(x + w * 0.5, y + h * 0.72, min(w, h) * 0.14)
        elif feat == "window":
            ic.cut().rect(x + w * 0.2, y + h * 0.3, w * 0.6, h * 0.4, 1.5)
        elif feat == "dial":
            r = min(w, h) * 0.2
            ic.acc().circle(x + w * 0.5, y + h * 0.5, r)
            ic.acc().line(x + w * 0.5, y + h * 0.5, x + w * 0.5 + r * 0.8, y + h * 0.5 - r * 0.5, AW)
        elif feat == "lever":
            ic.acc().line(x + w * 0.2, y + h * 0.72, x + w * 0.8, y + h * 0.42, AW + 0.6)
        elif feat == "vent":
            for i in range(3):
                ic.cut().rect(x + w * 0.22, y + h * (0.2 + i * 0.22), w * 0.56, 1.6, 0.6)
        elif feat == "bars":
            for i, f in enumerate((0.45, 0.8, 0.6)):
                ic.acc().rect(x + w * (0.22 + i * 0.22), y + h * (0.85 - f * 0.6), w * 0.14, h * f * 0.6, 0.5)
        elif feat == "keypad":
            for i in range(3):
                for j in range(2):
                    ic.cut().circle(x + w * (0.3 + i * 0.2), y + h * (0.35 + j * 0.28), min(w, h) * 0.05)
        elif feat == "terminal":
            ic.acc().rect(x + w * 0.3, y - 2.6, w * 0.4, 2.6, 0.5)
            ic.acc().rect(x + w * 0.3, y + h, w * 0.4, 2.0, 0.5)
    return fn


BODIES = [(w, h) for w in (16.0, 24.0, 32.0, 40.0) for h in (10.0, 16.0, 22.0, 28.0, 34.0, 40.0)]
FEATURES = ["screen", "knob", "window", "dial", "lever", "vent", "bars", "keypad"]
IDETAILS = ["plain", "feet", "hood"]
for _bi, (_w, _h) in enumerate(BODIES):
    for _f in FEATURES:
        for _d in IDETAILS:
            _reg("r22_instrument_%02d_%s_%s" % (_bi, _f, _d), "instruments",
                 "Instrument %d (%s, %s)" % (_bi + 1, _f, _d), ["instrument", _f, _d],
                 _instrument(_w, _h, _f, _d, None))

# ================================================================ 3. devices and cells (24 bodies x 6 x 2 = 288)
DEV_BODIES = [(w, h) for w in (14.0, 22.0, 30.0, 38.0) for h in (8.0, 12.0, 18.0, 24.0, 30.0, 36.0)]
DEV_FEATS = ["terminal", "window", "screen", "knob", "vent", "lever"]
DEV_DET = ["plain", "feet"]
for _bi, (_w, _h) in enumerate(DEV_BODIES):
    for _f in DEV_FEATS:
        for _d in DEV_DET:
            _reg("r22_device_%02d_%s_%s" % (_bi, _f, _d), "energy",
                 "Device %d (%s, %s)" % (_bi + 1, _f, _d), ["device", _f, _d],
                 _instrument(_w, _h, _f, _d, None))


# ================================================================ 4. material morphologies (21 shapes x 8 decorations x 3 sizes = 504)
def _ngon(n, r=22.0, rot=-90.0):
    return [(32 + r * math.cos(math.radians(rot + 360.0 * i / n)), 32 + r * math.sin(math.radians(rot + 360.0 * i / n))) for i in range(n)]


def _star(n, r1=22.0, r2=9.0):
    pts = []
    for i in range(2 * n):
        r = r1 if i % 2 == 0 else r2
        a = math.radians(-90 + 180.0 * i / n)
        pts.append((32 + r * math.cos(a), 32 + r * math.sin(a)))
    return pts


MORPH = [
    ("tri", "poly", _ngon(3)), ("square", "poly", _ngon(4, 20, -45)), ("pent", "poly", _ngon(5)),
    ("hex", "poly", _ngon(6)), ("oct", "poly", _ngon(8, 22, -22.5)),
    ("rod", "poly", [(26, 8), (38, 8), (40, 56), (24, 56)]),
    ("octahedron", "poly", [(32, 8), (52, 32), (32, 56), (12, 32)]),
    ("bipyramid", "poly", [(32, 8), (44, 32), (32, 56), (20, 32)]),
    ("nanocube", "poly", [(16, 16), (48, 16), (48, 48), (16, 48)]),
    ("sheet", "poly", [(14, 22), (50, 14), (54, 42), (18, 50)]),
    ("triplate", "poly", [(32, 10), (52, 48), (12, 48)]),
    ("rhombus", "poly", [(32, 8), (50, 32), (32, 56), (14, 32)]),
    ("star5", "poly", _star(5)), ("star6", "poly", _star(6, 22, 12)),
    ("sphere_cluster", "circles", [(22, 24, 9), (42, 22, 8), (32, 42, 10), (46, 44, 6)]),
    ("dumbbell", "circles", [(22, 32, 11), (42, 32, 11)]),
    ("core_shell", "circles", [(32, 32, 20), (32, 32, 12)]),
    ("dendrite", "lines", [[(32, 56), (32, 30)], [(32, 40), (20, 28)], [(32, 40), (44, 28)], [(32, 48), (24, 42)], [(32, 48), (40, 42)], [(32, 30), (30, 12)], [(32, 30), (36, 12)]]),
    ("wire_bundle", "lines", [[(26, 8), (28, 56)], [(32, 8), (32, 56)], [(38, 8), (36, 56)]]),
    ("tetrapod", "lines", [[(32, 32), (32, 8)], [(32, 32), (52, 44)], [(32, 32), (12, 44)], [(32, 32), (32, 56)]]),
    ("cage", "lines", [[(14, 14), (50, 14), (50, 50), (14, 50), (14, 14)], [(14, 32), (50, 32)], [(32, 14), (32, 50)]]),
]
DECOR = ["plain", "core", "hollow", "pores", "facets", "dots", "defect", "layers"]
MSIZES = [0.72, 0.86, 1.0]


def _morph(kind, data, dec, s):
    def fn(ic):
        if kind == "poly":
            pts = _scale(data, s)
            ic.spath(_pd(pts), SW)
        elif kind == "circles":
            for cx, cy, r in data:
                ic.scircle(32 + (cx - 32) * s, 32 + (cy - 32) * s, r * s, SW)
        else:
            for poly in data:
                ic.spath(_pd(_scale(poly, s), close=False), SW)
        if dec == "core":
            ic.acc().circle(32, 32, 6 * s)
        elif dec == "hollow":
            ic.cut().circle(32, 32, 8 * s)
        elif dec == "pores":
            for px, py in ((-6, -4), (5, -6), (2, 6)):
                ic.cut().circle(32 + px * s, 32 + py * s, 2.2 * s)
        elif dec == "facets":
            ic.line(32, 32 - 18 * s, 32, 32 + 18 * s, AW)
            ic.line(32 - 18 * s, 32, 32 + 18 * s, 32, AW)
        elif dec == "dots":
            for px, py in ((-8, -6), (6, -4), (-3, 7), (8, 7)):
                ic.acc().circle(32 + px * s, 32 + py * s, 2.1 * s)
        elif dec == "defect":
            ic.acc().circle(38 * s + 32 * (1 - s), 30, 4 * s)
            ic.cut().circle(38 * s + 32 * (1 - s), 30, 1.8 * s)
        elif dec == "layers" and kind == "poly":
            ic.spath(_pd(_scale(data, s * 0.6)), AW)
    return fn


for _name, _kind, _data in MORPH:
    for _dec in DECOR:
        for _si, _s in enumerate(MSIZES):
            _reg("r22_morph_%s_%s_%d" % (_name, _dec, _si), "materials",
                 "%s, %s (size %d)" % (_name.replace("_", " ").title(), _dec, _si + 1), ["material", _name, _dec],
                 _morph(_kind, _data, _dec, _s))


# ================================================================ 5. spectra and analysis traces (40 families x 5 = 200)
SPEC_FAMILIES = [
    ("xrd_pattern", "peaks", "XRD pattern"), ("pxrd_difference", "peaks", "PXRD difference"), ("rietveld_residual", "dips", "Rietveld residual"),
    ("ir_spectrum", "dips", "IR spectrum"), ("raman_spectrum", "peaks", "Raman spectrum"), ("uv_vis_band", "peaks", "UV-Vis absorbance"),
    ("pl_emission", "peaks", "PL emission"), ("nmr_multiplet", "sticks", "1H NMR multiplet"), ("xps_survey", "peaks", "XPS survey"),
    ("gc_chromatogram", "peaks", "GC chromatogram"), ("dls_distribution", "peaks", "DLS distribution"), ("pore_size", "peaks", "Pore size distribution"),
    ("eds_spectrum", "sticks", "EDS spectrum"), ("sec_chromatogram", "peaks", "SEC chromatogram"), ("ms_spectrum", "sticks", "Mass spectrum"),
    ("tga_steps", "steps", "TGA mass steps"), ("dsc_thermogram", "peaks", "DSC thermogram"), ("mossbauer", "dips", "Mossbauer spectrum"),
    ("sorption_isotherm", "hyst", "Sorption isotherm"), ("decay_kinetics", "decay", "Decay kinetics"), ("epr_spectrum", "peaks", "EPR spectrum"),
    ("cv_voltammogram", "peaks", "Cyclic voltammogram"), ("saxs_curve", "decay", "SAXS intensity"), ("tauc_plot", "decay", "Tauc plot"),
    ("impedance_bode", "decay", "Bode impedance"), ("ftir_difference", "dips", "FTIR difference"), ("sans_scattering", "decay", "SANS scattering"),
    ("pulse_response", "decay", "Pulse response"), ("thermal_conductivity", "steps", "Thermal conductivity"), ("nmr_relaxation", "decay", "NMR relaxation"),
    ("lifetime_fit", "decay", "Lifetime fit"), ("emission_spectrum", "peaks", "Emission spectrum"), ("fluorescence_anisotropy", "decay", "Anisotropy decay"),
    ("ocd_spectrum", "dips", "Optical constants"), ("reflectance", "dips", "Reflectance spectrum"), ("sieve_distribution", "steps", "Sieve distribution"),
    ("particle_histogram", "sticks", "Particle size histogram"), ("breakthrough", "steps", "Breakthrough curve"), ("stress_strain", "hyst", "Stress-strain"),
    ("creep_curve", "decay", "Creep curve"), ("rheology_flow", "decay", "Rheology flow curve"),
]


def _spec_data(kind, rng):
    if kind in ("peaks", "dips"):
        return [(rng.uniform(10, 54), rng.uniform(8, 28), rng.uniform(0.6, 1.6)) for _ in range(rng.randint(2, 5))]
    if kind == "sticks":
        return [(rng.uniform(10, 52), rng.uniform(6, 36)) for _ in range(rng.randint(4, 8))]
    if kind == "steps":
        return [(rng.uniform(14, 50), rng.uniform(4, 12)) for _ in range(rng.randint(2, 4))]
    if kind == "decay":
        return [(rng.uniform(0.05, 0.3),)]
    return [(rng.uniform(0.25, 0.6),)]


def _poly_d(pts):
    return "M" + " L".join("%.2f %.2f" % p for p in pts)


def _spectrum(kind, data):
    def fn(ic):
        axes(ic, SW)
        if kind in ("peaks", "dips"):
            pts = []
            for i in range(61):
                x = 11 + 44 * i / 60.0
                sig = sum(h * math.exp(-((x - c) / w) ** 2) for c, h, w in data)
                if kind == "peaks":
                    pts.append((x, 52 - min(38.0, sig * 38.0 / 28.0)))
                else:
                    pts.append((x, 14 + min(38.0, sig * 38.0 / 28.0)))
            ic.spath(_poly_d(pts), 2.0)
        elif kind == "sticks":
            for c, h in data:
                x = 11 + 44 * (c - 10) / 44.0
                ic.line(x, 52, x, 52 - h, 2.2)
        elif kind == "steps":
            pts = [(11, 14)]
            level = 14.0
            for x, drop in sorted(data):
                pts.append((x, level))
                level += drop
                pts.append((x, level))
            pts.append((55, level))
            ic.spath(_poly_d(pts), 2.0)
        elif kind == "decay":
            k = data[0][0]
            pts = [(11 + 44 * i / 40.0, 14 + 38 * (1 - math.exp(-k * 9 * (1 - i / 40.0)))) for i in range(41)]
            ic.spath(_poly_d(pts), 2.0)
        else:
            a = data[0][0]
            ads = [(11 + 44 * i / 30.0, 52 - 38 * (i / 30.0) ** 2) for i in range(31)]
            des = [(11 + 44 * i / 30.0, 52 - 38 * ((i / 30.0) ** 0.5) * (1 - a) - 2) for i in range(31)]
            ic.spath(_poly_d(ads), 2.0)
            ic.spath(_poly_d(des), 1.6)
            ic.acc().line(11, 52, 55, 52, AW)
    return fn


_SPEC_RNG = random.Random(220)
for _fam, _kind, _label in SPEC_FAMILIES:
    for _k in range(5):
        _data = _spec_data(_kind, _SPEC_RNG)
        _reg("r22_spec_%s_%d" % (_fam, _k), "characterisation", "%s %d" % (_label, _k + 1),
             ["spectrum", _fam.replace("_", " "), _kind], _spectrum(_kind, _data))


# ================================================================ 6. well plates (12 layouts x 6 fills x 3 frames = 216)
PLATES = [(1, 6), (2, 3), (2, 4), (3, 4), (4, 6), (6, 8), (8, 12), (4, 4), (3, 6), (2, 8), (5, 5), (1, 12)]
PFILLS = ["empty", "row", "checker", "full", "diagonal", "dots"]
PFRAMES = ["sharp", "rounded", "chamfer"]


def _plate(rows, cols, fill, frame):
    x0, y0, w, h = 6.0, 14.0, 52.0, 36.0
    well = min(w / cols, h / rows) / 2 - 1.2

    def fn(ic):
        if frame == "rounded":
            ic.srect(x0, y0, w, h, 4, SW)
        elif frame == "chamfer":
            ic.spath(_pd([(x0 + 4, y0), (x0 + w - 4, y0), (x0 + w, y0 + 4), (x0 + w, y0 + h - 4),
                          (x0 + w - 4, y0 + h), (x0 + 4, y0 + h), (x0, y0 + h - 4), (x0, y0 + 4)]), SW)
        else:
            ic.srect(x0, y0, w, h, 0.5, SW)
        for r in range(rows):
            for c in range(cols):
                cx = x0 + (c + 0.5) * w / cols
                cy = y0 + (r + 0.5) * h / rows
                on = {"empty": False, "row": r == 0, "checker": (r + c) % 2 == 0, "full": True,
                      "diagonal": c == r, "dots": (r + c) % 3 == 0}[fill]
                if on:
                    ic.acc().circle(cx, cy, well)
                else:
                    ic.scircle(cx, cy, well, AW)
    return fn


for _pi, (_r, _c) in enumerate(PLATES):
    for _f in PFILLS:
        for _fr in PFRAMES:
            _reg("r22_plate_%02d_%s_%s" % (_pi, _f, _fr), "instruments",
                 "Plate %dx%d (%s, %s)" % (_r, _c, _f, _fr), ["plate", "well", _f],
                 _plate(_r, _c, _f, _fr))


# ================================================================ 7. pipes and valves (8 styles x 5 diameters x 2 ends x 2 fills = 160)
PIPE_STYLES = ["straight", "elbow", "tee", "cross", "gate_valve", "ball_valve", "check_valve", "flange"]
PIPE_D = [2.0, 3.0, 4.0, 5.0, 6.5]
PIPE_END = ["open", "capped"]
PIPE_FILL = ["empty", "fluid"]


def _pipe(style, d, end, fill):
    def fn(ic):
        y, x0, x1 = 32.0, 8.0, 56.0
        if style == "elbow":
            ic.spath(_pd([(x0, y - d), (32 + d, y - d), (32 + d, 56 - 0), (32 - d, 56), (32 - d, y + d), (x0, y + d)], close=False), SW)
        elif style in ("straight", "flange"):
            ic.srect(x0, y - d, x1 - x0, 2 * d, 0.6, SW)
        else:
            ic.srect(x0, y - d, x1 - x0, 2 * d, 0.6, SW)
        if style == "tee":
            ic.srect(32 - d, y + d, 2 * d, 22 - d, 0.6, SW)
        if style == "cross":
            ic.srect(32 - d, 8, 2 * d, 48, 0.6, SW)
        if style == "flange":
            ic.acc().rect(x0 - 2, y - d - 2.5, 2.2, 2 * d + 5, 0.5)
            ic.acc().rect(x1 - 0.2, y - d - 2.5, 2.2, 2 * d + 5, 0.5)
        if style == "gate_valve":
            ic.spath(_pd([(28, y - 7), (36, y + 7), (28, y + 7), (36, y - 7)], close=True), AW)
        if style == "ball_valve":
            ic.acc().circle(32, y, 4.2)
            ic.scircle(32, y, 4.2, AW)
        if style == "check_valve":
            ic.spath(_pd([(28, y - 6), (28, y + 6), (36, y)]), AW)
        if end == "capped":
            ic.acc().rect(x1 + 1.6, y - d - 1.6, 2.4, 2 * d + 3.2, 0.6)
        if fill == "fluid":
            ic.sec().rect(x0 + 1, y - d + 1.2, (x1 - x0) - 2, 2 * d - 2.4, 0.5)
    return fn


for _st in PIPE_STYLES:
    for _d in PIPE_D:
        for _e in PIPE_END:
            for _fl in PIPE_FILL:
                _reg("r22_pipe_%s_%s_%s_%s" % (_st, str(_d).replace(".", "p"), _e, _fl), "process",
                     "Pipe %s %s, %s, %s" % (_st.replace("_", " "), _d, _e, _fl), ["pipe", "valve", _st.replace("_", " ")],
                     _pipe(_st, _d, _e, _fl))


# ================================================================ 8. reactors (4 shapes x 8 options x 4 levels = 128)
REACTORS = ["cylinder", "round_flask", "jar", "tall"]
ROPTS = ["plain", "stirrer", "heater", "condenser", "gas_inlet", "thermometer", "sampler", "jacket"]
RLEVELS = [0.3, 0.5, 0.7, 0.9]


def _reactor(shape, opt, lv):
    def fn(ic):
        if shape == "cylinder":
            ic.srect(18, 16, 28, 40, 3, SW)
            ic.line(16, 16, 48, 16, SW)
        elif shape == "round_flask":
            ic.spath(_pd([(20, 14), (20, 30), (12, 46), (14, 54), (50, 54), (52, 46), (44, 30), (44, 14)], close=False), SW)
            ic.line(18, 14, 46, 14, SW)
        elif shape == "jar":
            ic.srect(14, 20, 36, 36, 9, SW)
            ic.line(16, 20, 48, 20, SW)
        else:
            ic.srect(22, 10, 20, 46, 3, SW)
            ic.line(20, 10, 44, 10, SW)
        top = 16 if shape in ("cylinder", "tall") else 22
        bottom = 56 if shape in ("round_flask",) else 52
        level = bottom - (bottom - top) * lv
        ic.sec().rect(22 if shape == "tall" else 20, level, 24 if shape == "tall" else 24, bottom - level - 2, 0.5)
        if opt == "stirrer":
            ic.line(32, 12, 32, level + 6, AW)
            ic.acc().line(26, level + 6, 38, level + 6, AW + 0.4)
        elif opt == "heater":
            ic.acc().spath(_pd([(16, 58), (20, 61), (24, 58), (28, 61), (32, 58), (36, 61), (40, 58), (44, 61), (48, 58)], close=False), AW)
        elif opt == "condenser":
            for i in range(4):
                ic.acc().circle(50 + (i % 2) * 0, 18 + i * 6, 1.2)
            ic.acc().line(46, 18, 60, 18, AW)
        elif opt == "gas_inlet":
            ic.acc().line(32, 4, 32, 14, AW)
        elif opt == "thermometer":
            ic.acc().line(46, 12, 46, level, AW)
            ic.acc().circle(46, level + 1, 1.3)
        elif opt == "sampler":
            ic.acc().rect(54, level - 1, 6, 3, 0.4)
        elif opt == "jacket":
            ic.scircle(32, 36, 22, AW)
    return fn


for _si, _sh in enumerate(REACTORS):
    for _o in ROPTS:
        for _li, _lv in enumerate(RLEVELS):
            _reg("r22_reactor_%s_%s_%d" % (_sh, _o, _li), "instruments",
                 "Reactor %s, %s (level %d)" % (_sh.replace("_", " "), _o, _li + 1), ["reactor", _sh.replace("_", " "), _o],
                 _reactor(_sh, _o, _lv))


# ================================================================ 9. crystal lattices (16 cells x 8 decorations = 128)
LATTICES = [(sa, sb, sc, sk) for sa in (1.0, 1.25) for sb in (1.0, 1.2) for sc in (1.0, 1.4) for sk in (0.0, 0.5)]
LDECOR = ["plain", "corners", "centre", "faces", "edges", "shaded", "cut", "dots"]


def _lattice(sa, sb, sc, sk, dec):
    o = (16.0, 50.0)
    a = (26.0 * sa, 0.0)
    b = (sk * 14.0 * sb, -14.0 * sb)
    c = (0.0, -22.0 * sc)

    def P(i, j, k):
        return (o[0] + i * a[0] + j * b[0] + k * c[0], o[1] + i * a[1] + j * b[1] + k * c[1])

    def fn(ic):
        front = [P(0, 0, 0), P(1, 0, 0), P(1, 0, 1), P(0, 0, 1)]
        back = [P(0, 1, 0), P(1, 1, 0), P(1, 1, 1), P(0, 1, 1)]
        ic.spath(_pd(front), SW)
        ic.spath(_pd(back), AW)
        for p, q in ((P(0, 0, 0), P(0, 1, 0)), (P(1, 0, 0), P(1, 1, 0)), (P(1, 0, 1), P(1, 1, 1)), (P(0, 0, 1), P(0, 1, 1))):
            ic.line(p[0], p[1], q[0], q[1], AW)
        if dec == "shaded":
            ic.sec().path(_pd(front))
        if dec == "corners":
            for i in (0, 1):
                for j in (0, 1):
                    for k in (0, 1):
                        p = P(i, j, k)
                        ic.acc().circle(p[0], p[1], 2.0)
        if dec == "centre":
            p = P(0.5, 0.5, 0.5)
            ic.acc().circle(p[0], p[1], 3.2)
        if dec == "faces":
            p = P(0.5, 0.0, 0.5)
            ic.acc().circle(p[0], p[1], 2.4)
            p = P(0.5, 1.0, 0.5)
            ic.acc().circle(p[0], p[1], 2.4)
        if dec == "edges":
            p, q = P(0, 0, 0.5), P(1, 0, 0.5)
            ic.acc().line(p[0], p[1], q[0], q[1], AW)
        if dec == "cut":
            p = P(0.5, 0.0, 0.5)
            ic.cut().circle(p[0], p[1], 2.6)
        if dec == "dots":
            for i, j, k in ((0.25, 0.25, 0.25), (0.75, 0.25, 0.75), (0.5, 0.75, 0.5)):
                p = P(i, j, k)
                ic.acc().circle(p[0], p[1], 1.6)
    return fn


for _li, (_sa, _sb, _sc, _sk) in enumerate(LATTICES):
    for _dec in LDECOR:
        _reg("r22_lattice_%02d_%s" % (_li, _dec), "materials", "Crystal lattice %d (%s)" % (_li + 1, _dec.replace("_", " ")),
             ["crystal", "lattice", "unit cell", _dec], _lattice(_sa, _sb, _sc, _sk, _dec))


# ================================================================ 10. gas cylinders and ports (4 caps x 3 bands x 5 sizes x 2 valves = 120)
CAPS = ["flat", "valve", "dome", "collar"]
BANDS = ["none", "stripe", "label"]
GSIZES = [10.0, 12.0, 14.0, 16.0, 18.0]
GVALVE = ["plain", "regulator"]


def _gas(cap, band, w, valve):
    x = 32.0 - w / 2

    def fn(ic):
        ic.srect(x, 20, w, 36, w / 2.6, SW)
        if cap == "valve":
            ic.srect(x + w * 0.3, 12, w * 0.4, 8, 1.0, SW)
        elif cap == "dome":
            ic.spath(_pd([(x + 1, 20), (x + 1, 16), (32, 10), (x + w - 1, 16), (x + w - 1, 20)], close=False), SW)
        elif cap == "collar":
            ic.srect(x - 1, 20, w + 2, 4, 0.6, SW)
        else:
            ic.line(x + w * 0.3, 20, x + w * 0.7, 20, SW)
        if band == "stripe":
            ic.acc().rect(x, 34, w, 3, 0.2)
        elif band == "label":
            ic.acc().rect(x + 1.6, 40, w - 3.2, 8, 0.8)
        if valve == "regulator":
            ic.line(32, 10, 32, 6, AW)
            ic.acc().line(26, 4, 38, 4, AW)
    return fn


for _cap in CAPS:
    for _band in BANDS:
        for _w in GSIZES:
            for _v in GVALVE:
                _reg("r22_gas_%s_%s_%d_%s" % (_cap, _band, int(_w), _v), "instruments",
                     "Gas cylinder %s, %s, %d mm, %s" % (_cap, _band, int(_w), _v), ["gas", "cylinder", _cap, _v],
                     _gas(_cap, _band, _w, _v))


# ================================================================ 11. process arrows and boxes (8 x 10 + 4 x 5 x 4 = 160)
ASTYLES = ["straight", "curved", "elbow", "double", "step", "recycle", "tee", "split"]
APROPS = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]


def _arrow(style, prop):
    L = 40.0 * prop
    x0, x1 = 32 - L / 2, 32 + L / 2

    def fn(ic):
        head = _pd([(x1 - 8, 26), (x1, 32), (x1 - 8, 38)], close=True)
        if style == "straight":
            ic.line(x0, 32, x1 - 6, 32, 2.8)
            ic.acc().spath(head, 0.6)
        elif style == "curved":
            ic.spath("M%.2f 44 Q32 8 %.2f 30" % (x0, x1 - 6), 2.8)
            ic.acc().spath(_pd([(x1 - 10, 24), (x1 - 3, 30), (x1 - 9, 36)], close=False), 1.8)
        elif style == "elbow":
            ic.spath(_pd([(x0, 42), (32, 42), (32, 22), (x1 - 6, 22)], close=False), 2.8)
            ic.acc().spath(_pd([(x1 - 8, 16), (x1, 22), (x1 - 8, 28)], close=False), 1.8)
        elif style == "double":
            ic.line(x0 + 6, 32, x1 - 6, 32, 2.6)
            ic.acc().spath(head, 0.6)
            ic.acc().spath(_pd([(x0 + 8, 26), (x0, 32), (x0 + 8, 38)], close=True), 0.6)
        elif style == "step":
            ic.spath(_pd([(x0, 44), (29, 44), (29, 32), (x1 - 6, 32)], close=False), 2.8)
            ic.acc().spath(_pd([(x1 - 8, 26), (x1, 32), (x1 - 8, 38)], close=False), 1.8)
        elif style == "recycle":
            ic.spath("M%.2f 40 Q%.2f 16 32 16 Q%.2f 16 %.2f 40" % (x0, x0, x1, x1), 2.6)
            ic.acc().spath(_pd([(x1 - 5, 34), (x1, 42), (x1 + 4, 34)], close=False), 1.8)
        elif style == "tee":
            ic.line(x0, 32, x1 - 6, 32, 2.8)
            ic.line(32, 32, 32, 52, 2.8)
            ic.acc().spath(head, 0.6)
        else:
            ic.line(x0, 32, 32, 32, 2.8)
            ic.line(32, 32, x1 - 6, 22, 2.8)
            ic.line(32, 32, x1 - 6, 42, 2.8)
            ic.acc().spath(head, 0.6)
    return fn


for _st in ASTYLES:
    for _p in APROPS:
        _reg("r22_arrow_%s_%02d" % (_st, int(_p * 10)), "process", "Process arrow %s %d" % (_st, int(_p * 100)),
             ["arrow", "process flow", _st], _arrow(_st, _p))

BLOCKS = [(s, z) for s in ("box", "rounded", "diamond", "circle") for z in (8.0, 11.0, 14.0, 17.0, 20.0)]
AFLOW = ["straight", "double", "curved", "branch"]


def _block(shape, z, flow):
    cx, cy = 22.0, 32.0

    def fn(ic):
        if shape == "box":
            ic.srect(cx - z, cy - z * 0.7, 2 * z, 1.4 * z, 1.0, SW)
        elif shape == "rounded":
            ic.srect(cx - z, cy - z * 0.7, 2 * z, 1.4 * z, z * 0.5, SW)
        elif shape == "diamond":
            ic.spath(_pd([(cx, cy - z), (cx + z, cy), (cx, cy + z), (cx - z, cy)]), SW)
        else:
            ic.scircle(cx, cy, z, SW)
        x = cx + z + 3
        if flow == "straight":
            ic.line(x, cy, 54, cy, 2.6)
        elif flow == "double":
            ic.line(x, cy - 2, 54, cy - 2, 2.2)
            ic.line(x, cy + 2, 54, cy + 2, 2.2)
        elif flow == "curved":
            ic.spath("M%.2f %.2f Q%.2f %.2f 54 %.2f" % (x, cy + 4, 40, cy + 16, cy + 4), 2.4)
        else:
            ic.line(x, cy, 46, cy, 2.4)
            ic.line(46, cy, 54, cy - 8, 2.4)
            ic.line(46, cy, 54, cy + 8, 2.4)
    return fn


for _sh, _z in BLOCKS:
    for _fl in AFLOW:
        _reg("r22_block_%s_%02d_%s" % (_sh, int(_z), _fl), "process", "Process block %s %d (%s)" % (_sh, int(_z), _fl),
             ["process block", _sh, _fl], _block(_sh, _z, _fl))


# ================================================================ 12. charts and data glyphs (6 kinds x 20 = 120 x 2 = 240)
CHART_KINDS = ["bars", "lines", "scatter", "box", "heat", "area"]
CHART_NAMES = ["bar chart", "grouped bars", "line plot", "multi line", "scatter plot", "bubble plot", "box plot", "violin",
               "heatmap", "dot grid", "area plot", "stacked area", "histogram", "lollipop", "error bars", "step plot",
               "log-log line", "dose response", "radar plot", "stacked bars"]


def _chart(kind, rng):
    def fn(ic):
        axes(ic, SW)
        if kind == "bars":
            for i in range(5):
                h = rng.uniform(8, 36)
                ic.acc().rect(13 + i * 8.4, 52 - h, 5.2, h, 0.4)
        elif kind == "lines":
            for amp in (1.0, 0.6):
                pts = [(12 + 42 * i / 30.0, 50 - amp * 30 * (1 - math.exp(-3 * i / 30.0)) - rng.uniform(0, 3)) for i in range(31)]
                ic.spath(_poly_d(pts), 2.0)
        elif kind == "scatter":
            for _ in range(12):
                ic.acc().circle(rng.uniform(14, 52), rng.uniform(12, 50), 1.6)
        elif kind == "box":
            for i in range(3):
                x = 18 + i * 14
                ic.srect(x - 4, 26, 8, 14, 0.5, AW)
                ic.line(x - 4, 33, x + 4, 33, AW)
                ic.line(x, 20, x, 26, AW)
                ic.line(x, 40, x, 50, AW)
        elif kind == "heat":
            for r in range(4):
                for c in range(5):
                    if rng.random() < 0.5:
                        ic.acc().rect(13 + c * 8, 12 + r * 9, 7, 8, 0.4)
                    else:
                        ic.srect(13 + c * 8, 12 + r * 9, 7, 8, 0.4, AW)
        else:
            pts = [(12 + 42 * i / 30.0, 50 - 34 * (1 - math.exp(-2.6 * i / 30.0))) for i in range(31)]
            ic.sec().path(_pd(pts + [(54, 50), (12, 50)]))
            ic.spath(_poly_d(pts), 2.0)
    return fn


for _ki, _kind in enumerate(CHART_KINDS):
    for _ni in range(20):
        _rng = random.Random(1000 * _ki + _ni)
        _reg("r22_chart_%s_%02d" % (_kind, _ni), "characterisation", "%s %d" % (CHART_NAMES[_ni], _ni + 1),
             ["chart", CHART_NAMES[_ni], _kind], _chart(_kind, _rng))


# ================================================================ 13. tube racks (5 x 3 rows x 8 fills = 120)
RACKS = [(n, r) for n in (2, 3, 4, 5, 6) for r in (1, 2, 3)]
TFILL = ["empty", "low", "mid", "full", "mixed", "marked", "labelled", "capped"]


def _rack(n, rows, fill):
    def fn(ic):
        ic.srect(6, 54, 52, 4, 1.0, SW)
        tw = 44.0 / n
        for r in range(rows):
            top = 12 + r * 19
            for i in range(n):
                x = 10 + i * tw + tw * 0.12
                w = tw * 0.76
                h = 26.0
                ic.srect(x, top, w, h, min(w / 2, 3.2), SW)
                lvl = {"low": top + 18, "mid": top + 12, "full": top + 5, "mixed": top + 12, "marked": top + 12,
                       "labelled": top + 12, "capped": top + 12}.get(fill)
                if fill != "empty" and lvl is not None:
                    ic.sec().rect(x + 1.2, lvl, w - 2.4, top + h - lvl - 1.2, 0.5)
                if fill == "mixed" and i % 2 == 0:
                    ic.acc().rect(x + 1.2, top + 12, w - 2.4, 4, 0.5)
                if fill == "marked":
                    for k in range(3):
                        ic.acc().line(x, top + 6 + k * 6, x + 2.4, top + 6 + k * 6, AW)
                if fill == "labelled":
                    ic.acc().rect(x + 1.2, top + 9, w - 2.4, 4.5, 0.5)
                if fill == "capped":
                    ic.acc().rect(x - 0.8, top - 2.4, w + 1.6, 2.6, 0.5)
    return fn


for _n, _r in RACKS:
    for _fl in TFILL:
        _reg("r22_rack_%d_%d_%s" % (_n, _r, _fl), "glassware", "Tube rack %d x %d (%s)" % (_n, _r, _fl),
             ["rack", "test tube", _fl], _rack(_n, _r, _fl))


# ================================================================ 14. process symbols and annotations (20 numbered + 10 brackets + 10 callouts + 144 combined = 184)
def _numbered(n):
    def fn(ic):
        ic.acc().circle(32, 32, 20)
        for i in range(n):
            ang = 2 * math.pi * i / max(n, 1)
            ic.cut().circle(32 + 12 * math.cos(ang), 32 + 12 * math.sin(ang), 1.8 if n > 12 else 2.4)
    return fn


for _n in range(1, 21):
    _reg("r22_step_%02d" % _n, "process", "Step marker %d" % _n, ["step", "marker", "annotation"], _numbered(_n))


def _bracket(side, h):
    def fn(ic):
        x = 26 if side == "left" else 38
        d = 4 if side == "left" else -4
        ic.spath(_pd([(x + d, 32 - h), (x, 32 - h), (x, 32 + h), (x + d, 32 + h)], close=False), SW)
    return fn


for _i in range(10):
    _reg("r22_bracket_%s_%02d" % ("left" if _i < 5 else "right", _i), "process", "Bracket %d" % (_i + 1),
         ["bracket", "annotation"], _bracket("left" if _i < 5 else "right", 8 + _i * 2.4))


def _callout(i):
    def fn(ic):
        w, h = 30 + i, 16 + (i % 3) * 2
        x, y = 32 - w / 2, 22 - h / 2
        ic.srect(x, y, w, h, 2.5, SW)
        ic.spath(_pd([(x + 10, y + h), (x + 6, y + h + 6), (x + 16, y + h)], close=False), SW)
    return fn


for _i in range(10):
    _reg("r22_callout_%02d" % _i, "process", "Callout %d" % (_i + 1), ["callout", "annotation"], _callout(_i))

BASES = ["circle", "square", "triangle", "hexagon", "diamond", "octagon", "pentagon", "star", "ring", "plus", "cross", "bar"]
OVERLAYS = ["none", "check", "cross", "warning", "arrow_up", "arrow_down", "dot", "ring", "bar_h", "bar_v", "tri", "crosshair"]


def _combined(base, ov):
    def fn(ic):
        if base == "circle":
            ic.scircle(32, 32, 20, SW)
        elif base == "square":
            ic.srect(12, 12, 40, 40, 2, SW)
        elif base == "triangle":
            ic.spath(_pd([(32, 10), (54, 52), (10, 52)]), SW)
        elif base == "hexagon":
            ic.spath(_pd(_ngon(6, 22)), SW)
        elif base == "diamond":
            ic.spath(_pd([(32, 10), (54, 32), (32, 54), (10, 32)]), SW)
        elif base == "octagon":
            ic.spath(_pd(_ngon(8, 22, -22.5)), SW)
        elif base == "pentagon":
            ic.spath(_pd(_ngon(5, 22)), SW)
        elif base == "star":
            ic.spath(_pd(_star(5, 22, 10)), SW)
        elif base == "ring":
            ic.scircle(32, 32, 20, 4.0)
        elif base == "plus":
            ic.line(32, 12, 32, 52, 4.0)
            ic.line(12, 32, 52, 32, 4.0)
        elif base == "cross":
            ic.line(14, 14, 50, 50, 4.0)
            ic.line(50, 14, 14, 50, 4.0)
        else:
            ic.srect(12, 26, 40, 12, 2, SW)
        if ov == "check":
            ic.acc().spath(_pd([(22, 33), (29, 40), (42, 26)], close=False), AW + 1)
        elif ov == "cross":
            ic.acc().line(24, 24, 40, 40, AW + 1)
            ic.acc().line(40, 24, 24, 40, AW + 1)
        elif ov == "warning":
            ic.acc().spath(_pd([(32, 22), (42, 40), (22, 40)]), AW)
        elif ov == "arrow_up":
            ic.acc().spath(_pd([(32, 22), (40, 32), (24, 32)]), AW)
        elif ov == "arrow_down":
            ic.acc().spath(_pd([(32, 42), (40, 32), (24, 32)]), AW)
        elif ov == "dot":
            ic.acc().circle(32, 32, 3.4)
        elif ov == "ring":
            ic.acc().scircle(32, 32, 6, AW)
        elif ov == "bar_h":
            ic.acc().line(22, 32, 42, 32, AW + 1)
        elif ov == "bar_v":
            ic.acc().line(32, 22, 32, 42, AW + 1)
        elif ov == "tri":
            ic.acc().spath(_pd([(32, 26), (38, 36), (26, 36)]), AW)
        elif ov == "crosshair":
            ic.acc().line(32, 20, 32, 44, AW)
            ic.acc().line(20, 32, 44, 32, AW)
    return fn


for _b in BASES:
    for _o in OVERLAYS:
        _reg("r22_symbol_%s_%s" % (_b, _o), "process", "Symbol %s + %s" % (_b, _o.replace("_", " ")),
             ["symbol", _b, _o], _combined(_b, _o))


# ================================================================ 15. centrifuge rotors (8 hole counts x 3 radii x 4 modes = 96)
ROTOR_N = [2, 3, 4, 5, 6, 8, 10, 12]
ROTOR_R = [14.0, 19.0, 24.0]
ROTOR_MODE = ["empty", "tube", "filled", "labelled"]


def _rotor(n, rad, mode):
    def fn(ic):
        ic.scircle(32, 32, 27, SW)
        ic.scircle(32, 32, 4, AW)
        hr = min(4.0, 2.0 + rad / 12.0)
        for i in range(n):
            a = 2 * math.pi * i / n
            x, y = 32 + rad * math.cos(a), 32 + rad * math.sin(a)
            if mode == "empty":
                ic.scircle(x, y, hr, AW)
            elif mode == "tube":
                ic.srect(x - hr * 0.6, y - hr, hr * 1.2, hr * 2, hr * 0.5, AW)
            elif mode == "filled":
                ic.acc().circle(x, y, hr * 0.7)
            else:
                ic.scircle(x, y, hr, AW)
                ic.acc().line(x - hr * 0.5, y, x + hr * 0.5, y, 0.9)
    return fn


for _n in ROTOR_N:
    for _r in ROTOR_R:
        for _m in ROTOR_MODE:
            _reg("r22_rotor_%02d_%02d_%s" % (_n, int(_r), _m), "instruments",
                 "Centrifuge rotor %d holes, r %d (%s)" % (_n, int(_r), _m), ["rotor", "centrifuge", _m],
                 _rotor(_n, _r, _m))


# ================================================================ 16. cuvettes and sample cells (8 widths x 8 fills x 4 caps = 256)
CUV_W = [8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0]
CUV_FILL = [0.15, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
CUV_CAP = ["open", "screw", "stopper", "septum"]


def _cuvette(w, fill, cap):
    x0, y1 = 32.0 - w / 2, 54.0
    h = 30.0
    top = y1 - h

    def fn(ic):
        ic.srect(x0, top, w, h, 0.8, SW)
        lvl = y1 - h * fill
        ic.sec().rect(x0 + 1.2, lvl, w - 2.4, y1 - lvl - 1.2, 0.3)
        if cap == "screw":
            ic.srect(x0 + 1, top - 5, w - 2, 5, 1.0, AW)
        elif cap == "stopper":
            ic.srect(x0 + w * 0.2, top - 6, w * 0.6, 6, 1.5, AW)
        elif cap == "septum":
            ic.acc().rect(x0 + w * 0.25, top - 4, w * 0.5, 3, 0.8)
    return fn


for _w in CUV_W:
    for _f in CUV_FILL:
        for _c in CUV_CAP:
            _reg("r22_cuvette_%02d_%02d_%s" % (int(_w), int(_f * 100), _c), "glassware",
                 "Cuvette %d mm, fill %d %% (%s)" % (int(_w), int(_f * 100), _c), ["cuvette", "sample cell", _c],
                 _cuvette(_w, _f, _c))


# ================================================================ count
NEW_COUNT = len(_NEW)
