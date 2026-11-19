"""
Reticula 2.18 icon set: 300 original vector icons (CC0), built from four parametric families.

  spectra        (90)  analysis traces: XRD, IR, Raman, UV-Vis, PL, NMR, XPS, GC, DLS, TGA, DSC, ...
  vessels        (60)  twelve lab-vessel silhouettes x five contents (empty, half, full, crystals, gel)
  instruments    (60)  ten instrument bodies x six front features (screen, knob, window, dial, lever, vent)
  morphologies   (60)  ten material shapes x six decorations (core, hollow, pores, facets, aggregate, defect)
  process        (30)  six flow / process arrow styles x five proportions

Every icon is registered with @icon-style builders, so the toolkit gives each one its five variants.
Names are checked against the built-in registry at import time: a clash raises instead of overwriting.
"""
import math

from reticon.library import icon, REGISTRY, axes

_TAKEN_AT_IMPORT = set(REGISTRY)   # built-ins are registered before user icons are loaded
_NEW = []


def _reg(name, category, label, tags, fn):
    if name in _TAKEN_AT_IMPORT or name in _NEW:
        raise ValueError(f"icon name already in use: {name}")
    _NEW.append(name)
    icon(name, category, label=label, tags=tags)(fn)


def _poly_d(pts):
    return "M" + " L".join(f"{x:.2f} {y:.2f}" for x, y in pts)


# ------------------------------------------------------------------ spectra ---------------------------------------
SPECTRA = {
    # family: (kind, label, tags, peak sets per variant)
    "xrd_pattern":        ("peaks", "XRD pattern", ["PXRD", "diffraction", "Bragg"],
                           [[(16, 9, 0.6), (21, 22, 0.6), (26, 13, 0.7), (33, 30, 0.7), (40, 8, 0.8), (46, 11, 0.8)],
                            [(14, 14, 0.5), (19, 28, 0.6), (30, 10, 0.7), (37, 24, 0.7), (44, 6, 0.9)],
                            [(18, 18, 0.6), (24, 8, 0.6), (28, 26, 0.7), (35, 12, 0.8), (42, 20, 0.7), (49, 7, 0.9)],
                            [(12, 6, 0.5), (22, 30, 0.6), (31, 16, 0.7), (38, 10, 0.8)],
                            [(15, 20, 0.5), (20, 12, 0.6), (27, 28, 0.6), (34, 9, 0.7), (45, 15, 0.8)]]),
    "ir_spectrum":        ("dips", "IR spectrum", ["FTIR", "vibrational", "ATR"],
                           [[(20, 14, 1.2), (36, 26, 1.5), (46, 9, 1.0)],
                            [(14, 10, 1.4), (28, 20, 1.0), (40, 30, 1.6), (52, 8, 1.2)],
                            [(24, 28, 1.2), (34, 12, 1.0)],
                            [(18, 8, 2.0), (32, 22, 1.8), (48, 16, 1.4)],
                            [(12, 18, 1.1), (26, 10, 1.3), (38, 24, 1.2), (50, 12, 1.0)]]),
    "raman_spectrum":     ("peaks", "Raman spectrum", ["Raman", "vibrational", "Stokes"],
                           [[(18, 10, 0.7), (30, 30, 0.8), (42, 16, 1.0)],
                            [(16, 22, 0.6), (26, 12, 0.7), (38, 28, 0.9), (50, 8, 0.8)],
                            [(22, 26, 1.0), (34, 14, 0.7)],
                            [(14, 14, 0.8), (24, 30, 0.7), (32, 9, 0.9), (48, 20, 1.1)],
                            [(20, 18, 0.9), (36, 26, 1.2), (46, 7, 0.7)]]),
    "uv_vis_band":        ("peaks", "UV-Vis absorbance", ["UV-Vis", "absorbance", "band gap"],
                           [[(24, 30, 4.0)], [(20, 18, 3.5), (40, 26, 4.5)], [(30, 22, 6.0)],
                            [(18, 28, 3.0), (34, 14, 3.8)], [(26, 24, 5.0), (44, 10, 4.0)]]),
    "pl_emission":        ("peaks", "PL emission", ["photoluminescence", "emission", "phosphor"],
                           [[(34, 30, 3.0)], [(28, 24, 2.5), (42, 10, 2.2)], [(38, 28, 4.0)],
                            [(24, 14, 2.0), (34, 30, 3.5)], [(40, 22, 2.8)]]),
    "nmr_multiplet":      ("peaks", "NMR multiplet", ["1H NMR", "multiplet", "J coupling"],
                           [[(26, 26, 0.6), (30, 26, 0.6)], [(22, 14, 0.5), (26, 30, 0.5), (30, 30, 0.5), (34, 14, 0.5)],
                            [(20, 10, 0.4), (28, 30, 0.4), (36, 10, 0.4)], [(30, 30, 0.7)],
                            [(18, 20, 0.5), (24, 28, 0.5), (30, 20, 0.5)]]),
    "xps_survey":         ("peaks", "XPS survey", ["XPS", "photoelectron", "core level"],
                           [[(16, 20, 1.5), (30, 14, 1.8), (44, 26, 2.0)], [(20, 26, 1.4), (38, 12, 1.6)],
                            [(14, 12, 1.2), (28, 28, 1.7), (42, 18, 1.5), (50, 8, 1.4)], [(24, 22, 2.2)],
                            [(18, 16, 1.3), (32, 24, 1.6), (46, 14, 1.8)]]),
    "gc_chromatogram":    ("peaks", "GC chromatogram", ["GC", "retention time", "chromatography"],
                           [[(16, 10, 0.5), (26, 28, 0.6), (38, 16, 0.7), (48, 22, 0.8)],
                            [(20, 24, 0.5), (30, 12, 0.6), (44, 26, 0.7)],
                            [(14, 18, 0.5), (24, 8, 0.5), (34, 30, 0.6), (42, 10, 0.6), (52, 14, 0.7)],
                            [(22, 14, 0.5), (36, 28, 0.9)], [(18, 26, 0.6), (30, 18, 0.6), (46, 12, 0.7)]]),
    "dls_distribution":   ("peaks", "DLS size distribution", ["DLS", "hydrodynamic size", "PDI"],
                           [[(30, 30, 5.0)], [(22, 20, 3.0), (42, 26, 4.0)], [(28, 14, 2.5)],
                            [(34, 24, 6.0), (20, 8, 2.0)], [(26, 28, 4.5)]]),
    "pore_size_distribution": ("peaks", "Pore size distribution", ["pore size", "BJH", "mesopore"],
                               [[(22, 26, 3.0)], [(16, 14, 2.5), (34, 28, 3.5)], [(28, 22, 5.0)],
                                [(20, 18, 2.0), (30, 12, 2.2), (42, 10, 2.5)], [(36, 26, 4.0)]]),
    "eds_spectrum":       ("peaks", "EDS spectrum", ["EDX", "X-ray energy", "elemental map"],
                           [[(14, 26, 0.5), (26, 14, 0.6), (40, 28, 0.6)], [(18, 20, 0.5), (32, 30, 0.7)],
                            [(12, 12, 0.5), (22, 24, 0.5), (36, 16, 0.6), (48, 30, 0.7)],
                            [(26, 30, 0.6), (38, 12, 0.5)], [(16, 16, 0.5), (30, 22, 0.8), (44, 10, 0.5)]]),
    "sec_chromatogram":   ("peaks", "SEC chromatogram", ["SEC", "GPC", "molar mass"],
                           [[(20, 26, 3.0), (38, 16, 2.5)], [(24, 14, 2.5), (42, 28, 3.5)],
                            [(16, 20, 2.0), (30, 26, 2.8), (46, 10, 2.2)], [(34, 30, 4.5)], [(22, 24, 2.6), (40, 12, 2.0)]]),
    "ms_spectrum":        ("sticks", "Mass spectrum", ["MS", "mass spectrum", "m/z"],
                           [[(10, 14), (18, 30), (26, 10), (34, 22), (42, 40), (50, 16)],
                            [(12, 24), (22, 8), (30, 36), (46, 18)],
                            [(14, 10), (20, 28), (28, 14), (38, 32), (54, 12)],
                            [(16, 20), (24, 36), (32, 14), (40, 26), (48, 10)],
                            [(10, 30), (26, 16), (36, 40), (44, 10), (52, 22)]]),
    "tga_steps":          ("steps", "TGA mass steps", ["TGA", "thermogravimetric", "mass loss"],
                           [[(24, 10), (42, 22)], [(20, 6), (34, 12), (48, 24)], [(30, 30)],
                            [(18, 8), (40, 18)], [(26, 14), (38, 10), (50, 16)]]),
    "dsc_thermogram":     ("peaks", "DSC thermogram", ["DSC", "heat flow", "endotherm", "exotherm"],
                           [[(24, 14, 2.5), (40, 22, 3.0)], [(30, 26, 2.0)], [(18, 12, 2.0), (34, 18, 2.2), (46, 10, 2.4)],
                            [(26, 24, 3.5)], [(20, 16, 2.6), (36, 26, 2.2)]]),
    "mossbauer_spectrum": ("dips", "Mössbauer spectrum", ["Mossbauer", "isomer shift", "hyperfine"],
                           [[(24, 22, 1.5), (40, 22, 1.5)], [(20, 16, 1.2), (32, 28, 1.4), (44, 12, 1.2)],
                            [(30, 30, 2.0)], [(18, 10, 1.6), (36, 26, 1.8)], [(26, 18, 1.3), (42, 14, 1.4)]]),
    "sorption_hysteresis": ("hyst", "Sorption hysteresis", ["isotherm", "hysteresis", "BET"],
                            [[(0.35,)], [(0.45,)], [(0.25,)], [(0.55,)], [(0.4,)]]),
    "decay_kinetics":     ("decay", "Decay kinetics", ["kinetics", "lifetime", "decay"],
                           [[(0.14,)], [(0.09,)], [(0.22,)], [(0.06,)], [(0.3,)]]),
}


def _trace_d(peaks, x0=11, x1=55, base=52, top=44, scale=1.0):
    pts = []
    n = 60
    for i in range(n + 1):
        x = x0 + (x1 - x0) * i / n
        s = 0.0
        for c, h, w in peaks:
            s += h * math.exp(-((x - c) / (w * scale)) ** 2)
        pts.append((x, base - min(top, s * top / 30.0)))
    return _poly_d(pts)


def _sticks_d(sticks, x0=11, x1=55, base=52, top=44):
    segs = []
    for c, h in sticks:
        x = x0 + (x1 - x0) * c / 60.0
        segs.append(f"M{x:.2f} {base:.2f} V{base - top * h / 40.0:.2f}")
    return " ".join(segs)


def _steps_d(steps, x0=11, x1=55, top=12, bottom=50):
    pts = [(x0, top)]
    level = top
    for x, drop in steps:
        pts.append((x0 + (x - 0) * (x1 - x0) / 60.0, level))
        level += drop
        pts.append((x0 + (x - 0) * (x1 - x0) / 60.0, level))
    pts.append((x1, level))
    return _poly_d(pts)


def _spectrum_fn(kind, data):
    def fn(ic):
        axes(ic, 3.0)
        if kind == "peaks" or kind == "dips":
            d = _trace_d(data if kind == "peaks" else [(c, -h, w) for c, h, w in data])
            if kind == "dips":
                d = _trace_d([(c, -h, w) for c, h, w in data], base=20, top=44)
            ic.spath(d, 2.2)
            ic.acc().rect(12, 48, 40, 3, 0)
        elif kind == "sticks":
            ic.spath(_sticks_d(data), 2.6)
        elif kind == "steps":
            ic.spath(_steps_d(data), 2.4)
        elif kind == "hyst":
            ads = [(12 + 42 * i / 30, 50 - 38 * (i / 30) ** 2) for i in range(31)]
            des = [(12 + 42 * i / 30, 50 - 38 * ((i / 30) ** 0.5) * 0.9 - 2) for i in range(31)]
            ic.spath(_poly_d(ads), 2.2)
            ic.spath(_poly_d(des), 2.0)
        elif kind == "decay":
            k = data[0][0]
            pts = [(12 + 42 * i / 40, 14 + 36 * (1 - math.exp(-k * 9 * (1 - i / 40)))) for i in range(41)]
            ic.spath(_poly_d(pts), 2.2)
    return fn


def _register_spectra():
    for fam, (kind, label, tags, variants) in SPECTRA.items():
        for k, data in enumerate(variants, 1):
            name = f"{fam}_{k}"
            _reg(name, "characterisation", f"{label} {k}", tags, _spectrum_fn(kind, data))


# ------------------------------------------------------------------ vessels ---------------------------------------
# silhouette: (neck half-width, body half-width, neck bottom y, body bottom y, corner radius)
VESSELS = {
    "tapered_tube":   (4.0, 6.5, 22, 56, 3),
    "pear_flask":     (3.5, 16.0, 30, 56, 14),
    "kjeldahl_flask": (3.0, 12.0, 26, 56, 10),
    "bottle_round":   (5.0, 14.0, 20, 56, 10),
    "tall_beaker":    (14.0, 14.0, 10, 56, 2),
    "ampoule":        (3.0, 9.0, 30, 56, 9),
    "crimp_vial":     (9.0, 12.0, 14, 56, 4),
    "microtube":      (6.0, 6.0, 26, 56, 3),
    "bell_jar":       (16.0, 22.0, 14, 56, 12),
    "cylinder_tube":  (9.0, 9.0, 8, 56, 3),
    "conical_bottom": (9.0, 16.0, 12, 56, 4),
    "wide_dish":      (18.0, 26.0, 30, 52, 6),
}
CONTENTS = ["empty", "half", "full", "crystals", "gel"]


def _vessel_fn(shape, content):
    nw, bw, neck_bottom, body_bottom, r = VESSELS[shape]
    cx = 32.0
    top_y = 7.0
    body_top = neck_bottom
    def fn(ic):
        ic.line(cx - nw - 2, top_y, cx + nw + 2, top_y, 3.2)             # rim
        left = [(cx - nw, top_y), (cx - nw, body_top)]
        right = [(cx + nw, top_y), (cx + nw, body_top)]
        bl = (cx - bw, body_top + 4)
        br = (cx + bw, body_top + 4)
        wall = (f"M{cx - nw:.2f} {top_y} V{body_top:.2f} "
                f"L{bl[0]:.2f} {bl[1]:.2f} V{body_bottom - r:.2f} Q{bl[0]:.2f} {body_bottom} {bl[0] + r:.2f} {body_bottom} "
                f"H{br[0] - r:.2f} Q{br[0]:.2f} {body_bottom} {br[0]:.2f} {body_bottom - r:.2f} V{br[1]:.2f} "
                f"L{cx + nw:.2f} {body_top:.2f} V{top_y}")
        ic.spath(wall, 3.2, glass=True)
        if content == "empty":
            return
        level = {"half": body_bottom - 12, "full": body_top + 8, "crystals": body_bottom - 10,
                 "gel": body_bottom - 14}[content]
        inner = bw - 2.5
        liquid = (f"M{cx - inner:.2f} {level:.2f} H{cx + inner:.2f} V{body_bottom - 3:.2f} "
                  f"H{cx - inner:.2f} Z")
        ic.sec().path(liquid)
        if content == "crystals":
            for i, (dx, dy) in enumerate([(-6, 6), (2, 4), (-2, 14), (6, 13), (-8, 22)]):
                ic.acc().rect(cx + dx, body_bottom - 12 + dy * 0.5, 3.5, 3.5, 0.8)
        if content == "gel":
            ic.acc().path(f"M{cx - inner + 1:.2f} {level + 2:.2f} Q{cx:.2f} {level - 2:.2f} {cx + inner - 1:.2f} {level + 2:.2f}")
    return fn


def _register_vessels():
    for shape in VESSELS:
        for content in CONTENTS:
            name = f"{shape}_{content}"
            label = f"{shape.replace('_', ' ').title()} ({content})"
            _reg(name, "glassware", label, [shape.replace("_", " "), content, "lab glass"], _vessel_fn(shape, content))


# ------------------------------------------------------------------ instruments -----------------------------------
INSTRUMENTS = {
    "bench_box":      (14, 10, 36, 28),
    "tall_cabinet":   (18, 4, 28, 48),
    "rotary_head":    (14, 18, 36, 18),
    "round_chamber":  (16, 12, 32, 32),
    "stirrer_plate":  (10, 36, 44, 14),
    "balance_pan":    (12, 36, 40, 10),
    "pump_head":      (16, 24, 28, 22),
    "dryer_drum":     (12, 14, 40, 36),
    "mixer_frame":    (20, 20, 24, 30),
    "tube_furnace":   (6, 22, 52, 20),
}
FEATURES = ["screen", "knob", "window", "dial", "lever", "vent"]


def _instrument_fn(body, feat):
    x, y, w, h = INSTRUMENTS[body]
    def fn(ic):
        if body == "round_chamber":
            ic.scircle(x + w / 2, y + h / 2, w / 2, 3.2)
        elif body == "rotary_head":
            ic.srect(x, y, w, h, 6, 3.2)
            ic.line(x + w / 2, y, x + w / 2, y - 6, 3.0)
        elif body == "tube_furnace":
            ic.srect(x, y, w, h, 8, 3.2)
            ic.line(x + 4, y + h / 2, x + w - 4, y + h / 2, 1.4)
        else:
            ic.srect(x, y, w, h, 3, 3.2)
        if feat == "screen":
            ic.cut().rect(x + w * 0.22, y + h * 0.18, w * 0.56, h * 0.26, 1.5)
            ic.acc().line(x + w * 0.28, y + h * 0.31, x + w * 0.7, y + h * 0.31, 1.4)
        elif feat == "knob":
            ic.acc().circle(x + w * 0.5, y + h * 0.66, min(w, h) * 0.12)
        elif feat == "window":
            ic.cut().rect(x + w * 0.2, y + h * 0.3, w * 0.6, h * 0.4, 2)
        elif feat == "dial":
            ic.acc().circle(x + w * 0.5, y + h * 0.5, min(w, h) * 0.2)
            ic.line(x + w * 0.5, y + h * 0.5, x + w * 0.5 + min(w, h) * 0.16, y + h * 0.5 - min(w, h) * 0.1, 1.6)
        elif feat == "lever":
            ic.acc().line(x + w * 0.2, y + h * 0.7, x + w * 0.8, y + h * 0.4, 2.2)
        elif feat == "vent":
            for i in range(3):
                ic.cut().rect(x + w * 0.22, y + h * (0.2 + i * 0.22), w * 0.56, 2, 1)
    return fn


def _register_instruments():
    for body in INSTRUMENTS:
        for feat in FEATURES:
            name = f"{body}_{feat}"
            _reg(name, "instruments", f"{body.replace('_', ' ').title()} with {feat}", [body.replace("_", " "), feat, "instrument"],
                 _instrument_fn(body, feat))


# ------------------------------------------------------------------ morphologies ----------------------------------
MORPH = ["octahedron", "truncated_cube", "hex_plate", "rod", "sheet", "sphere_cluster", "dendrite", "cage", "wire_bundle", "star"]
DECOR = ["core", "hollow", "pores", "facets", "aggregate", "defect"]


def _shape_pts(name):
    if name == "octahedron":
        return [(32, 12), (50, 32), (32, 52), (14, 32)]
    if name == "truncated_cube":
        return [(20, 12), (44, 12), (52, 20), (52, 44), (44, 52), (20, 52), (12, 44), (12, 20)]
    if name == "hex_plate":
        return [(32 + 20 * math.cos(math.radians(a)), 32 + 14 * math.sin(math.radians(a))) for a in range(0, 360, 60)]
    if name == "rod":
        return [(24, 10), (40, 10), (42, 54), (22, 54)]
    if name == "sheet":
        return [(12, 18), (50, 12), (54, 44), (16, 52)]
    if name == "star":
        pts = []
        for i in range(10):
            r = 22 if i % 2 == 0 else 9
            a = math.radians(-90 + i * 36)
            pts.append((32 + r * math.cos(a), 32 + r * math.sin(a)))
        return pts
    if name == "nanocube":
        return [(18, 18), (46, 18), (46, 46), (18, 46)]
    if name == "bipyramid":
        return [(32, 8), (46, 32), (32, 56), (18, 32)]
    if name == "triangular_plate":
        return [(32, 12), (52, 48), (12, 48)]
    if name == "nanoflower":
        return [(32 + (20 if i % 2 == 0 else 10) * math.cos(math.radians(i * 30)),
                 32 + (20 if i % 2 == 0 else 10) * math.sin(math.radians(i * 30))) for i in range(12)]
    if name in ("mesoporous_sphere", "janus_particle"):
        return [(32 + 20 * math.cos(math.radians(a)), 32 + 20 * math.sin(math.radians(a))) for a in range(0, 360, 20)]
    return None


def _morph_fn(shape, decor):
    def fn(ic):
        if shape in ("sphere_cluster", "dendrite", "cage", "wire_bundle"):
            if shape == "sphere_cluster":
                for cx, cy, r in [(22, 24, 9), (42, 22, 8), (32, 42, 10), (46, 44, 6)]:
                    ic.scircle(cx, cy, r, 2.6)
            elif shape == "dendrite":
                ic.spath("M32 52 V30 M32 36 L20 26 M32 36 L44 26 M32 44 L24 40 M32 44 L40 40 M32 30 L30 14 M32 30 L36 14", 2.8)
            elif shape == "cage":
                ic.srect(14, 14, 36, 36, 6, 3.0)
                ic.spath("M14 32 H50 M32 14 V50", 2.0)
            else:
                for dx in (-8, 0, 8):
                    ic.spath(f"M{32 + dx} 10 Q{28 + dx} 32 {32 + dx} 54", 2.4)
        else:
            pts = _shape_pts(shape)
            if pts:
                ic.spath(_poly_d(pts) + " Z", 3.0)
        if decor == "core":
            ic.acc().circle(32, 32, 7)
        elif decor == "hollow":
            ic.cut().circle(32, 32, 9)
        elif decor == "pores":
            for px, py in [(24, 26), (36, 24), (30, 38), (42, 40)]:
                ic.cut().circle(px, py, 2.4)
        elif decor == "facets":
            ic.spath("M32 12 L32 52 M14 32 L50 32", 1.4)
        elif decor == "aggregate":
            for px, py in [(20, 40), (30, 30), (40, 42), (46, 26)]:
                ic.acc().circle(px, py, 3.2)
        elif decor == "defect":
            ic.acc().circle(38, 30, 4)
            ic.cut().circle(38, 30, 1.8)
    return fn


def _register_morph():
    for shape in MORPH:
        for decor in DECOR:
            name = f"{shape}_{decor}"
            _reg(name, "materials", f"{shape.replace('_', ' ').title()} ({decor})", [shape.replace("_", " "), decor, "morphology"],
                 _morph_fn(shape, decor))


# ------------------------------------------------------------------ process arrows --------------------------------
ARROW_STYLES = ["straight", "curved", "elbow", "double", "step", "recycle"]
PROPORTIONS = [0.7, 0.85, 1.0, 1.15, 1.3]


def _arrow_fn(style, prop):
    L = 40 * prop
    x0 = 32 - L / 2
    x1 = 32 + L / 2
    def fn(ic):
        if style == "straight":
            ic.spath(f"M{x0:.2f} 32 H{x1 - 6:.2f}", 3.6)
            ic.acc().path(f"M{x1 - 9:.2f} 22 L{x1:.2f} 32 L{x1 - 9:.2f} 42 Z")
        elif style == "curved":
            ic.spath(f"M{x0:.2f} 44 Q32 8 {x1 - 6:.2f} 30", 3.6)
            ic.acc().path(f"M{x1 - 10:.2f} 24 L{x1 - 3:.2f} 30 L{x1 - 9:.2f} 36 Z")
        elif style == "elbow":
            ic.spath(f"M{x0:.2f} 42 H32 V22 H{x1 - 6:.2f}", 3.6)
            ic.acc().path(f"M{x1 - 9:.2f} 14 L{x1:.2f} 22 L{x1 - 9:.2f} 30 Z")
        elif style == "double":
            ic.spath(f"M{x0 + 6:.2f} 32 H{x1 - 6:.2f}", 3.4)
            ic.acc().path(f"M{x1 - 9:.2f} 22 L{x1:.2f} 32 L{x1 - 9:.2f} 42 Z")
            ic.acc().path(f"M{x0 + 9:.2f} 22 L{x0:.2f} 32 L{x0 + 9:.2f} 42 Z")
        elif style == "step":
            ic.spath(f"M{x0:.2f} 44 H{32 - 3:.2f} V32 H{x1 - 6:.2f}", 3.6)
            ic.acc().path(f"M{x1 - 9:.2f} 24 L{x1:.2f} 32 L{x1 - 9:.2f} 40 Z")
        elif style == "recycle":
            ic.spath(f"M{x0:.2f} 40 Q{x0:.2f} 16 32 16 Q{x1:.2f} 16 {x1:.2f} 40", 3.4)
            ic.acc().path(f"M{x1 - 5:.2f} 34 L{x1:.2f} 42 L{x1 + 4:.2f} 34 Z")
    return fn


def _register_arrows():
    for style in ARROW_STYLES:
        for p in (PROPORTIONS if style in ("straight", "curved", "elbow", "double", "step", "recycle") else PROPORTIONS_EXTRA):
            name = f"process_{style}_{int(p * 100)}"
            _reg(name, "process", f"Process arrow {style} {int(p * 100)} %", [style, "arrow", "process flow"], _arrow_fn(style, p))



# ---- second batch (200 more, brings the total to 500 new icons = 3000 shapes with variants) ----------------------
import random as _random

_EXTRA_SPECTRA = {
    "epr_spectrum":        ("peaks", "EPR spectrum", ["EPR", "g-factor", "radical"]),
    "cv_voltammogram":     ("peaks", "Cyclic voltammogram", ["CV", "redox", "voltammetry"]),
    "saxs_curve":          ("decay", "SAXS intensity", ["SAXS", "scattering", "q-range"]),
    "tauc_plot":           ("decay", "Tauc plot", ["Tauc", "band gap", "UV-Vis"]),
    "impedance_bode":      ("decay", "Bode impedance", ["EIS", "impedance", "Bode"]),
    "gcms_tic":            ("peaks", "GC-MS total ion current", ["GC-MS", "TIC", "mass spectrometry"]),
    "ftir_difference":     ("dips", "FTIR difference", ["FTIR", "difference spectrum", "adsorption"]),
    "xrd_difference":      ("peaks", "XRD difference", ["PXRD", "difference", "phase change"]),
    "sans_scattering":     ("decay", "SANS scattering", ["SANS", "neutron", "scattering"]),
    "pulse_response":      ("decay", "Pulse response", ["transient", "breakthrough", "response"]),
    "thermal_conductivity": ("steps", "Thermal conductivity", ["thermal", "conductivity", "transient"]),
    "nmr_relaxation":      ("decay", "NMR relaxation", ["T1", "T2", "relaxation"]),
}
for _fam, (_kind, _label, _tags) in _EXTRA_SPECTRA.items():
    _rng = _random.Random(_fam)
    _vs = []
    for _k in range(5):
        if _kind == "peaks" or _kind == "dips":
            _vs.append([(_rng.uniform(10, 54), _rng.uniform(6, 30), _rng.uniform(0.6, 1.6)) for _ in range(_rng.randint(2, 5))])
        elif _kind == "decay":
            _vs.append([(_rng.uniform(0.05, 0.3),)])
        else:
            _vs.append([(_rng.uniform(14, 48), _rng.uniform(4, 12)) for _ in range(_rng.randint(1, 3))])
    SPECTRA[_fam] = (_kind, _label, _tags, _vs)

VESSELS.update({
    "florence_flask":  (4.5, 17.0, 26, 56, 12),
    "schlenk_tube":    (3.5, 5.5, 18, 56, 3),
    "gas_washing":     (3.0, 11.0, 22, 56, 8),
    "dropping_funnel": (3.0, 12.0, 14, 56, 6),
    "separatory_flask": (3.5, 15.0, 18, 56, 10),
    "burette_tube":    (2.5, 3.0, 10, 56, 2),
    "screwcap_vial":   (8.0, 11.0, 18, 56, 4),
    "square_bottle":   (6.0, 13.0, 18, 56, 2),
    "pipette_tube":    (2.0, 2.5, 8, 56, 2),
    "reactor_jacket":  (14.0, 14.0, 12, 56, 3),
    "cuvette_square":  (13.0, 13.0, 28, 56, 1),
    "round_test_tube": (7.0, 7.0, 10, 56, 6),
})
INSTRUMENTS.update({
    "glovebox":        (8, 14, 48, 36),
    "fume_hood":       (6, 8, 52, 44),
    "sieve_shaker":    (14, 20, 36, 16),
    "ultrasonic_bath": (10, 30, 44, 18),
    "spin_coater":     (16, 18, 32, 20),
    "thermocycler":    (14, 16, 36, 30),
})
MORPH.extend(["nanocube", "bipyramid", "triangular_plate", "nanoflower", "mesoporous_sphere", "janus_particle"])
ARROW_STYLES.extend(["tee", "split"])
PROPORTIONS_EXTRA = [0.8, 1.0, 1.2, 1.4]

_register_spectra()
_register_vessels()
_register_instruments()
_register_morph()
_register_arrows()
