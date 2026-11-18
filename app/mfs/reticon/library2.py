"""
Extended set: the gaps no free icon library covers for materials science,
heterogeneous catalysis, MOFs and proteomics.
"""
from __future__ import annotations

import math

from .core import polar, smooth_path, lorentz_trace, _f, Op
from .library import icon, axes


def gauss_trace(x0, x1, base, peaks, n=160, sign=-1):
    pts = []
    for i in range(n + 1):
        x = x0 + (x1 - x0) * i / n
        y = base
        for c, h, w in peaks:
            y += sign * h * math.exp(-((x - c) ** 2) / (2 * w * w))
        pts.append((x, y))
    return pts


def rot_rect(ic, cx, cy, w, h, r, deg, role="m"):
    ic.ops.append(Op(role, "rect", {"x": cx - w / 2, "y": cy - h / 2, "width": w, "height": h, "rx": r,
                                    "transform": f"rotate({_f(deg)} {_f(cx)} {_f(cy)})"}))


# ==========================================================================
# DATA / PLOT ICONS
# ==========================================================================

@icon("isotherm_type_iv", "data", label="Isotherm (Type IV)", tags=["BET", "mesoporous", "hysteresis", "N2 sorption"])
def _(ic):
    axes(ic)
    ic.acc().spath("M12 50 C13 34 17 32 28 31 C38 30 40 27 44 18 C47 13 51 12 56 11", 3)
    ic.sec().spath("M56 11 C50 12 49 14 47.5 22 C46.5 27 44 29 36 30", 2.4)


@icon("pore_size_distribution", "data", label="Pore Size Distribution", tags=["PSD", "NLDFT", "BJH", "porosity"])
def _(ic):
    axes(ic)
    for x, h in ((14, 4), (19, 12), (24, 36), (29, 22), (34, 8), (39, 4), (44, 10), (49, 17), (54, 6)):
        ic.acc().rect(x - 2, 55 - h, 4, h - 1.75, 0.8)


@icon("tauc_plot", "data", label="Tauc Plot", tags=["band gap", "UV-Vis DRS", "optical"])
def _(ic):
    axes(ic)
    ic.acc().spath("M12 51 C24 51 30 49 34 42 C38 33 42 20 47 10", 3)
    ic.sec().line(26, 52, 47.5, 9, 2)
    ic.circle(26, 52, 2.6)


@icon("dls_distribution", "data", label="DLS Size", tags=["DLS", "hydrodynamic", "particle size", "PDI"])
def _(ic):
    axes(ic)
    ic.acc().spoly(gauss_trace(12, 56, 52, [(32, 38, 5)], n=120), 3)
    ic.sec().spoly(gauss_trace(12, 56, 52, [(44, 14, 3.5)], n=120), 2.2)


@icon("cyclic_voltammetry", "data", label="Cyclic Voltammetry", tags=["CV", "electrochemistry", "redox"])
def _(ic):
    ic.line(8, 32, 58, 32, 2.2)
    ic.line(32, 6, 32, 58, 2.2)
    ic.acc().spath("M10 38 C22 38 26 37 31 30 C34 26 36 14 40 14 C44 14 45 22 55 22 "
                   "C50 28 44 30 40 34 C37 37 36 50 31 50 C26 50 25 40 10 38 Z", 3)


@icon("dsc", "data", label="DSC", tags=["calorimetry", "melting", "glass transition", "thermal"])
def _(ic):
    axes(ic)
    ic.acc().spoly(gauss_trace(12, 56, 30, [(26, -14, 2.2), (44, 10, 2.6)], n=160, sign=-1), 3)


@icon("michaelis_menten", "data", label="Enzyme Kinetics", tags=["Michaelis-Menten", "Km", "Vmax", "nanozyme", "kinetics"])
def _(ic):
    axes(ic)
    ic.sec().spath("M12 16 H57", 2)
    ic.acc().spath("M12 54 C16 30 26 20 56 17", 3)
    for x, y in ((15, 40), (20, 29), (28, 23), (38, 19.5), (50, 18)):
        ic.circle(x, y, 2.3)


@icon("arrhenius_plot", "data", label="Arrhenius Plot", tags=["activation energy", "Ea", "kinetics"])
def _(ic):
    axes(ic)
    ic.acc().line(13, 13, 55, 48, 3)
    for x in (16, 24, 32, 40, 48):
        ic.circle(x, 13 + (x - 13) * 35 / 42, 2.6)


@icon("calibration_curve", "data", label="Calibration Curve", tags=["quantification", "linear fit", "standard curve"])
def _(ic):
    axes(ic)
    ic.sec().line(12, 52, 56, 12, 2.4)
    for x, dy in ((16, 1), (24, -1.5), (32, 1), (40, -0.8), (48, 1.4)):
        ic.circle(x, 52 - (x - 12) * 40 / 44 + dy, 2.6)
    ic.cut().circle(36, 30.5, 4.2)
    ic.acc().circle(36, 30.5, 2.6)


@icon("xps_spectrum", "data", label="XPS Peak Fit", tags=["XPS", "deconvolution", "Zr 3d", "Hf 4f"])
def _(ic):
    axes(ic)
    a = gauss_trace(12, 56, 50, [(27, 30, 3)], n=140)
    b = gauss_trace(12, 56, 50, [(38, 20, 3)], n=140)
    env = gauss_trace(12, 56, 50, [(27, 30, 3), (38, 20, 3)], n=140)
    ic.sec().spoly(a, 1.8)
    ic.sec().spoly(b, 1.8)
    ic.acc().spoly(env, 3)


@icon("edx_spectrum", "data", label="EDX Spectrum", tags=["EDS", "EDX", "elemental analysis"])
def _(ic):
    axes(ic)
    ic.acc().spoly(lorentz_trace(12, 56, 50, [(16, 26, .5), (22, 12, .5), (30, 36, .5), (35, 8, .5),
                                              (44, 18, .5), (50, 6, .5)], n=240), 2.4)


@icon("conversion_time", "data", label="Conversion vs Time", tags=["kinetics", "yield", "hydrolysis", "catalysis"])
def _(ic):
    axes(ic)
    ic.acc().spath("M12 52 C20 30 30 18 56 14", 3)
    ic.sec().spath("M12 52 C24 46 36 40 56 32", 2.6)
    for x, y in ((18, 36), (26, 25), (36, 19), (48, 15.5)):
        ic.circle(x, y, 2.2)


@icon("uv_vis_spectrum", "data", label="Absorbance Spectrum", tags=["UV-Vis", "absorbance", "colorimetric"])
def _(ic):
    axes(ic)
    ic.acc().spoly(gauss_trace(12, 56, 50, [(16, 34, 4), (34, 22, 6)], n=160), 3)


@icon("pl_spectrum", "data", label="Photoluminescence", tags=["PL", "fluorescence", "emission"])
def _(ic):
    axes(ic)
    ic.acc().spoly(gauss_trace(12, 56, 52, [(34, 40, 4.5)], n=140), 3)
    ic.sec().spoly(gauss_trace(12, 56, 52, [(34, 20, 4.5)], n=140), 2.2)


@icon("zeta_potential", "data", label="Zeta Potential", tags=["surface charge", "colloid", "stability"])
def _(ic):
    ic.line(8, 55, 58, 55, 3.5)
    ic.line(33, 8, 33, 55, 1.8)
    ic.acc().spoly(gauss_trace(10, 56, 54, [(22, 36, 4)], n=140), 3)
    ic.sec().rect(43, 11, 10, 3, 1.5)
    ic.sec().rect(46.5, 7.5, 3, 10, 1.5)
    ic.sec().rect(13, 11, 10, 3, 1.5)


@icon("sds_page", "proteomics", label="SDS-PAGE", tags=["gel", "electrophoresis", "protein bands", "western"])
def _(ic):
    ic.rect(9, 6, 46, 52, 4)
    for x in (14, 25, 36, 47):
        ic.cut().rect(x - 2, 9, 9, 4, 1)
    bands = {0: (18, 24, 31, 40, 47), 1: (22, 33, 47), 2: (26, 38, 44, 50), 3: (20, 30, 42)}
    for lane, ys in bands.items():
        x = 12 + lane * 11
        for y in ys:
            ic.cut().rect(x, y, 9, 2.2 if lane else 1.6, 1)


@icon("chromatogram_fragments", "proteomics", label="Peptide Map", tags=["peptide mapping", "LC-MS", "digest", "TIC"])
def _(ic):
    axes(ic)
    ic.acc().spoly(lorentz_trace(12, 56, 52, [(15, 10, .7), (19, 24, .6), (23, 14, .6), (27, 34, .6),
                                              (31, 9, .6), (35, 22, .6), (39, 30, .6), (44, 12, .6),
                                              (48, 18, .6), (53, 8, .7)], n=260), 2.2)


# ==========================================================================
# TECHNIQUE SCHEMATICS
# ==========================================================================

@icon("bragg_diffraction", "characterisation", label="Bragg's Law", tags=["XRD", "lattice planes", "d-spacing"])
def _(ic):
    for y in (38, 48, 58):
        for x in range(8, 60, 8):
            ic.circle(x, y, 2.6)
    ic.acc().line(6, 14, 28, 38, 2.6)
    ic.acc().line(28, 38, 50, 14, 2.6)
    ic.sec().line(14, 14, 36, 48, 2.2)
    ic.sec().line(36, 48, 58, 14, 2.2)


@icon("afm", "characterisation", label="AFM", tags=["atomic force microscopy", "topography", "roughness"])
def _(ic):
    ic.rect(4, 12, 30, 6, 2)
    ic.poly([(30, 15), (40, 15), (37, 32), (35, 32)])
    ic.acc().line(46, 6, 37, 15, 2.2)
    ic.acc().line(37, 15, 58, 10, 2.2)
    ic.acc().circle(46, 6, 2.6)
    ic.sec().path("M4 58 V46 C8 46 9 40 13 40 C17 40 17 46 21 46 C25 46 26 36 31 36 C36 36 36 44 40 44 "
                  "C44 44 45 39 49 39 C53 39 54 46 60 46 V58 Z")


@icon("contact_angle", "characterisation", label="Contact Angle", tags=["wettability", "hydrophobic", "surface"])
def _(ic):
    ic.rect(4, 44, 56, 8, 2)
    ic.sec().path("M14 44 C14 22 50 22 50 44 Z")
    ic.acc().spath("M25 44 A11 11 0 0 0 14.6 33", 2.6)


@icon("icp", "characterisation", label="ICP-OES / MS", tags=["ICP", "elemental", "metal leaching", "digest"])
def _(ic):
    ic.rect(22, 40, 20, 18, 3)
    for y in (31, 35):
        ic.rect(21, y, 22, 2.5, 1.2)
    ic.acc().path("M32 4 C40 12 42 20 40 26 C38 30 26 30 24 26 C22 20 24 12 32 4 Z")
    ic.cut().path("M32 13 C36 18 36 22 34.5 24.5 C33.5 26 30.5 26 29.5 24.5 C28 22 28 18 32 13 Z")
    ic.cut().line(32, 43, 32, 54, 2)


@icon("dls_instrument", "characterisation", label="Light Scattering", tags=["DLS", "laser", "scattering"])
def _(ic):
    ic.rect(3, 26, 14, 12, 2.5)
    ic.acc().line(17, 32, 27, 32, 3)
    ic.srect(27.75, 20.75, 12.5, 22.5, 2, 3, glass=True)
    ic.sec().rect(31, 28, 6, 12, 1)
    for a in (-40, -15, 15, 40):
        x2, y2 = polar(34, 32, 22, a)
        x1, y1 = polar(34, 32, 11, a)
        ic.acc().line(x1, y1, x2, y2, 2)
    ic.rect(52, 6, 9, 10, 2)


@icon("esi_source", "proteomics", label="ESI Source", tags=["electrospray", "ionisation", "mass spectrometry"])
def _(ic):
    ic.rect(2, 28, 18, 8, 2)
    ic.poly([(20, 28), (28, 31), (28, 33), (20, 36)])
    ic.acc().poly([(28, 32), (46, 22), (46, 42)])
    ic.cut().poly([(32, 32), (46, 25), (46, 39)])
    for x, y in ((36, 30), (39, 34), (42, 28.5), (43, 36), (40, 31)):
        ic.sec().circle(x, y, 1.4)
    ic.path("M50 18 H62 V46 H50 V36 H54 V28 H50 Z")


@icon("maldi_plate", "proteomics", label="MALDI Target", tags=["MALDI", "target plate", "spots"])
def _(ic):
    ic.rect(4, 12, 56, 40, 4)
    for r in range(4):
        for c in range(6):
            ic.cut().circle(12 + c * 8, 20 + r * 8, 2.8)
    ic.sec().circle(20, 28, 2.8)
    ic.sec().circle(36, 36, 2.8)
    ic.sec().circle(44, 20, 2.8)


@icon("lc_column", "proteomics", label="LC Column", tags=["HPLC", "C18", "column", "chromatography"])
def _(ic):
    rot_rect(ic, 32, 32, 12, 46, 3, 45)
    rot_rect(ic, 14, 14, 10, 8, 2, 45)
    rot_rect(ic, 50, 50, 10, 8, 2, 45)
    for k in (-10, 0, 10):
        cx, cy = 32 + k * .707, 32 + k * .707
        ic.ops.append(Op("x", "line", {"x1": cx - 3.2, "y1": cy + 3.2, "x2": cx + 3.2, "y2": cy - 3.2}, 1.8))
    ic.acc().line(6, 6, 10.5, 10.5, 2.4)
    ic.acc().line(53.5, 53.5, 58, 58, 2.4)


# ==========================================================================
# SYNTHESIS & PREP EQUIPMENT
# ==========================================================================

@icon("reflux_setup", "synthesis", label="Reflux", tags=["condenser", "reflux", "heating"])
def _(ic):
    ic.srect(28, 4, 8, 30, 2, 3, glass=True)
    ic.sec().line(32, 8, 32, 30, 2.6)
    ic.line(36, 10, 42, 10, 2.6)
    ic.line(28, 27, 22, 27, 2.6)
    ic.spath("M28.5 34 V40 A11 11 0 1 0 35.5 40 V34", 3, glass=True)
    ic.sec().path("M23.5 49 A8 8 0 0 0 40.5 49 Z")
    ic.acc().path("M18 62 H46 L42 58 H22 Z")


@icon("rotary_evaporator", "synthesis", label="Rotary Evaporator", tags=["rotavap", "solvent removal", "evaporation"])
def _(ic):
    # condenser (vertical coil) on the right
    ic.srect(42, 4, 12, 28, 3, 3, glass=True)
    ic.spath("M45 9 l6 3 l-6 3 l6 3 l-6 3 l6 3", 2)
    # motor head + tilted vapour duct
    ic.rect(36, 30, 18, 9, 2.5)
    rot_rect(ic, 30, 41, 16, 5, 2.5, 35)
    # receiving flask under condenser
    ic.line(50, 39, 50, 44, 3)
    ic.scircle(50, 50, 5.5, 3, glass=True)
    # evaporation flask in bath
    ic.scircle(19, 46, 8, 3, glass=True)
    ic.sec().path("M6 46 H36 V56 Q36 60 32 60 H10 Q6 60 6 56 Z")
    ic.cut().circle(19, 46, 9.75)
    ic.scircle(19, 46, 8, 3, glass=True)
    ic.sec().path("M12.6 48.5 A7 7 0 0 0 25.4 48.5 Z")


@icon("schlenk_line", "synthesis", label="Schlenk Line", tags=["air-free", "manifold", "vacuum", "inert gas"])
def _(ic):
    ic.rect(4, 10, 56, 7, 3.5)
    ic.rect(4, 20, 56, 7, 3.5)
    for x in (14, 28, 42):
        ic.rect(x - 3, 27, 6, 9, 1.5)
        ic.cut().circle(x, 31.5, 1.6)
        ic.line(x, 36, x, 42, 2.6)
    ic.spath("M39 42 V46 A7 7 0 1 0 45 46 V42", 2.6, glass=True)
    ic.sec().path("M37.2 51 A5 5 0 0 0 46.8 51 Z")
    ic.acc().circle(56, 13.5, 2)


@icon("ball_mill", "synthesis", label="Ball Mill", tags=["mechanochemistry", "milling", "grinding", "planetary"])
def _(ic):
    rot_rect(ic, 32, 32, 26, 44, 8, 30)
    ic.ops.append(Op("x", "rect", {"x": 23, "y": 14, "width": 18, "height": 36, "rx": 5,
                                   "transform": "rotate(30 32 32)"}))
    for x, y, r in ((28, 36, 4.2), (36, 30, 4.2), (33, 41, 3.2), (27, 26.5, 3), (38, 38.5, 2.6)):
        ic.sec().circle(x, y, r)
    ic.acc().spath("M8 18 A26 26 0 0 1 22 6", 2.6)
    ic.acc().poly([(20, 2.5), (26, 5), (21, 9.5)])


@icon("spin_coater", "synthesis", label="Spin Coater", tags=["thin film", "coating", "deposition"])
def _(ic):
    ic.rect(6, 40, 52, 18, 4)
    ic.cut().rect(12, 46, 14, 6, 1.5)
    ic.rect(29, 32, 6, 8)
    ic.sec().rect(16, 28, 32, 4, 2)
    ic.acc().spath("M10 22 A24 9 0 0 1 54 22", 2.4)
    ic.acc().poly([(53, 17), (58, 22), (51, 24)])
    ic.acc().path("M32 4 C35 9 35 12 32 12 C29 12 29 9 32 4 Z")


@icon("freeze_dryer", "synthesis", label="Freeze Dryer", tags=["lyophilisation", "lyophilizer", "drying"])
def _(ic):
    ic.rect(10, 38, 44, 20, 4)
    ic.cut().rect(16, 44, 12, 6, 1.5)
    ic.path("M14 38 V20 A18 18 0 0 1 50 20 V38 Z")
    ic.cut().path("M19 35 V21 A13 13 0 0 1 45 21 V35 Z")
    for x in (25, 32, 39):
        ic.sec().rect(x - 2.5, 24, 5, 9, 1.5)
    ic.acc().circle(44, 47, 2.2)


@icon("vacuum_desiccator", "synthesis", label="Desiccator", tags=["drying", "vacuum", "storage"])
def _(ic):
    ic.path("M8 26 Q32 4 56 26 Z")
    ic.rect(29, 6, 6, 8, 2)
    ic.spath("M10 30 V48 Q10 56 18 56 H46 Q54 56 54 48 V30", 3.5, glass=True)
    ic.line(8, 28, 56, 28, 3)
    ic.sec().rect(14, 38, 36, 3, 1.5)
    for x in (22, 32, 42):
        ic.sec().circle(x, 47, 3)


@icon("vacuum_pump", "synthesis", label="Vacuum Pump", tags=["vacuum", "pump", "oil pump"])
def _(ic):
    ic.rect(6, 22, 40, 28, 6)
    ic.cut().circle(26, 36, 9)
    ic.circle(26, 36, 5.5)
    ic.cut().circle(26, 36, 2)
    ic.rect(46, 30, 12, 6, 2)
    ic.rect(10, 50, 8, 6, 1.5)
    ic.rect(34, 50, 8, 6, 1.5)
    ic.sec().rect(14, 14, 6, 8, 1.5)
    ic.acc().spath("M52 26 q3 -4 0 -8 q-3 -4 0 -8", 2.4)


@icon("syringe_pump", "synthesis", label="Syringe Pump", tags=["dosing", "flow", "addition"])
def _(ic):
    ic.rect(4, 34, 56, 22, 4)
    ic.cut().rect(10, 40, 14, 8, 1.5)
    ic.srect(20, 22, 30, 10, 2, 2.6, glass=True)
    ic.sec().rect(36, 25, 12, 4, 1)
    ic.rect(8, 25, 12, 4, 1)
    ic.rect(6, 20, 4, 14, 1.5)
    ic.line(50, 27, 60, 27, 2.4)
    ic.acc().circle(48, 44, 3)


@icon("microwave_reactor", "synthesis", label="Microwave Reactor", tags=["microwave synthesis", "MW", "rapid heating"])
def _(ic):
    ic.rect(4, 12, 56, 40, 5)
    ic.cut().rect(10, 18, 32, 28, 3)
    for k in range(3):
        ic.acc().spath(f"M{14 + k * 9} 22 q3 4 0 8 t0 8 t0 8", 2.2)
    ic.cut().rect(47, 20, 8, 5, 1.2)
    ic.cut().circle(51, 36, 3.5)
    ic.rect(10, 52, 6, 4, 1)
    ic.rect(48, 52, 6, 4, 1)


@icon("crucible", "synthesis", label="Crucible", tags=["calcination", "TGA pan", "alumina", "ashing"])
def _(ic):
    ic.path("M12 22 H52 L46 54 Q45 58 41 58 H23 Q19 58 18 54 Z")
    ic.cut().path("M18 27 H46 L45 31 H19 Z")
    ic.sec().path("M21 40 H43 L41.6 50 H22.4 Z")
    ic.rect(8, 16, 48, 5, 2.5)
    for x in (22, 32, 42):
        ic.acc().spath(f"M{x} 12 q-3 -3 0 -5 q3 -2 0 -5", 2.2)


@icon("pellet_press", "synthesis", label="Pellet Press", tags=["KBr pellet", "hydraulic press", "die"])
def _(ic):
    ic.rect(8, 4, 48, 8, 2)
    ic.rect(8, 52, 48, 8, 2)
    ic.rect(10, 12, 5, 40)
    ic.rect(49, 12, 5, 40)
    ic.rect(28, 12, 8, 14)
    ic.rect(22, 26, 20, 4, 1.5)
    ic.sec().rect(24, 31, 16, 3, 1.5)
    ic.rect(20, 35, 24, 17, 2)
    ic.acc().poly([(32, 14), (28, 20), (36, 20)])


@icon("sieve", "synthesis", label="Sieve", tags=["sieving", "particle size", "mesh"])
def _(ic):
    ic.path("M6 18 H58 V26 Q58 30 54 30 H10 Q6 30 6 26 Z")
    ic.ellipse(32, 18, 26, 5)
    ic.cut().ellipse(32, 18, 22, 3.4)
    for x in range(12, 54, 6):
        ic.cut().line(x, 22, x, 28, 1.4)
    for x, y in ((20, 38), (32, 42), (44, 37), (26, 48), (39, 51), (31, 57)):
        ic.sec().circle(x, y, 2.2)


@icon("oil_bath", "synthesis", label="Oil Bath", tags=["heating", "bath", "temperature control"])
def _(ic):
    ic.spath("M10 22 V48 Q10 54 16 54 H48 Q54 54 54 48 V22", 3.5, glass=True)
    ic.sec().rect(13.5, 30, 37, 20.5, 2.5)
    ic.cut().circle(32, 34, 10)
    ic.spath("M29 6 V24 A8 8 0 1 0 35 24 V6", 2.8, glass=True)
    ic.rect(6, 56, 52, 4, 2)
    ic.acc().line(46, 4, 46, 36, 2)


@icon("ice_bath", "synthesis", label="Ice Bath", tags=["cooling", "0 C", "ice"])
def _(ic):
    ic.spath("M10 22 V48 Q10 54 16 54 H48 Q54 54 54 48 V22", 3.5, glass=True)
    for x, y, a in ((18, 32, 15), (44, 31, -12), (20, 44, -20), (44, 45, 10), (32, 48, 0)):
        rot_rect(ic, x, y, 9, 9, 2, a, "s")
    ic.cut().circle(32, 32, 10)
    ic.spath("M29 6 V24 A8 8 0 1 0 35 24 V6", 2.8, glass=True)


@icon("dialysis", "proteomics", label="Dialysis", tags=["desalting", "membrane", "buffer exchange"])
def _(ic):
    ic.spath("M10 14 V50 Q10 56 16 56 H48 Q54 56 54 50 V14", 3.5, glass=True)
    ic.sec().rect(13.5, 24, 37, 28.5, 2.5)
    ic.cut().rect(23, 4, 18, 46, 6)
    ic.rect(24, 4, 16, 5, 1.5)
    ic.srect(26, 10, 12, 36, 5, 2.6, glass=True)
    for x, y in ((30, 22), (34, 30), (30, 38)):
        ic.acc().circle(x, y, 2.4)
    ic.sec().circle(18, 30, 1.6)
    ic.sec().circle(46, 40, 1.6)


@icon("centrifugal_filter", "proteomics", label="Spin Filter", tags=["Amicon", "ultrafiltration", "MWCO", "concentrator"])
def _(ic):
    ic.rect(20, 3, 24, 6, 2)
    ic.spath("M22.5 12 V44 L29 60 H35 L41.5 44 V12", 3.5, glass=True)
    ic.srect(26, 12, 12, 26, 2, 2.4, glass=True)
    ic.sec().rect(28.5, 26, 7, 10, 1)
    ic.acc().line(26, 40, 38, 40, 2.4)
    ic.sec().path("M28 48 H36 L33.6 55 H30.4 Z")


@icon("spe_cartridge", "proteomics", label="SPE / ZipTip", tags=["solid-phase extraction", "desalting", "C18"])
def _(ic):
    ic.rect(22, 3, 20, 5, 1.5)
    ic.spath("M24 10 V38 L30 58 H34 L40 38 V10", 3.5, glass=True)
    ic.sec().rect(27.5, 30, 9, 8, 1)
    ic.acc().path("M32 44 C34 47 34 49 32 49 C30 49 30 47 32 44 Z")


@icon("thermomixer", "proteomics", label="Thermomixer", tags=["incubation", "shaking", "digestion", "37 C"])
def _(ic):
    ic.rect(6, 34, 52, 22, 4)
    ic.cut().rect(12, 42, 16, 8, 1.5)
    ic.cut().circle(46, 46, 4)
    ic.sec().rect(8, 28, 48, 6, 2)
    for x in (16, 26, 36, 46):
        ic.path(f"M{x - 4} 10 H{x + 4} V22 L{x + 1} 27 H{x - 1} L{x - 4} 22 Z")
    ic.acc().spath("M2 18 q2 -3 0 -6", 2)
    ic.acc().spath("M62 18 q-2 -3 0 -6", 2)


# ==========================================================================
# MATERIALS & STRUCTURE
# ==========================================================================

def _cube(ic, nodes_extra=()):
    f = [(12, 24), (40, 24), (40, 52), (12, 52)]
    b = [(x + 12, y - 12) for x, y in f]
    for i in range(4):
        for p, q in ((f[i], f[(i + 1) % 4]), (b[i], b[(i + 1) % 4]), (f[i], b[i])):
            ic.sec().line(*p, *q, 2.2)
    return f, b


@icon("unit_cell_fcc", "materials", label="FCC Unit Cell", tags=["fcc", "crystal structure", "lattice"])
def _(ic):
    f, b = _cube(ic)
    faces = [((f[0][0] + f[2][0]) / 2, (f[0][1] + f[2][1]) / 2),
             ((b[0][0] + b[2][0]) / 2, (b[0][1] + b[2][1]) / 2),
             ((f[0][0] + b[1][0]) / 2, (f[0][1] + b[1][1]) / 2),
             ((f[3][0] + b[2][0]) / 2, (f[3][1] + b[2][1]) / 2),
             ((f[0][0] + b[3][0]) / 2, (f[0][1] + b[3][1]) / 2),
             ((f[1][0] + b[2][0]) / 2, (f[1][1] + b[2][1]) / 2)]
    for p in f + b:
        ic.circle(*p, 4.2)
    for p in faces:
        ic.acc().circle(*p, 3.6)


@icon("unit_cell_bcc", "materials", label="BCC Unit Cell", tags=["bcc", "crystal structure", "lattice"])
def _(ic):
    f, b = _cube(ic)
    for p in f + b:
        ic.circle(*p, 4.2)
    ic.acc().circle(32, 32, 5.5)


@icon("graphene_lattice", "materials", label="Hexagonal Lattice", tags=["graphene", "2D material", "honeycomb", "carbon"])
def _(ic):
    r = 8.5
    centres = [(32, 32), (32 - 1.732 * r, 32), (32 + 1.732 * r, 32),
               (32 - .866 * r, 32 - 1.5 * r), (32 + .866 * r, 32 - 1.5 * r),
               (32 - .866 * r, 32 + 1.5 * r), (32 + .866 * r, 32 + 1.5 * r)]
    pts = set()
    for cx, cy in centres:
        hp = [polar(cx, cy, r, 30 + 60 * k) for k in range(6)]
        ic.sec().spoly(hp + [hp[0]], 2.2)
        pts.update((round(x, 1), round(y, 1)) for x, y in hp)
    for x, y in pts:
        ic.circle(x, y, 2.4)


@icon("layered_material", "materials", label="Layered / 2D Stack", tags=["MXene", "LDH", "2D", "exfoliation", "nanosheet"])
def _(ic):
    for k, y in enumerate((14, 28, 42)):
        ic.path(f"M6 {y + 8} L22 {y} H58 L42 {y + 8} Z", "s" if k == 1 else "m")
    ic.acc().spath("M52 26 v6", 2.2)
    ic.acc().spath("M52 40 v6", 2.2)


@icon("zeolite_cage", "materials", label="Zeolite / Cage", tags=["zeolite", "sodalite", "cage", "microporous"])
def _(ic):
    outer = [polar(32, 32, 26, 22.5 + 45 * k) for k in range(8)]
    inner = [polar(32, 32, 13, 45 * k) for k in range(8)]
    ic.sec().spoly(outer + [outer[0]], 2.2)
    ic.sec().spoly(inner + [inner[0]], 2.2)
    for i in range(8):
        ic.sec().line(*outer[i], *inner[i], 2.2)
        ic.sec().line(*outer[i], *inner[(i + 1) % 8], 2.2)
    for p in outer + inner:
        ic.circle(*p, 2.8)


@icon("nanorod", "materials", label="Nanorods / Wires", tags=["nanorod", "nanowire", "1D", "anisotropic"])
def _(ic):
    for (cx, cy, a, l) in ((22, 22, -35, 34), (40, 30, 20, 36), (26, 44, 60, 26), (46, 50, -10, 22)):
        rot_rect(ic, cx, cy, l, 8, 4, a, "m")
        rot_rect(ic, cx, cy, l - 8, 2, 1, a, "x")


@icon("nanosheet", "materials", label="Nanosheet", tags=["2D", "nanosheet", "monolayer", "exfoliated"])
def _(ic):
    ic.path("M4 36 C14 24 22 30 32 24 C42 18 50 20 60 14 L60 26 C50 32 42 30 32 36 C22 42 14 36 4 48 Z")
    ic.cut().spath("M8 40 C16 32 23 36 32 30.5 C41 25 49 26 56 21", 1.4)
    for x, y in ((14, 50), (30, 48), (46, 42), (54, 36)):
        ic.acc().circle(x, y, 2.6)


@icon("quantum_dot", "materials", label="Quantum Dots", tags=["QD", "luminescent", "nanocrystal"])
def _(ic):
    for cx, cy, r in ((20, 24, 9), (42, 22, 7), (30, 44, 8), (50, 44, 5.5)):
        ic.circle(cx, cy, r)
        ic.cut().circle(cx - r * .35, cy - r * .35, r * .25)
    for k in range(8):
        x1, y1 = polar(20, 24, 12, k * 45)
        x2, y2 = polar(20, 24, 15.5, k * 45)
        ic.acc().line(x1, y1, x2, y2, 2)


@icon("polymer_chain", "materials", label="Polymer", tags=["polymer", "chain", "macromolecule", "PEG"])
def _(ic):
    pts = [(6, 40), (12, 28), (20, 36), (24, 22), (32, 30), (36, 16), (44, 24), (48, 38), (40, 46), (48, 54), (58, 46)]
    ic.sec().spath(smooth_path(pts), 3)
    for i, p in enumerate(pts):
        ic.circle(*p, 3.6) if i % 2 == 0 else ic.acc().circle(*p, 3)


@icon("micelle", "materials", label="Micelle", tags=["surfactant", "self-assembly", "vesicle", "emulsion"])
def _(ic):
    for k in range(14):
        a = k * 360 / 14
        x1, y1 = polar(32, 32, 10, a)
        x2, y2 = polar(32, 32, 21, a)
        ic.sec().line(x1, y1, x2, y2, 2.2)
        hx, hy = polar(32, 32, 24.5, a)
        ic.circle(hx, hy, 4)
    ic.acc().circle(32, 32, 6)


@icon("dendrimer", "materials", label="Dendrimer", tags=["dendrimer", "branched", "PAMAM"])
def _(ic):
    ic.circle(32, 32, 5)
    for k in range(3):
        a = 90 + k * 120
        p1 = polar(32, 32, 13, a)
        ic.sec().line(32, 32, *p1, 2.4)
        ic.circle(*p1, 3.4)
        for d in (-35, 35):
            p2 = polar(*p1, 12, a + d)
            ic.sec().line(*p1, *p2, 2.2)
            ic.acc().circle(*p2, 3)


@icon("guest_in_pore", "materials", label="Host-Guest / Adsorption", tags=["adsorption", "guest molecule", "pore", "encapsulation", "CO2 capture"])
def _(ic):
    ic.rect(6, 6, 52, 52, 6)
    ic.cut().rect(14, 14, 36, 36, 4)
    for x, y in ((6, 6), (58, 6), (6, 58), (58, 58)):
        ic.acc().circle(x, y, 5)
    ic.sec().circle(26, 32, 5)
    ic.sec().circle(38, 32, 5)
    ic.sec().line(26, 32, 38, 32, 3)
    ic.ops.append(Op("x", "circle", {"cx": 32, "cy": 32, "r": 1.6}))


@icon("missing_linker_defect", "materials", label="Missing-Linker Defect", tags=["defect", "MOF", "open metal site", "Lewis acid"])
def _(ic):
    nodes = [(12, 12), (32, 12), (52, 12), (12, 32), (32, 32), (52, 32), (12, 52), (32, 52), (52, 52)]
    edges = [(0, 1), (1, 2), (3, 4), (6, 7), (7, 8), (0, 3), (3, 6), (1, 4), (4, 7), (2, 5), (5, 8)]
    for i, j in edges:
        ic.sec().line(*nodes[i], *nodes[j], 3)
    ic.acc().spath("M38 32 H46", 2.6)
    for p in nodes:
        ic.circle(*p, 5)
    ic.acc().circle(52, 32, 2.6)
    ic.acc().circle(32, 32, 2.6)


@icon("sulfonate_group", "materials", label="Sulfonate / SO3H", tags=["sulfonic acid", "Brønsted acid", "functionalisation", "SO3H"])
def _(ic):
    ic.circle(32, 34, 9)
    ic.cut().spath("M35.5 29.5 C33.5 28 28.5 28 28.5 31.5 C28.5 35 35.5 33 35.5 36.5 C35.5 40 30 40 28 38.5", 2.4)
    for a in (-90, 30, 150):
        x1, y1 = polar(32, 34, 9, a)
        x2, y2 = polar(32, 34, 19, a)
        ic.sec().line(x1, y1, x2, y2, 3)
        ic.sec().circle(x2, y2, 5)
    hx, hy = polar(32, 34, 19, 30)
    ic.acc().circle(hx + 7, hy - 4, 3)


@icon("metal_oxide_cluster", "materials", label="Oxide Cluster on Support", tags=["ZrO2", "HfO2", "nanocluster", "supported", "catalyst"])
def _(ic):
    ic.rect(4, 44, 56, 14, 3)
    for x in range(10, 58, 8):
        ic.cut().circle(x, 51, 1.6)
    for x, y, r, role in ((22, 36, 6, "s"), (32, 30, 7, "a"), (42, 36, 6, "s"), (27, 24, 4.5, "s"), (38, 23, 4.5, "s")):
        ic.circle(x, y, r, role=role)


@icon("active_site", "materials", label="Active Site", tags=["catalysis", "active site", "binding", "substrate"])
def _(ic):
    ic.path("M4 58 V30 C14 30 18 44 32 44 C46 44 50 30 60 30 V58 Z")
    ic.acc().circle(32, 44, 4)
    ic.sec().path("M22 18 C22 10 42 10 42 18 C42 28 32 32 32 32 C32 32 22 28 22 18 Z")
    ic.acc().spath("M48 10 l4 -4 m-4 0 l4 4", 2)


@icon("catalytic_cycle", "materials", label="Catalytic Cycle", tags=["catalysis", "turnover", "mechanism", "cycle"])
def _(ic):
    for a0 in (0, 120, 240):
        p1 = polar(32, 32, 20, a0 + 12)
        p2 = polar(32, 32, 20, a0 + 100)
        ic.sec().spath(f"M{_f(p1[0])} {_f(p1[1])} A20 20 0 0 1 {_f(p2[0])} {_f(p2[1])}", 3)
        tip = polar(32, 32, 20, a0 + 108)
        l = polar(32, 32, 15.5, a0 + 98)
        r = polar(32, 32, 24.5, a0 + 98)
        ic.sec().poly([tip, l, r])
        ic.circle(*polar(32, 32, 20, a0 + 4), 5)
    ic.acc().circle(32, 32, 5)


@icon("hydrolysis", "proteomics", label="Hydrolysis", tags=["hydrolysis", "water", "bond cleavage", "amide"])
def _(ic):
    ic.circle(12, 36, 6)
    ic.circle(52, 36, 6)
    ic.line(18, 36, 26, 36, 3)
    ic.line(38, 36, 46, 36, 3)
    ic.acc().poly([(31, 26), (35, 26), (33.5, 34), (37, 34), (30, 46), (31.5, 38), (28, 38)])
    ic.sec().circle(32, 12, 5)
    ic.sec().circle(24.5, 7, 3)
    ic.sec().circle(39.5, 7, 3)


@icon("trypsin_digest", "proteomics", label="Proteolytic Digest", tags=["trypsin", "digestion", "protease", "cleavage site"])
def _(ic):
    xs = [6, 14, 22, 30, 38, 46, 54]
    for i in range(len(xs) - 1):
        if i == 3:
            continue
        ic.sec().line(xs[i], 44, xs[i + 1], 44, 2.6)
    for i, x in enumerate(xs):
        ic.circle(x, 44, 3.6) if i not in (3,) else ic.acc().circle(x, 44, 3.6)
    # scissors
    ic.scircle(26, 14, 4, 2.6)
    ic.scircle(38, 14, 4, 2.6)
    ic.line(28.5, 17, 36, 33, 2.8)
    ic.line(35.5, 17, 28, 33, 2.8)


@icon("peptide_fragmentation", "proteomics", label="MS/MS Fragmentation", tags=["b ions", "y ions", "CID", "HCD", "sequencing"])
def _(ic):
    xs = [8, 20, 32, 44, 56]
    for i in range(4):
        ic.sec().line(xs[i], 36, xs[i + 1], 36, 2.6)
    for x in xs:
        ic.circle(x, 36, 4.2)
    for x in (14, 26, 38, 50):
        ic.acc().spath(f"M{x + 3} 22 H{x} V50 H{x - 3}", 2.4)


@icon("protein_ribbon", "proteomics", label="Protein Structure", tags=["helix", "secondary structure", "enzyme", "fold"])
def _(ic):
    pts = []
    for i in range(0, 61):
        t = i / 60
        x = 10 + 44 * t
        y = 22 + 8 * math.sin(t * math.pi * 5)
        pts.append((x, y))
    ic.spath(smooth_path(pts[::3]), 5)
    ic.acc().poly([(12, 44), (40, 44), (40, 40), (50, 47), (40, 54), (40, 50), (12, 50)])
