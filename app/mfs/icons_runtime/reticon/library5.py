"""
v5 batch 1: characterisation techniques, energy & electrochemistry, adsorption / environment,
synthesis methods. All original geometry (CC0).
"""
from __future__ import annotations

import math

from .core import Op, polar, lorentz_trace
from .library import icon, axes
from .library2 import gauss_trace
from .library4 import arrowhead


def bubbles(ic, pts, r=1.4, role="a"):
    for x, y in pts:
        ic.circle(x, y, r, role=role)


def zig(ic, x0, y0, x1, y1, n=6, amp=2.5, w=2, role="a"):
    pts = []
    for i in range(n * 2 + 1):
        t = i / (n * 2)
        x, y = x0 + (x1 - x0) * t, y0 + (y1 - y0) * t
        nx, ny = -(y1 - y0), (x1 - x0)
        L = math.hypot(nx, ny) or 1
        o = amp * (1 if i % 2 else -1) if 0 < i < n * 2 else 0
        pts.append((x + nx / L * o, y + ny / L * o))
    ic.spoly(pts, w, role=role)


# ==========================================================================
# 1. CHARACTERISATION
# ==========================================================================

@icon("scxrd", "characterisation", label="Single-Crystal XRD", tags=["SCXRD", "single crystal", "goniometer", "structure determination"])
def _(ic):
    ic.rect(4, 26, 14, 10, 2)
    ic.acc().line(18, 31, 29, 31, 2.6)
    ic.line(32, 58, 32, 40, 4)
    ic.rect(24, 54, 16, 6, 2)
    ic.circle(32, 38, 4)
    ic.line(32, 34, 32, 31, 1.6)
    ic.sec().poly([(32, 27), (35, 30), (32, 34), (29, 30)])
    for a in (-50, -20, 15, 45):
        x2, y2 = polar(32, 31, 22, a)
        ic.acc().line(35, 31, x2, y2, 1.4)
    ic.rect(55, 12, 6, 38, 2)


@icon("crystal_on_loop", "characterisation", label="Crystal on Loop", tags=["single crystal", "mounting loop", "cryoloop", "SCXRD"])
def _(ic):
    ic.rect(29, 40, 6, 22, 1.5)
    ic.line(32, 40, 32, 24, 1.6)
    ic.scircle(32, 17, 8, 1.6)
    ic.sec().poly([(32, 11), (37, 15), (35.5, 21.5), (28.5, 21.5), (27, 15)])
    ic.cut().line(32, 11, 32, 21.5, 0.9)


@icon("saxs", "characterisation", label="SAXS / WAXS", tags=["SAXS", "WAXS", "small-angle scattering", "nanostructure"])
def _(ic):
    ic.rect(2, 27, 12, 10, 2)
    ic.acc().line(14, 32, 26, 32, 2.6)
    ic.sec().circle(28, 32, 3)
    ic.ops.append(Op("a", "path", {"d": "M31 32 L54 18 L54 46 Z", "opacity": 0.35}))
    ic.rect(54, 8, 8, 48, 2)
    for r in (5, 9, 13):
        ic.ops.append(Op("a", "ellipse", {"cx": 58, "cy": 32, "rx": 1.8, "ry": r}, 1.2))


@icon("saxs_pattern", "data", label="SAXS Curve", tags=["SAXS", "I(q)", "log-log", "scattering"])
def _(ic):
    axes(ic)
    pts = [(12 + i, 12 + 40 * (1 - math.exp(-i / 6)) * 0.85 + 4 * math.sin(i / 2.2) * math.exp(-i / 14)) for i in range(0, 44)]
    ic.acc().spoly(pts, 2.6)


@icon("synchrotron", "characterisation", label="Synchrotron", tags=["synchrotron", "beamline", "XAS", "storage ring"])
def _(ic):
    ic.scircle(26, 32, 20, 5)
    ic.cut().scircle(26, 32, 20, 1.2)
    ic.acc().circle(26, 12, 2.6)
    ic.acc().line(40, 18, 58, 8, 2.4)
    ic.rect(52, 4, 10, 8, 2)
    ic.acc().line(44, 24, 60, 22, 2)
    ic.rect(56, 18, 6, 8, 1.5)


@icon("xas_spectrum", "data", label="XAS / XANES", tags=["XAS", "XANES", "EXAFS", "absorption edge"])
def _(ic):
    axes(ic)
    pts = []
    for i in range(0, 45):
        x = 12 + i
        e = 1 / (1 + math.exp(-(x - 24) / 1.2))
        y = 50 - 28 * e - 10 * math.exp(-((x - 26) ** 2) / 4) + 3 * math.sin((x - 28) / 2.5) * e * math.exp(-(x - 28) / 20)
        pts.append((x, y))
    ic.acc().spoly(pts, 2.6)


@icon("neutron_diffraction", "characterisation", label="Neutron Diffraction", tags=["neutron", "diffraction", "spallation", "deuterated"])
def _(ic):
    ic.rect(2, 22, 16, 20, 3)
    ic.cut().circle(10, 32, 3.5)
    for x in (22, 28, 34):
        ic.acc().circle(x, 32, 2.2)
        ic.acc().line(x - 3, 32, x - 5.5, 32, 1.4)
    ic.sec().rect(38, 26, 6, 12, 1)
    ic.spath("M48 12 A26 26 0 0 1 48 52", 4)


@icon("epr_spectrometer", "characterisation", label="EPR Spectrometer", tags=["EPR", "ESR", "radicals", "electron spin"])
def _(ic):
    ic.rect(4, 16, 14, 34, 3)
    ic.rect(46, 16, 14, 34, 3)
    ic.rect(18, 24, 28, 18, 2)
    ic.cut().rect(28, 24, 8, 18)
    ic.srect(29.5, 4, 5, 40, 2, 2, glass=True)
    ic.sec().rect(30.5, 30, 3, 8, 1)
    ic.acc().spath("M20 54 q3 -6 6 0 t6 0 t6 0 t6 0", 2)


@icon("epr_spectrum", "data", label="EPR Spectrum", tags=["EPR", "ESR", "derivative", "g-value"])
def _(ic):
    ic.line(8, 32, 58, 32, 1.6)
    ic.line(9, 6, 9, 58, 3)
    pts = []
    for i in range(0, 46):
        x = 12 + i
        y = 32 - 220 * (x - 34) * math.exp(-((x - 34) ** 2) / 18) / 30
        pts.append((x, y))
    ic.acc().spoly(pts, 2.6)


@icon("mossbauer_spectrum", "data", label="Mössbauer Spectrum", tags=["Mössbauer", "57Fe", "quadrupole doublet", "iron"])
def _(ic):
    axes(ic)
    ic.acc().spoly(lorentz_trace(12, 56, 14, [(28, 24, 1.6), (40, 22, 1.6)], n=160, sign=1), 2.6)
    for x, y in lorentz_trace(12, 56, 14, [(28, 24, 1.6), (40, 22, 1.6)], n=28, sign=1)[1::2]:
        ic.circle(x, y, 1)


@icon("tof_sims", "characterisation", label="ToF-SIMS", tags=["SIMS", "secondary ions", "surface", "time of flight"])
def _(ic):
    ic.rect(4, 48, 30, 8, 2)
    ic.acc().line(6, 10, 18, 44, 2.6)
    ic.acc().circle(18, 46, 2.2)
    for k, (x, y) in enumerate(((26, 40), (30, 34), (34, 28))):
        ic.sec().circle(x, y, 2 - k * 0.3)
    ic.spath("M36 26 H58 V12 H40", 3)
    ic.rect(36, 6, 8, 10, 2)


@icon("chn_analyzer", "instruments", label="CHN / Elemental Analyser", tags=["CHN", "CHNS", "elemental analysis", "combustion"])
def _(ic):
    ic.rect(6, 18, 52, 40, 4)
    ic.cut().rect(12, 24, 14, 26, 2)
    ic.acc().path("M19 30 C23 35 23 40 21 43 C20 44.5 18 44.5 17 43 C15 40 15 35 19 30 Z")
    ic.cut().rect(32, 24, 20, 9, 1.5)
    for k, role in enumerate(("m", "s", "a")):
        ic.rect(33 + k * 7, 40, 5, 10, 1, role)
    ic.sec().circle(32, 12, 4)
    ic.line(32, 16, 32, 18, 2)


@icon("karl_fischer", "instruments", label="Karl Fischer Titrator", tags=["Karl Fischer", "water content", "titration", "moisture"])
def _(ic):
    ic.rect(4, 36, 30, 22, 3)
    ic.cut().rect(8, 41, 16, 8, 1.2)
    ic.spath("M38 26 V52 Q38 58 44 58 H54 Q60 58 60 52 V26", 2.6, glass=True)
    ic.sec().rect(41.5, 38, 15, 17.5, 2)
    ic.line(46, 8, 46, 48, 2.4)
    ic.line(52, 12, 52, 48, 2.4)
    ic.acc().rect(44, 4, 10, 6, 1.5)
    ic.line(34, 46, 46, 46, 1.6)


@icon("chemisorption_analyzer", "instruments", label="Chemisorption (TPD / TPR)", tags=["TPD", "TPR", "TPO", "chemisorption", "acid sites"])
def _(ic):
    ic.rect(4, 8, 30, 50, 4)
    ic.cut().rect(9, 14, 20, 30, 2)
    ic.spath("M17 10 V34 A2 2 0 0 0 21 34 V10", 2, glass=True)
    ic.sec().rect(17.5, 28, 3, 6, 1)
    for y in (18, 24, 30, 36):
        ic.acc().line(10, y, 15, y, 1.4)
    ic.rect(38, 22, 22, 30, 3)
    ic.cut().rect(42, 27, 14, 8, 1.2)
    ic.acc().line(34, 30, 38, 30, 2)


@icon("tpd_profile", "data", label="TPD / TPR Profile", tags=["TPD", "TPR", "desorption", "acid strength"])
def _(ic):
    axes(ic)
    ic.acc().spoly(gauss_trace(12, 56, 52, [(24, 24, 3.5), (42, 16, 5)], n=160), 2.6)
    ic.sec().spoly(gauss_trace(12, 56, 52, [(24, 24, 3.5)], n=100), 1.6)


@icon("hg_porosimeter", "instruments", label="Hg Porosimeter", tags=["mercury intrusion", "porosimetry", "macropores"])
def _(ic):
    ic.rect(6, 30, 52, 28, 4)
    ic.cut().rect(12, 36, 16, 8, 1.2)
    ic.spath("M38 4 V22 A6 6 0 0 0 50 22 V4", 2.4, glass=True)
    ic.sec().circle(44, 22, 4)
    ic.ops[-1].attrs["_fill"] = "#AAB4BF"
    ic.line(44, 28, 44, 30, 2)
    ic.acc().circle(46, 47, 3)


@icon("pycnometer", "instruments", label="Gas Pycnometer", tags=["He pycnometry", "skeletal density", "true density"])
def _(ic):
    ic.rect(4, 20, 56, 36, 4)
    ic.circle(22, 38, 11)
    ic.cut().circle(22, 38, 7.5)
    ic.sec().circle(22, 40, 4.5)
    ic.cut().rect(38, 26, 16, 9, 1.5)
    ic.acc().circle(46, 46, 4)
    ic.rect(18, 12, 8, 8, 2)


@icon("rheometer", "instruments", label="Rheometer", tags=["rheology", "viscosity", "gel", "storage modulus"])
def _(ic):
    ic.rect(8, 50, 48, 8, 2)
    ic.rect(10, 6, 8, 44, 2)
    ic.rect(10, 6, 34, 10, 3)
    ic.rect(29, 16, 6, 18)
    ic.rect(20, 34, 24, 3, 1.5)
    ic.sec().rect(22, 37, 20, 4, 1)
    ic.rect(18, 41, 28, 4, 1)
    ic.acc().spath("M20 28 A12 4 0 0 0 44 28", 1.8)
    arrowhead(ic, 45, 27, -20, 3.5)


@icon("tensile_tester", "instruments", label="Mechanical Testing", tags=["tensile", "compression", "stress-strain", "Young's modulus"])
def _(ic):
    ic.rect(6, 56, 52, 5, 2)
    ic.rect(8, 4, 6, 52, 1.5)
    ic.rect(50, 4, 6, 52, 1.5)
    ic.rect(8, 4, 48, 7, 2)
    ic.rect(14, 18, 36, 6, 1.5)
    ic.rect(28, 24, 8, 6, 1)
    ic.sec().path("M29 30 H35 V33 H33.5 V41 H35 V44 H29 V41 H30.5 V33 H29 Z")
    ic.rect(28, 44, 8, 6, 1)
    arrowhead(ic, 32, 13.5, -90, 3.5)


@icon("stress_strain", "data", label="Stress–Strain Curve", tags=["stress-strain", "mechanical", "modulus", "toughness"])
def _(ic):
    axes(ic)
    ic.acc().spath("M10 54 L22 26 C26 18 32 16 40 16 C46 16 50 20 54 30", 2.8)
    ic.acc().circle(54, 30, 2.4)


@icon("conductivity_cell", "instruments", label="Conductivity Cell", tags=["conductivity", "impedance", "pellet", "proton conduction"])
def _(ic):
    ic.rect(8, 8, 20, 8, 2)
    ic.rect(8, 48, 20, 8, 2)
    ic.rect(14, 16, 8, 10, 1)
    ic.rect(14, 38, 8, 10, 1)
    ic.sec().rect(10, 26, 16, 12, 2)
    ic.acc().spath("M28 12 H44 V24", 1.8)
    ic.acc().spath("M28 52 H44 V40", 1.8)
    ic.rect(36, 24, 22, 16, 3)
    ic.cut().rect(40, 28, 14, 8, 1)


# ==========================================================================
# 2. ENERGY & ELECTROCHEMISTRY
# ==========================================================================

@icon("coin_cell", "energy", label="Coin Cell", tags=["coin cell", "CR2032", "battery", "electrode testing"])
def _(ic):
    ic.ellipse(32, 38, 26, 12)
    ic.rect(6, 28, 52, 10)
    ic.sec().ellipse(32, 28, 26, 12)
    ic.cut().ellipse(32, 28, 18, 7)
    ic.ellipse(32, 28, 15, 5.5)


@icon("pouch_cell", "energy", label="Pouch Cell", tags=["pouch cell", "battery", "Li-ion", "laminate"])
def _(ic):
    ic.rect(20, 4, 7, 12, 1)
    ic.ops[-1].attrs["_fill"] = "#C8CDD3"
    ic.acc().rect(37, 4, 7, 12, 1)
    ic.sec().rect(10, 14, 44, 46, 5)
    ic.cut().rect(16, 20, 32, 34, 3)
    ic.sec().rect(18, 22, 28, 30, 2)


@icon("battery_cell", "energy", label="Battery", tags=["battery", "cylindrical cell", "18650", "energy storage"])
def _(ic):
    ic.rect(26, 4, 12, 5, 1.5)
    ic.rect(18, 8, 28, 52, 6)
    ic.sec().rect(18, 30, 28, 30, 6)
    ic.cut().rect(18, 28, 28, 2.5)
    ic.cut().line(32, 14, 32, 22, 2.4)
    ic.cut().line(28, 18, 36, 18, 2.4)


@icon("supercapacitor", "energy", label="Supercapacitor", tags=["supercapacitor", "EDLC", "double layer", "capacitor"])
def _(ic):
    ic.rect(6, 8, 8, 48, 2)
    ic.rect(50, 8, 8, 48, 2)
    ic.sec().rect(14, 10, 14, 44, 1)
    ic.acc().rect(36, 10, 14, 44, 1)
    ic.cut().rect(30, 8, 4, 48)
    for y in (16, 26, 36, 46):
        ic.acc().circle(17, y, 2)
        ic.sec().circle(47, y, 2)


@icon("fuel_cell", "energy", label="Fuel Cell", tags=["fuel cell", "PEMFC", "MEA", "hydrogen"])
def _(ic):
    ic.rect(4, 10, 10, 44, 2)
    ic.rect(50, 10, 10, 44, 2)
    ic.sec().rect(16, 10, 8, 44, 1)
    ic.rect(26, 10, 12, 44, 1)
    ic.ops[-1].attrs["_fill"] = "#F2E6B5"
    ic.sec().rect(40, 10, 8, 44, 1)
    for y in (18, 30, 42):
        ic.acc().line(18, y, 44, y, 1.6)
        arrowhead(ic, 46, y, 0, 3)
    ic.spath("M9 10 V4 H55 V10", 2)


@icon("electrolyzer", "energy", label="Electrolyser", tags=["water splitting", "electrolysis", "OER", "HER", "hydrogen"])
def _(ic):
    ic.spath("M6 22 V52 Q6 58 12 58 H52 Q58 58 58 52 V22", 3, glass=True)
    ic.sec().rect(9.5, 30, 45, 25.5, 2)
    ic.rect(18, 10, 5, 40, 1.5)
    ic.rect(41, 10, 5, 40, 1.5)
    bubbles(ic, [(15, 36), (16, 44), (14.5, 28), (26, 40), (25, 32)], 1.6)
    bubbles(ic, [(49, 38), (50, 30), (38, 44)], 2.2, "m")
    ic.line(20.5, 10, 20.5, 4, 2)
    ic.line(43.5, 10, 43.5, 4, 2)
    ic.line(20.5, 4, 43.5, 4, 2)


@icon("working_electrode", "energy", label="Working Electrode", tags=["working electrode", "glassy carbon", "GC electrode", "RDE"])
def _(ic):
    ic.rect(26, 4, 12, 44, 3)
    ic.circle(32, 50, 7)
    ic.sec().circle(32, 50, 4)
    ic.acc().line(32, 4, 32, 1, 2)


@icon("counter_electrode", "energy", label="Counter Electrode (Pt)", tags=["counter electrode", "Pt wire", "auxiliary electrode"])
def _(ic):
    ic.rect(27, 2, 10, 20, 3)
    ic.spath("M32 22 V30 C24 32 24 36 32 38 C40 40 40 44 32 46 C24 48 24 52 32 54 V60", 2.2)


@icon("reference_electrode", "energy", label="Reference Electrode", tags=["reference electrode", "Ag/AgCl", "SCE", "frit"])
def _(ic):
    ic.rect(26, 2, 12, 8, 2)
    ic.srect(27, 10, 10, 46, 4, 2.6, glass=True)
    ic.sec().rect(29.5, 26, 5, 28, 2)
    ic.line(32, 10, 32, 44, 1.4)
    ic.acc().rect(30.5, 40, 3, 8, 1)
    ic.rect(28, 54, 8, 4, 1)


@icon("gas_diffusion_electrode", "energy", label="Gas Diffusion Electrode", tags=["GDE", "CO2 reduction", "electrocatalysis", "flow cell"])
def _(ic):
    ic.rect(22, 6, 6, 52, 1)
    ic.sec().rect(28, 6, 5, 52, 1)
    for y in range(9, 58, 6):
        ic.acc().circle(35, y, 1.8)
    for y in (16, 32, 48):
        ic.line(4, y, 16, y, 2)
        arrowhead(ic, 19, y, 0, 3.5, "m")
    bubbles(ic, [(44, 14), (50, 26), (44, 38), (52, 48)], 2.4, "s")


@icon("solar_cell", "energy", label="Solar Cell", tags=["photovoltaic", "PV", "perovskite solar cell", "solar panel"])
def _(ic):
    ic.poly([(6, 44), (18, 20), (58, 20), (46, 44)])
    for t in (0.33, 0.66):
        ic.cut().line(6 + (18 - 6) * t, 44 - 24 * t, 46 + 12 * t, 44 - 24 * t, 1.2)
    for t in (0.25, 0.5, 0.75):
        ic.cut().line(6 + 40 * t, 44, 18 + 40 * t, 20, 1.2)
    ic.rect(28, 44, 6, 12, 1)
    ic.rect(20, 56, 22, 4, 1.5)
    ic.acc().circle(10, 9, 5)


@icon("photoelectrode", "energy", label="Photoelectrode", tags=["photoelectrochemistry", "PEC", "photoanode", "water splitting"])
def _(ic):
    ic.spath("M14 24 V52 Q14 58 20 58 H52 Q58 58 58 52 V24", 2.8, glass=True)
    ic.sec().rect(17.5, 32, 37, 23.5, 2)
    ic.rect(28, 8, 8, 44, 1.5)
    ic.acc().rect(28, 32, 8, 20, 1)
    for k in range(3):
        ic.acc().spath(f"M{2 + k * 2} {14 + k * 8} q3 3 6 0 t6 0", 1.8)
    bubbles(ic, [(40, 40), (42, 34), (39, 46)], 1.6, "m")


@icon("h2_tank", "energy", label="H₂ Storage Tank", tags=["hydrogen storage", "tank", "H2", "gas storage"])
def _(ic):
    ic.rect(6, 18, 52, 30, 15)
    ic.cut().rect(12, 30, 40, 6, 3)
    ic.acc().rect(14, 31.5, 26, 3, 1.5)
    ic.rect(28, 10, 8, 8, 2)
    ic.rect(14, 48, 6, 8, 1)
    ic.rect(44, 48, 6, 8, 1)


@icon("nyquist_plot", "data", label="Nyquist Plot", tags=["EIS", "impedance", "Nyquist", "charge transfer"])
def _(ic):
    axes(ic)
    ic.acc().spath("M14 53 C14 38 36 38 36 53 L54 24", 2.6)
    for t in (0.2, 0.5, 0.8):
        ic.circle(14 + 22 * t, 53 - 14 * math.sin(math.pi * t), 1.8)


@icon("tafel_plot", "data", label="Tafel Plot", tags=["Tafel slope", "overpotential", "kinetics", "electrocatalysis"])
def _(ic):
    axes(ic)
    ic.acc().spath("M12 50 C18 46 22 40 26 34 L52 14", 2.6)
    ic.sec().line(22, 44, 56, 12, 1.6)
    ic.ops[-1].attrs["stroke-dasharray"] = "3 2"


@icon("lsv_curve", "data", label="LSV / Polarisation", tags=["LSV", "polarisation curve", "onset", "current density"])
def _(ic):
    ic.line(56, 10, 56, 56, 3)
    ic.line(8, 56, 56, 56, 3)
    ic.acc().spath("M54 54 C40 54 34 52 28 44 C22 34 18 20 14 8", 2.6)
    ic.sec().spath("M54 54 C44 54 40 52 34 46 C28 38 24 26 20 14", 2)


@icon("charge_discharge", "data", label="Charge–Discharge", tags=["galvanostatic", "GCD", "capacity", "battery"])
def _(ic):
    axes(ic)
    ic.acc().spoly([(12, 50), (24, 16), (36, 50), (48, 16), (56, 36)], 2.6)


@icon("cycling_stability", "data", label="Cycling Stability", tags=["cycle life", "capacity retention", "coulombic efficiency"])
def _(ic):
    axes(ic)
    for k in range(13):
        x = 13 + k * 3.4
        ic.acc().circle(x, 18 + k * 0.6 + (0.8 if k % 2 else 0), 1.6)
        ic.sec().circle(x, 12, 1.2)


# ==========================================================================
# 3. ADSORPTION, SEPARATION & ENVIRONMENT
# ==========================================================================

@icon("breakthrough_column", "characterisation", label="Breakthrough Column", tags=["breakthrough", "fixed bed", "gas separation", "column"])
def _(ic):
    ic.srect(24, 8, 16, 48, 4, 2.6, glass=True)
    ic.sec().rect(27, 16, 10, 34, 2)
    for y in range(20, 48, 6):
        ic.acc().circle(30, y, 1.3)
        ic.circle(34, y + 3, 1.3)
    ic.line(32, 2, 32, 8, 2.4)
    ic.line(32, 56, 32, 62, 2.4)
    arrowhead(ic, 32, 13, 90, 3.5)


@icon("breakthrough_curve", "data", label="Breakthrough Curve", tags=["breakthrough", "C/C0", "selectivity", "dynamic adsorption"])
def _(ic):
    axes(ic)
    ic.sec().spath("M12 52 H18 C22 52 22 14 28 14 H56", 2.6)
    ic.acc().spath("M12 52 H32 C38 52 38 14 44 14 H56", 2.6)


@icon("co2_capture", "materials", label="CO₂ Capture", tags=["CO2", "carbon capture", "adsorption", "DAC"])
def _(ic):
    ic.rect(30, 8, 28, 48, 5)
    ic.cut().rect(35, 13, 18, 38, 3)
    for x, y in ((40, 20), (47, 30), (40, 40)):
        ic.acc().circle(x, y, 2.8)
        ic.sec().circle(x - 4.5, y, 2)
        ic.sec().circle(x + 4.5, y, 2)
    for y in (18, 32, 46):
        ic.acc().circle(12, y, 2.4)
        ic.sec().circle(8, y, 1.6)
        ic.sec().circle(16, y, 1.6)
    ic.line(20, 32, 26, 32, 2)
    arrowhead(ic, 29, 32, 0, 3.5, "m")


@icon("water_harvesting", "materials", label="Water Harvesting", tags=["atmospheric water", "MOF-801", "humidity", "water adsorption"])
def _(ic):
    ic.acc().circle(12, 12, 6)
    for a in range(0, 360, 45):
        x1, y1 = polar(12, 12, 8, a)
        x2, y2 = polar(12, 12, 10.5, a)
        ic.acc().line(x1, y1, x2, y2, 1.4)
    ic.rect(20, 26, 36, 14, 3)
    for x in range(24, 54, 6):
        ic.cut().circle(x, 33, 1.8)
    for x, y in ((28, 48), (38, 52), (48, 47)):
        ic.sec().path(f"M{x} {y - 5} C{x + 3} {y - 1} {x + 3} {y + 2} {x} {y + 2} C{x - 3} {y + 2} {x - 3} {y - 1} {x} {y - 5} Z")
    for x, y in ((26, 16), (36, 12), (46, 18)):
        ic.sec().circle(x, y, 1.6)


@icon("gas_bag", "glassware", label="Gas Sampling Bag", tags=["gas bag", "Tedlar", "gas sampling", "headspace"])
def _(ic):
    ic.sec().rect(10, 14, 44, 44, 10)
    ic.cut().rect(16, 22, 32, 10, 3)
    ic.rect(28, 6, 8, 10, 2)
    ic.acc().rect(25, 3, 14, 4, 2)


@icon("water_treatment", "materials", label="Water Treatment", tags=["water purification", "adsorbent", "remediation", "pollutants"])
def _(ic):
    ic.sec().path("M8 26 C14 22 20 30 26 26 C32 22 38 30 44 26 C50 22 56 30 58 26 V56 H8 Z")
    for x, y in ((16, 38), (30, 44), (46, 36), (22, 50), (40, 50)):
        ic.acc().circle(x, y, 2)
    ic.rect(26, 4, 12, 18, 3)
    ic.cut().rect(29, 8, 6, 10, 1.5)
    ic.line(32, 22, 32, 28, 2.4)


@icon("dye_adsorption", "materials", label="Dye Adsorption", tags=["dye removal", "methylene blue", "decolourisation", "adsorption"])
def _(ic):
    for x, frac in ((18, 1.0), (46, 0.15)):
        ic.spath(f"M{x - 10} 16 V52 Q{x - 10} 58 {x - 4} 58 H{x + 4} Q{x + 10} 58 {x + 10} 52 V16", 2.6, glass=True)
        ic.sec().rect(x - 7, 26, 14, 29.5, 2)
        ic.ops[-1].attrs["_fill"] = "#3B6FD8" if frac > 0.5 else "#CFE0FA"
    ic.line(30, 34, 34, 34, 2.4)
    arrowhead(ic, 37, 34, 0, 3.5, "m")
    for x, y in ((44, 50), (48, 52), (46, 47)):
        ic.acc().circle(x, y, 1.6)


@icon("heavy_metal_removal", "materials", label="Heavy-Metal Removal", tags=["heavy metals", "Pb2+", "Hg2+", "sorbent", "ion capture"])
def _(ic):
    ic.circle(40, 36, 18)
    ic.cut().circle(40, 36, 12)
    ic.sec().circle(40, 36, 9)
    for x, y in ((10, 14), (8, 34), (12, 52)):
        ic.acc().circle(x, y, 4)
        ic.cut().line(x - 1.6, y, x + 1.6, y, 1.2)
        ic.cut().line(x, y - 1.6, x, y + 1.6, 1.2)
        ic.line(x + 5, y + (36 - y) * 0.2, 21, 36 + (y - 36) * 0.3, 1.4)
    ic.acc().circle(38, 34, 2.6)
    ic.acc().circle(43, 39, 2.2)


@icon("langmuir_fit", "data", label="Langmuir / Freundlich Fit", tags=["Langmuir", "Freundlich", "isotherm fit", "qmax"])
def _(ic):
    axes(ic)
    ic.sec().spath("M12 52 C16 26 26 20 56 18", 2)
    for x, y in ((15, 38), (19, 30), (25, 25), (33, 22), (43, 20), (53, 19)):
        ic.acc().circle(x, y, 2.2)
    ic.line(12, 16, 56, 16, 1.2)
    ic.ops[-1].attrs["stroke-dasharray"] = "3 2"


@icon("selectivity_plot", "data", label="Selectivity (IAST)", tags=["IAST", "selectivity", "CO2/N2", "separation"])
def _(ic):
    axes(ic)
    ic.acc().spath("M12 20 C22 22 34 30 56 34", 2.6)
    ic.sec().spath("M12 44 C22 44 34 46 56 48", 2.4)


# ==========================================================================
# 4. SYNTHESIS METHODS
# ==========================================================================

@icon("sol_gel", "synthesis", label="Sol–Gel", tags=["sol-gel", "hydrolysis condensation", "gelation", "silica"])
def _(ic):
    ic.srect(4, 18, 18, 30, 3, 2.4, glass=True)
    ic.sec().rect(7, 28, 12, 17, 1.5)
    for x, y in ((10, 34), (15, 38), (12, 41)):
        ic.acc().circle(x, y, 1)
    ic.line(25, 33, 31, 33, 2.2)
    arrowhead(ic, 34, 33, 0, 3.5, "m")
    ic.srect(38, 18, 22, 30, 3, 2.4, glass=True)
    ic.sec().rect(41, 22, 16, 23, 1.5)
    pts = [(44, 26), (50, 30), (54, 25), (46, 36), (53, 40), (44, 42)]
    for i in range(len(pts) - 1):
        ic.acc().line(*pts[i], *pts[i + 1], 1.2)
    for p in pts:
        ic.acc().circle(*p, 1.3)


@icon("cvd_reactor", "synthesis", label="CVD / ALD Reactor", tags=["CVD", "ALD", "vapour deposition", "thin film growth"])
def _(ic):
    ic.rect(10, 18, 44, 28, 4)
    ic.cut().rect(1, 26, 62, 12, 6)
    ic.srect(2.5, 27.5, 59, 9, 4.5, 2, glass=True)
    ic.sec().rect(26, 33, 12, 3, 0.8)
    for x in (8, 14):
        ic.acc().circle(x, 31, 1.4)
    for x in (22, 32, 42):
        ic.acc().line(x, 14, x, 18, 1.6)
    ic.rect(15, 46, 5, 10, 1.5)
    ic.rect(44, 46, 5, 10, 1.5)


@icon("spray_dryer", "synthesis", label="Spray Dryer", tags=["spray drying", "atomisation", "microspheres", "powder"])
def _(ic):
    ic.rect(26, 2, 12, 8, 2)
    ic.spath("M14 12 H50 V36 L38 52 H26 L14 36 Z", 2.6, glass=True)
    ic.ops.append(Op("a", "path", {"d": "M32 12 L22 34 H42 Z", "opacity": 0.35}))
    for x, y in ((28, 26), (34, 22), (32, 32), (38, 30), (26, 34), (30, 42), (34, 46)):
        ic.sec().circle(x, y, 1.3)
    ic.rect(28, 52, 8, 4)
    ic.sec().rect(24, 56, 16, 6, 2)


@icon("electrospinning", "synthesis", label="Electrospinning", tags=["electrospinning", "nanofibres", "fibre mat", "high voltage"])
def _(ic):
    ic.srect(4, 10, 22, 8, 2, 2, glass=True)
    ic.sec().rect(8, 12, 14, 4, 1)
    ic.rect(1, 12, 4, 4, 1)
    ic.line(26, 14, 30, 14, 1.6)
    ic.acc().poly([(30, 12), (36, 14), (30, 16)])
    ic.acc().spath("M36 14 C42 10 40 22 46 18 C50 16 48 30 54 26", 1.4)
    ic.rect(54, 6, 6, 52, 1.5)
    ic.acc().spath("M8 40 l4 -6 l-3 0 l4 -6", 2)
    ic.line(30, 50, 57, 50, 1.4)


@icon("layering_crystallisation", "synthesis", label="Layering / Slow Diffusion", tags=["crystal growth", "layering", "slow diffusion", "single crystals"])
def _(ic):
    ic.spath("M24 6 V54 A8 8 0 0 0 40 54 V6", 3, glass=True)
    ic.sec().path("M27.5 42 V54 A4.5 4.5 0 0 0 36.5 54 V42 Z")
    ic.acc().rect(27.5, 32, 9, 10, 0.5)
    ic.ops[-1].attrs["opacity"] = 0.55
    ic.sec().rect(27.5, 18, 9, 14, 0.5)
    ic.ops[-1].attrs["_fill"] = "#CFE6F7"
    ic.acc().poly([(31, 38), (33, 36), (35, 38), (33, 40)])
    ic.acc().poly([(29, 44), (31, 42), (33, 44), (31, 46)])


@icon("lbl_growth", "synthesis", label="Layer-by-Layer Growth", tags=["LbL", "SURMOF", "liquid-phase epitaxy", "thin film"])
def _(ic):
    ic.rect(6, 50, 52, 8, 1.5)
    for k, role in enumerate(("s", "a", "s", "a")):
        ic.rect(6, 42 - k * 7, 52, 5, 1, role)
    for k, x in enumerate((14, 32, 50)):
        ic.spath(f"M{x - 5} 4 V12 A5 5 0 0 0 {x + 5} 12 V4", 2, glass=True)
        ic.sec().path(f"M{x - 3} 9 V12 A3 3 0 0 0 {x + 3} 12 V9 Z")


@icon("dip_coater", "synthesis", label="Dip Coating", tags=["dip coating", "withdrawal", "thin film", "coating"])
def _(ic):
    ic.rect(8, 4, 48, 6, 2)
    ic.line(32, 10, 32, 18, 1.6)
    ic.rect(26, 18, 12, 30, 1)
    ic.ops[-1].attrs["_fill"] = "#DCEBF5"
    ic.sec().rect(26, 36, 12, 12, 0.5)
    ic.spath("M14 32 V54 Q14 58 18 58 H46 Q50 58 50 54 V32", 2.6, glass=True)
    ic.sec().rect(17.5, 40, 29, 15.5, 2)
    arrowhead(ic, 44, 14, -90, 4)
    ic.acc().line(44, 14, 44, 26, 2)


@icon("soxhlet", "synthesis", label="Soxhlet Extractor", tags=["Soxhlet", "extraction", "washing", "purification"])
def _(ic):
    ic.srect(28, 2, 8, 14, 2, 2.4, glass=True)
    ic.srect(22, 16, 20, 22, 3, 2.6, glass=True)
    ic.sec().rect(25, 22, 14, 13, 1.5)
    ic.ops[-1].attrs["opacity"] = 0.6
    ic.rect(27, 20, 10, 14, 1)
    ic.ops[-1].attrs["_fill"] = "#E8D9B8"
    ic.spath("M42 20 H47 V40 H38", 1.8, glass=True)
    ic.spath("M29 38 V42 A11 11 0 1 0 35 42 V38", 2.6, glass=True)
    ic.sec().path("M24 52 A8 8 0 0 0 40 52 Z")


@icon("dean_stark", "synthesis", label="Dean–Stark Trap", tags=["Dean-Stark", "azeotropic distillation", "water removal"])
def _(ic):
    ic.srect(40, 2, 8, 18, 2, 2.4, glass=True)
    ic.spath("M44 20 V44", 3, glass=True)
    ic.srect(40, 30, 8, 16, 2, 2.4, glass=True)
    ic.sec().rect(41.5, 38, 5, 7, 1)
    ic.spath("M44 24 H22 V34", 2.6, glass=True)
    ic.spath("M19 34 V38 A11 11 0 1 0 25 38 V34", 2.6, glass=True)
    ic.sec().path("M14 48 A8 8 0 0 0 30 48 Z")
    ic.acc().path("M10 62 H34 L30 58 H14 Z")


@icon("distillation", "synthesis", label="Distillation", tags=["distillation", "condenser", "receiver", "purification"])
def _(ic):
    ic.spath("M11 28 V32 A9 9 0 1 0 17 32 V28", 2.6, glass=True)
    ic.sec().path("M7 42 A7 7 0 0 0 21 42 Z")
    ic.line(14, 28, 14, 22, 2.4)
    ic.ops.append(Op("m", "rect", {"x": 14, "y": 18, "width": 34, "height": 8, "rx": 3, "transform": "rotate(22 14 22)"}, 2.4, True))
    ic.ops.append(Op("s", "rect", {"x": 17, "y": 20, "width": 28, "height": 4, "rx": 2, "transform": "rotate(22 14 22)"}))
    ic.spath("M47 37 V46 L41 58 H59 L53 46 V37", 2.4, glass=True)
    ic.sec().path("M43.6 53 H56.4 L57.6 56 H42.4 Z")
    ic.acc().path("M6 62 H22 L19 58 H9 Z")


@icon("kugelrohr", "synthesis", label="Kugelrohr", tags=["Kugelrohr", "bulb-to-bulb", "short-path distillation", "sublimation"])
def _(ic):
    for k, cx in enumerate((12, 28, 44)):
        ic.scircle(cx, 36, 7, 2.4, glass=True)
        if k < 2:
            ic.line(cx + 7, 36, cx + 9, 36, 2.4)
    ic.sec().path("M7 39 A6 6 0 0 0 17 39 Z")
    ic.sec().circle(44, 39, 2)
    ic.line(51, 36, 60, 36, 2.4)
    ic.rect(2, 28, 22, 16, 3)
    ic.ops[-1].attrs["opacity"] = 0.35
    ic.acc().spath("M58 30 v-6 M58 42 v6", 1.8)
