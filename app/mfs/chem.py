"""Atomic masses, formula parser and the species library used by the defect calculator."""
from __future__ import annotations
import re

# IUPAC 2021 standard atomic weights (abridged, conventional values)
AW = {"H": 1.008, "D": 2.014, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Na": 22.990,
      "P": 30.974, "S": 32.06, "Cl": 35.45, "K": 39.098, "Br": 79.904, "Zr": 91.224, "Ce": 140.116,
      "Hf": 178.49, "Ti": 47.867, "Fe": 55.845, "Zn": 65.38, "Al": 26.982, "Si": 28.085}
# v2.11: the rest of the metals/elements met in MOFs (node library + custom nodes/linkers)
AW.update({"Li": 6.94, "Be": 9.0122, "B": 10.81, "Mg": 24.305, "Ca": 40.078, "Sc": 44.956, "V": 50.942, "Cr": 51.996,
           "Mn": 54.938, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971,
           "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Nb": 92.906, "Mo": 95.95, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42,
           "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71, "Sb": 121.76, "I": 126.90, "Cs": 132.91, "Ba": 137.33,
           "La": 138.91, "Pr": 140.91, "Nd": 144.24, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93, "Dy": 162.50,
           "Ho": 164.93, "Er": 167.26, "Yb": 173.05, "Lu": 174.97, "Ta": 180.95, "W": 183.84, "Re": 186.21, "Ir": 192.22,
           "Pt": 195.08, "Au": 196.97, "Hg": 200.59, "Pb": 207.2, "Bi": 208.98, "Th": 232.04, "U": 238.03})


def molar_mass(formula: str) -> float:
    """Molar mass of a formula such as 'Zr6O4(OH)4(C8H4O4)6' or 'C4H4O7S'. Decimal subscripts allowed."""
    f = formula.replace(" ", "").replace("·", "").replace("μ3-", "").replace("µ3-", "")
    tokens = re.findall(r"([A-Z][a-z]?|\(|\)|\[|\]|\d+(?:\.\d+)?|\.\d+)", f)
    stack = [{}]
    i = 0
    while i < len(tokens):
        t = tokens[i]
        if t in "([":
            stack.append({}); i += 1; continue
        if t in ")]":
            grp = stack.pop(); i += 1
            mult = 1.0
            if i < len(tokens) and re.fullmatch(r"\d+(?:\.\d+)?|\.\d+", tokens[i]):
                mult = float(tokens[i]); i += 1
            for k, v in grp.items():
                stack[-1][k] = stack[-1].get(k, 0) + v * mult
            continue
        if t in AW:
            i += 1
            n = 1.0
            if i < len(tokens) and re.fullmatch(r"\d+(?:\.\d+)?|\.\d+", tokens[i]):
                n = float(tokens[i]); i += 1
            stack[-1][t] = stack[-1].get(t, 0) + n
            continue
        raise ValueError(f"Cannot parse '{t}' in formula '{formula}'")
    comp = stack[0]
    return sum(AW[k] * v for k, v in comp.items())


M_O, M_H = AW["O"], AW["H"]
M_OH = M_O + M_H
M_H2O = M_O + 2 * M_H
M_CL = AW["Cl"]


def metal_mass(f_hf: float = 0.0, base: str = "Zr") -> float:
    """Average metal molar mass for a Zr/Hf mixed node; f_hf = Hf/(Hf+Zr) molar fraction."""
    if base == "Ce":
        return AW["Ce"]
    return (1 - f_hf) * AW["Zr"] + f_hf * AW["Hf"]


def oxide_mass(f_hf: float = 0.0, base: str = "Zr") -> float:
    """Molar mass of MO2 for the (mixed) metal."""
    return metal_mass(f_hf, base) + 2 * M_O


# Species library: coordinated (deprotonated) form, charge z = number of carboxylates bound,
# nH and shift give the default 1H NMR signal used for quantification after digestion.
SPECIES = {
    # name: (formula of bound form, z, role, nH, typical δ / ppm and comment)
    "BDC (terephthalate)": ("C8H4O4", 2, "linker", 4, "7.8-8.0 (s, 4H), base or HF digest"),
    "Fumarate": ("C4H2O4", 2, "linker", 2, "6.5-6.6 (s, 2H)"),
    "TBAPy (NU-1000)": ("C44H22O8", 4, "linker", 22, "aromatic 7.5-8.3, use all 22 H or a resolved 8H doublet"),
    "BTC (trimesate)": ("C9H3O6", 3, "linker", 3, "8.4-8.6 (s, 3H)"),
    "BDC-SO3H (2-sulfoterephthalate, SO3- form)": ("C8H3O7S", 2, "linker", 3, "7.9-8.6 (3H); sulfonate counted as non-coordinating"),
    "2-Sulfosuccinate (dicarboxylate, SO3H kept)": ("C4H4O7S", 2, "linker", 3, "CH 3.6-3.8 (dd, 1H), CH2 2.8-3.1 (2H)"),
    "2-Sulfosuccinate (monodentate cap)": ("C4H5O7S", 1, "cap", 3, "as above; use when bound through one COO-"),
    "Formate": ("CHO2", 1, "cap", 1, "8.4-8.5 (s, 1H); volatile as HCOOH in HF/evaporation workups"),
    "Acetate": ("C2H3O2", 1, "cap", 3, "1.9 (s, 3H)"),
    "Benzoate": ("C7H5O2", 1, "cap", 5, "7.4-7.9 (5H)"),
    "Trifluoroacetate": ("C2F3O2", 1, "cap", 0, "no 1H: quantify by 19F NMR"),
    "Propionate": ("C3H5O2", 1, "cap", 3, "1.05 (t, 3H)"),
    "DMF (neutral, in pores)": ("C3H7NO", 0, "neutral", 1, "7.9 (s, 1H, CHO) or 2.9/2.7 (2x3H)"),
    "Dimethylammonium (DMA+)": ("C2H8N", -1, "neutral", 6, "2.7 (s, 6H); cation - enter as charge -1 (adds +1)"),
    "Methanol (in pores)": ("CH4O", 0, "neutral", 3, "3.34 (s, 3H)"),
    # v2.11: linkers/guests of the non-Zr node library
    "NH2-BDC (2-aminoterephthalate)": ("C8H5NO4", 2, "linker", 3, "7.0-7.8 (3H)"),
    "Isophthalate (CAU-10)": ("C8H4O4", 2, "linker", 4, "8.5 (1H), 8.1 (2H), 7.5 (1H)"),
    "BPDC (biphenyldicarboxylate, UiO-67)": ("C14H8O4", 2, "linker", 8, "7.8-8.1 (8H)"),
    "dobdc (MOF-74, carboxylate + phenolate)": ("C8H2O6", 4, "linker", 2, "7.3 (s, 2H)"),
    "2-Methylimidazolate (ZIF-8/67)": ("C4H5N2", 1, "linker", 3, "2.6 (s, 3H, CH3) or 7.3 (2H)"),
    "Imidazolate": ("C3H3N2", 1, "linker", 3, "7.7 (1H), 7.1 (2H)"),
    "DABCO (pillar, neutral)": ("C6H12N2", 0, "neutral", 12, "3.0 (s, 12H)"),
    "Fluoride (terminal)": ("F", 1, "cap", 0, "no 1H: 19F NMR / ion chromatography"),
    "Hydroxide (terminal, explicit)": ("OH", 1, "cap", 0, "not seen by 1H NMR"),
}

# Framework definitions. CN = ideal node connectivity (carboxylates per node, charge units); node/metal from mfs.nodes.
FRAMEWORKS = {
    "UiO-66": dict(CN=12, linker="BDC (terephthalate)", L_ideal=6, a=20.75, ideal="M6O4(OH)4(BDC)6"),
    "MOF-801": dict(CN=12, linker="Fumarate", L_ideal=6, a=17.83, ideal="M6O4(OH)4(fum)6"),
    "NU-1000": dict(CN=8, linker="TBAPy (NU-1000)", L_ideal=2, a=39.6, ideal="M6O4(OH)4(OH)4(H2O)4(TBAPy)2"),
    "MOF-808": dict(CN=6, linker="BTC (trimesate)", L_ideal=2, a=35.1, ideal="M6O4(OH)4(BTC)2(HCOO)6"),
    "UiO-66-SO3H": dict(CN=12, linker="BDC-SO3H (2-sulfoterephthalate, SO3- form)", L_ideal=6, a=20.75, ideal="M6O4(OH)4(BDC-SO3H)6"),
    # v2.11 node library: 'node' = key in mfs.nodes.NODES (default the M6O4(OH)4 node), 'metal' = default metal
    "UiO-66-NH2": dict(CN=12, linker="NH2-BDC (2-aminoterephthalate)", L_ideal=6, a=20.75, ideal="M6O4(OH)4(NH2-BDC)6"),
    "UiO-67": dict(CN=12, linker="BPDC (biphenyldicarboxylate, UiO-67)", L_ideal=6, a=26.8, ideal="M6O4(OH)4(BPDC)6"),
    "Ce-UiO-66": dict(CN=12, linker="BDC (terephthalate)", L_ideal=6, a=21.5, ideal="Ce6O4(OH)4(BDC)6", metal="Ce"),
    "Th-UiO-66": dict(CN=12, linker="BDC (terephthalate)", L_ideal=6, a=21.8, ideal="Th6O4(OH)4(BDC)6", metal="Th"),
    "HKUST-1": dict(CN=4, linker="BTC (trimesate)", L_ideal=4 / 3, a=26.34, ideal="Cu3(BTC)2 = Cu2(BTC)4/3 per paddlewheel",
                    node="Cu2 paddlewheel", metal="Cu"),
    "DMOF-1": dict(CN=4, linker="BDC (terephthalate)", L_ideal=2, a=None, ideal="Zn2(BDC)2(dabco)", node="Zn2 paddlewheel (pillared)", metal="Zn"),
    "MOF-5": dict(CN=6, linker="BDC (terephthalate)", L_ideal=3, a=25.83, ideal="Zn4O(BDC)3", node="Zn4O", metal="Zn"),
    "MIL-53(Al)": dict(CN=2, linker="BDC (terephthalate)", L_ideal=1, a=None, ideal="Al(OH)(BDC)", node="M(OH) chain (Al-OH, MIL-53 / CAU-10)", metal="Al"),
    "CAU-10(Al)": dict(CN=2, linker="Isophthalate (CAU-10)", L_ideal=1, a=None, ideal="Al(OH)(ipa)", node="M(OH) chain (Al-OH, MIL-53 / CAU-10)", metal="Al"),
    "Al-fumarate (A520)": dict(CN=2, linker="Fumarate", L_ideal=1, a=None, ideal="Al(OH)(fum)", node="M(OH) chain (Al-OH, MIL-53 / CAU-10)", metal="Al"),
    "MIL-100(Fe)": dict(CN=6, linker="BTC (trimesate)", L_ideal=2, a=73.3, ideal="Fe3O(F/OH)(H2O)2(BTC)2", node="M3O trimer (MIL-100 / MIL-101 / MIL-88)", metal="Fe"),
    "MIL-101(Cr)": dict(CN=6, linker="BDC (terephthalate)", L_ideal=3, a=88.9, ideal="Cr3O(F/OH)(H2O)2(BDC)3", node="M3O trimer (MIL-100 / MIL-101 / MIL-88)", metal="Cr"),
    "MIL-88B(Fe)": dict(CN=6, linker="BDC (terephthalate)", L_ideal=3, a=None, ideal="Fe3O(OH)(H2O)2(BDC)3", node="M3O trimer (MIL-100 / MIL-101 / MIL-88)", metal="Fe"),
    "Mg-MOF-74": dict(CN=4, linker="dobdc (MOF-74, carboxylate + phenolate)", L_ideal=1, a=None, ideal="Mg2(dobdc)", node="M2 MOF-74 / CPO-27", metal="Mg"),
    "Co-MOF-74": dict(CN=4, linker="dobdc (MOF-74, carboxylate + phenolate)", L_ideal=1, a=None, ideal="Co2(dobdc)", node="M2 MOF-74 / CPO-27", metal="Co"),
    "Ni-MOF-74": dict(CN=4, linker="dobdc (MOF-74, carboxylate + phenolate)", L_ideal=1, a=None, ideal="Ni2(dobdc)", node="M2 MOF-74 / CPO-27", metal="Ni"),
    "MIL-125(Ti)": dict(CN=12, linker="BDC (terephthalate)", L_ideal=6, a=None, ideal="Ti8O8(OH)4(BDC)6", node="Ti8O8(OH)4 (MIL-125)", metal="Ti"),
    "NH2-MIL-125(Ti)": dict(CN=12, linker="NH2-BDC (2-aminoterephthalate)", L_ideal=6, a=None, ideal="Ti8O8(OH)4(NH2-BDC)6", node="Ti8O8(OH)4 (MIL-125)", metal="Ti"),
    "ZIF-8": dict(CN=2, linker="2-Methylimidazolate (ZIF-8/67)", L_ideal=2, a=16.99, ideal="Zn(mIm)2 (CN in charge units: 4 N donors = 2 mIm-)", node="M (ZIF, tetrahedral)", metal="Zn"),
    "ZIF-67": dict(CN=2, linker="2-Methylimidazolate (ZIF-8/67)", L_ideal=2, a=16.95, ideal="Co(mIm)2", node="M (ZIF, tetrahedral)", metal="Co"),
    "custom framework": dict(CN=12, linker="BDC (terephthalate)", L_ideal=6, a=None, ideal="(set node + linker)"),
}


def species_mass(name: str, custom: dict | None = None) -> float:
    if custom and name in custom:
        return molar_mass(custom[name][0])
    return molar_mass(SPECIES[name][0])
