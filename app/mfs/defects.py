r"""Defect-engineering calculator for M6-oxo-cluster MOFs (M = Zr, Hf or Zr/Hf mixed).

One `DefectInput` is pushed through every applicable method in parallel
("shotgun"). Each method returns a `MethodResult` that carries its numbers AND
its full explanation: rendered equations, symbol definitions, assumptions, a
step-by-step substitution with the user's own numbers, pitfalls and references.

All counts are per M6 node:
    L   linkers per node         (ideal L0 = CN / t; 6 for UiO-66/MOF-801, 2 for NU-1000)
    x   missing linkers per node  x = L0 - L
    m   monocarboxylate caps per node (modulator)
    Cl  chloride per node
    OS  open (OH/OH2-terminated) sites per node
        OS = (12 - CN) + t*x - m - Cl      (reduces to 2x - m - Cl for 12-c UiO-66/MOF-801)
"""
from __future__ import annotations
from dataclasses import dataclass, field, asdict
import math
import numpy as np
from .chem import (molar_mass, SPECIES, FRAMEWORKS, AW, M_O, M_H, M_OH, M_H2O, M_CL,
                   metal_mass, oxide_mass)

NA = 6.02214076e23

REFS = {
    "valenzano": ("Valenzano, L. et al. Disclosing the Complex Structure of UiO-66 Metal Organic Framework: "
                  "A Synergic Combination of Experiment and Theory. Chem. Mater. 2011, 23, 1700-1718.", "10.1021/cm1022882"),
    "shearer2014": ("Shearer, G. C. et al. Tuned to Perfection: Ironing Out the Defects in Metal-Organic Framework "
                    "UiO-66. Chem. Mater. 2014, 26, 4068-4071.", "10.1021/cm501859p"),
    "shearer2016": ("Shearer, G. C. et al. Defect Engineering: Tuning the Porosity and Composition of the "
                    "Metal-Organic Framework UiO-66 via Modulated Synthesis. Chem. Mater. 2016, 28, 3749-3761.",
                    "10.1021/acs.chemmater.6b00602"),
    "oien2014": ("Øien, S. et al. Detailed Structure Analysis of Atomic Positions and Defects in Zirconium "
                 "Metal-Organic Frameworks. Cryst. Growth Des. 2014, 14, 5370-5372.", "10.1021/cg501386j"),
    "sannes2023": ("Sannes, D. K.; Øien-Ødegaard, S.; Aunan, E.; Nova, A.; Olsbye, U. Quantification of Linker "
                   "Defects in UiO-Type Metal-Organic Frameworks. Chem. Mater. 2023, 35, 3793-3800.",
                   "10.1021/acs.chemmater.2c03744"),
    "cliffe2014": ("Cliffe, M. J. et al. Correlated defect nanoregions in a metal-organic framework. "
                   "Nat. Commun. 2014, 5, 4176.", "10.1038/ncomms5176"),
    "deroo2025": ("Pulparayil Mathew, J.; Simms, C.; Salazar Marcano, D. E.; Dhaene, E.; Parac-Vogt, T. N.; "
                  "De Roo, J. Determination of the Experimental Minimal Formula of Metal-Organic Frameworks. "
                  "Adv. Sci. 2025, 12, e04713.", "10.1002/advs.202504713"),
    "klet2016": ("Klet, R. C.; Liu, Y.; Wang, T. C.; Hupp, J. T.; Farha, O. K. Evaluation of Brønsted acidity and "
                 "proton topology in Zr- and Hf-based metal-organic frameworks using potentiometric acid-base "
                 "titration. J. Mater. Chem. A 2016, 4, 1479-1485.", "10.1039/C5TA07687K"),
    "wu2013": ("Wu, H. et al. Unusual and Highly Tunable Missing-Linker Defects in Zirconium Metal-Organic "
               "Framework UiO-66 and Their Important Effects on Gas Adsorption. J. Am. Chem. Soc. 2013, 135, "
               "10525-10532.", "10.1021/ja404514r"),
    "planas2014": ("Planas, N. et al. Defining the Proton Topology of the Zr6-Based Metal-Organic Framework "
                   "NU-1000. J. Phys. Chem. Lett. 2014, 5, 3716-3723.", "10.1021/jz501521x"),
    "atzori2017": ("Atzori, C. et al. Effect of Benzoic Acid as a Modulator in the Structure of UiO-66: An "
                   "Experimental and Computational Study. J. Phys. Chem. C 2017, 121, 9312-9324.", "10.1021/acs.jpcc.7b00483"),
    "vermoortele2013": ("Vermoortele, F. et al. Synthesis Modulation as a Tool To Increase the Catalytic Activity "
                        "of Metal-Organic Frameworks: The Unique Case of UiO-66(Zr). J. Am. Chem. Soc. 2013, 135, "
                        "11465-11468.", "10.1021/ja405078u"),
    "furukawa2014": ("Furukawa, H. et al. Water Adsorption in Porous Metal-Organic Frameworks and Related Materials. "
                     "J. Am. Chem. Soc. 2014, 136, 4369-4381.", "10.1021/ja500330a"),
}


# --------------------------------------------------------------------------- input model
@dataclass
class Species:
    name: str
    formula: str          # bound (deprotonated) form, e.g. C8H4O4
    z: int                # carboxylate charge / number of node sites occupied (0 for neutral guests)
    role: str             # 'linker' | 'cap' | 'neutral'
    nH: int = 1           # protons under the integrated 1H signal
    integral: float | None = None   # integral of that signal (any common scale)
    at_plateau: bool = False        # still present at the TGA dehydroxylated plateau (~400 °C)?
    main: bool = False              # the framework's main linker

    @property
    def M(self):
        return molar_mass(self.formula)


@dataclass
class DefectInput:
    sample: str = "sample"
    framework: str = "MOF-801"
    CN: int | None = None           # ideal connectivity (12/12/8)
    metal: str = "Hf"               # 'Zr' | 'Hf' | 'mixed'
    f_hf: float | None = None       # Hf/(Hf+Zr) molar fraction for 'mixed' (from ICP); None -> from ICP fields
    extraframework_Zr: bool = False # True: all Zr is extra-framework ZrO2 (e.g. ZrO2@Hf-MOF-801 composites)
    species: list = field(default_factory=list)
    # TGA (all % of the initial mass at room temperature; any consistent scale works)
    tga_W_initial: float | None = 100.0
    tga_W_200: float | None = None       # mass at ~200 °C (dry, node still hydroxylated) - Sannes M-200
    tga_W_plateau: float | None = None   # mass at the dehydroxylated plateau (~350-450 °C) - M-400
    tga_W_end: float | None = None       # residue (MO2) after full combustion (≥ 600-800 °C, air)
    # qNMR with internal standard (De Roo protocol)
    nmr_std_conc_mM: float | None = None   # [TMSP-d4] in the NMR tube
    nmr_std_integral: float | None = 9.0   # integral of the standard (set to 9 for 9H)
    nmr_std_nH: int = 9
    nmr_sample_mg: float | None = None     # MOF mass digested
    nmr_volume_mL: float | None = None     # final NMR solution volume
    # Chloride
    cl_per_node: float | None = None       # from XRF/EDX (6*Cl/M) - direct
    cl_conc_mM: float | None = None        # De Roo UV-vis (Spectroquant) on the same digest as qNMR
    cl_wt_pct: float | None = None         # Cl wt% of the sample (combustion / XRF quantitative)
    # ICP-OES (wt% of metal in the sample, sample state = state at which it was weighed)
    icp_Zr_wt: float | None = None
    icp_Hf_wt: float | None = None
    icp_state: str = "as-measured (RT)"    # or 'dried 200 °C'
    linker_umol_per_mg: float | None = None  # absolute qNMR linker content (µmol mg-1) for ICP/qNMR ratio
    # XRF / EDX ratios (atomic)
    xrf_S_per_M: float | None = None       # S/(framework metal) atomic ratio -> sulfonated species per node
    # Porosity
    V_meas: float | None = None            # total (or micro) pore volume, cm3 g-1
    V_ideal: float | None = None           # ideal (defect-free) value from the crystal structure
    S_meas: float | None = None            # BET m2 g-1
    S_ideal: float | None = None
    v_linker_A3: float | None = None       # van der Waals volume freed per missing linker (Å3)
    # Potentiometric titration (Klet et al. 2016)
    tit_mass_mg: float | None = None
    tit_c_NaOH: float | None = 0.1         # mol L-1
    tit_V_eq: list | None = None           # equivalence volumes, mL (2 or 3 values, ascending)
    tit_HCl_mmol: float | None = None      # HCl added before titration (to subtract from eq. 1)
    tit_M_node: float | None = None        # molar mass per node used to convert mass to mol (None -> consensus/ideal)
    # PXRD missing-cluster check
    reo_peaks_seen: bool | None = None
    # v2.11 node library (any metal / node): '' -> the framework's node (default M6O4(OH)4), see mfs.nodes
    node: str = ""
    node_custom: dict | None = None        # custom node: mfs.nodes.custom_node(...) kwargs or a full definition
    icp_M_wt: float | None = None          # ICP framework-metal wt% for non-Zr/Hf nodes (Zr/Hf keep icp_Zr_wt/icp_Hf_wt)

    # ---- derived helpers -------------------------------------------------
    def fw(self):
        return FRAMEWORKS.get(self.framework, FRAMEWORKS["UiO-66"])

    def cn(self):
        return self.CN or self.fw()["CN"]

    def main_linker(self) -> Species:
        for s in self.species:
            if s.main:
                return s
        for s in self.species:
            if s.role == "linker":
                return s
        nm = self.fw()["linker"]
        f, z, r, nH, _ = SPECIES[nm]
        return Species(nm, f, z, "linker", nH, None, True, True)

    def t(self):
        return self.main_linker().z

    def L0(self):
        return self.cn() / self.t()

    # ---- v2.11 node helpers (M6O4(OH)4 Zr/Hf give exactly the v2.10 numbers) ----
    def node_def(self):
        from . import nodes as N
        return N.get(self.node or self.fw().get("node") or N.DEFAULT_NODE, self.node_custom)

    def nM(self):
        """Metal atoms per node (6 for M6O4(OH)4)."""
        return self.node_def()["nM"]

    def Q(self):
        """Charge available for linkers/caps/Cl/terminal OH per node (12 for M6O4(OH)4)."""
        from . import nodes as N
        return N.sites(self.node_def())

    def is_zrhf(self):
        return self.metal in ("Zr", "Hf", "mixed")

    def M_clust(self):
        """Hydroxylated node core, e.g. M6O4(OH)4 = 6M + 8O + 4H."""
        nd = self.node_def()
        return self.nM() * self.MM() + (nd.get("core_O", 0) + nd.get("core_OH", 0)) * M_O + nd.get("core_OH", 0) * M_H

    def M_core(self):
        nd = self.node_def()
        return self.nM() * self.MM() + nd.get("core_O", 0) * M_O + nd.get("core_OH", 0) * M_OH

    def M_oxidic_node(self):
        """Fully oxide-compensated node M_nO_(n·val/2) (= n MO2 for Zr/Hf)."""
        return self.nM() * self.MM() + self.nM() * self.node_def()["val"] / 2 * M_O

    def M_plateau_base(self):
        """Node mass at the TGA plateau without linkers, every ligand charge compensated by O2- (method a convention).
        Dehydroxylating nodes (M6O4(OH)4 -> M6O6): = M_oxidic_node. Nodes whose µ-OH survives the plateau (Al-OH chains):
        M_core + Q/2 O."""
        if self.node_def().get("plateau_dehydroxylated", True):
            return self.M_oxidic_node()
        return self.M_core() + self.Q() / 2 * M_O

    def residue_is_oxidic(self):
        """True if the plateau node base equals the residue mass (ZrO2, CuO, ZnO ...; False for Co3O4 or Al-OH chains)."""
        return abs(self.nM() * self.MMO2() - self.M_plateau_base()) < 1e-6

    def residue_label(self):
        from . import nodes as N
        if self.is_zrhf():
            return "MO2"
        r = self.node_def().get("residue")
        return r[0] if r else N.residue_label(self.metal)

    def plateau_excess(self, ratio):
        """W_plat/W_end ratio -> (organic + Cl mass at the plateau) per node, relative to the oxide-compensated node."""
        if self.residue_is_oxidic():
            return (ratio - 1) * self.nM() * self.MMO2()
        return ratio * self.nM() * self.MMO2() - self.M_plateau_base()

    def n_fw_metal(self):
        """mol framework metal per g from ICP."""
        if self.is_zrhf():
            nZ = (self.icp_Zr_wt or 0) / 100 / AW["Zr"]; nH_ = (self.icp_Hf_wt or 0) / 100 / AW["Hf"]
            return nH_ if self.extraframework_Zr else nZ + nH_
        return (self.icp_M_wt or 0) / 100 / AW[self.metal]

    def w_fw_metal(self):
        if self.is_zrhf():
            return (self.icp_Hf_wt or 0) + (0 if self.extraframework_Zr else (self.icp_Zr_wt or 0))
        return self.icp_M_wt or 0

    def fhf(self):
        """Framework Hf fraction."""
        if not self.is_zrhf():
            return 0.0
        if self.metal == "Zr":
            return 0.0
        if self.metal == "Hf":
            return 1.0
        if self.f_hf is not None:
            return self.f_hf
        if self.icp_Zr_wt and self.icp_Hf_wt:
            nz, nh = self.icp_Zr_wt / AW["Zr"], self.icp_Hf_wt / AW["Hf"]
            return nh / (nz + nh)
        return 0.0

    def MM(self):
        return metal_mass(self.fhf()) if self.is_zrhf() else AW[self.metal]

    def MMO2(self):
        """TGA residue per mol of framework metal (MO2 for Zr/Hf; CuO, Al2O3/2, Co3O4/3 ... otherwise)."""
        if self.is_zrhf():
            return oxide_mass(self.fhf())
        from . import nodes as N
        r = self.node_def().get("residue")
        return molar_mass(r[0]) / r[1] if r else N.residue_per_metal(self.metal)

    def M_symbol(self):
        if not self.is_zrhf():
            return self.metal
        f = self.fhf()
        if f <= 1e-6: return "Zr"
        if f >= 1 - 1e-6: return "Hf"
        return f"(Zr$_{{{1-f:.2f}}}$Hf$_{{{f:.2f}}}$)"

    def framework_oxide_fraction(self):
        """Fraction of the TGA residue that is framework MO2 (corrects for extra-framework ZrO2 in composites)."""
        if self.extraframework_Zr and self.icp_Zr_wt and self.icp_Hf_wt:
            m_zro2 = self.icp_Zr_wt / AW["Zr"] * oxide_mass(0.0)
            m_hfo2 = self.icp_Hf_wt / AW["Hf"] * oxide_mass(1.0)
            return m_hfo2 / (m_hfo2 + m_zro2)
        return 1.0

    def nmr_ratios(self):
        """Relative molar amounts a_i = A_i / H_i for all species with an integral."""
        out = {}
        for s in self.species:
            if s.integral is not None and s.nH:
                out[s.name] = s.integral / s.nH
        return out

    def cl_node(self, M_node=None):
        if self.cl_per_node is not None:
            return self.cl_per_node, f"XRF/EDX Cl per node ({self.nM()}·Cl/M)"
        if self.cl_wt_pct is not None and M_node:
            return self.cl_wt_pct / 100 * M_node / M_CL, f"Cl wt% × M_node / M_Cl (M_node = {M_node:.1f})"
        return 0.0, "not measured (assumed 0)"


@dataclass
class MethodResult:
    key: str
    title: str
    status: str = "n/a"            # 'ok' | 'n/a' | 'warn'
    semi_quant: bool = False
    L: float | None = None
    x: float | None = None
    m: float | None = None
    Cl: float | None = None
    OS: float | None = None
    formula: str | None = None
    formula_tex: str | None = None
    M_formula: float | None = None
    extra: dict = field(default_factory=dict)
    equations: list = field(default_factory=list)
    symbols: dict = field(default_factory=dict)
    assumptions: list = field(default_factory=list)
    steps: list = field(default_factory=list)
    pitfalls: list = field(default_factory=list)
    refs: list = field(default_factory=list)
    message: str = ""

    def row(self):
        return {"method": self.title, "key": self.key, "status": self.status,
                "semi-quantitative": self.semi_quant,
                "linkers/node L": _r(self.L), "missing linkers/node x": _r(self.x),
                "caps/node m": _r(self.m), "Cl/node": _r(self.Cl), "open sites/node": _r(self.OS),
                "M (g/mol per node)": _r(self.M_formula, 1), "formula": self.formula or "", "note": self.message}


def _r(v, n=2):
    return None if v is None or (isinstance(v, float) and math.isnan(v)) else round(float(v), n)


def _f(v, n=3):
    if v is None:
        return "-"
    return f"{v:.{n}g}" if abs(v) < 1e-3 or abs(v) >= 1e5 else f"{v:.{n}f}".rstrip("0").rstrip(".")


def open_sites(inp: DefectInput, L, m, cl, other_sites=None):
    """OS = (12 - CN) + t*x - (sites held by caps and secondary linkers) - Cl.
    other_sites = sum(z_i * n_i) over all non-main anionic species; defaults to m (monocarboxylate caps only)."""
    CN, t = inp.cn(), inp.t()
    occ = other_sites if other_sites is not None else (m or 0)
    return (inp.Q() - CN) + t * (inp.L0() - L) - occ - (cl or 0)


def _other_sites(inp, others: dict):
    sp = {s.name: s for s in inp.species}
    tot = 0.0
    for n, v in others.items():
        s = sp.get(n)
        if s is not None and s.z > 0 and v:
            tot += s.z * v
    return tot


# --------------------------------------------------------------------------- formula builder
def build_formula(inp: DefectInput, L, others: dict, cl: float, water_coord: bool = True, n_ads: float = 0.0):
    """Charge-balanced minimal formula on the M6O4(µ3-OH)4 node.
    others: {species_name: count per node} (caps, extra linkers, neutral guests).
    OH- fills the remaining charge (De Roo 2025 / Sannes 2023 convention), one H2O per terminal OH/Cl."""
    ml = inp.main_linker()
    sp = {s.name: s for s in inp.species}
    charge = ml.z * L + cl
    mass = inp.M_core() + L * ml.M + cl * M_CL
    parts = [(ml.name.split(" (")[0], L)]
    for nm, n in others.items():
        if n is None or n <= 1e-4:
            continue
        s = sp.get(nm)
        if s is None:
            f, z, *_ = SPECIES[nm]; M = molar_mass(f)
        else:
            z, M = s.z, s.M
        charge += z * n
        mass += n * M
        parts.append((nm.split(" (")[0], n))
    w = inp.Q() - charge
    mass += w * M_OH
    h2o = (max(w, 0) + cl) if water_coord else 0.0
    mass += (h2o + n_ads) * M_H2O
    Ms = inp.M_symbol()
    msym_plain = Ms.replace("$_{", "").replace("}$", "").replace("{", "").replace("}", "")
    plain = inp.node_def()["core"].format(M=msym_plain) + (f"(OH){w:.2f}" if abs(w) > 0.005 else "")
    for nm, n in parts:
        plain += f"({nm}){n:.2f}"
    if cl > 0.005:
        plain += f"(Cl){cl:.2f}"
    if h2o + n_ads > 0.005:
        plain += f"(H2O){h2o + n_ads:.2f}"
    return plain, mass, w, h2o


# --------------------------------------------------------------------------- method (a)
def method_a(inp: DefectInput) -> MethodResult:
    r = MethodResult("a", "(a) TGA plateau, M-400 (Valenzano 2011 / Shearer 2014)")
    r.refs = [REFS["valenzano"], REFS["shearer2014"], REFS["shearer2016"]]
    ml = inp.main_linker(); t = ml.z
    r.equations = [
        r"\mathrm{M_6O_{12-tL/2}(Linker)_L}\ \xrightarrow{\ \mathrm{air}\ }\ 6\,\mathrm{MO_2}",
        r"\frac{W_\mathrm{plat}}{W_\mathrm{end}} = \frac{6M_\mathrm{MO_2} + L\,\left(M_\mathrm{Lk} - \frac{t}{2}M_\mathrm{O}\right)}{6M_\mathrm{MO_2}}",
        r"L = \left(\frac{W_\mathrm{plat}}{W_\mathrm{end}} - 1\right)\frac{6M_\mathrm{MO_2}}{M_\mathrm{Lk} - \frac{t}{2}M_\mathrm{O}},\qquad x = L_0 - L",
        r"W_\mathrm{plat}^\mathrm{ideal} = 100\%\times\left(1 + \frac{L_0(M_\mathrm{Lk} - \frac{t}{2}M_\mathrm{O})}{6M_\mathrm{MO_2}}\right)\quad(\mathrm{UiO\text{-}66: 220.2\%})",
    ]
    r.symbols = {"W_plat": "mass at the dehydroxylated plateau (≈350–450 °C), any scale",
                 "W_end": "residual mass after combustion (MO2), same scale",
                 "M_MO2": "molar mass of ZrO2 (123.22), HfO2 (210.49) or the ICP-weighted mixture",
                 "M_Lk": "molar mass of the bound (deprotonated) linker", "t": "linker topicity (carboxylates per linker)",
                 "L": "linkers per M6 node", "L0": "ideal linkers per node (CN/t)", "x": "missing linkers per node"}
    r.assumptions = ["Node is fully dehydroxylated at the plateau: M6O4(OH)4 → M6O6 (Valenzano 2011).",
                     "Every charge left by a missing linker is compensated by O2-/OH- (no modulator survives the plateau).",
                     "End product is pure monoclinic MO2; no carbon residue, no extra-framework oxide.",
                     "No adsorbed water or solvent remains at the plateau temperature.",
                     f"Connectivity CN = {inp.cn()}, topicity t = {t}, ideal L0 = {inp.L0():g}."]
    r.pitfalls = ["Benzoate/TFA/sulfonated linkers may still be present at 400 °C → mass counted as linker → L overestimated (Sannes 2023).",
                  "Extra-framework oxide (e.g. ZrO2 nanoparticles on the MOF) inflates W_end → L underestimated. Use the ICP-corrected framework fraction.",
                  "No clear plateau (De Roo 2025 for un-calcined UiO-66): read W at the DTG minimum and report the temperature.",
                  "Chloride bound to the node adds mass at the plateau.",
                  "Under N2 instead of air, carbon residue raises W_end."]
    if not (inp.tga_W_plateau and inp.tga_W_end):
        r.message = "needs TGA plateau and end masses"; return r
    frac = inp.framework_oxide_fraction()
    Wend_fw = inp.tga_W_end * frac
    Wpl_fw = inp.tga_W_plateau - inp.tga_W_end * (1 - frac)   # extra-framework oxide mass is inert
    Mox = inp.MMO2(); dM = ml.M - t / 2 * M_O; nM = inp.nM()
    ratio = Wpl_fw / Wend_fw
    L = inp.plateau_excess(ratio) / dM
    Wid = (100 * (1 + inp.L0() * dM / (nM * Mox)) if inp.residue_is_oxidic()
           else 100 * (inp.M_plateau_base() + inp.L0() * dM) / (nM * Mox))
    r.L, r.x = L, inp.L0() - L
    cl, cls = inp.cl_node()
    r.m, r.Cl = 0.0, 0.0
    r.OS = open_sites(inp, L, 0, 0)
    f, M, w, h2o = build_formula(inp, L, {}, 0.0)
    r.formula, r.M_formula = f, M
    r.status = "ok"
    r.steps = [
        f"M_{{MO_2}} = {Mox:.2f}\\ \\mathrm{{g\\,mol^{{-1}}}}\\ (f_\\mathrm{{Hf}} = {inp.fhf():.3f}),\\quad M_\\mathrm{{Lk}} - \\tfrac{{{t}}}{{2}}M_\\mathrm{{O}} = {ml.M:.2f} - {t/2*M_O:.2f} = {dM:.2f}",
    ]
    if frac < 1:
        r.steps.append(f"W_\\mathrm{{end,fw}} = W_\\mathrm{{end}}\\times f_\\mathrm{{fw}} = {inp.tga_W_end:.2f}\\times{frac:.3f} = {Wend_fw:.2f}\\ (\\mathrm{{extra\\text{{-}}framework\\ ZrO_2\\ removed}})")
    r.steps += [
        f"W_\\mathrm{{plat}}/W_\\mathrm{{end}} = {Wpl_fw:.2f}/{Wend_fw:.2f} = {ratio:.4f}\\ \\ (\\mathrm{{ideal}}\\ {Wid/100:.4f})",
        (f"L = ({ratio:.4f} - 1)\\times\\frac{{{nM}\\times{Mox:.2f}}}{{{dM:.2f}}} = {L:.3f}" if inp.residue_is_oxidic() else
         f"L = \\frac{{{ratio:.4f}\\times{nM}\\times{Mox:.2f} - {inp.M_plateau_base():.2f}}}{{{dM:.2f}}} = {L:.3f}\\ (\\mathrm{{residue\\ {inp.residue_label()}}})"),
        f"x = {inp.L0():g} - {L:.3f} = {r.x:.3f}\\ \\mathrm{{missing\\ linkers\\ per\\ node}}",
    ]
    r.extra = {"W_plat_ideal_%": Wid, "ratio": ratio}
    return r


# --------------------------------------------------------------------------- method (c) NMR-corrected plateau
def method_c(inp: DefectInput) -> MethodResult:
    r = MethodResult("c", "(c) TGA plateau + 1H NMR modulator ratio (Shearer 2016; Atzori 2017)")
    r.refs = [REFS["shearer2016"], REFS["atzori2017"], REFS["sannes2023"]]
    ml = inp.main_linker(); t = ml.z
    r.equations = [
        r"r_i = \frac{n_i}{n_\mathrm{Lk}} = \frac{A_i/H_i}{A_\mathrm{Lk}/H_\mathrm{Lk}}",
        r"\frac{W_\mathrm{plat}}{W_\mathrm{end}}\,6M_\mathrm{MO_2} = 6M_\mathrm{MO_2} + L\left[\left(M_\mathrm{Lk}-\frac{t}{2}M_\mathrm{O}\right) + \sum_{i\in\mathrm{plat}} r_i\left(M_i - \frac{z_i}{2}M_\mathrm{O}\right)\right] + n_\mathrm{Cl}\left(M_\mathrm{Cl}-\frac{1}{2}M_\mathrm{O}\right)",
        r"m = L\sum_{i\in\mathrm{caps}} r_i",
    ]
    r.symbols = {"A_i, H_i": "1H integral of species i and the number of protons under it",
                 "r_i": "molar ratio species : main linker from the digested-MOF 1H NMR",
                 "z_i": "charge (carboxylates bound) of species i", "i ∈ plat": "species flagged as still present at the plateau",
                 "n_Cl": "Cl per node (if measured)", "m": "modulator caps per node"}
    r.assumptions = ["The NMR ratio (measured on the as-synthesised/activated MOF) still holds at the plateau for the species flagged 'at plateau'.",
                     "Species not flagged (formate, acetate, DMF) are gone before the plateau; their per-node amount is obtained from L × r_i.",
                     "Dehydroxylated node, charge compensated by O2- (as in method a).",
                     "Digestion does not lose or create any species (no volatile acids lost, no linker degradation)."]
    r.pitfalls = ["HF digestion followed by evaporation to dryness (as in the Paper 1 draft) removes HCOOH/AcOH → caps underestimated. Use 1 M NaOH/D2O, NH4HCO3 (De Roo 2025) or D2SO4/DMSO-d6 without evaporation.",
                  "Relaxation delay must be ≥ 5×T1 of the slowest signal (aromatic linker protons: use d1 ≥ 30–60 s), otherwise ratios are biased (De Roo 2025).",
                  "Sulfosuccinate CH/CH2 signals (2.8–3.8 ppm) can overlap with residual solvent; integrate cleanly resolved multiplets.",
                  "If a modulator survives 400 °C but is not flagged, L is overestimated."]
    ratios = inp.nmr_ratios()
    if not (inp.tga_W_plateau and inp.tga_W_end) or ml.name not in ratios:
        r.message = "needs TGA plateau/end masses and an NMR integral for the main linker"; return r
    aL = ratios[ml.name]
    sp = {s.name: s for s in inp.species}
    Mox = inp.MMO2(); frac = inp.framework_oxide_fraction()
    Wend = inp.tga_W_end * frac
    ratio = (inp.tga_W_plateau - inp.tga_W_end * (1 - frac)) / Wend
    denom = ml.M - t / 2 * M_O
    terms = []
    for nm, a in ratios.items():
        if nm == ml.name: continue
        s = sp[nm]
        if s.at_plateau and s.z >= 0:
            ri = a / aL
            denom += ri * (s.M - s.z / 2 * M_O)
            terms.append(f"{ri:.3f}\\times({s.M:.2f}-{s.z/2*M_O:.2f})")
    cl, cls = inp.cl_node()
    num = inp.plateau_excess(ratio) - cl * (M_CL - M_O / 2)
    L = num / denom
    others = {nm: L * a / aL for nm, a in ratios.items() if nm != ml.name}
    m = sum(v for nm, v in others.items() if sp[nm].role == "cap")
    r.L, r.x, r.m, r.Cl = L, inp.L0() - L, m, cl
    r.OS = open_sites(inp, L, m, cl, _other_sites(inp, others))
    f, M, w, _ = build_formula(inp, L, others, cl)
    r.formula, r.M_formula = f, M
    r.status = "ok" if w >= -0.05 else "warn"
    if w < -0.05:
        r.message = f"charge balance gives OH = {w:.2f} < 0: more anions than the node can hold - check integrals/flags"
    r.steps = [f"r_i\\ \\mathrm{{from\\ NMR}}:\\ " + ",\\ ".join(f"\\mathrm{{{nm.split(' (')[0]}}}: {a/aL:.3f}" for nm, a in ratios.items() if nm != ml.name),
               f"\\mathrm{{denominator}} = {ml.M - t/2*M_O:.2f}" + ("".join(" + " + s for s in terms)) + f" = {denom:.2f}",
               f"\\mathrm{{numerator}} = ({ratio:.4f}-1)\\times {inp.nM()}\\times{Mox:.2f} - {cl:.2f}\\times({M_CL:.2f}-{M_O/2:.2f}) = {num:.2f}",
               f"L = {num:.2f}/{denom:.2f} = {L:.3f},\\quad m = L\\sum r_\\mathrm{{cap}} = {m:.3f},\\quad x = {r.x:.3f}"]
    return r


# --------------------------------------------------------------------------- method (b) Sannes M-200 + missing cluster
def method_b(inp: DefectInput, L_ref: float | None = None) -> list[MethodResult]:
    out = []
    r = MethodResult("b1", "(b) M-200: TGA at 200 °C + NMR + inorganic counter-ions (Sannes/Øien-Ødegaard/Olsbye 2023)")
    r.refs = [REFS["sannes2023"], REFS["oien2014"], REFS["shearer2016"]]
    ml = inp.main_linker()
    r.equations = [
        r"M_\mathrm{MOF} = \frac{m_{200}}{m_{800}}\,6M_\mathrm{MO_2}",
        r"a_i = \frac{A_i}{H_i},\qquad b_j = 6\,\frac{B_j}{B_\mathrm{M}}",
        r"M_\mathrm{MOF} = M_\mathrm{clust} + \alpha\sum_i a_iM_i + \sum_j b_jM_j + \beta\,(M_\mathrm{OH^-} + M_\mathrm{H_2O})",
        r"12 = \alpha\sum_i a_iz_i + \sum_j b_jz_j + \beta",
        r"\alpha = \frac{M_\mathrm{MOF} - M_\mathrm{clust} - \sum_j b_jM_j - (12 - \sum_j b_jz_j)M_\mathrm{pair}}{\sum_i a_iM_i - M_\mathrm{pair}\sum_i a_iz_i},\qquad L = \alpha\,a_\mathrm{Lk},\quad x = L_0 - L",
    ]
    r.symbols = {"m_200, m_800": "TGA masses at 200 °C and after combustion", "M_clust": "M6O4(OH)4 node (12+), 6M_M + 8M_O + 4M_H",
                 "a_i": "relative molar amount of organic species i from NMR", "b_j": "inorganic counter-ion per node (Cl from EDX/XRF)",
                 "α": "scaling factor NMR → per node", "β": "number of OH-/H2O pairs per node", "M_pair": "M(OH-) + M(H2O) = 35.02 g/mol"}
    r.assumptions = ["At 200 °C the node is still M6O4(OH)4 (not dehydroxylated) and the pores are empty of physisorbed water/solvent "
                     "(or the remaining guests are quantified by NMR and included as neutral species, z = 0).",
                     "Every charge not balanced by carboxylates/Cl is balanced by one terminal OH- with one coordinated H2O (an OH-/H2O pair).",
                     "No missing-cluster defects (condition 2 of Sannes 2023).", "Full conversion to MO2 at 800 °C."]
    r.pitfalls = ["Requires the 200 °C mass on the same TGA trace (hold 30–60 min at 200 °C or read after the water step).",
                  "Physisorbed water not removed at 200 °C is counted as OH/H2O pairs (β too high → L slightly low).",
                  "Cl that is not quantified ends up as OH/H2O pairs (mass 35.0 vs 35.5 - small effect on mass, large on site counting)."]
    ratios = inp.nmr_ratios()
    if not (inp.tga_W_200 and inp.tga_W_end) or ml.name not in ratios:
        r.message = "needs TGA 200 °C + end masses and NMR integrals"
    else:
        Mox = inp.MMO2(); frac = inp.framework_oxide_fraction()
        # extra-framework ZrO2 (composites) is inert: remove its mass from both numerator and residue
        MMOF = (inp.tga_W_200 - inp.tga_W_end * (1 - frac)) / (inp.tga_W_end * frac) * inp.nM() * Mox
        Mclust = inp.M_clust(); Q = inp.Q()
        sp = {s.name: s for s in inp.species}
        cl, cls = inp.cl_node(MMOF)
        Mpair = M_OH + M_H2O
        SM = sum(a * sp[n].M for n, a in ratios.items())
        Sz = sum(a * sp[n].z for n, a in ratios.items())
        alpha = (MMOF - Mclust - cl * M_CL - (Q - cl) * Mpair) / (SM - Mpair * Sz)
        beta = Q - alpha * Sz - cl
        L = alpha * ratios[ml.name]
        others = {n: alpha * a for n, a in ratios.items() if n != ml.name}
        m = sum(v for n, v in others.items() if sp[n].role == "cap")
        r.L, r.x, r.m, r.Cl = L, inp.L0() - L, m, cl
        r.OS = open_sites(inp, L, m, cl, _other_sites(inp, others))
        f, M, w, _ = build_formula(inp, L, others, cl)
        r.formula, r.M_formula = f, M
        r.extra = {"M_MOF": MMOF, "alpha": alpha, "beta": beta}
        r.status = "ok" if beta >= -0.1 else "warn"
        if beta < -0.1:
            r.message = f"β = {beta:.2f} < 0 - inconsistent (pores not empty at 200 °C or species missing)"
        r.steps = [f"M_\\mathrm{{MOF}} = \\frac{{{inp.tga_W_200:.2f}}}{{{inp.tga_W_end*frac:.2f}}}\\times {inp.nM()}\\times{Mox:.2f} = {MMOF:.1f}\\ \\mathrm{{g\\,mol^{{-1}}}}",
                   f"M_\\mathrm{{clust}} = {Mclust:.2f},\\quad b_\\mathrm{{Cl}} = {cl:.2f}\\ ({cls})".replace("%", "\\%"),
                   f"\\sum a_iM_i = {SM:.2f},\\quad \\sum a_iz_i = {Sz:.3f}",
                   f"\\alpha = \\frac{{{MMOF:.1f} - {Mclust:.1f} - {cl*M_CL:.1f} - {Q-cl:.2f}\\times{Mpair:.2f}}}{{{SM:.2f} - {Mpair:.2f}\\times{Sz:.3f}}} = {alpha:.4f}",
                   f"\\beta = {Q} - {alpha*Sz:.3f} - {cl:.2f} = {beta:.3f}\\ \\mathrm{{OH^-/H_2O\\ pairs}}",
                   f"L = \\alpha a_\\mathrm{{Lk}} = {L:.3f},\\ m = {m:.3f},\\ x = {r.x:.3f}"]
    out.append(r)

    # ---- missing cluster interpretation
    q = MethodResult("b2", "(b') Missing-cluster interpretation (Cliffe 2014; Shearer 2016)")
    q.refs = [REFS["cliffe2014"], REFS["shearer2016"], REFS["oien2014"]]
    q.semi_quant = True
    q.equations = [r"L = L_0\,(1-\varphi)\quad\Rightarrow\quad \varphi = 1 - \frac{L}{L_0}",
                   r"m_\mathrm{pred} = CN\,\varphi\quad(\mathrm{one\ monocarboxylate\ cap\ per\ lost\ connection})",
                   r"\mathrm{reo\ check:}\ d_{100}=a,\ d_{110}=a/\sqrt{2},\quad 2\theta = 2\arcsin\left(\frac{\lambda}{2d}\right)"]
    q.symbols = {"φ": "fraction of missing nodes (randomly distributed nanoregions)", "m_pred": "caps per remaining node predicted by the model",
                 "a": "cubic cell parameter", "λ": "1.5406 Å (Cu Kα1)"}
    q.assumptions = ["Same TGA/NMR data as (a)/(b)/(c); composition per remaining node is identical to the capped missing-linker model, so composition alone cannot distinguish the two.",
                     "Discrimination needs diffraction: reo nanoregions give symmetry-forbidden (100)/(110) reflections at low angle (Cliffe 2014), plus pair-distribution function analysis or low-angle diffuse scattering.",
                     "A cap ratio m_meas/m_pred ≈ 1 is consistent with missing clusters or capped missing linkers; ≈ 0 indicates OH/H2O-compensated missing linkers."]
    q.pitfalls = ["For MOF-801 (Pn-3) and for small, broad crystallites the low-angle check is weak; do not over-interpret.",
                  "Missing clusters raise the pore volume much more per missing linker than missing linkers do (method f)."]
    Lr = L_ref if L_ref is not None else next((x.L for x in out if x.L is not None), None)
    if Lr is not None:
        phi = 1 - Lr / inp.L0()
        mpred = inp.cn() * phi
        a = inp.fw().get("a")
        tth = [2 * math.degrees(math.asin(1.5406 / (2 * d))) for d in (a, a / math.sqrt(2))] if a else [float("nan")] * 2
        mmeas = out[0].m if out[0].m is not None else None
        q.status = "ok"; q.L = Lr; q.x = inp.L0() - Lr; q.m = mpred; q.Cl = 0
        q.OS = open_sites(inp, Lr, mpred, 0)
        q.extra = {"phi_missing_cluster": phi, "m_pred": mpred, "m_meas": mmeas, "reo_2theta": tth}
        q.steps = [f"\\varphi = 1 - {Lr:.3f}/{inp.L0():g} = {phi:.3f}\\ ({100*phi:.1f}\\%\\ \\mathrm{{of\\ nodes\\ missing\\ if\\ all\\ defects\\ were\\ clusters}})",
                   f"m_\\mathrm{{pred}} = {inp.cn()}\\times{phi:.3f} = {mpred:.2f}" + (f",\\quad m_\\mathrm{{meas}} = {mmeas:.2f}\\ (\\mathrm{{ratio}}\\ {mmeas/mpred:.2f})" if mmeas is not None and mpred > 0.01 else ""),
                   (f"\\mathrm{{no\\ cubic\\ cell\\ for\\ {inp.framework}:\\ reo\\ check\\ n/a}}" if not a else f"\\mathrm{{reo\\ reflections\\ (Cu\\,K\\alpha_1,\\ a = {a}\\ \\AA)}}:\\ (100)\\ {tth[0]:.2f}^\\circ,\\ (110)\\ {tth[1]:.2f}^\\circ")]
        if inp.reo_peaks_seen is not None:
            q.message = "reo peaks observed → missing-cluster nanoregions likely" if inp.reo_peaks_seen else "no reo peaks → missing-linker picture preferred"
        q.formula, q.M_formula = None, None
    else:
        q.message = "needs a linker count from (a)/(b)/(c)"
    out.append(q)
    return out


# --------------------------------------------------------------------------- method (d) De Roo
def method_d(inp: DefectInput) -> MethodResult:
    r = MethodResult("d", "(d) Experimental minimal formula: qNMR + TGA + Cl (Pulparayil Mathew/De Roo 2025)")
    r.refs = [REFS["deroo2025"]]
    r.equations = [
        r"M_\mathrm{exp} = \frac{W_i}{W_\mathrm{MO_2}}\times 6M_\mathrm{MO_2}\qquad(\mathrm{Zr:}\ 6M_\mathrm{ZrO_2} = 739.3\ \mathrm{g\,mol^{-1}})",
        r"[\mathrm{MOF}] = \frac{m_\mathrm{dig}/M_\mathrm{exp}}{V},\qquad [i] = \frac{A_i/H_i}{A_\mathrm{std}/H_\mathrm{std}}\,[\mathrm{std}]",
        r"x_i = \frac{[i]}{[\mathrm{MOF}]},\qquad i_\mathrm{Cl} = \frac{[\mathrm{Cl^-}]}{[\mathrm{MOF}]}",
        r"24 = 4\cdot 2 + 4 + w + \sum_i z_ix_i + i_\mathrm{Cl}\quad\Rightarrow\quad w = 12 - \sum_i z_ix_i - i_\mathrm{Cl}",
        r"n_\mathrm{H_2O} = \frac{M_\mathrm{exp} - M_\mathrm{anhy}}{18.0},\qquad n_\mathrm{H_2O,coord} = w + i_\mathrm{Cl}",
    ]
    r.symbols = {"W_i, W_MO2": "initial TGA mass (RT, as weighed) and residue", "M_exp": "experimental molar mass per M6 at lab conditions (includes water)",
                 "m_dig, V": "MOF mass digested and NMR volume", "[std]": "TMSP-d4 concentration (9H)", "x_i": "species i per node",
                 "w": "terminal OH- per node from charge balance", "M_anhy": "molar mass of the formula without water", "n_H2O": "total water per formula"}
    r.assumptions = ["The TGA starts from the same (humid, as-weighed) state as the NMR sample - M_exp includes adsorbed water by design.",
                     "Residue is pure MO2 (checked by PXRD of the residue in the original paper).",
                     "Remaining charge is balanced by OH-, not O2- (oxide would leave under-coordinated M).",
                     "One coordinated H2O per monodentate OH-/Cl- ligand; all other water is physisorbed.",
                     "Quantitative NMR conditions (d1 ≥ 5 T1, no loss on digestion; 1 M NH4HCO3 or NaOH/D2O digest)."]
    r.pitfalls = ["Using the ideal molar mass to calculate [MOF] instead of M_exp gives wrong x_i (water uptake changes M by 10–20 %).",
                  "Chloride is invisible in NMR/IR; without the Cl test it appears as extra OH-.",
                  "HF digestion evaporated to dryness loses formate/acetate.",
                  "TGA and NMR must be done on the same batch, stored under the same conditions."]
    ratios = inp.nmr_ratios(); ml = inp.main_linker()
    have = inp.tga_W_end and inp.nmr_std_conc_mM and inp.nmr_sample_mg and inp.nmr_volume_mL and ml.name in ratios
    if not have:
        r.message = "needs TGA residue, internal-standard concentration, digested mass/volume and integrals"; return r
    Mox = inp.MMO2(); frac = inp.framework_oxide_fraction()
    Wi = inp.tga_W_initial or 100.0
    M_total = Wi / (inp.tga_W_end * frac) * inp.nM() * Mox          # weighed mass per mol of nodes (incl. any extra-framework ZrO2)
    Mexp = (Wi - inp.tga_W_end * (1 - frac)) / (inp.tga_W_end * frac) * inp.nM() * Mox   # framework only (= M_total if no ZrO2)
    # mg / (g mol-1) = mmol*1e-3 ... -> concentration in mM
    cMOF = inp.nmr_sample_mg / M_total / inp.nmr_volume_mL * 1000
    std_per_H = inp.nmr_std_integral / inp.nmr_std_nH
    conc = {n: (a / std_per_H) * inp.nmr_std_conc_mM for n, a in ratios.items()}
    xs = {n: c / cMOF for n, c in conc.items()}
    if inp.cl_conc_mM is not None:
        cl = inp.cl_conc_mM / cMOF; cls = "UV-vis [Cl-]/[MOF]"
    else:
        cl, cls = inp.cl_node(Mexp)
    sp = {s.name: s for s in inp.species}
    L = xs[ml.name]
    others = {n: v for n, v in xs.items() if n != ml.name}
    m = sum(v for n, v in others.items() if sp[n].role == "cap")
    f_anhy, Manhy, w, h2o = build_formula(inp, L, others, cl, water_coord=False)
    nw = (Mexp - Manhy) / M_H2O
    f_full, Mfull, _, _ = build_formula(inp, L, others, cl, water_coord=False, n_ads=nw)
    r.L, r.x, r.m, r.Cl = L, inp.L0() - L, m, cl
    r.OS = open_sites(inp, L, m, cl, _other_sites(inp, others))
    r.formula, r.M_formula = f_full, Mfull
    r.extra = {"M_exp": Mexp, "M_anhy": Manhy, "n_H2O_total": nw, "n_H2O_coord": max(w, 0) + cl, "w_OH": w, "[MOF]_mM": cMOF,
               "formula_anhydrous": f_anhy}
    r.status = "ok" if w >= -0.05 and nw >= -0.5 else "warn"
    if w < -0.05: r.message = f"w = {w:.2f} < 0: too many anions for {inp.Q():g}+ node (check [std], mass, integrals)"
    if nw < -0.5: r.message += f" negative water ({nw:.1f}): M_exp smaller than anhydrous formula"
    r.steps = [f"M_\\mathrm{{exp}} = \\frac{{{Wi:.1f}}}{{{inp.tga_W_end*frac:.2f}}}\\times {inp.nM()}\\times{Mox:.2f} = {Mexp:.1f}\\ \\mathrm{{g\\,mol^{{-1}}}}",
               f"[\\mathrm{{MOF}}] = \\frac{{{inp.nmr_sample_mg:.2f}\\ \\mathrm{{mg}}/{M_total:.1f}}}{{{inp.nmr_volume_mL:.3f}\\ \\mathrm{{mL}}}} = {cMOF:.3f}\\ \\mathrm{{mM}}",
               ",\\ ".join(f"[\\mathrm{{{n.split(' (')[0]}}}] = {conc[n]:.3f}\\ \\mathrm{{mM}} \\Rightarrow x = {xs[n]:.3f}" for n in xs),
               f"i_\\mathrm{{Cl}} = {cl:.3f}\\ ({cls})",
               f"w = {inp.Q()} - {sum(sp[n].z*v for n, v in xs.items()):.3f} - {cl:.3f} = {w:.3f}",
               f"M_\\mathrm{{anhy}} = {Manhy:.1f},\\quad n_\\mathrm{{H_2O}} = ({Mexp:.1f} - {Manhy:.1f})/18.02 = {nw:.2f}\\ (\\mathrm{{coord.}}\\ {max(w,0)+cl:.2f})"]
    return r


# --------------------------------------------------------------------------- method (e) ICP
def method_e(inp: DefectInput) -> list[MethodResult]:
    res = []
    r = MethodResult("e1", "(e) ICP-OES metal : qNMR linker molar ratio")
    r.refs = [REFS["deroo2025"], REFS["sannes2023"]]
    r.equations = [r"n_\mathrm{M} = \frac{w_\mathrm{Zr}}{100\,M_\mathrm{Zr}} + \frac{w_\mathrm{Hf}}{100\,M_\mathrm{Hf}}\ \ (\mathrm{mol\,g^{-1}}),\qquad f_\mathrm{Hf} = \frac{n_\mathrm{Hf}}{n_\mathrm{Zr}+n_\mathrm{Hf}}",
                   r"L = 6\,\frac{n_\mathrm{Lk}}{n_\mathrm{M}}\qquad(n_\mathrm{Lk}\ \mathrm{from\ qNMR\ in\ \mu mol\ per\ mg\ of\ the\ same\ sample})"]
    r.symbols = {"w_Zr, w_Hf": "metal mass fractions (wt%) from ICP-OES", "n_M": "mol framework metal per g", "n_Lk": "mol linker per g (qNMR, internal standard)"}
    r.assumptions = ["ICP and qNMR performed on the same batch in the same hydration state (or both corrected to dry mass).",
                     "Complete digestion for ICP (HF or microwave HNO3/HF); Zr/Hf emission lines free of spectral interference (Zr 343.823/339.198, Hf 277.336/264.141 nm).",
                     "For ZrO2@Hf-MOF composites only Hf is framework metal when 'extra-framework Zr' is ticked."]
    r.pitfalls = ["Zr and Hf are chemically almost identical - incomplete Hf/Zr separation in commercial HfCl4 (typically 0.2–2 % Zr) shifts f_Hf.",
                  "Sample weighing state (humid vs dried) does not matter for the ratio, but does for method e2."]
    n_fw = inp.n_fw_metal()
    if inp.linker_umol_per_mg and n_fw > 0:
        nL = inp.linker_umol_per_mg * 1e-6 * 1e3  # mol per g
        L = inp.nM() * nL / n_fw
        r.L, r.x = L, inp.L0() - L
        r.Cl = inp.cl_node()[0]
        r.status = "ok"
        r.OS = open_sites(inp, L, 0, r.Cl)
        r.message = "caps not determined by this method (OS assumes m = 0)"
        r.steps = [f"n_\\mathrm{{M}} = {n_fw*1e3:.4f}\\ \\mathrm{{mmol\\,g^{{-1}}}}\\ (f_\\mathrm{{Hf}} = {inp.fhf():.3f})",
                   f"n_\\mathrm{{Lk}} = {nL*1e3:.4f}\\ \\mathrm{{mmol\\,g^{{-1}}}}",
                   f"L = {inp.nM()}\\times{nL*1e3:.4f}/{n_fw*1e3:.4f} = {L:.3f}"]
    else:
        r.message = "needs ICP metal wt% and absolute qNMR linker content (µmol/mg)"
    res.append(r)

    q = MethodResult("e2", "(e') ICP mass balance: molar mass per node from metal wt% + NMR ratios")
    q.refs = [REFS["sannes2023"], REFS["deroo2025"]]
    q.equations = [r"M_\mathrm{node} = \frac{6\,\bar M_\mathrm{M}}{w_\mathrm{M}/100}",
                   r"M_\mathrm{node} = M_\mathrm{clust} + \alpha\sum_i a_iM_i + n_\mathrm{Cl}M_\mathrm{Cl} + \beta M_\mathrm{pair}\ (+\,n_\mathrm{w}M_\mathrm{H_2O}),\quad 12 = \alpha\sum_i a_iz_i + n_\mathrm{Cl} + \beta",
                   r"L = \alpha\,a_\mathrm{Lk}",
                   r"w_\mathrm{M,dry} = w_\mathrm{M,RT}\,\frac{W_i}{W_{200}}\quad(\mathrm{if\ weighed\ as\text{-}measured})"]
    q.symbols = {"w_M": "framework metal wt% (ICP)", "M̄_M": "ICP-weighted average metal mass", "W_i/W_200": "TGA mass ratio RT/200 °C (removes physisorbed water)"}
    q.assumptions = ["Same as M-200 (b) but the molar mass comes from ICP instead of TGA: an independent check of the TGA residue.",
                     "Sample must be weighed for ICP in the stated state: 'dried 200 °C' → no physisorbed water; 'as-measured' → water term unknown (use method d instead)."]
    q.pitfalls = ["Weighing errors on 2–5 mg ICP samples dominate (±2–3 % in w_M → ±0.2 linkers per node).",
                  "Physisorbed water in an 'as-measured' sample is interpreted as OH/H2O pairs."]
    ratios = inp.nmr_ratios(); ml = inp.main_linker()
    wM = inp.w_fw_metal()
    dry_corr = 1.0
    if "RT" in inp.icp_state and inp.tga_W_200 and inp.tga_W_initial:
        dry_corr = inp.tga_W_initial / inp.tga_W_200     # as-measured wt% -> dry (200 °C) basis via the TGA water loss
    if wM > 0 and ml.name in ratios:
        wM_dry = wM * dry_corr
        Mnode = inp.nM() * inp.MM() / (wM_dry / 100)
        Mclust = inp.M_clust(); Q = inp.Q()
        if inp.extraframework_Zr and inp.icp_Zr_wt:
            # remove ZrO2 mass carried per node
            zro2_per_g = inp.icp_Zr_wt / 100 / AW["Zr"] * oxide_mass(0.0)
            Mnode = Mnode * (1 - zro2_per_g)
        sp = {s.name: s for s in inp.species}
        cl, cls = inp.cl_node(Mnode)
        Mpair = M_OH + M_H2O
        SM = sum(a * sp[n].M for n, a in ratios.items()); Sz = sum(a * sp[n].z for n, a in ratios.items())
        alpha = (Mnode - Mclust - cl * M_CL - (Q - cl) * Mpair) / (SM - Mpair * Sz)
        beta = Q - alpha * Sz - cl
        L = alpha * ratios[ml.name]
        others = {n: alpha * a for n, a in ratios.items() if n != ml.name}
        m = sum(v for n, v in others.items() if sp[n].role == "cap")
        q.L, q.x, q.m, q.Cl = L, inp.L0() - L, m, cl
        q.OS = open_sites(inp, L, m, cl, _other_sites(inp, others))
        f, M, w, _ = build_formula(inp, L, others, cl)
        q.formula, q.M_formula = f, M
        q.status = "ok" if beta > -0.1 else "warn"
        if "RT" in inp.icp_state:
            if dry_corr != 1.0:
                q.message = f"as-measured ICP wt% converted to the 200 °C dry basis with the TGA water loss (×{dry_corr:.4f})"
            else:
                q.status = "warn"; q.message = "ICP sample weighed as-measured and no TGA 200 °C mass: water counted as OH/H2O pairs (unreliable)"
        q.extra = {"M_node": Mnode, "alpha": alpha, "beta": beta}
        q.steps = [f"w_\\mathrm{{M,dry}} = {wM:.3f}\\times{dry_corr:.4f} = {wM_dry:.3f}\\ \\mathrm{{wt\\%}}",
                   f"M_\\mathrm{{node}} = {inp.nM()}\\times{inp.MM():.2f}/({wM_dry:.3f}/100) = {Mnode:.1f}\\ \\mathrm{{g\\,mol^{{-1}}}}",
                   f"\\alpha = {alpha:.4f},\\ \\beta = {beta:.3f}",
                   f"L = {L:.3f},\\ m = {m:.3f},\\ x = {q.x:.3f}"]
    else:
        q.message = "needs ICP metal wt% and NMR integrals"
    res.append(q)
    return res


# --------------------------------------------------------------------------- method (f) porosity
def method_f(inp: DefectInput) -> MethodResult:
    r = MethodResult("f", "(f) Porosity-based estimate (pore volume / BET vs ideal) - semi-quantitative")
    r.semi_quant = True
    r.refs = [REFS["wu2013"], REFS["shearer2016"], REFS["vermoortele2013"]]
    ml = inp.main_linker()
    r.equations = [r"V(x) = \frac{V_0M_0 + x\,v_\mathrm{Lk}N_\mathrm{A}}{M_0 - x\,\Delta M},\qquad \Delta M = M_\mathrm{Lk} - t\,(M_\mathrm{OH^-}+M_\mathrm{H_2O})",
                   r"x = \frac{M_0\,(V - V_0)}{v_\mathrm{Lk}N_\mathrm{A} + V\,\Delta M}\qquad(\mathrm{same\ form\ with}\ S/S_0\ \mathrm{if\ only\ BET\ is\ known})"]
    r.symbols = {"V, V0": "measured and ideal (defect-free) pore volume, cm3 g-1", "M0": "ideal molar mass per node",
                 "v_Lk": "van der Waals volume freed per missing linker (Å3 → cm3 per mol via ×0.6022)",
                 "ΔM": "mass lost per missing linker replaced by t OH-/H2O pairs"}
    r.assumptions = ["Each missing linker frees its own van der Waals volume and the replacement OH/H2O occupy negligible volume.",
                     "Crystallinity and particle size identical to the ideal reference (no mesopores, no amorphous phase).",
                     "V0 from a geometric (probe-occupiable, N2 probe 1.82 Å) calculation on the defect-free crystal structure (Porosity helper: UiO-66(Zr) 0.42, MOF-801(Zr) 0.25, MOF-801(Hf) 0.18 cm3 g-1 from the bundled COD CIFs), or a defect-free reference measured on the same instrument.",
                     "For ZrO2 composites V is rescaled to the MOF mass using the ICP Zr content."]
    r.pitfalls = ["Missing clusters and interparticle/mesoporosity raise V far more than missing linkers → x overestimated.",
                  "Pore blocking by residual modulator/solvent or partial collapse lowers V → x underestimated or negative.",
                  "BET areas of microporous MOFs depend on the pressure range chosen - always apply the Rouquerol criteria (Sorption tab).",
                  "Treat as a trend indicator within one sample series only."]
    v = inp.v_linker_A3 or {"BDC (terephthalate)": 135.0, "Fumarate": 80.0, "TBAPy (NU-1000)": 560.0, "BTC (trimesate)": 150.0,
                             "NH2-BDC (2-aminoterephthalate)": 145.0, "BPDC (biphenyldicarboxylate, UiO-67)": 215.0,
                             "Isophthalate (CAU-10)": 135.0, "2-Methylimidazolate (ZIF-8/67)": 70.0,
                             "dobdc (MOF-74, carboxylate + phenolate)": 145.0}.get(ml.name, 120.0)
    # ideal molar mass per node (hydrated, OH/H2O-free ideal)
    M0 = inp.M_core() + inp.L0() * ml.M + (inp.Q() - inp.cn()) * (M_OH + M_H2O)
    dM = ml.M - ml.z * (M_OH + M_H2O)
    vmol = v * 0.6022  # cm3 per mol
    use = None
    if inp.V_meas and inp.V_ideal:
        use = ("V", inp.V_meas, inp.V_ideal)
    elif inp.S_meas and inp.S_ideal:
        use = ("S", inp.S_meas / inp.S_ideal * 1.0, 1.0)
    if use is None:
        r.message = "needs measured + ideal pore volume (or BET)"; return r
    kind, V, V0 = use
    zcorr = 1.0
    if inp.extraframework_Zr and inp.icp_Zr_wt:
        zcorr = 1 / (1 - inp.icp_Zr_wt / 100 / AW["Zr"] * oxide_mass(0.0))   # per g of MOF instead of per g of composite
        V = V * zcorr
    if kind == "S":
        # scale: treat S/S0 like V/V0 with V0 = 1 arbitrary units; need V0 in cm3/g for consistency -> use typical V0
        V0 = inp.V_ideal or 0.40
        V = use[1] * V0
    x = M0 * (V - V0) / (vmol + V * dM)
    L = inp.L0() - x
    r.L, r.x = L, x; r.m, r.Cl = 0.0, 0.0
    r.OS = open_sites(inp, L, 0, 0)
    r.status = "ok"
    r.steps = ([f"V_\\mathrm{{MOF}} = V_\\mathrm{{meas}}/(1-w_\\mathrm{{ZrO_2}}) = {V:.3f}\\ (\\times{zcorr:.3f})"] if zcorr != 1 else []) + [f"M_0 = {M0:.1f},\\ v_\\mathrm{{Lk}} = {v:.0f}\\ \\AA^3 = {vmol:.1f}\\ \\mathrm{{cm^3\\,mol^{{-1}}}},\\ \\Delta M = {dM:.2f}",
               f"x = \\frac{{{M0:.1f}\\times({V:.3f} - {V0:.3f})}}{{{vmol:.1f} + {V:.3f}\\times{dM:.2f}}} = {x:.3f}" + ("\\ (\\mathrm{from}\\ S/S_0)" if kind == "S" else "")]
    return r


# --------------------------------------------------------------------------- method (g) titration
def method_g(inp: DefectInput, M_node_fallback: float | None = None) -> MethodResult:
    r = MethodResult("g", "(g) Potentiometric acid-base titration: µ3-OH, M-OH2, M-OH (Klet 2016)")
    r.refs = [REFS["klet2016"], REFS["planas2014"]]
    r.equations = [r"n_\mathrm{\mu_3\text{-}OH} = c\,V_\mathrm{eq1} - n_\mathrm{HCl},\quad n_\mathrm{M\text{-}OH_2} = c\,(V_\mathrm{eq2}-V_\mathrm{eq1}),\quad n_\mathrm{M\text{-}OH} = c\,(V_\mathrm{eq3}-V_\mathrm{eq2})",
                   r"n_\mathrm{node} = \frac{m}{M_\mathrm{node}},\qquad p = \frac{n_\mathrm{M\text{-}OH_2}+n_\mathrm{M\text{-}OH}}{2\,n_\mathrm{node}}\ \ (\mathrm{OH/OH_2\ pairs\ per\ node})",
                   r"\mathrm{missing\ connections} = p - (12 - CN),\qquad L = \frac{12 - p - m - n_\mathrm{Cl}}{t}",
                   r"\mathrm{p}K_a = \mathrm{pH}\ \mathrm{at}\ V_\mathrm{eq}/2\ \ (\mathrm{UiO\text{-}66(Zr):}\ 3.52,\ 6.79,\ 8.30;\ \mathrm{Hf:}\ 3.37,\ 6.67,\ 8.14;\ \mathrm{NU\text{-}1000:}\ 3.59,\ 5.75,\ 8.2)"]
    r.symbols = {"c": "NaOH concentration (0.1 M)", "V_eq1..3": "equivalence volumes (first-derivative maxima)",
                 "n_HCl": "HCl added to bring the suspension to pH 3 (subtract from eq. 1)", "p": "terminal OH/OH2 pairs per node",
                 "m, n_Cl": "caps and Cl per node (they occupy sites but are not titrated)"}
    r.assumptions = ["Each lost carboxylate is replaced by one M-OH2 + one M-OH (Klet 2016): two titratable protons per open site.",
                     "Klet's 'missing linkers of ideally 12' counts missing connections per node; here x = missing connections / t.",
                     "Equivalence points well separated (fit first derivative with Lorentzians); sample in 0.01 M NaNO3, equilibrated 18 h.",
                     "No other titratable protons (free -COOH of dangling linkers, -SO3H, ZrO2 surface OH) - otherwise they add to p."]
    r.pitfalls = ["Sulfonated samples and ZrO2-decorated composites contain extra acidic protons (SO3H pKa < 2, ZrO2 surface OH) → p overestimated.",
                  "Linker dissolution above pH 9–10 releases protons/consumes base; stop at pH 10.5.",
                  "For NU-1000 the pKa2/pKa3 steps overlap - only the total base consumption is reliable (Klet 2016).",
                  "Converting mass → mol needs a realistic M_node (use M_exp from method d, not the ideal)."]
    V = inp.tit_V_eq
    if not (inp.tit_mass_mg and inp.tit_c_NaOH and V and len(V) >= 2):
        r.message = "needs titrated mass, NaOH concentration and ≥ 2 equivalence volumes"; return r
    Mn = inp.tit_M_node or M_node_fallback
    if Mn is None:
        ml = inp.main_linker()
        Mn = inp.M_core() + inp.L0() * ml.M
    n_node = inp.tit_mass_mg / Mn  # mmol
    c = inp.tit_c_NaOH
    V = sorted(V)
    if len(V) == 2:
        n_def = c * (V[1] - V[0]); n_oh2 = n_oh = n_def / 2
    else:
        n_oh2, n_oh = c * (V[1] - V[0]), c * (V[2] - V[1]); n_def = n_oh2 + n_oh
    n_mu3 = (c * V[0] - inp.tit_HCl_mmol) if inp.tit_HCl_mmol is not None else None
    p = n_def / (2 * n_node)
    cl = inp.cl_node(Mn)[0]
    ratios = inp.nmr_ratios(); ml = inp.main_linker()
    m = 0.0
    L = (inp.Q() - p - m - cl) / ml.z
    r.L, r.x, r.m, r.Cl = L, inp.L0() - L, m, cl
    r.OS = p
    r.status = "ok"
    r.extra = {"n_node_mmol": n_node, "OH2_per_node": n_oh2 / n_node, "OH_per_node": n_oh / n_node,
               "mu3OH_per_node": (n_mu3 / n_node) if n_mu3 is not None else None,
               "missing_connections_Klet": p - (inp.Q() - inp.cn())}
    r.message = "open sites measured directly (= p); caps not seen by titration (m = 0 assumed for L)"
    r.steps = [f"n_\\mathrm{{node}} = {inp.tit_mass_mg:.2f}\\ \\mathrm{{mg}}/{Mn:.1f} = {n_node*1e3:.3f}\\ \\mathrm{{\\mu mol}}",
               f"n_\\mathrm{{M\\text{{-}}OH_2}} = {c}\\times({V[1]:.3f}-{V[0]:.3f}) = {n_oh2*1e3:.3f}\\ \\mathrm{{\\mu mol}},\\quad n_\\mathrm{{M\\text{{-}}OH}} = {n_oh*1e3:.3f}\\ \\mathrm{{\\mu mol}}",
               (f"n_\\mathrm{{\\mu_3\\text{{-}}OH}}/n_\\mathrm{{node}} = {n_mu3/n_node:.2f}\\ (\\mathrm{{ideal}}\\ 4)" if n_mu3 is not None else "\\mu_3\\text{-}OH\\ \\mathrm{not\\ quantified\\ (enter\\ HCl\\ added)}"),
               f"p = \\frac{{{n_def*1e3:.3f}}}{{2\\times{n_node*1e3:.3f}}} = {p:.3f}\\ \\mathrm{{pairs/node;\\ missing\\ connections}} = {p-(inp.Q()-inp.cn()):.3f}",
               f"L = ({inp.Q()} - {p:.3f} - {m:.2f} - {cl:.2f})/{ml.z} = {L:.3f}"]
    return r


# --------------------------------------------------------------------------- method (h) and consensus
def method_h_text():
    return {
        "title": "(h) Charge-balance open-site rule (applied to every method)",
        "equations": [r"OS = (12 - CN) + t\,x - m - n_\mathrm{Cl}\qquad\left(\mathrm{general:}\ m \rightarrow \sum_{i\neq\mathrm{Lk}} z_in_i\right)",
                      r"\mathrm{UiO\text{-}66,\ MOF\text{-}801}\ (CN = 12,\ t = 2):\quad OS = 2x - m - n_\mathrm{Cl}",
                      r"\mathrm{NU\text{-}1000}\ (CN = 8,\ t = 4):\quad OS = 4 + 4x - m - n_\mathrm{Cl}"],
        "symbols": {"OS": "open M sites per node terminated by OH-/OH2 (catalytically relevant Lewis/Brønsted pairs)",
                    "x": "missing linkers per node", "m": "monocarboxylate caps per node", "n_Cl": "Cl per node"},
        "assumptions": ["Each missing ditopic linker removes one carboxylate from each of two nodes → t vacated sites per node per missing linker (per-node formula basis).",
                        "Each cap or Cl- occupies exactly one site; everything else is an OH-/OH2 pair (open site).",
                        "Dangling (monodentate) linkers count as caps if they are entered with role 'cap'.",
                        "Secondary linkers / multidentate modulators (e.g. 2-sulfosuccinate bound through both carboxylates) occupy z sites each: the generalised rule subtracts Σ z_i n_i."],
        "refs": [REFS["klet2016"], REFS["sannes2023"], REFS["deroo2025"]],
    }


SUMMARY_PITFALLS = [
    ("(a) ≫ (c)/(b)", "Modulator or sulfonated species survive 400 °C: (a) counts their mass as linker. Trust (b)/(c)/(d)."),
    ("(a), (b) ≪ (e1)", "Extra-framework oxide in the TGA residue (ZrO2 NPs!) or carbon residue under N2: tick 'extra-framework Zr' and use air."),
    ("(d) ≠ (b)", "Different hydration states of the TGA and NMR samples, or [MOF] computed with the wrong molar mass."),
    ("(f) ≫ others", "Missing clusters, mesopores or interparticle porosity; (f) is semi-quantitative only."),
    ("(f) < 0", "Pore blocking (residual DMF/modulator) or partial amorphisation - check PXRD and activation."),
    ("(g) OS ≫ (h) OS", "Extra titratable protons (SO3H, ZrO2 surface OH, free COOH) or base-induced dissolution."),
    ("caps in (c) ≈ 0 but (b) β small", "Formate lost during HF digestion/evaporation - switch digestion method."),
]


def run_all(inp: DefectInput) -> dict:
    a = method_a(inp)
    c = method_c(inp)
    b = method_b(inp)
    ref_L = next((x.L for x in (b[0], c, a) if x.L is not None), None)
    b = method_b(inp, ref_L)
    d = method_d(inp)
    e = method_e(inp)
    f = method_f(inp)
    Mn = d.extra.get("M_exp") if d.status != "n/a" else (b[0].extra.get("M_MOF") if b[0].status != "n/a" else None)
    g = method_g(inp, Mn)
    methods = [a, b[0], b[1], c, d, e[0], e[1], f, g]
    nd = inp.node_def()
    if nd["core"] != "{M}6O4(OH)4" or not inp.is_zrhf():
        from . import nodes as N
        note = (f"Node {nd['core'].format(M=inp.M_symbol())}: equations are written for M6O4(OH)4; here 6 → n_M = {inp.nM()}, "
                f"12 → node charge for ligands Q = {inp.Q():g}, MO2 → TGA residue {inp.residue_label()} "
                f"({inp.MMO2():.2f} g/mol per metal), M6O4(OH)4 → {nd['core'].format(M=inp.M_symbol())}.")
        for m_ in methods:
            m_.assumptions = [note] + list(m_.assumptions)
    quant = [m for m in methods if m.status in ("ok", "warn") and not m.semi_quant and m.L is not None]
    cons = {}
    if quant:
        Ls = np.array([m.L for m in quant])
        cons["L"] = float(np.median(Ls)); cons["L_min"] = float(Ls.min()); cons["L_max"] = float(Ls.max())
        cons["L_sd"] = float(Ls.std(ddof=1)) if len(Ls) > 1 else 0.0
        cons["n_methods"] = len(Ls)
        cons["x"] = inp.L0() - cons["L"]
        ms = [m.m for m in quant if m.m is not None and m.key in ("b1", "c", "d", "e2")]
        cons["m"] = float(np.median(ms)) if ms else 0.0
        cls = [m.Cl for m in quant if m.Cl is not None]
        cons["Cl"] = float(np.median(cls)) if cls else 0.0
        ratios_ = inp.nmr_ratios(); ml_ = inp.main_linker()
        oth_ = {n: cons["L"] * a_ / ratios_[ml_.name] for n, a_ in ratios_.items() if n != ml_.name} if ml_.name in ratios_ else {}
        cons["other_sites"] = _other_sites(inp, oth_) if oth_ else cons["m"]
        cons["OS"] = open_sites(inp, cons["L"], cons["m"], cons["Cl"], cons["other_sites"])
        OSs = [m.OS for m in methods if m.OS is not None and m.status in ("ok", "warn")]
        cons["OS_min"], cons["OS_max"] = (min(OSs), max(OSs)) if OSs else (None, None)
        # consensus composition: use cap/guest ratios from NMR
        ratios = inp.nmr_ratios(); ml = inp.main_linker()
        others = {}
        if ml.name in ratios:
            for n, a_ in ratios.items():
                if n != ml.name:
                    others[n] = cons["L"] * a_ / ratios[ml.name]
        f_, M_, w_, h_ = build_formula(inp, cons["L"], others, cons["Cl"])
        cons["formula"], cons["M"] = f_, M_
        f2, M2, _, _ = build_formula(inp, cons["L"], others, cons["Cl"], water_coord=False)
        cons["formula_anhydrous"], cons["M_anhydrous"] = f2, M2
        cons["metal_wt_pct"] = 100 * inp.nM() * inp.MM() / M_
    return {"methods": methods, "consensus": cons, "h": method_h_text(), "input": inp}


def results_table(res) -> "pd.DataFrame":
    import pandas as pd
    rows = [m.row() for m in res["methods"]]
    c = res["consensus"]
    if c:
        rows.append({"method": f"CONSENSUS (median of {c['n_methods']} quantitative methods)", "key": "consensus", "status": "ok",
                     "semi-quantitative": False, "linkers/node L": _r(c["L"]), "missing linkers/node x": _r(c["x"]),
                     "caps/node m": _r(c["m"]), "Cl/node": _r(c["Cl"]), "open sites/node": _r(c["OS"]),
                     "M (g/mol per node)": _r(c["M"], 1), "formula": c["formula"],
                     "note": f"L range {c['L_min']:.2f}–{c['L_max']:.2f}, SD {c['L_sd']:.2f}"})
    return pd.DataFrame(rows)
