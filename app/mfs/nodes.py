"""v2.11 node (SBU) library for the defect calculator and formula tools.

A node is described by
    nM      metal atoms per node (the formula basis: every count in the calculator is "per node")
    val     oxidation state of the metal in the node
    core_O  µ-O2- in the node core,  core_OH  µ-OH- in the node core
    sites   = nM*val - 2*core_O - core_OH   charge left for linkers/caps/Cl/terminal OH (12 for M6O4(OH)4)
    CN      ideal connectivity in charge units (= donor groups for carboxylates; 2 for a ZIF Zn: 4 N donors, 2 imidazolates):
            sites - CN terminal ligands remain in the ideal framework (MIL-100: 7 - 6 = one F-/OH- per trimer)
    plateau_dehydroxylated  (default True) the core µ-OH is lost before the TGA plateau (M6O4(OH)4 -> M6O6); False for
            Al-OH chains whose µ-OH survives to linker combustion
The TGA residue in air is looked up per metal (RESIDUE): ZrO2, CeO2, CuO, ZnO, Al2O3, Fe2O3, Co3O4, ...
Everything is plain data so a custom node is just another dict (see `custom_node`)."""
from __future__ import annotations
from .chem import AW, molar_mass, M_O, M_OH

# residue after combustion in air: oxide formula and metal atoms per formula unit
RESIDUE = {"Zr": ("ZrO2", 1), "Hf": ("HfO2", 1), "Ce": ("CeO2", 1), "Th": ("ThO2", 1), "U": ("U3O8", 3), "Ti": ("TiO2", 1),
           "Cu": ("CuO", 1), "Zn": ("ZnO", 1), "Co": ("Co3O4", 3), "Ni": ("NiO", 1), "Mg": ("MgO", 1), "Mn": ("Mn2O3", 2),
           "Fe": ("Fe2O3", 2), "Cr": ("Cr2O3", 2), "Al": ("Al2O3", 2), "Sc": ("Sc2O3", 2), "In": ("In2O3", 2), "Ga": ("Ga2O3", 2),
           "V": ("V2O5", 2), "Y": ("Y2O3", 2), "La": ("La2O3", 2), "Eu": ("Eu2O3", 2), "Tb": ("Tb4O7", 4), "Ca": ("CaO", 1),
           "Sr": ("SrO", 1), "Ba": ("BaCO3", 1), "Cd": ("CdO", 1), "Ag": ("Ag", 1), "Bi": ("Bi2O3", 2), "Mo": ("MoO3", 1)}


def residue_per_metal(metal: str) -> float:
    """g of air-TGA residue per mol of metal (ZrO2 -> 123.22, Fe -> Fe2O3/2 = 79.84, Co -> Co3O4/3 = 80.31)."""
    f, n = RESIDUE.get(metal, (f"{metal}O", 1))
    return molar_mass(f) / n


def residue_label(metal: str) -> str:
    return RESIDUE.get(metal, (f"{metal}O", 1))[0]


# key -> definition. 'core' is the plain node formula with {M} for the metal symbol.
NODES = {
    "M6O4(OH)4 (Zr6 / Hf6 / Ce6 / Th6)": dict(nM=6, val=4, core_O=4, core_OH=4, CN=12, metals=["Zr", "Hf", "mixed", "Ce", "Th"],
                                            core="{M}6O4(OH)4", examples="UiO-66/67, MOF-801, NU-1000 (CN 8), MOF-808 (CN 6), DUT-67"),
    "Cu2 paddlewheel": dict(nM=2, val=2, core_O=0, core_OH=0, CN=4, metals=["Cu", "Zn", "Co", "Ni", "Mo"], core="{M}2",
                            examples="HKUST-1 (BTC), MOF-505, PCN-14; axial sites hold H2O"),
    "Zn2 paddlewheel (pillared)": dict(nM=2, val=2, core_O=0, core_OH=0, CN=4, metals=["Zn", "Cu", "Co", "Ni"], core="{M}2",
                                       examples="DMOF-1 Zn2(BDC)2(dabco): enter dabco as a neutral species"),
    "Zn4O": dict(nM=4, val=2, core_O=1, core_OH=0, CN=6, metals=["Zn", "Co", "Be"], core="{M}4O", examples="MOF-5/IRMOF-n, MOF-177"),
    "M(OH) chain (Al-OH, MIL-53 / CAU-10)": dict(nM=1, val=3, core_O=0, core_OH=1, CN=2, metals=["Al", "Fe", "Cr", "Ga", "In", "Sc", "V"],
                                              core="{M}(OH)", plateau_dehydroxylated=False, examples="MIL-53(Al), CAU-10, A520 Al-fumarate, MIL-118 (per metal)"),
    "M3O trimer (MIL-100 / MIL-101 / MIL-88)": dict(nM=3, val=3, core_O=1, core_OH=0, CN=6, metals=["Fe", "Cr", "Al", "Sc", "V"],
                                                 core="{M}3O", examples="MIL-100 (BTC), MIL-101 (BDC), MIL-88; 1 terminal F-/OH- per trimer"),
    "M2 MOF-74 / CPO-27": dict(nM=2, val=2, core_O=0, core_OH=0, CN=4, metals=["Mg", "Co", "Ni", "Zn", "Mn", "Fe", "Cu"], core="{M}2",
                               examples="M2(dobdc): dobdc4- binds through 2 carboxylates + 2 phenolates (z = 4); open M sites"),
    "Ti8O8(OH)4 (MIL-125)": dict(nM=8, val=4, core_O=8, core_OH=4, CN=12, metals=["Ti"], core="{M}8O8(OH)4", examples="MIL-125, NH2-MIL-125"),
    "M (ZIF, tetrahedral)": dict(nM=1, val=2, core_O=0, core_OH=0, CN=2, metals=["Zn", "Co", "Cd"], core="{M}",
                                 examples="ZIF-8 (Zn, 2-mIm), ZIF-67 (Co); CN counted in charge units: 4 N donors = 2 imidazolates (z = 1), L0 = 2"),
}
DEFAULT_NODE = "M6O4(OH)4 (Zr6 / Hf6 / Ce6 / Th6)"
CUSTOM = "custom node…"


def sites(nd: dict) -> int | float:
    return nd["nM"] * nd["val"] - 2 * nd.get("core_O", 0) - nd.get("core_OH", 0)


def custom_node(nM=1, val=2, core_O=0, core_OH=0, CN=None, metal="Zn", core=None, residue=None):
    nd = dict(nM=int(nM), val=val, core_O=core_O, core_OH=core_OH, metals=[metal], custom=True, examples="custom",
              core=core or ("{M}" + (str(int(nM)) if int(nM) > 1 else "") + (f"O{core_O:g}" if core_O else "") + (f"(OH){core_OH:g}" if core_OH else "")))
    nd["CN"] = CN if CN is not None else sites(nd)
    if residue:
        nd["residue"] = residue           # (formula, metal atoms per formula)
    return nd


def get(key: str | None, custom: dict | None = None) -> dict:
    if key == CUSTOM and custom:
        return custom if "nM" in custom and "core" in custom else custom_node(**custom)
    return NODES.get(key or DEFAULT_NODE, NODES[DEFAULT_NODE])


def core_mass(nd: dict, MM: float) -> float:
    """Molar mass of the node core (M6O4(OH)4 = 6M + 4O + 4OH)."""
    return nd["nM"] * MM + nd.get("core_O", 0) * M_O + nd.get("core_OH", 0) * M_OH


def table_rows():
    rows = []
    for k, nd in NODES.items():
        ms = [m for m in nd["metals"] if m != "mixed"]
        rows.append({"node": k, "metals": ", ".join(nd["metals"]), "metals/node": nd["nM"], "oxidation state": nd["val"],
                     "core": nd["core"].format(M="M"), "charge for linkers (sites)": sites(nd), "ideal CN": nd["CN"],
                     "TGA residue": ", ".join(f"{m}→{residue_label(m)}" for m in ms), "examples": nd["examples"],
                     f"M core, {ms[0]} (g/mol)": round(core_mass(nd, AW[ms[0]]), 2)})
    return rows
