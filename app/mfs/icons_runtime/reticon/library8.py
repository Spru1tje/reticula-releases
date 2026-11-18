"""
v5 batch 4 (generated sets): periodic-table tiles (118 elements, family colours, hand-made variants),
small molecules (solvents, modulators, gases, linkers) and topology nets. All original (CC0).
"""
from __future__ import annotations

import copy
import math

from .core import Icon, Op, polar
from .library import icon, REGISTRY
from .generators import icon_from_molecule
from . import variants as V

# ----------------------------------------------------------------- periodic table
ELEMENTS = """1 H Hydrogen nonmetal|2 He Helium noble|3 Li Lithium alkali|4 Be Beryllium alkaline|5 B Boron metalloid|
6 C Carbon nonmetal|7 N Nitrogen nonmetal|8 O Oxygen nonmetal|9 F Fluorine halogen|10 Ne Neon noble|11 Na Sodium alkali|
12 Mg Magnesium alkaline|13 Al Aluminium post|14 Si Silicon metalloid|15 P Phosphorus nonmetal|16 S Sulfur nonmetal|
17 Cl Chlorine halogen|18 Ar Argon noble|19 K Potassium alkali|20 Ca Calcium alkaline|21 Sc Scandium transition|
22 Ti Titanium transition|23 V Vanadium transition|24 Cr Chromium transition|25 Mn Manganese transition|26 Fe Iron transition|
27 Co Cobalt transition|28 Ni Nickel transition|29 Cu Copper transition|30 Zn Zinc transition|31 Ga Gallium post|
32 Ge Germanium metalloid|33 As Arsenic metalloid|34 Se Selenium nonmetal|35 Br Bromine halogen|36 Kr Krypton noble|
37 Rb Rubidium alkali|38 Sr Strontium alkaline|39 Y Yttrium transition|40 Zr Zirconium transition|41 Nb Niobium transition|
42 Mo Molybdenum transition|43 Tc Technetium transition|44 Ru Ruthenium transition|45 Rh Rhodium transition|
46 Pd Palladium transition|47 Ag Silver transition|48 Cd Cadmium transition|49 In Indium post|50 Sn Tin post|
51 Sb Antimony metalloid|52 Te Tellurium metalloid|53 I Iodine halogen|54 Xe Xenon noble|55 Cs Caesium alkali|
56 Ba Barium alkaline|57 La Lanthanum lanthanide|58 Ce Cerium lanthanide|59 Pr Praseodymium lanthanide|
60 Nd Neodymium lanthanide|61 Pm Promethium lanthanide|62 Sm Samarium lanthanide|63 Eu Europium lanthanide|
64 Gd Gadolinium lanthanide|65 Tb Terbium lanthanide|66 Dy Dysprosium lanthanide|67 Ho Holmium lanthanide|
68 Er Erbium lanthanide|69 Tm Thulium lanthanide|70 Yb Ytterbium lanthanide|71 Lu Lutetium lanthanide|
72 Hf Hafnium transition|73 Ta Tantalum transition|74 W Tungsten transition|75 Re Rhenium transition|76 Os Osmium transition|
77 Ir Iridium transition|78 Pt Platinum transition|79 Au Gold transition|80 Hg Mercury transition|81 Tl Thallium post|
82 Pb Lead post|83 Bi Bismuth post|84 Po Polonium post|85 At Astatine halogen|86 Rn Radon noble|87 Fr Francium alkali|
88 Ra Radium alkaline|89 Ac Actinium actinide|90 Th Thorium actinide|91 Pa Protactinium actinide|92 U Uranium actinide|
93 Np Neptunium actinide|94 Pu Plutonium actinide|95 Am Americium actinide|96 Cm Curium actinide|97 Bk Berkelium actinide|
98 Cf Californium actinide|99 Es Einsteinium actinide|100 Fm Fermium actinide|101 Md Mendelevium actinide|
102 No Nobelium actinide|103 Lr Lawrencium actinide|104 Rf Rutherfordium transition|105 Db Dubnium transition|
106 Sg Seaborgium transition|107 Bh Bohrium transition|108 Hs Hassium transition|109 Mt Meitnerium unknown|
110 Ds Darmstadtium unknown|111 Rg Roentgenium unknown|112 Cn Copernicium unknown|113 Nh Nihonium unknown|
114 Fl Flerovium unknown|115 Mc Moscovium unknown|116 Lv Livermorium unknown|117 Ts Tennessine unknown|118 Og Oganesson unknown"""

FAMILY = {"alkali": ("#F7A072", "alkali metal"), "alkaline": ("#F9D27A", "alkaline earth metal"),
          "transition": ("#9CC3E6", "transition metal"), "post": ("#B7D7A8", "post-transition metal"),
          "metalloid": ("#A8DCCF", "metalloid"), "nonmetal": ("#C9E8A0", "nonmetal"), "halogen": ("#F5E28A", "halogen"),
          "noble": ("#C8B6E8", "noble gas"), "lanthanide": ("#F4B6C8", "lanthanide"), "actinide": ("#EBA9A9", "actinide"),
          "unknown": ("#DADFE5", "superheavy element")}
INK = "#1F2328"


def _parse():
    out = []
    for item in ELEMENTS.replace("\n", "").split("|"):
        n, sym, name, fam = item.split()
        out.append((int(n), sym, name, fam))
    return out


def periodic_position(n: int) -> tuple[int, float]:
    """(group 1-18, row) for the picker grid. Lanthanides/actinides sit in rows 9.4 / 10.4 under the table."""
    if 57 <= n <= 71:
        return n - 57 + 3, 9.4
    if 89 <= n <= 103:
        return n - 89 + 3, 10.4
    if n == 1:
        return 1, 1
    if n == 2:
        return 18, 1
    for lo, hi, p in ((3, 10, 2), (11, 18, 3), (19, 36, 4), (37, 54, 5), (55, 86, 6), (87, 118, 7)):
        if lo <= n <= hi:
            i = n - lo
            if p in (2, 3):
                return (i + 1 if i < 2 else i + 11), p
            if p in (4, 5):
                return i + 1, p
            return (i + 1 if i < 2 else i - 13), p      # 72 (Hf) / 104 (Rf) -> group 4
    raise ValueError(n)


def _tile(ic, n, sym, fam, number=True, outline=False):
    col = FAMILY[fam][0]
    if outline:
        ic.rect(6, 6, 52, 52, 7).fixed(col)
        ic.cut().rect(10, 10, 44, 44, 4)
    else:
        ic.rect(6, 6, 52, 52, 7).fixed(col)
    size = 26 if len(sym) == 1 else 22
    role = "m" if outline else "x"          # on a hollow tile the letters must be body ink, not a cut
    ic.text(32, 44 if number else 41, sym, size, "bold", role=role, fill=INK)
    if outline:
        ic.fixed(INK)
    if number:
        ic.text(14 if outline else 13, 21 if outline else 20, str(n), 9, "bold", anchor="start", role=role, fill=INK)
        if outline:
            ic.fixed(INK)


def _make(n, sym, name, fam):
    nm = f"element_{sym.lower()}"

    @icon(nm, "elements", label=f"{sym} – {name}", tags=[name.lower(), sym, f"element {n}", FAMILY[fam][1], "periodic table",
                f"group {periodic_position(n)[0]}", f"row {periodic_position(n)[1]}"])
    def _(ic):
        _tile(ic, n, sym, fam)

    def plain(base, n=n, sym=sym, fam=fam):
        ic = Icon(base.name, f"{base.label} (symbol only)", base.category, base.tags)
        _tile(ic, n, sym, fam, number=False)
        return ic

    def outline(base, n=n, sym=sym, fam=fam):
        ic = Icon(base.name, f"{base.label} (outline)", base.category, base.tags)
        _tile(ic, n, sym, fam, outline=True)
        return ic
    def dark(base, n=n, sym=sym, fam=fam):
        ic = Icon(base.name, f"{base.label} (dark)", base.category, base.tags)
        ic.rect(6, 6, 52, 52, 7).fixed(INK)
        ic.acc().rect(6, 50, 52, 8, 0).fixed(FAMILY[fam][0])          # family colour strip
        ic.cut().rect(6, 49, 52, 1.2)
        ic.text(32, 40, sym, 26 if len(sym) == 1 else 22, "bold", role="x", fill="#FFFFFF")
        ic.text(13, 19, str(n), 9, "bold", anchor="start", role="x", fill="#FFFFFF")
        return ic

    def named(base, n=n, sym=sym, fam=fam, name=name):
        ic = Icon(base.name, f"{base.label} (with name)", base.category, base.tags)
        ic.rect(6, 6, 52, 52, 7).fixed(FAMILY[fam][0])
        ic.text(32, 37, sym, 22 if len(sym) == 1 else 19, "bold", role="x", fill=INK)
        ic.text(13, 19, str(n), 8.5, "bold", anchor="start", role="x", fill=INK)
        fs = 7.2 if len(name) <= 9 else (6.2 if len(name) <= 11 else 5.3)
        ic.text(32, 50, name, fs, "normal", role="x", fill=INK)
        return ic

    def atom(base, sym=sym, fam=fam):
        ic = Icon(base.name, f"{base.label} (atom)", base.category, base.tags)
        ic.circle(32, 32, 25).fixed(FAMILY[fam][0])
        ic.acc().ellipse(23.5, 21.5, 8, 5.2).fixed("#FFFFFF")
        ic.ops[-1].attrs["opacity"] = 0.55
        ic.ops[-1].attrs["transform"] = "rotate(-35 23.5 21.5)"
        ic.text(32, 41, sym, 25 if len(sym) == 1 else 21, "bold", role="x", fill=INK)
        return ic
    V.OVERRIDES[nm] = (("symbol-only", plain), ("outline", outline), ("dark", dark), ("named", named), ("atom", atom))


for _e in _parse():
    _make(*_e)


# ----------------------------------------------------------------- molecules (hand-placed 2D coordinates)
def _chain(n, x0=0.0, y0=0.0, d=1.0):
    return [(x0 + i * d * 0.866, y0 + (0.5 * d if i % 2 else 0)) for i in range(n)]


def _ring(cx, cy, r=1.0, n=6, rot=90):
    return [(cx + r * math.cos(math.radians(rot + 360 * k / n)), cy + r * math.sin(math.radians(rot + 360 * k / n))) for k in range(n)]


MOLECULES = {}


def mol(name, label, tags):
    def deco(fn):
        atoms, bonds = fn()
        MOLECULES[name] = (label, tags, atoms, bonds)
        return fn
    return deco


@mol("mol_water", "Water (H₂O)", ["water", "H2O", "solvent"])
def _():
    return [("O", 0, 0), ("H", -0.8, -0.6), ("H", 0.8, -0.6)], [(0, 1, 1), (0, 2, 1)]


@mol("mol_methanol", "Methanol", ["methanol", "MeOH", "solvent"])
def _():
    return [("C", 0, 0), ("O", 1, 0)], [(0, 1, 1)]


@mol("mol_ethanol", "Ethanol", ["ethanol", "EtOH", "solvent"])
def _():
    p = _chain(3)
    return [("C", *p[0]), ("C", *p[1]), ("O", *p[2])], [(0, 1, 1), (1, 2, 1)]


@mol("mol_acetone", "Acetone", ["acetone", "propanone", "solvent"])
def _():
    return [("C", -0.87, -0.5), ("C", 0, 0), ("C", 0.87, -0.5), ("O", 0, 1)], [(0, 1, 1), (1, 2, 1), (1, 3, 2)]


@mol("mol_dmf", "DMF", ["DMF", "dimethylformamide", "solvent", "MOF synthesis"])
def _():
    return [("O", -0.87, 1.0), ("C", 0, 0.5), ("N", 0.87, 1.0), ("C", 0.87, 2.0), ("C", 1.73, 0.5)], \
        [(0, 1, 2), (1, 2, 1), (2, 3, 1), (2, 4, 1)]


@mol("mol_dmso", "DMSO", ["DMSO", "dimethyl sulfoxide", "solvent"])
def _():
    return [("C", -0.87, -0.5), ("S", 0, 0), ("C", 0.87, -0.5), ("O", 0, 1)], [(0, 1, 1), (1, 2, 1), (1, 3, 2)]


@mol("mol_acetonitrile", "Acetonitrile", ["acetonitrile", "MeCN", "solvent", "HPLC"])
def _():
    return [("C", 0, 0), ("C", 1, 0), ("N", 2, 0)], [(0, 1, 1), (1, 2, 2)]


@mol("mol_thf", "THF", ["THF", "tetrahydrofuran", "solvent"])
def _():
    r = _ring(0, 0, 1, 5, -90)
    return [("O", *r[0])] + [("C", *p) for p in r[1:]], [(i, (i + 1) % 5, 1) for i in range(5)]


@mol("mol_toluene", "Toluene", ["toluene", "aromatic", "solvent"])
def _():
    r = _ring(0, 0, 1)
    return [("C", *p) for p in r] + [("C", 0, 2)], [(i, (i + 1) % 6, 2 if i % 2 else 1) for i in range(6)] + [(0, 6, 1)]


@mol("mol_chloroform", "Chloroform", ["chloroform", "CHCl3", "CDCl3", "solvent"])
def _():
    return [("C", 0, 0), ("Cl", 0, 1), ("Cl", -0.87, -0.5), ("Cl", 0.87, -0.5)], [(0, 1, 1), (0, 2, 1), (0, 3, 1)]


@mol("mol_dcm", "Dichloromethane", ["DCM", "dichloromethane", "CH2Cl2", "solvent"])
def _():
    return [("Cl", -0.87, -0.5), ("C", 0, 0), ("Cl", 0.87, -0.5)], [(0, 1, 1), (1, 2, 1)]


@mol("mol_formic_acid", "Formic Acid", ["formic acid", "modulator", "HCOOH"])
def _():
    return [("C", 0, 0), ("O", 0.87, 0.5), ("O", -0.87, 0.5)], [(0, 1, 2), (0, 2, 1)]


@mol("mol_acetic_acid", "Acetic Acid", ["acetic acid", "modulator", "AcOH"])
def _():
    return [("C", 0, 0), ("C", 0.87, 0.5), ("O", 1.73, 0), ("O", 0.87, 1.5)], [(0, 1, 1), (1, 2, 1), (1, 3, 2)]


@mol("mol_benzoic_acid", "Benzoic Acid", ["benzoic acid", "modulator", "monocarboxylic acid"])
def _():
    r = _ring(0, 0, 1)
    return [("C", *p) for p in r] + [("C", 0, -2), ("O", -0.87, -2.5), ("O", 0.87, -2.5)], \
        [(i, (i + 1) % 6, 2 if i % 2 else 1) for i in range(6)] + [(3, 6, 1), (6, 7, 2), (6, 8, 1)]


@mol("mol_2_methylimidazole", "2-Methylimidazole", ["2-methylimidazole", "Hmim", "ZIF-8 linker", "imidazolate"])
def _():
    r = _ring(0, 0, 1, 5, -90)
    return [("C", *r[0]), ("N", *r[1]), ("C", *r[2]), ("C", *r[3]), ("N", *r[4]), ("C", 0, -2)], \
        [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 0, 1), (0, 5, 1)]


@mol("mol_terephthalic_acid", "Terephthalic Acid (H₂BDC)", ["terephthalic acid", "BDC", "UiO-66 linker", "MOF-5 linker"])
def _():
    r = _ring(0, 0, 1)
    a = [("C", *p) for p in r] + [("C", 0, 2), ("O", -0.87, 2.5), ("O", 0.87, 2.5), ("C", 0, -2), ("O", -0.87, -2.5), ("O", 0.87, -2.5)]
    return a, [(i, (i + 1) % 6, 2 if i % 2 else 1) for i in range(6)] + [(0, 6, 1), (6, 7, 2), (6, 8, 1), (3, 9, 1), (9, 10, 2), (9, 11, 1)]


@mol("mol_glycylglycine", "Glycylglycine", ["GlyGly", "dipeptide", "peptide bond", "hydrolysis substrate"])
def _():
    p = _chain(7)
    els = ["N", "C", "C", "N", "C", "C", "O"]
    atoms = [(e, *xy) for e, xy in zip(els, p)] + [("O", p[2][0], p[2][1] + 1.0), ("O", p[5][0], p[5][1] - 1.0)]
    return atoms, [(i, i + 1, 1) for i in range(6)] + [(2, 7, 2), (5, 8, 2)]


@mol("mol_co2", "CO₂", ["CO2", "carbon dioxide", "gas"])
def _():
    return [("O", -1, 0), ("C", 0, 0), ("O", 1, 0)], [(0, 1, 2), (1, 2, 2)]


@mol("mol_n2", "N₂", ["N2", "nitrogen", "gas"])
def _():
    return [("N", 0, 0), ("N", 1, 0)], [(0, 1, 2)]


@mol("mol_ch4", "Methane", ["methane", "CH4", "gas"])
def _():
    return [("C", 0, 0), ("H", 0, 1), ("H", -0.95, -0.4), ("H", 0.95, -0.4), ("H", 0.3, -0.95)], \
        [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)]


@mol("mol_nh3", "Ammonia", ["ammonia", "NH3", "gas"])
def _():
    return [("N", 0, 0), ("H", 0, -1), ("H", -0.87, 0.5), ("H", 0.87, 0.5)], [(0, 1, 1), (0, 2, 1), (0, 3, 1)]


@mol("mol_h2o2", "Hydrogen Peroxide", ["H2O2", "peroxide", "oxidant"])
def _():
    return [("H", -1.4, 0.6), ("O", -0.5, 0), ("O", 0.5, 0), ("H", 1.4, -0.6)], [(0, 1, 1), (1, 2, 1), (2, 3, 1)]


MOL_STYLES = (("ball-and-stick", "ball"), ("space-filling", "space"), ("labelled", "label"))

for _n, (_lab, _tags, _atoms, _bonds) in MOLECULES.items():
    def _build(atoms=_atoms, bonds=_bonds, n=_n, lab=_lab, tags=_tags):
        show_h = any(a[0] == "H" for a in atoms) and len(atoms) <= 5
        ic = icon_from_molecule(atoms, bonds, n, lab, show_h=show_h)
        ic.tags = tags + ["molecule"]
        return ic
    REGISTRY[_n] = _build

    def _mk(mode, atoms=_atoms, bonds=_bonds, n=_n, lab=_lab, tags=_tags, title=""):
        def f(base):
            ic = icon_from_molecule(atoms, bonds, n, f"{lab} ({title})", show_h=True, mode=mode)
            ic.tags = tags + ["molecule"]
            return ic
        return f
    V.OVERRIDES[_n] = tuple((key, _mk(mode, title=key.replace("-", " "))) for key, mode in MOL_STYLES)


# ----------------------------------------------------------------- topology nets (2D projections)
@icon("net_sql", "materials", label="sql Net", tags=["sql", "square lattice", "2D net", "topology", "layer"])
def _(ic):
    xs = (10, 24, 40, 54)
    for x in xs:
        ic.sec().line(x, 8, x, 56, 2.2)
        ic.sec().line(8, x, 56, x, 2.2)
    for x in xs:
        for y in xs:
            ic.circle(x, y, 3.6)


@icon("net_kgm", "materials", label="kgm Net (Kagome)", tags=["kagome", "kgm", "topology", "2D net"])
def _(ic):
    # one kagome patch: central hexagon with six corner-sharing triangles
    hexp = [polar(32, 32, 15, 30 + 60 * k) for k in range(6)]
    for k in range(6):
        p, q = hexp[k], hexp[(k + 1) % 6]
        tip = polar(32, 32, 26, 60 + 60 * k)
        for a, b in ((p, q), (q, tip), (tip, p)):
            ic.sec().line(*a, *b, 2.2)
    for k in range(6):
        ic.circle(*hexp[k], 3.4)
        ic.circle(*polar(32, 32, 26, 60 + 60 * k), 3.0)


@icon("net_hxl", "materials", label="hxl Net (Triangular)", tags=["hxl", "triangular lattice", "6-connected", "topology"])
def _(ic):
    pts = []
    for j in range(-2, 3):
        for i in range(-2, 3):
            x = 32 + i * 13 + (6.5 if j % 2 else 0)
            y = 32 + j * 11.3
            if 4 < x < 60:
                pts.append((x, y))
    for a in pts:
        for b in pts:
            if a < b and abs(math.dist(a, b) - 13) < 0.6:
                ic.sec().line(*a, *b, 1.8)
    for p in pts:
        ic.circle(*p, 3.2)


@icon("net_bcu", "materials", label="bcu Net", tags=["bcu", "body-centred cubic", "8-connected", "topology"])
def _(ic):
    f = [(12, 24), (40, 24), (40, 52), (12, 52)]
    b = [(x + 12, y - 12) for x, y in f]
    for i in range(4):
        for p, q in ((f[i], f[(i + 1) % 4]), (b[i], b[(i + 1) % 4]), (f[i], b[i])):
            ic.line(*p, *q, 1.2)
    for p in f + b:
        ic.sec().line(32, 32, *p, 2.2)
        ic.circle(*p, 3.4)
    ic.acc().circle(32, 32, 5)


@icon("net_srs", "materials", label="srs Net (Chiral)", tags=["srs", "chiral net", "3-connected", "gyroid", "topology"])
def _(ic):
    pts = [(32, 8), (48, 18), (48, 36), (32, 46), (16, 36), (16, 18), (32, 27), (40, 56), (24, 56), (56, 46), (8, 46)]
    for a, b in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (6, 0), (6, 2), (6, 4), (3, 7), (3, 8), (2, 9), (4, 10)):
        ic.sec().line(*pts[a], *pts[b], 2.2)
    for p in pts:
        ic.circle(*p, 3.2)
    ic.acc().circle(32, 27, 3.6)
