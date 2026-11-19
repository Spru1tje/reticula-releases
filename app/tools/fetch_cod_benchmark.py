"""Developer tool: downloads a benchmark set of crystal structures (CIF) from the Crystallography Open Database (COD, public domain) for the indexing and
structure-solution benchmarks (tests/data/cod_bench). Usage:  python tools/fetch_cod_benchmark.py [target_dir]
Each entry: (label, formula query, crystal system wanted, max cell volume (Å3)). The first COD entry that fits the filter is kept."""
import os
import re
import sys
import time
import urllib.parse
import urllib.request

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tests", "data", "cod_bench")
UA = {"User-Agent": "reticula-benchmark/1.0"}

WANT = [  # label, formula, system, max volume
    ("NaCl", "Cl Na", "cubic", 400), ("MgO", "Mg O", "cubic", 200), ("CaF2", "Ca F2", "cubic", 400), ("ZnS_sphalerite", "S Zn", "cubic", 400), ("Si", "Si", "cubic", 400),
    ("Cu", "Cu", "cubic", 100), ("CeO2", "Ce O2", "cubic", 400), ("SrTiO3", "O3 Sr Ti", "cubic", 100), ("Y2O3", "O3 Y2", "cubic", 1500),
    ("anatase", "O2 Ti", "tetragonal", 200), ("rutile", "O2 Ti", "tetragonal", 100), ("SnO2", "O2 Sn", "tetragonal", 100), ("ZrO2_tetragonal", "O2 Zr", "tetragonal", 200),
    ("ZnO", "O Zn", "hexagonal", 100), ("quartz", "O2 Si", "hexagonal", 200), ("corundum", "Al2 O3", "trigonal", 400), ("graphite", "C", "hexagonal", 100),
    ("Mg", "Mg", "hexagonal", 100), ("GaN", "Ga N", "hexagonal", 100),
    ("aragonite", "C Ca O3", "orthorhombic", 500), ("barite", "Ba O4 S", "orthorhombic", 1000), ("forsterite", "Mg2 O4 Si", "orthorhombic", 1000),
    ("brookite", "O2 Ti", "orthorhombic", 500), ("celestine", "O4 S Sr", "orthorhombic", 1000), ("PbSO4_anglesite", "O4 Pb S", "orthorhombic", 1000),
    ("baddeleyite", "O2 Zr", "monoclinic", 400), ("glycine_gamma", "C2 H5 N O2", "trigonal", 1500), ("alpha_glycine", "C2 H5 N O2", "monoclinic", 1000),
    ("gypsum", "Ca H4 O6 S", "monoclinic", 1500), ("calcite", "C Ca O3", "trigonal", 1000), ("wollastonite", "Ca O3 Si", "triclinic", 1500),
]


def get(url, timeout=40):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def cell_of(cif):
    def val(tag):
        m = re.search(rf"^{tag}\s+([0-9.]+)", cif, re.M)
        return float(m.group(1)) if m else None
    return [val(t) for t in ("_cell_length_a", "_cell_length_b", "_cell_length_c", "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma")]


def system_of(cell, cif):
    a, b, c, al, be, ga = cell
    close = lambda x, y, t=0.02: abs(x - y) / max(x, y) < t
    sg = re.search(r"_symmetry_space_group_name_H-M\s+'?([^'\n]+)", cif) or re.search(r"_space_group_name_H-M_alt\s+'?([^'\n]+)", cif)
    if close(a, b) and close(b, c) and abs(al - 90) < 0.1 and abs(be - 90) < 0.1 and abs(ga - 90) < 0.1:
        return "cubic"
    if close(a, b) and abs(al - 90) < 0.1 and abs(be - 90) < 0.1 and abs(ga - 90) < 0.1:
        return "tetragonal"
    if close(a, b) and abs(al - 90) < 0.1 and abs(be - 90) < 0.1 and abs(ga - 120) < 0.1:
        return "hexagonal"
    if close(a, b) and close(b, c) and close(al, be, 0.001) and close(be, ga, 0.001):
        return "trigonal"
    if abs(al - 90) < 0.1 and abs(be - 90) < 0.1 and abs(ga - 90) < 0.1:
        return "orthorhombic"
    if abs(al - 90) < 0.1 and abs(ga - 90) < 0.1:
        return "monoclinic"
    return "triclinic"


def main():
    os.makedirs(OUT, exist_ok=True)
    for label, formula, system, vmax in WANT:
        dst = os.path.join(OUT, f"{label}.cif")
        if os.path.exists(dst):
            continue
        try:
            urls = get("https://www.crystallography.net/cod/result?" + urllib.parse.urlencode({"formula": formula, "format": "urls", "max_r": 1}) if False else
                       "https://www.crystallography.net/cod/result?" + urllib.parse.urlencode({"formula": formula, "format": "urls"})).split()
        except Exception as e:
            print(label, "search failed", e); continue
        got = False
        for u in urls[:60]:
            try:
                cif = get(u)
            except Exception:
                continue
            cell = cell_of(cif)
            if None in cell:
                continue
            vol = cell[0] * cell[1] * cell[2]
            ok = system_of(cell, cif)
            if ok == system and vol <= vmax * 1.2 and ("_atom_site_fract_x" in cif) and ("_atom_site_occupancy" in cif or True):
                open(dst, "w", encoding="utf-8").write(cif)
                print(label, "->", u.split("/")[-1], system, [round(v, 3) for v in cell]); got = True
                break
            time.sleep(0.05)
        if not got:
            print(label, "no suitable entry")


if __name__ == "__main__":
    main()
