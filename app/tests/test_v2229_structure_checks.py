"""v2.22.9 (S-Structuur-1..3). Only repo reference files and textbook absence rules; no private data.
  (1) the cell in references/library/index.json equals the cell of the CIF file it points to;
  (2) space-group absences of two non-I groups through the cell-only route (Fm-3m, R-3c, hexagonal setting);
  (3) every reflection of the structure route is in the cell-only reflection list (same cell, CeO2).
Textbook rules used (International Tables): Fm-3m: h,k,l all even or all odd (F-centring); R-3c hexagonal setting:
-h+k+l = 3n (R-centring) and 000l: l = 6n (c-glide/3-fold screw), hhl: l = 2n (c-glide)."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest
from mfs import rietveld as RV

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIB = os.path.join(APP, "references", "library")
INDEX = json.load(open(os.path.join(LIB, "index.json"), encoding="utf-8"))


def _cif_cell_text(cell, symbol):
    a, b, c, al, be, ga = cell
    return (f"data_cellonly\n_cell_length_a {a}\n_cell_length_b {b}\n_cell_length_c {c}\n"
            f"_cell_angle_alpha {al}\n_cell_angle_beta {be}\n_cell_angle_gamma {ga}\n"
            f"_symmetry_space_group_name_H-M '{symbol}'\n")


def _cell_only_phase(tmp_path, cell, symbol):
    p = tmp_path / "cell_only.cif"
    p.write_text(_cif_cell_text(cell, symbol), encoding="utf-8")
    return RV.Phase.from_cif(str(p))


# ---------------------------------------------------------------- (1) index cell equals the CIF cell
@pytest.mark.parametrize("entry", [e for e in INDEX if os.path.exists(os.path.normpath(os.path.join(LIB, e["file"])))],
                         ids=lambda e: e["file"])
def test_index_cell_equals_the_cif_cell(entry):
    pytest.importorskip("pymatgen.core")
    from pymatgen.core import Structure
    s = Structure.from_file(os.path.normpath(os.path.join(LIB, entry["file"])))
    L = s.lattice
    cif_cell = [L.a, L.b, L.c, L.alpha, L.beta, L.gamma]
    assert entry["cell"] == pytest.approx(cif_cell, abs=1e-3)


# ---------------------------------------------------------------- (2) absences of non-I groups via the cell-only route
def test_fm3m_absences_through_the_cell_only_route(tmp_path):
    pytest.importorskip("gemmi")
    ph = _cell_only_phase(tmp_path, (5.4113, 5.4113, 5.4113, 90, 90, 90), "F m -3 m")
    ph.build(1.0)
    present = {tuple(int(v) for v in r) for r in np.asarray(ph.m_hkl)}
    assert (1, 0, 0) not in present and (1, 1, 0) not in present     # mixed parity: absent
    assert (2, 0, 0) in present and (1, 1, 1) in present and (2, 2, 0) in present   # all even / all odd: allowed
    assert all((h % 2 == k % 2 == l % 2) for h, k, l in present)     # F-centring rule for every reflection


def test_r3c_hexagonal_absences_through_the_cell_only_route(tmp_path):
    pytest.importorskip("gemmi")
    # corundum-type cell (repo reference cell, hexagonal setting)
    ph = _cell_only_phase(tmp_path, (4.7606, 4.7606, 12.994, 90, 90, 120), "R -3 c")
    ph.build(1.0)
    present = {tuple(int(v) for v in r) for r in np.asarray(ph.m_hkl)}
    # R-centring: -h+k+l = 3n; 0 0 l: l = 6n; h h l: l = 2n
    assert (1, 0, 0) not in present          # -1+0+0 not 3n
    assert (1, 1, 1) not in present          # -1+1+1 = 1 not 3n
    assert (0, 0, 3) not in present          # 000l needs l = 6n
    assert (0, 0, 6) in present
    # (1 0 1) is absent in corundum (glide rule), even though -h+k+l = 0; (1 0 4) is observed in corundum (strong 104)
    assert (1, 0, 1) not in present
    assert (1, 0, 4) in present
    assert all(((-h + k + l) % 3 == 0) for h, k, l in present)


# ---------------------------------------------------------------- (3) cell-only list contains the structure-route list (CeO2)
def test_cell_only_reflections_contain_the_structure_route_reflections(tmp_path):
    pytest.importorskip("gemmi")
    struct_ph = RV.Phase.from_cif(os.path.join(LIB, "CeO2.cif")).build(1.2)
    cell = (5.4113, 5.4113, 5.4113, 90, 90, 90)
    cell_ph = _cell_only_phase(tmp_path, cell, "F m -3 m").build(1.2)
    s_set = {tuple(int(v) for v in r) for r in np.asarray(struct_ph.hkl)}
    c_set = {tuple(int(v) for v in r) for r in np.asarray(cell_ph.hkl)}
    # the structure route drops reflections whose intensity is numerically zero, so it is a subset of the cell list
    # up to the sign of the half-space; compare |h|,|k|,|l| sorted
    key = lambda t: tuple(sorted(abs(v) for v in t))
    s_keys = {key(t) for t in s_set}
    c_keys = {key(t) for t in c_set}
    assert s_keys <= c_keys
