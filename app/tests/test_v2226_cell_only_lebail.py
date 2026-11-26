"""v2.22.6: cell-only reserve for Le Bail. A CIF with a cell and a space group but no readable atom sites gives a
phase without an atom model: reflections are generated from the cell and the space-group absences (gemmi).
Only Le Bail can use it; Rietveld needs atoms.

Limits of these tests (stated on purpose):
  * the no-atoms CIF is a SYNTHETIC file written in tmp_path, not a real CIF. It is rejected by real pymatgen.
  * the fallback order is tested with a monkeypatched reader only where noted (pymatgen structure reader and the
    gemmi structure reader raised ValueError). Other exception types are not caught and propagate.
  * the absence rule is gemmi's Operations.is_systematically_absent; it is checked on I-centred cubic reflections.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest
from mfs import rietveld as RV

CELL_ONLY_I432 = """data_cell_only
_cell_length_a 31.006
_cell_length_b 31.006
_cell_length_c 31.006
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'I 4 3 2'
_symmetry_Int_Tables_number 211
"""


@pytest.fixture
def cell_cif(tmp_path):
    p = tmp_path / "cell_only.cif"
    p.write_text(CELL_ONLY_I432, encoding="utf-8")
    return str(p)


def test_pymatgen_rejects_the_cell_only_file_and_the_reserve_gives_a_phase(cell_cif):
    pytest.importorskip("gemmi")
    ph = RV.Phase.from_cif(cell_cif)  # real pymatgen rejects this file (no atoms); the cell-only reserve is used
    assert getattr(ph, "cell_only", False) is True
    assert ph.sg_symbol == "I 4 3 2"
    assert ph.cell0[0] == pytest.approx(31.006, abs=1e-6)


def test_reflections_follow_the_i432_absence_rules(cell_cif):
    """Checks the I-centring rule (h+k+l even) and a set of explicit reflections of I432.
    Only the I-centring absences are tested here; the other I432 absences (e.g. the 4-fold screw and the 3-fold
    rules of the symmetry operations) are NOT tested one by one: they come from gemmi and are not reproduced here."""
    pytest.importorskip("gemmi")
    ph = RV.Phase.from_cif(cell_cif).build(1.0)
    hkl = np.asarray(ph.m_hkl)
    assert hkl.shape[0] > 0
    assert np.all((hkl.sum(1) % 2) == 0), "I-centred: h+k+l must be even"
    present = {tuple(int(v) for v in r) for r in hkl}
    assert (2, 0, 0) in present and (2, 2, 0) in present and (4, 0, 0) in present   # allowed
    assert (1, 0, 0) not in present and (1, 1, 1) not in present and (2, 1, 0) not in present  # absent


def test_rietveld_mode_refuses_a_cell_only_phase(cell_cif):
    pytest.importorskip("gemmi")
    ph = RV.Phase.from_cif(cell_cif)
    with pytest.raises(ValueError, match="cell-only"):
        RV.Refinement(np.linspace(5, 60, 200), np.ones(200) * 100.0, [ph], "Cu Kα1", 6, "rietveld")


def test_lebail_mode_accepts_a_cell_only_phase(cell_cif):
    pytest.importorskip("gemmi")
    ph = RV.Phase.from_cif(cell_cif)
    R = RV.Refinement(np.linspace(5, 60, 200), np.ones(200) * 100.0, [ph], "Cu Kα1", 6, "lebail")
    assert R.mode == "lebail"


def test_space_group_absences_remove_reflections_compared_with_plain_cell(cell_cif):
    pytest.importorskip("gemmi")
    with_sg = RV.Phase.from_cif(cell_cif).build(1.0)
    plain = RV.Phase.from_cell(with_sg.cell0, "plain").build(1.0)
    assert with_sg.n_ref < plain.n_ref


def test_monkeypatched_structure_readers_both_failing_use_cell_only(cell_cif, monkeypatch):
    pytest.importorskip("gemmi")
    from pymatgen.core import Structure

    def pm_fail(*a, **k):
        raise ValueError("pymatgen rejects this file")

    def gemmi_fail(path):
        raise ValueError("gemmi finds no atom sites")

    monkeypatch.setattr(Structure, "from_file", staticmethod(pm_fail))
    monkeypatch.setattr(RV, "_structure_from_cif_gemmi", gemmi_fail)
    ph = RV.Phase.from_cif(cell_cif)
    assert ph.cell_only is True
    assert ph.sg_symbol == "I 4 3 2"


def test_a_file_without_a_cell_or_space_group_still_fails(tmp_path):
    pytest.importorskip("gemmi")
    bad = tmp_path / "no_cell.cif"
    bad.write_text("data_x\n_chemical_name_common 'nothing useful'\n", encoding="utf-8")
    with pytest.raises(ValueError):
        RV.Phase.from_cif(str(bad))
