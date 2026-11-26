"""v2.22.5: gemmi is a reserve CIF reader for the phase loader. pymatgen stays the first choice.
Checks use only repo reference CIFs (no private data).

Limits of these tests (stated on purpose):
  * the FALLBACK PATH is tested only with a monkeypatched pymatgen (Structure.from_file raises ValueError).
    No real CIF that pymatgen rejects and gemmi reads is in the repo, so the real-file path is NOT tested.
  * only ValueError from pymatgen triggers the fallback. Other exception types (e.g. OSError, TypeError,
    KeyError from pymatgen internals) are NOT caught and propagate unchanged."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pytest
from mfs import rietveld as RV

REF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "references", "library")
CEO2 = os.path.join(REF, "CeO2.cif")


def test_gemmi_reader_matches_pymatgen_on_a_reference_cif():
    pytest.importorskip("gemmi")
    from pymatgen.core import Structure
    s_gemmi = RV._structure_from_cif_gemmi(CEO2)
    s_pmg = Structure.from_file(CEO2)
    assert s_gemmi.lattice.abc == pytest.approx(s_pmg.lattice.abc, abs=1e-3)
    assert len(s_gemmi) == len(s_pmg)
    assert s_gemmi.composition.reduced_formula == s_pmg.composition.reduced_formula


def test_fallback_is_used_only_when_pymatgen_fails(monkeypatch):
    pytest.importorskip("gemmi")
    from pymatgen.core import Structure
    calls = {"gemmi": 0}
    real = RV._structure_from_cif_gemmi

    def spy(path):
        calls["gemmi"] += 1
        return real(path)

    monkeypatch.setattr(RV, "_structure_from_cif_gemmi", spy)
    ph = RV.Phase.from_cif(CEO2)  # pymatgen reads this file: no fallback
    assert calls["gemmi"] == 0
    a_ref = ph._structure.lattice.a

    def boom(*a, **k):
        raise ValueError("pymatgen rejects this file")

    monkeypatch.setattr(Structure, "from_file", staticmethod(boom))
    ph2 = RV.Phase.from_cif(CEO2)  # pymatgen 'fails' -> gemmi reserve
    assert calls["gemmi"] == 1
    assert ph2._structure.lattice.a == pytest.approx(a_ref, abs=1e-3)


def test_reserve_reports_a_clear_error_when_gemmi_cannot_read(tmp_path):
    bad = tmp_path / "bad.cif"
    bad.write_text("this is not a cif file\n", encoding="utf-8")
    with pytest.raises(ValueError):
        RV._structure_from_cif_gemmi(str(bad))
