"""v2.22.8: (S-Structuur-1) deuterium in the MIL-53_Cr reference is treated like hydrogen in the structure factors;
(S-Structuur-2) the reserve reader (gemmi) and pymatgen agree on the repo reference CIFs.
Only repo reference files are used (no private data). These are checks of the code path, not physical validation."""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest
from mfs import rietveld as RV

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIB = os.path.join(APP, "references", "library")
MIL53 = os.path.join(LIB, "MIL-53_Cr.cif")


# ---------------------------------------------------------------- S-Structuur-1: D is treated as H
def test_deuterium_scatters_like_hydrogen_in_mil53_structure_factors():
    from pymatgen.core import Element
    assert Element("D").Z == Element("H").Z == 1
    ph_d = RV.Phase.from_cif(MIL53)                              # D as in the CIF
    ph_h = RV.Phase.from_cif(MIL53, substitute={"D": "H"})       # the same atoms, D replaced by H
    # the loader already maps the CIF's deuterium to hydrogen: no 'D' species reaches the structure factors
    assert "D" not in set(ph_d.elems)
    assert "H" in set(ph_d.elems)
    ph_d.build(1.0)
    ph_h.build(1.0)
    assert ph_d.m_F2.shape == ph_h.m_F2.shape
    assert np.allclose(ph_d.m_F2, ph_h.m_F2, rtol=1e-9, atol=1e-9)


# ---------------------------------------------------------------- S-Structuur-2: reader parity on repo references
def _reference_cifs():
    files = sorted(glob.glob(os.path.join(LIB, "*.cif")))
    files += sorted(glob.glob(os.path.join(APP, "references", "*.cif")))
    return files


@pytest.mark.parametrize("path", _reference_cifs(), ids=lambda p: os.path.basename(p))
def test_gemmi_reserve_and_pymatgen_agree_on_reference_cif(path):
    pytest.importorskip("gemmi")
    from pymatgen.core import Structure
    s_pm = Structure.from_file(path)
    s_gm = RV._structure_from_cif_gemmi(path)
    assert len(s_gm) == len(s_pm)
    assert s_gm.lattice.abc == pytest.approx(s_pm.lattice.abc, abs=2e-3)
    assert s_gm.lattice.angles == pytest.approx(s_pm.lattice.angles, abs=2e-2)
    assert s_gm.composition.reduced_formula == s_pm.composition.reduced_formula
