"""v2.18.4: outputs of any refinement (mfs.refout): CIF with the refined cell, 3D payload, .xy files, reflection list, ZIP."""
import io
import os
import sys
import warnings
import zipfile

import numpy as np
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
CIF = os.path.join(HERE, "references", "MOF-801_Zr_COD4121459.cif")


@pytest.fixture(scope="module")
def refined():
    warnings.filterwarnings("ignore")
    from mfs import analysis as A, rietveld as RV
    sim, sticks, formula = A.simulate_pxrd(CIF, (4.0, 40.0), fwhm=0.12)
    rng = np.random.default_rng(2)
    x = sim.x.values
    y = 1000 * sim.y.values / sim.y.values.max() + 30 + rng.normal(0, 2, len(x))
    ph = [RV.Phase.from_cif(CIF, "MOF-801")]
    res = RV.quick_refine(x, y, ph, "Cu Kα1", "rietveld", tt_range=(4.0, 40.0))
    return res, ph


def test_cif_has_the_refined_cell_and_all_atoms(refined):
    from mfs import refout as RO
    from pymatgen.core import Structure
    res, ph = refined
    txt = RO.cif_for_phase(res, 0, ph[0], "test MOF-801", "x")
    assert "_cell_length_a" in txt and "NOT refined: the atom positions" in txt
    s = Structure.from_str(txt, fmt="cif")
    a_ref = res["phases"][0]["cell"][0]
    assert s.lattice.a == pytest.approx(a_ref, abs=2e-4)                       # esd digits may round the last decimal
    assert len(s) == len(ph[0]._structure)                                    # the symmetry-expanded cell content is complete


def test_3d_payload_in_the_refined_cell(refined):
    from mfs import refout as RO
    res, ph = refined
    d = RO.payload_for_phase(res, 0, ph[0])
    assert d and len(d["atoms"]) >= len(ph[0]._structure) and d["bonds"]
    assert d["cell"][0][0] == pytest.approx(res["phases"][0]["cell"][0], abs=1e-6)


def test_reflection_list_and_hkl_file(refined):
    from mfs import refout as RO
    res, ph = refined
    t = RO.reflection_table(res, 0, ph[0])
    assert len(t) > 20 and t["d (Å)"].is_monotonic_decreasing and t["intensity"].max() == pytest.approx(100.0)
    assert abs(t["2theta (°)"].iloc[0] - 2 * np.degrees(np.arcsin(1.5406 / (2 * t["d (Å)"].iloc[0])))) < 0.05
    h = RO.hkl_file(t).splitlines()
    assert len(h) == len(t) + 1 and h[-1].split() == ["0", "0", "0", "0.00", "0.00"]


def test_bundle_contains_every_file_to_continue_with(refined):
    from mfs import refout as RO
    res, ph = refined
    z = zipfile.ZipFile(io.BytesIO(RO.bundle(res, ph, "sample 1", "x")))
    names = z.namelist()
    for ext in (".cif", "_reflections.csv", ".hkl", "_observed.xy", "_calculated.xy", "_difference.xy", "_background.xy", "_pattern.csv", "_summary.json"):
        assert any(n.endswith(ext) for n in names), (ext, names)
    obs = [n for n in names if n.endswith("_observed.xy")][0]
    first = [l for l in z.read(obs).decode().splitlines() if not l.startswith("#")][0].split()
    assert float(first[0]) == pytest.approx(res["x"][0], abs=1e-4)


def test_page_shows_3d_cif_and_files_after_a_plain_rietveld_run():
    from streamlit.testing.v1 import AppTest
    warnings.filterwarnings("ignore")
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600)
    at.session_state["_lv"] = {"global": "Standard"}
    at.run()
    at.button(key="nv_" + "PXRD").click().run()
    [t for t in at.toggle if "Rietveld" in t.label][0].set_value(True).run()
    ms = at.multiselect(key="rv_phases"); ms.set_value([o for o in ms.options if "801" in o][:1]).run()
    [b for b in at.button if b.key == "rv_go"][0].click().run()
    assert not at.exception, [e.value for e in at.exception][:2]
    assert [t for t in at.tabs if "3D structure" in t.label], [t.label for t in at.tabs]
    labels = [getattr(d, "label", "") for d in at.get("download_button")]
    assert any("CIF (refined cell)" in l for l in labels) and any("Everything to continue with" in l for l in labels), labels


def test_cell_only_result_gets_a_cell_cif():
    from mfs import refout as RO, rietveld as RV
    warnings.filterwarnings("ignore")
    x = np.linspace(8, 60, 2000)
    cub = RV.Phase.from_cell((8.0, 8.0, 8.0, 90, 90, 90), "cell")
    r0 = RV.Refinement(x, np.full_like(x, 100.0), [cub], "Cu Kα1", 6, "lebail", 0.0, instrument={"U": 0.002, "V": -0.001, "W": 0.006, "X": 0.03, "Y": 0.01, "fixed": True})
    res = RV.quick_refine(x, r0._solve_linear(r0._values())[2], [RV.Phase.from_cell((8.02, 8.02, 8.02, 90, 90, 90), "cell")], mode="lebail")
    txt = RO.cif_for_phase(res, 0, None, "unknown", "x")
    assert "_cell_length_a" in txt and "No atoms yet" in txt and "_atom_site" not in txt
