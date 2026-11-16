"""v2.11: defect calculator for any node/metal. Each framework's ideal composition is written as an independent formula
string; the TGA ratio computed from it must give back the ideal linker count with methods (a), (b) M-200 and (d)."""
import os, sys, json, glob, copy
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import pytest
from mfs import defects as D, nodes as N
from mfs.chem import molar_mass, SPECIES, FRAMEWORKS

# framework, plateau formula (per node, dry, as assumed by method a), residue formula per node, L0 per node
CASES = [
    ("HKUST-1", "Cu2(C9H3O6)1.3333333333", "Cu2O2", 4 / 3),
    ("MOF-5", "Zn4O(C8H4O4)3", "Zn4O4", 3),
    ("MIL-53(Al)", "Al(OH)(C8H4O4)", "AlO1.5", 1),
    ("CAU-10(Al)", "Al(OH)(C8H4O4)", "AlO1.5", 1),
    ("MIL-100(Fe)", "Fe3O1.5(C9H3O6)2", "Fe3O4.5", 2),
    ("MIL-101(Cr)", "Cr3O1.5(C8H4O4)3", "Cr3O4.5", 3),
    ("Mg-MOF-74", "Mg2(C8H2O6)", "Mg2O2", 1),
    ("Co-MOF-74", "Co2(C8H2O6)", "Co2O2.6666666667", 1),     # residue Co3O4 in air
    ("Ni-MOF-74", "Ni2(C8H2O6)", "Ni2O2", 1),
    ("MIL-125(Ti)", "Ti8O10(C8H4O4)6", "Ti8O16", 6),
    ("ZIF-8", "Zn(C4H5N2)2", "ZnO", 2),
    ("ZIF-67", "Co(C4H5N2)2", "CoO1.3333333333", 2),
    ("Ce-UiO-66", "Ce6O6(C8H4O4)6", "Ce6O12", 6),
    ("Th-UiO-66", "Th6O6(C8H4O4)6", "Th6O12", 6),
    ("DMOF-1", "Zn2(C8H4O4)2", "Zn2O2", 2),
]


def _inp(fw, **kw):
    d = FRAMEWORKS[fw]; lk = d["linker"]; f, z, role, nH, _ = SPECIES[lk]
    return D.DefectInput(sample=fw, framework=fw, metal=d.get("metal", "Zr"), species=[D.Species(lk, f, z, "linker", nH, 1.0 * nH, False, True)], **kw)


@pytest.mark.parametrize("fw,plat,res,L0", CASES)
def test_method_a_recovers_ideal_linkers_any_node(fw, plat, res, L0):
    inp = _inp(fw, tga_W_plateau=100 * molar_mass(plat) / molar_mass(res), tga_W_end=100.0)
    assert inp.L0() == pytest.approx(L0)
    r = D.method_a(inp)
    assert r.status == "ok" and r.L == pytest.approx(L0, abs=2e-4), (fw, r.L, r.steps)
    assert r.x == pytest.approx(0, abs=2e-4) and r.OS == pytest.approx(inp.Q() - inp.cn(), abs=1e-3)
    assert r.extra["W_plat_ideal_%"] == pytest.approx(100 * molar_mass(plat) / molar_mass(res), rel=1e-5)
    assert any("Node " in a for a in r.assumptions)


@pytest.mark.parametrize("fw,plat,res,L0", [c for c in CASES if c[0] in ("HKUST-1", "MOF-5", "Co-MOF-74", "ZIF-8", "MIL-100(Fe)")])
def test_m200_and_deroo_any_node(fw, plat, res, L0):
    """Ideal dry framework at 200 °C: M-200 (b) gives L0; De Roo (d) with qNMR on a weighed digest gives L0."""
    inp = _inp(fw)
    Mdry = D.build_formula(inp, L0, {}, 0.0)[1]          # ideal formula incl. one H2O per terminal OH (MIL-100)
    inp.tga_W_200 = inp.tga_W_initial = 100 * Mdry / molar_mass(res); inp.tga_W_end = 100.0
    b = D.method_b(inp)[0]
    assert b.status == "ok" and b.L == pytest.approx(L0, abs=1e-3), b.steps
    # De Roo: 10 mg digested in 0.6 mL, TMSP 5 mM (9 H); linker integral set from the ideal content
    ml = inp.main_linker(); c_mof = 10 / Mdry / 0.6 * 1000
    inp.nmr_std_conc_mM, inp.nmr_std_integral, inp.nmr_sample_mg, inp.nmr_volume_mL = 5.0, 9.0, 10.0, 0.6
    ml.integral = L0 * c_mof / 5.0 * ml.nH
    d = D.method_d(inp)
    assert d.L == pytest.approx(L0, rel=1e-6) and abs(d.extra["n_H2O_total"] - (inp.Q() - inp.cn())) < 1e-3


def test_cu_paddlewheel_formula_and_residue():
    inp = _inp("HKUST-1")
    assert inp.nM() == 2 and inp.Q() == 4 and inp.MMO2() == pytest.approx(molar_mass("CuO"))
    f, M, w, h2o = D.build_formula(inp, 4 / 3, {}, 0.0, water_coord=False)
    assert f.startswith("Cu2(BTC)1.33") and abs(w) < 1e-9 and M == pytest.approx(molar_mass("Cu2(C9H3O6)1.3333333333"), rel=1e-6)
    co = _inp("Co-MOF-74"); assert co.residue_label() == "Co3O4" and not co.residue_is_oxidic()


def test_icp_ratio_any_metal():
    inp = _inp("MOF-5")
    M = molar_mass("Zn4O(C8H4O4)3")
    inp.icp_M_wt = 100 * 4 * 65.38 / M; inp.linker_umol_per_mg = 3 / M * 1e3
    e1, e2 = D.method_e(inp)
    assert e1.L == pytest.approx(3, rel=1e-4)
    inp.icp_state = "dried 200 °C"; e1, e2 = D.method_e(inp)
    assert e2.L == pytest.approx(3, rel=1e-3)


def test_custom_node_equals_library_node():
    lib = _inp("HKUST-1", tga_W_plateau=250.0, tga_W_end=100.0)
    cus = copy.deepcopy(lib); cus.node = N.CUSTOM; cus.node_custom = dict(nM=2, val=2, core_O=0, core_OH=0, CN=4, metal="Cu")
    assert D.method_a(cus).L == pytest.approx(D.method_a(lib).L) and cus.Q() == 4
    cus.node_custom["residue"] = ("Cu2O", 2)          # e.g. under N2
    assert D.method_a(cus).L != pytest.approx(D.method_a(lib).L)


def test_zr_hf_unchanged_against_v210_snapshot():
    """Byte-identical to the v2.10 engine on every bundled input x metal x framework (snapshot recorded before the refactor)."""
    snap = os.path.join(ROOT, "tests", "data", "defects_v210_snapshot.json")
    if not os.path.exists(snap):
        pytest.skip("snapshot missing")
    sys.path.insert(0, os.path.join(ROOT, "tests", "data"))
    import defects_snapshot as G
    assert G.cases() == json.load(open(snap))


def test_node_library_table():
    rows = N.table_rows()
    assert len(rows) >= 9 and {r["charge for linkers (sites)"] for r in rows} >= {12, 4, 6, 2, 7}
    assert all(r["TGA residue"] for r in rows)


def test_defect_page_any_framework_apptest(tmp_path):
    from streamlit.testing.v1 import AppTest
    os.environ.update(MFS_STORE=str(tmp_path / "s.json"), MFS_DB=str(tmp_path / "r.db"), MFS_EXPORT_DIR=str(tmp_path))
    at = AppTest.from_file(os.path.join(ROOT, "app.py"), default_timeout=300); at.run()
    at.sidebar.radio[0].set_value("Defect calculator").run()
    assert at.selectbox(key="def_fw").value == "UiO-66"            # neutral default
    plat = 100 * molar_mass("Cu2(C9H3O6)1.3333333333") / molar_mass("Cu2O2")
    at.selectbox(key="def_fw").set_value("HKUST-1").run()
    at.number_input(key="dn_wpl").set_value(plat).run(); at.number_input(key="dn_wend").set_value(100.0).run()
    assert not at.exception, [e.value for e in at.exception][:1]
    m = {x.label: x.value for x in at.metric}
    assert m["Linkers / node (median)"] == "1.33", m
    assert any(n.label == "W residue (CuO)" for n in at.number_input)
    assert any("Node Cu2" in str(c.value) for c in at.caption)
    k = [s.key for s in at.selectbox if s.key and s.key.startswith("def_node")][0]
    at.selectbox(key=k).set_value(N.CUSTOM).run()
    assert not at.exception, [e.value for e in at.exception][:1]
    assert any("charges for ligands" in str(c.value) for c in at.caption)


def test_guess_framework_from_sample_name():
    sys.path.insert(0, ROOT)
    import importlib
    src = open(os.path.join(ROOT, "mfs", "pages", "defects.py"), encoding="utf-8").read()     # v2.13: page code lives in mfs/pages
    ns = {}
    import re
    code = src[src.index("def guess_framework"):src.index("def default_defect_input")]
    exec("from mfs.chem import FRAMEWORKS\n" + code, ns)
    g = ns["guess_framework"]
    assert g("MOF-801(Hf)") == ("MOF-801", "Hf")
    assert g("UiO-66-NH2 Zr batch 3") == ("UiO-66-NH2", "Zr")
    assert g("HKUST-1_activated")[0] == "HKUST-1"
    assert g("my sample") == (None, None)
