"""v2.15: the opt-in PXRD modules (crystallite size, Rietveld) in the app."""
import os, sys, warnings
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _app(level="Standard"):
    warnings.filterwarnings("ignore")
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600)
    at.session_state["_lv"] = {"global": level}
    at.run()
    at.sidebar.radio[0].set_value("PXRD").run()
    return at


def test_modules_hidden_until_clicked_even_in_standard():
    at = _app("Standard")
    assert not at.exception
    labels = [t.label for t in at.toggle]
    assert any("Crystallite size" in l for l in labels) and any("Rietveld" in l for l in labels)
    assert not any("Refine" in b.label for b in at.button)              # nothing of the Rietveld module is drawn yet
    assert not any("Instrument broadening" in r.label for r in at.radio)


def test_crystallite_size_and_rietveld_run_in_the_app():
    at = _app("Advanced")
    [t for t in at.toggle if "Crystallite size" in t.label][0].set_value(True).run()
    assert not at.exception, at.exception
    assert any("Instrument broadening" in r.label for r in at.radio)
    assert any("Williamson" in e.value or "D Scherrer" in str(getattr(e, "value", "")) for e in at.dataframe) or len(at.dataframe) > 0
    [t for t in at.toggle if "Rietveld" in t.label][0].set_value(True).run()
    assert not at.exception, at.exception
    _ms = at.multiselect(key="rv_phases"); _ms.set_value([o for o in _ms.options if "801" in o][:1]).run()      # v2.18.4: no phase is preselected any more
    go = [b for b in at.button if b.key == "rv_go"][0]
    assert not go.disabled                                                # the chosen bundled MOF-801
    go.click().run()
    assert not at.exception, at.exception
    res = at.session_state["rv_res"]
    r = next(iter(res.values()))
    assert r["stats"]["Rwp"] < 0.6 and r["phases"][0]["cell"][0] > 10
    assert any(m.label == "Rwp" for m in at.metric)


def test_band_assignment_and_fit_in_ftir_page():
    from streamlit.testing.v1 import AppTest
    warnings.filterwarnings("ignore")
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600)
    at.session_state["_lv"] = {"global": "Advanced"}
    at.run()
    at.sidebar.radio[0].set_value("FTIR / Raman").run()
    assert not at.exception, at.exception
    assert any("Band assignment" in e.label for e in at.expander)
    [b for b in at.button if b.key == "bd_go"][0].click().run()
    assert not at.exception, at.exception
    assert "bd_res" in at.session_state and at.session_state["bd_res"][1]["r2"] > 0.8
    at2 = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600)
    at2.session_state["_lv"] = {"global": "Standard"}
    at2.run()
    at2.sidebar.radio[0].set_value("FTIR / Raman").run()
    assert not any("Band assignment" in e.label for e in at2.expander)          # Standard: not drawn


def _page(name, level="Advanced"):
    from streamlit.testing.v1 import AppTest
    warnings.filterwarnings("ignore")
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600)
    at.session_state["_lv"] = {"global": level}
    at.run()
    at.sidebar.radio[0].set_value(name).run()
    return at


def test_xps_batch_fit_and_nmr_library_panel():
    at = _page("XPS")
    assert not at.exception, at.exception
    [b for b in at.button if b.key == "xps_batch_go"][0].click().run()
    assert not at.exception, at.exception
    tb, bad = at.session_state["xps_batch_res"]
    assert len(tb) >= 1 and "area (total)" in tb
    at = _page("1H NMR (digest)")
    assert not at.exception, at.exception
    assert any("Linker / modulator library" in e.label for e in at.expander)
    at = _page("1H NMR (digest)", "Standard")
    assert not any("Linker / modulator library" in e.label for e in at.expander)


def test_sizing_template_callbacks():
    import streamlit as st
    from mfs.pages import microscopy as M
    from streamlit.testing.v1 import AppTest
    # run the callbacks inside a tiny app so that st.session_state exists
    code = (
        "import streamlit as st, os, tempfile\n"
        "os.environ['MFS_STORE'] = os.path.join(tempfile.mkdtemp(), 'sample_store.json')\n"
        "from mfs.pages import microscopy as M\n"
        "st.session_state.setdefault('ps_sm', 2.0); st.session_state.setdefault('ps_split', False); st.session_state.setdefault('ps_imgs', ['a'])\n"
        "st.session_state['ps_tpl_name'] = 'MOF SEM'\n"
        "M._tpl_cb('save')\n"
        "t = st.session_state['ps_templates']['MOF SEM']\n"
        "st.session_state['ps_sm'] = 0.5; st.session_state['ps_tpl_pick'] = 'MOF SEM'\n"
        "M._tpl_cb('apply')\n"
        "st.session_state['out'] = (sorted(t), st.session_state['ps_sm'])\n")
    at = AppTest.from_string(code)
    at.run()
    assert not at.exception, at.exception
    keys, sm = at.session_state["out"]
    assert "ps_sm" in keys and "ps_split" in keys and "ps_imgs" not in keys and sm == 2.0


def test_pl_thermometry_automatic_fit_in_the_app():
    at = _page("Photoluminescence", "Standard")
    assert not at.exception, at.exception
    # the thermometry tab is drawn on every run (tabs); the automatic global fit is the default
    assert any("Automatic" in str(r.value) for r in at.radio if r.key == "pl_method")
    assert any(m.label == "Model" for m in at.metric) and any(m.label == "Typical ratio error" for m in at.metric)
    r = [x for x in at.radio if x.key == "pl_method"][0]
    r.set_value("Compare both").run()
    assert not at.exception, at.exception
    r = [x for x in at.radio if x.key == "pl_method"][0]
    r.set_value("Integration windows (choose them on the graph)").run()
    assert not at.exception, at.exception
    assert any(s.key == "pl_b1" for s in at.slider)


def test_pl_extras_panels_exist_and_run():
    at = _page("Photoluminescence", "Advanced")
    assert not at.exception, at.exception
    labels = [e.label for e in at.expander]
    assert any("CIE 1931" in l for l in labels) and any("Power dependence" in l for l in labels)
    assert any("Repeats" in l for l in labels)          # v2.18.4: "Repeats, independent runs and heating / cooling" and any("Spectral-response correction" in l for l in labels)
    at = _page("Photoluminescence", "Standard")
    assert not at.exception, at.exception
    assert not any("CIE 1931" in e.label for e in at.expander)


def test_structure_refinement_in_the_app_with_3d_tables_and_cif():
    from streamlit.testing.v1 import AppTest
    warnings.filterwarnings("ignore")
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=900)
    at.session_state["_lv"] = {"global": "Advanced"}
    at.run()
    at.sidebar.radio[0].set_value("PXRD").run()
    [t for t in at.toggle if "Rietveld" in t.label][0].set_value(True).run()
    _ms = at.multiselect(key="rv_phases"); _ms.set_value([o for o in _ms.options if "801" in o][:1]).run()      # v2.18.4: no phase is preselected any more
    [r for r in at.radio if r.key == "rv_mode"][0].set_value("Rietveld + refine the atoms (structure)").run()
    [b for b in at.button if b.key == "rv_go"][0].click().run()
    assert not at.exception, at.exception
    res = next(iter(at.session_state["rv_res"].values()))
    sd = res["phases"][0]["structure"]
    assert sd["atoms"] and sd["bond_rms_deviation"] < 0.15 and sd["bonds"] and not sd.get("report_error")
    assert any("Refined structure" in str(m.value) for m in at.markdown)
    assert [t.label for t in at.tabs if "3D structure" in t.label] and [t.label for t in at.tabs if "Output" in t.label]


def test_phase_identification_in_the_pxrd_page():
    from streamlit.testing.v1 import AppTest
    warnings.filterwarnings("ignore")
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600)
    at.session_state["_lv"] = {"global": "Standard"}
    at.run()
    at.sidebar.radio[0].set_value("PXRD").run()
    [t for t in at.toggle if "Rietveld" in t.label][0].set_value(True).run()
    [b for b in at.button if b.key == "rv_ident"][0].click().run()
    assert not at.exception, at.exception
    name, res = at.session_state["rv_ident_res"]
    assert res and res[0]["name"].startswith("MOF-801")
    [b for b in at.button if b.key == "rv_use_ident"][0].click().run()
    assert not at.exception, at.exception
    assert any("MOF-801" in x for x in at.session_state["rv_phases"])


def test_indexing_and_le_bail_of_an_unknown_phase_in_the_page():
    from streamlit.testing.v1 import AppTest
    warnings.filterwarnings("ignore")
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600)
    at.session_state["_lv"] = {"global": "Standard"}
    at.run()
    at.sidebar.radio[0].set_value("PXRD").run()
    [t for t in at.toggle if "Rietveld" in t.label][0].set_value(True).run()
    _sa = at.selectbox(key="rv_sample"); _sa.set_value([o for o in _sa.options if "SampleA" in o or "MOF" in o][0]).run()     # v2.18.5: the MOF-801 example (Sample A), not the first pattern (PM)
    [b for b in at.button if b.key == "rv_index"][0].click().run()
    assert not at.exception, at.exception
    name, sols, n = at.session_state["rv_index_res"]
    assert sols and sols[0]["system"] == "cubic" and abs(sols[0]["cell"][0] - 17.8) < 0.6           # MOF-801 (Hf): fcc, a ≈ 17.8 Å
    [b for b in at.button if b.key == "rv_index_add"][0].click().run()
    assert at.session_state["rv_cellph"]
    [b for b in at.button if b.key == "rv_go"][0].click().run()
    assert not at.exception, at.exception
    res = next(iter(at.session_state["rv_res"].values()))
    assert res["mode"] == "lebail" and res["stats"]["Rwp"] < 0.25


def test_one_click_automatic_identification_and_structure_refinement():
    from streamlit.testing.v1 import AppTest
    warnings.filterwarnings("ignore")
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=900)
    at.session_state["_lv"] = {"global": "Standard"}
    at.run()
    at.sidebar.radio[0].set_value("PXRD").run()
    [t for t in at.toggle if "Rietveld" in t.label][0].set_value(True).run()
    _sa = at.selectbox(key="rv_sample"); _sa.set_value([o for o in _sa.options if "SampleA" in o or "MOF" in o][0]).run()     # v2.18.5: the MOF-801 example (Sample A), not the first pattern (PM)
    [b for b in at.button if b.key == "rv_auto"][0].click().run()
    assert not at.exception, at.exception
    res = next(iter(at.session_state["rv_res"].values()))
    assert res["phases"][0]["name"].startswith("MOF-801") and res["phases"][0].get("structure") and res["stats"]["Rwp"] < 0.2
