"""Batch A: shared sample store, SI export, N2 computed PSD / t-plot, kinetics CI / Arrhenius / Eyring / TOF / pD,
PXRD Scherrer with instrument correction, crystallinity index, peak shifts."""
import io, os, glob, json
import numpy as np, pandas as pd, pytest
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX = os.path.join(HERE, "example_data")


# ---------------- modules ----------------
def test_samplestore_roundtrip(tmp_path):
    from mfs import samplestore as SS
    st = {}
    SS.put(st, "Sample C", "tga", {"W_end": np.float64(41.2), "bad": float("nan")}, "TGA")
    SS.put(st, "samplec", "icp", {"Zr_wt": 30.1}, "ICP")
    assert len(SS.ids(st)) == 1
    r = SS.get(st, "Sample C")
    assert set(r["sections"]) == {"tga", "icp"} and r["sections"]["tga"]["values"]["bad"] is None
    p = tmp_path / "s.json"; SS.save(st, str(p))
    st2 = {}; SS.load(st2, str(p)); assert SS.get(st2, "samplec", "icp")["Zr_wt"] == 30.1
    assert SS.summary_rows(st2)[0]["tga"]


def test_si_export_xlsx_docx():
    from mfs import samplestore as SS, si_export as SI
    st = {}
    for nm, L in (("Sample A", 5.2), ("Sample B", 4.6)):
        SS.put(st, nm, "defects", {"framework": "MOF-801", "metal": "Hf", "L": L, "x": 6 - L, "m": 0.4, "OS": 2 * (6 - L) - 0.4,
                                   "M": 1650.0, "formula": "Hf6O4(OH)4(fum)5", "methods": "TGA; NMR"})
        SS.put(st, nm, "n2", {"S_BET": 800.0, "V_total": 0.4})
    df = SI.defect_rows(st["sample_store"] if "sample_store" in st else st)
    assert len(df) == 2 and any("S_BET" in c for c in df.columns)
    x = SI.to_xlsx(df, notes="n"); assert x[:2] == b"PK"
    pd.read_excel(io.BytesIO(x))
    d = SI.to_docx(df, notes="n"); assert d[:2] == b"PK"
    docx = pytest.importorskip("docx")      # python-docx is optional: the app writes DOCX without it (the Windows build does not ship it)
    doc = docx.Document(io.BytesIO(d)); assert len(doc.tables[0].rows) == 3


def test_kinetics_ci_and_activation():
    from mfs import kinetics as K
    rng = np.random.default_rng(1)
    t = np.repeat([0, 2, 4, 8, 12, 24, 36, 48.0], 3)
    c = 2 * np.exp(-0.05 * t) * (1 + rng.normal(0, 0.01, len(t)))
    f = K.first_order_ci(t, c)
    assert f["k_lo"] < 0.05 < f["k_hi"] and f["k_ci"] > f["k_se"]
    lo, hi = f["band"](np.array([10.0])); assert lo[0] < f["model"](10.0) < hi[0]
    ir = K.initial_rate(t, c, 0.2); assert abs(ir["r0"] - 0.1) < 0.02
    T = np.array([40, 50, 60, 70.0]); Ea = 72e3; k = np.exp(23) * np.exp(-Ea / 8.314 / (T + 273.15)) / 3600
    a = K.arrhenius(T, k); assert abs(a["Ea_kJ"] - 72) < 0.1
    e = K.eyring(T, k); assert abs(e["dH_kJ"] - (72 - 8.314e-3 * 328)) < 0.6
    to = K.tof(0.05 / 3600, 2.0, 0.6, 2.0, 1.5, 1650.0)
    assert abs(to["TOF_h"] - 0.05 * 2e-3 * 0.6e-3 / (2e-3 / 1650 * 1.5)) < 1e-9
    pD = np.linspace(4, 10, 11); kk = 0.05 / (1 + 10 ** (6.2 - pD))
    fr = K.fit_pd_profile(pD, kk); assert abs(fr["params"]["pKa"] - 6.2) < 0.01
    kb = K.pd_profile_bell(pD, 0.08, 5.5, 8.5); fb = K.fit_pd_profile(pD, kb, "bell")
    assert abs(fb["params"]["pKa1"] - 5.5) < 0.05 and abs(fb["params"]["pKa2"] - 8.5) < 0.05


def test_xrd_tools():
    from mfs import xrd as X
    x = np.linspace(5, 15, 4001)
    g = lambda c, w, a=100: a * np.exp(-4 * np.log(2) * (x - c) ** 2 / w ** 2)
    y1 = g(7.40, 0.15) + g(8.55, 0.15) + 5; y2 = g(7.35, 0.15, 50) + g(8.50, 0.15, 50) + 5
    r = X.fit_peak(x, y1, 7.4); assert abs(r["fwhm"] - 0.15) < 0.005 and abs(r["center"] - 7.4) < 1e-3
    s = X.scherrer_size(r["center"], r["fwhm"], 0.05); s0 = X.scherrer_size(r["center"], r["fwhm"], 0.0)
    assert s["size_nm"] > s0["size_nm"] > 40
    assert np.isnan(X.scherrer_size(7.4, 0.04, 0.05)["size_nm"])
    assert abs(X.instrument_correction(0.15, 0.05, "lorentzian") - 0.10) < 1e-12
    ci = X.crystallinity_index(x, y2, x, y1, [(7.0, 7.8), (8.1, 9.0)]); assert abs(ci["CI_pct"] - 50) < 2
    tr = X.track_shifts([("ref", x, y1), ("s", x, y2)], [7.4, 8.55], hkl=[(1, 1, 1), (2, 0, 0)])
    d = [r_ for r_ in tr if r_["sample"] == "s"]
    assert abs(d[0]["Δ2θ vs ref (°)"] + 0.05) < 2e-3 and "a cubic (Å)" in d[0]


def test_psd_methods_on_example_isotherm():
    from mfs import psd as PS, readers as R
    d = R.read_sorption(os.path.join(EX, "Sorption", "SampleA_3Flex_MicroActive.xlsx"), "SampleA_3Flex_MicroActive.xlsx")
    iso = d["iso"]; a_ = iso[iso.branch.astype(str).str.startswith("ads")]
    p, q = a_.p_rel.values.astype(float), a_.q.values.astype(float)
    r = PS.dft_fit(p, q)
    assert r["residual_rel"] < 0.15 and 0.5 < r["w"][np.argmax(r["dVdw"])] < 2.5
    h = PS.hk_psd(p, q, "cylinder", "AlSi oxide ion (SF)"); assert np.all(np.isfinite(h["w"]))
    assert abs(PS.harkins_jura(np.array([0.5]))[0] - 0.646) < 0.005   # t = sqrt(13.99/(0.034 − log p)) Å


# ---------------- app ----------------
def _goto(at, page):
    if page == "1H NMR (digest)":
        at.session_state["nmr_mode"] = "⚗️ Digestion / composition (ratios, qNMR)"   # v2.13
    at.sidebar.radio[0].set_value(page).run()
    assert not at.exception, [e.value for e in at.exception][:2]
    return at


@pytest.fixture
def at(tmp_path, monkeypatch):
    monkeypatch.setenv("MFS_STORE", str(tmp_path / "store.json"))
    from streamlit.testing.v1 import AppTest
    a = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=300); a.run()
    return a


def test_app_store_flow_to_defects(at, tmp_path):
    for pg in ("TGA", "1H NMR (digest)", "ICP-OES / XRF", "N2 sorption (BET/PSD)"):
        _goto(at, pg)
    store = at.session_state["sample_store"]
    secs = {s for v in store.values() for s in v["sections"]}
    assert {"tga", "nmr", "n2"} <= secs, secs
    assert os.path.exists(tmp_path / "store.json")
    _goto(at, "Defect calculator")
    opts = at.selectbox(key="def_pull").options
    assert len(opts) > 1
    with_tga = [k for k, v in store.items() if "tga" in v["sections"]]
    at.selectbox(key="def_pull").set_value(store[with_tga[0]]["name"] if False else with_tga[0]).run()
    assert not at.exception, [e.value for e in at.exception][:2]
    assert any("Pulled for" in s.value for s in at.success)
    assert any("SI table" in e.label for e in at.expander)


def test_app_n2_psd_tplot(at):
    _goto(at, "N2 sorption (BET/PSD)")
    assert any("Rouquerol" in (c.value if hasattr(c, "value") else "") for c in at.markdown) or len(at.dataframe) >= 2
    for m in at.selectbox(key="sorp_psdm").options:
        at.selectbox(key="sorp_psdm").set_value(m).run()
        assert not at.exception, (m, [e.value for e in at.exception][:2])
    at.checkbox(key="sorp_sbet_leg").check().run(); assert not at.exception


def test_app_kinetics_tseries(at):
    _goto(at, "Kinetics & adsorption")
    assert not at.error
    at.checkbox(key="kin_tser").check().run()
    assert not at.exception, [e.value for e in at.exception][:2]
    md = " ".join(m.value for m in at.markdown)
    assert "Ea = **7" in md and "ΔH‡" in md and "**pKa = 6." in md, md[:600]


def test_app_pxrd_tools(at):
    _goto(at, "PXRD")
    assert not at.exception and not at.error
    cap = [d for d in at.dataframe]
    assert len(cap) >= 3
