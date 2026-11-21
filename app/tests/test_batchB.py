"""Batch B: FTIR/Raman band normalisation, peak labels, insets, difference spectra; XPS Shirley/Tougaard + doublets +
atomic %; TGA steps, plateau → L, TGA–MS; DSC cycle-to-cycle comparison and cp (sapphire); NMR ratio string."""
import os
import numpy as np, pandas as pd, pytest
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX = os.path.join(HERE, "example_data")


def test_spectra_tools():
    from mfs import spectra as SP
    x = np.linspace(400, 2000, 3201)
    g = lambda c, w, a: a * np.exp(-0.5 * ((x - c) / w) ** 2)
    y1 = g(1580, 12, 2.0) + g(1400, 10, 1.0) + 0.1
    y2 = g(1580, 12, 1.0) + g(1400, 10, 0.5) + g(1050, 8, 0.4) + 0.1
    n1, n2 = SP.normalise_to_band(x, y1, 1540, 1620), SP.normalise_to_band(x, y2, 1540, 1620)
    assert abs(SP.band_value(x, n1, 1540, 1620) - 1) < 1e-9 and abs(SP.band_value(x, n2, 1540, 1620) - 1) < 1e-9
    pk = SP.pick_peaks(x, y2, 0.05, 5)
    assert all(any(abs(p - c) < 2 for p in pk) for c in (1050, 1400, 1580))
    f = SP.auto_factor((x, y1), (x, y2), 1500, 1650); assert abs(f - 0.5) < 0.01
    gx, d = SP.difference((x, y1), (x, y2), f)
    assert abs(d[np.argmin(abs(gx - 1050))] - 0.45) < 0.01 and abs(d[np.argmin(abs(gx - 1580))] - 0.05) < 0.02   # 0.1 − 0.5·0.1 offset
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from mfs import plots as P
    fig, ax = plt.subplots()
    P.draw_stacked(ax, [("a", x, n1), ("b", x, n2)], "ftir_A", scale="common")
    SP.label_peaks(ax, x, n2, pk); ia = SP.add_inset(ax, lambda a_: P.draw_stacked(a_, [("a", x, n1)], "ftir_A"), (1100, 1000))
    assert len([t for t in ax.texts if t.get_gid() == "mfs_peaklabel"]) == len(pk) and ia.get_gid() == "mfs_inset"
    info = {e[0]: e for e in ax._mfs_stack}; assert info["a"][2] == info["b"][2]      # common scale


def test_xps_backgrounds_doublets_atomic():
    from mfs import xps as X, readers as R
    x = np.arange(160, 174, 0.05)
    y = X.gl(x, 168.6, 1.4, 0.3, 1000) + X.gl(x, 169.78, 1.4, 0.3, 500) + 200
    f = X.fit(x, y, [{"name": "S", "region": "S 2p", "center": 168.4}], "linear")
    r = f["rows"][0]
    assert abs(r["BE major (eV)"] - 168.6) < 0.01 and abs(r["area (total)"] - 1500) < 30 and abs(r["BE minor (eV)"] - r["BE major (eV)"] - 1.18) < 1e-9
    d = R.read_xps(os.path.join(EX, "XPS", "SampleC_S2p_CasaXPS.csv"), "S2p.csv")
    assert X.guess_region(d.x.values) == "S 2p"
    fs = X.fit(d.x.values, d.y.values, [{"name": "S", "region": "S 2p", "center": 168.5}], "shirley")
    assert abs(fs["rows"][0]["area (total)"] - 2250) / 2250 < 0.03
    ft = X.fit(d.x.values, d.y.values, [{"name": "S", "region": "S 2p", "center": 168.5}], "tougaard")
    assert np.all(np.isfinite(ft["bg"])) and ft["rows"][0]["area (total)"] > 0
    b = X.shirley(d.x.values, d.y.values); assert abs(b[0] - d.y.values[:3].mean()) < 1 or abs(b[-1] - d.y.values[-3:].mean()) < 1
    at = X.atomic_percent({"Zr 3d": 7.04, "O 1s": 2 * 2.93})
    assert abs(at["O 1s"] - 200 / 3) < 1e-9
    for k, v in {"Zr 3d": (2.43, 2 / 3), "Hf 4f": (1.71, 0.75), "S 2p": (1.18, 0.5)}.items():
        assert X.DOUBLETS[k]["split"] == v[0] and abs(X.DOUBLETS[k]["ratio"] - v[1]) < 1e-12


def test_tga_steps_ms():
    from mfs import readers as R, analysis as A, tga as TG
    t = R.read_tga(os.path.join(EX, "TGA", "SampleA_Netzsch.txt"), "SampleA_Netzsch.txt"); a = A.tga_analysis(t["T"], t["w"], 21)
    st = TG.steps(a, 1.0)
    assert len(st) == 3 and abs(sum(s["loss"] for s in st) - (a["W_initial"] - a["W_end"])) < 0.5
    assert 450 < st[-1]["T_peak"] < 550
    ms = TG.read_ms(os.path.join(EX, "TGA", "SampleA_QMS_ion_currents_EXAMPLE.csv"))
    assert {"T", "m/z 18", "m/z 44", "m/z 28"} <= set(ms.columns)
    assert abs(ms["T"][ms["m/z 44"].idxmax()] - 504) < 5


def test_dsc_cp_and_cycle_compare():
    from mfs import dsc
    Dp = os.path.join(EX, "DSC", "Cp_sapphire")
    rd = lambda f: dsc.read_dsc(os.path.join(Dp, f))
    b, s, x = rd("Cp_blank_EXAMPLE.csv"), rd("Cp_sapphire_EXAMPLE.csv"), rd("Cp_sample_MOF-801_EXAMPLE.csv")
    assert s["meta"]["mass_mg"] == 25.30 and x["meta"]["mass_mg"] == 9.80
    r = dsc.cp_from_sapphire(x["curves"][0], s["curves"][0], b["curves"][0], 9.80, 25.30)
    tr = 0.62 + 1.6e-3 * (r["T"] - 25)
    assert np.max(np.abs(r["Cp"] - tr)) < 0.015 and np.median(np.abs(r["Cp_direct"] - tr)) < 0.01
    assert abs(dsc.cp_sapphire(26.85) - 0.7792) < 1e-4
    d = dsc.read_dsc(os.path.join(EX, "DSC", "full_cycle", "SG_DSC_heat-cool-heat_EXAMPLE.csv"))
    peaks = {}
    for i, c in enumerate(d["curves"]):
        c["uid"] = str(i); c["sample"] = "x"; peaks[str(i)] = dsc.evaluate_all(dict(c, id=i), "linear", d["meta"]["mass_mg"])
    rows = dsc.cycle_compare(d["curves"], peaks)
    h2 = [r_ for r_ in rows if r_["kind"] == "heating"][1]
    assert abs(rows[0]["ΔH (J/g)"] + 105.6) < 0.1 and abs(h2["ΔH / ΔH(cycle 1)"] - 0.93) < 0.01 and abs(h2["ΔT onset vs cycle 1 (K)"] + 0.8) < 0.05


def _at(monkeypatch, tmp_path):
    monkeypatch.setenv("MFS_STORE", str(tmp_path / "s.json"))
    from streamlit.testing.v1 import AppTest
    a = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=300); a.run(); return a


def _ok(a):
    assert not a.exception, [str(e.value)[:400] for e in a.exception]


def test_app_ftir_tools(monkeypatch, tmp_path):
    a = _at(monkeypatch, tmp_path); a.button(key="nv_" + "FTIR / Raman").click().run(); _ok(a)
    for k in ("ft_nb", "ft_pl", "ft_ins"):
        a.checkbox(key=k).check().run(); _ok(a)
    assert any("Difference" in e.label for e in a.expander)
    a.radio(key="ft_dfm").set_value("least squares in a window").run(); _ok(a)
    assert any("peaks (cm⁻¹)" in d.value.columns for d in a.dataframe)


def test_app_xps_multi(monkeypatch, tmp_path):
    a = _at(monkeypatch, tmp_path); a.button(key="nv_" + "XPS").click().run(); _ok(a)
    a.checkbox(key="xps_multi").check().run(); _ok(a)
    at = [d.value for d in a.dataframe if "S 2p at.%" in d.value.columns]
    assert at and 0 < at[0]["S/(Zr+Hf)"].iloc[0] < 2
    a.radio(key="xps_bg").set_value("Tougaard").run(); _ok(a)
    a.selectbox(key="xps_fit_pick").set_value(a.selectbox(key="xps_fit_pick").options[-1]).run(); _ok(a)


def test_app_tga_dsc(monkeypatch, tmp_path):
    a = _at(monkeypatch, tmp_path); a.button(key="nv_" + "TGA").click().run(); _ok(a)
    assert any("linkers / node L" in d.value.columns for d in a.dataframe)
    a.button(key="nv_" + "DSC").click().run(); _ok(a)
    a.checkbox(key="dsc_hch_ex").check().run(); _ok(a)
    assert any("ΔH / ΔH(cycle 1)" in d.value.columns for d in a.dataframe)
    cp = [d.value for d in a.dataframe if "cp 100 °C" in d.value.columns]
    assert cp and abs(cp[0]["cp 100 °C"].iloc[0] - 0.74) < 0.015


def test_app_nmr_ratio_string(monkeypatch, tmp_path):
    a = _at(monkeypatch, tmp_path); a.session_state["nmr_mode"] = "⚗️ Digestion / composition (ratios, qNMR)"   # v2.13
    a.button(key="nv_" + "1H NMR (digest)").click().run(); _ok(a)
    assert any("composition" in d.value.columns for d in a.dataframe)
