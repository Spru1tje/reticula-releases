"""v2.18.4: PL quality, uncertainty, manifold fit, decays, colour, transitions, header metadata - known-answer tests."""
import math
import os
import sys
import warnings

import numpy as np
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from mfs import pl as PL, plqa as QA, plratio as PR, pldecay as DC, plextras as PX, pllines as LN  # noqa: E402

KB = 8.617333262e-5
CM = 8065.54429


def _series(dE_cm1=750.0, alpha=8.0, T=np.arange(293.0, 354.0, 10.0), peak=4e4, seed=0, lines=False):
    """Two bands (or two groups of overlapping lines) whose area ratio follows the Boltzmann law exactly; Poisson noise on a pedestal."""
    rng = np.random.default_rng(seed)
    x = np.linspace(500, 580, 401)
    g = lambda c, w: np.exp(-0.5 * ((x - c) / w) ** 2) / (w * math.sqrt(2 * math.pi))
    out = []
    for t in T:
        r = alpha * math.exp(-dE_cm1 / CM / (KB * t))
        a2 = peak * 6.0
        a1 = r * a2
        if lines:
            y = a1 * (0.5 * g(520, 1.6) + 0.3 * g(524, 1.6) + 0.2 * g(527, 1.6)) + a2 * (0.45 * g(545, 1.6) + 0.35 * g(549, 1.6) + 0.2 * g(553, 1.6))
        else:
            y = a1 * g(525, 4.0) + a2 * g(548, 4.0)
        raw = rng.poisson(y + 300.0).astype(float)
        out.append({"T": float(t), "x": x, "y": raw - 300.0, "var": raw, "truth": r})
    return out


# ------------------------------------------------------------------------------------------------------------------------------ ratio and Boltzmann
def test_window_ratio_error_matches_its_scatter():
    sp = _series()
    reps = [QA.window_ratio(_series(seed=k), (510, 536), (536, 562))[2] for k in range(30)]
    _, _, D, sD = QA.window_ratio(sp, (510, 536), (536, 562))
    assert np.all(sD > 0) and abs(np.std(np.array(reps), axis=0).mean() / sD.mean() - 1) < 0.35


def test_manifold_fit_sums_the_lines_and_recovers_the_ratio():
    sp = _series(lines=True)
    f = PR.fit_manifolds(sp, (515, 532), (540, 560), True, 505, 570, 1)
    assert f is not None and f["groups"] and len(f["groups"][0]) >= 2 and len(f["groups"][1]) >= 2
    D = np.array([r["ratio"] for r in f["rows"]])
    truth = np.array([s["truth"] for s in sp])
    assert np.allclose(D, truth, rtol=0.03)
    assert PR.flag_share(f) <= 0.3


def test_information_and_checks_flag_a_flat_ratio():
    T = np.arange(293.0, 354.0, 10.0)
    D = 0.8 + np.random.default_rng(1).normal(0, 0.002, len(T))
    sD = np.full(len(T), 0.002)
    bf = PL.fit_boltzmann(T, np.abs(D), sigma=sD)
    ch = QA.boltzmann_checks(T, D, sD, bf)
    assert QA.verdict(ch) == "low" and any(c["name"] == "Thermometric information" and c["status"] == "fail" for c in ch)


def test_good_thermometer_gets_high_confidence_and_the_gap_back():
    sp = _series(seed=3)
    _, _, D, sD = QA.window_ratio(sp, (505, 536.5), (536.5, 575))
    T = np.array([s["T"] for s in sp])
    bf = PL.fit_boltzmann(T, D, sigma=sD)
    assert abs(bf["dE_cm1"] - 750) < 4 * bf["dE_err_cm1"] + 5
    gap = QA.spectroscopic_gap_cm1(525, 548)
    ch = QA.boltzmann_checks(T, D, sD, bf, gap, 0, len(T))
    assert QA.verdict(ch) in ("high", "moderate"), [(c["name"], c["status"], c["text"]) for c in ch]
    assert not any(c["status"] == "fail" for c in ch)


def test_monte_carlo_budget_agrees_with_the_fit_error_and_adds_the_sensor():
    sp = _series(seed=4)
    T = np.array([s["T"] for s in sp])
    _, _, D, sD = QA.window_ratio(sp, (505, 536.5), (536.5, 575))
    bf = PL.fit_boltzmann(T, D, sigma=sD)
    B0 = QA.budget(T, D, sD, bf)
    assert abs(B0["rows"][0]["u_dE_cm1"] / bf["dE_err_cm1"] - 1) < 0.25          # Monte Carlo = linear propagation for small errors
    B1 = QA.budget(T, D, sD, bf, sT_rand=1.0, sT_sys=2.0)
    assert B1["u_dE_cm1"] > B0["u_dE_cm1"] and np.all(B1["T_cal"] > B0["T_cal"])
    # resolution: (sD/D)/Sr
    Sr = 100 * bf["dE_eV"] / (KB * T ** 2)
    assert np.allclose(B0["T_res"], (sD / D) / (Sr / 100))
    # a processing alternative that shifts the gap shows up in the budget
    D2 = D * np.exp(-30 / CM / KB * (1 / T - 1 / T.mean()))
    B2 = QA.budget(T, D, sD, bf, alt_sets=[("test", D2)])
    row = [r for r in B2["rows"] if r["component"] == "processing: test"][0]
    assert abs(row["u_dE_cm1"] - 30) < 3


def test_predict_T_inverts_the_calibration():
    T = np.array([300.0, 320.0])
    D = 8 * np.exp(-0.09 / (KB * T))
    assert np.allclose(QA.predict_T(D, 8.0, 0.09), T)


def test_hysteresis_and_laser_heating():
    dE = 0.093
    Th = np.array([293, 303, 313, 323, 333.0])
    Tc = Th[::-1][1:]
    T = np.concatenate([Th, Tc])
    Dtrue = 8 * np.exp(-dE / (KB * T))
    D = Dtrue.copy()
    D[len(Th):] = 8 * np.exp(-dE / (KB * (Tc + 2.0)))                       # the cooling branch is 2 K warmer
    h = QA.hysteresis(T, D, np.arange(len(T)), dE)
    assert h["status"] == "ok" and abs(h["max_dT"] - 2.0) < 0.3
    P = np.array([0.1, 0.2, 0.4, 0.8])
    Dp = 8 * np.exp(-dE / (KB * (300 + 5.0 * P)))                            # 5 K per W
    lh = QA.laser_heating(P, Dp, 8.0, dE)
    assert lh["status"] == "ok" and abs(lh["slope"] - 5.0) < 0.05 and abs(lh["T0"] - 300) < 0.05


def test_saturation_detection():
    x = np.linspace(500, 600, 200)
    y = 60000 * np.exp(-0.5 * ((x - 550) / 5) ** 2)
    clip = np.minimum(y, 65535 * 0.9)
    clip[95:100] = 52000.0
    ok = QA.saturation([{"name": "fine", "I": y * 0.5, "wl": x, "meta": {"detector": "CCD"}}])[0]
    flat = QA.saturation([{"name": "flat", "I": np.minimum(y, 52000.0), "wl": x, "meta": {"detector": "CCD"}}])[0]
    adc = QA.saturation([{"name": "adc", "I": np.minimum(y * 2, 65535.0), "wl": x, "meta": {"detector": "CCD"}}])[0]
    pmt = QA.saturation([{"name": "pmt", "I": y * 1.5, "wl": x, "meta": {"detector": "NIR PMT"}}])[0]
    assert ok["status"] == "ok" and flat["status"] == "fail" and adc["status"] == "fail" and pmt["status"] == "ok"


def test_settings_that_differ_are_listed():
    x = np.linspace(500, 600, 50)
    rows, differ = QA.settings_table([{"name": "a", "wl": x, "I": x, "meta": {"dwell_s": 0.1}}, {"name": "b", "wl": x, "I": x, "meta": {"dwell_s": 0.2}}])
    assert "dwell / integration (s)" in differ and len(rows) == 2


def test_runs_test_finds_a_trend():
    z, _ = QA.runs_test(np.r_[np.ones(8), -np.ones(8)])
    assert z < -2
    z2, _ = QA.runs_test(np.tile([1, -1], 8))
    assert z2 > 2


def test_self_check_is_unbiased_on_a_clean_series():
    sp = _series(seed=6)
    T = np.array([s["T"] for s in sp])
    f = PR.fit_global(sp, [(525.0, 9.4), (548.0, 9.4)], (505, 575), 1, "none", True, 60)
    f["bands"], f["extra"] = [(525.0, 9.4), (548.0, 9.4)], []
    D = np.array([r["ratio"] for r in f["rows"]])
    sD = np.array([r["ratio_err"] for r in f["rows"]])
    bf = PL.fit_boltzmann(T, D, sigma=sD)
    r = QA.self_check(f, T, bf["dE_eV"], bf["alpha"], True, [300.0] * len(T), n=4)
    assert r["status"] == "ok" and abs(r["bias_pct"]) < 5 and r["coverage"] >= 0.5


def test_analysis_record_is_json():
    import json
    rec = json.loads(QA.analysis_record("2.18.4", [{"file": "a.FS", "md5": "x"}], {"m": np.float64(1.5)}, {"T": np.array([1.0, np.nan])}, [QA.check("a", "ok", "b")]))
    assert rec["program"] == "Reticula 2.18.4" and rec["results"]["T"] == [1.0, None] and rec["quality_checks"][0]["status"] == "ok"


# ------------------------------------------------------------------------------------------------------------------------------ decays
def test_decay_fits_recover_known_lifetimes():
    rng = np.random.default_rng(2)
    t = np.arange(0, 2000, 2.0)
    tt = np.clip(t - 50, 0, None)
    y = rng.poisson((3000 * np.exp(-tt / 40) + 2000 * np.exp(-tt / 300)) * (t >= 50) + 15).astype(float)
    b, _ = DC.fit_auto(t, y)
    assert b["n_exp"] == 2 and abs(b["taus"][0] - 40) < 2 and abs(b["taus"][1] - 300) < 6 and b["flag"] == "good"
    assert abs(b["frac_amp"][0] - 0.6) < 0.04
    y1 = rng.poisson(5000 * np.exp(-tt / 150) * (t >= 50) + 20).astype(float)
    b1, _ = DC.fit_auto(t, y1)
    assert b1["n_exp"] == 1 and abs(b1["taus"][0] - 150) < 2
    yr = rng.poisson(4000 * (np.exp(-tt / 400) - np.exp(-tt / 30)) * (t >= 50) + 10).astype(float)
    br, _ = DC.fit_auto(t, yr, model="rise")
    assert abs(br["taus"][0] - 400) < 8 and abs(br["rise"] - 30) < 2


def test_reconvolution_with_an_irf():
    rng = np.random.default_rng(5)
    t = np.arange(0, 50, 0.02)
    irf = np.exp(-0.5 * ((t - 5) / 0.64) ** 2)
    true = np.convolve(irf / irf.sum(), np.exp(-t / 2.0))[: len(t)]
    y = rng.poisson(true / true.max() * 1e4 + 5).astype(float)
    b, _ = DC.fit_auto(t, y, model="reconv", irf=irf, max_exp=2)
    assert b["n_exp"] == 1 and abs(b["taus"][0] - 2.0) < 0.03


def test_read_decay_units_and_columns():
    txt = "Time (us),300K,320K\n" + "\n".join(f"{i * 2.0},{1000 * math.exp(-i * 2 / 100):.1f},{1000 * math.exp(-i * 2 / 80):.1f}" for i in range(200))
    c = DC.read_decay(txt.encode(), "d.csv")
    assert len(c) == 2 and c[0]["unit"] == "µs" and c[0]["name"] == "300K" and len(c[0]["t"]) == 200


def test_quench_fit_has_errors():
    T = np.arange(280.0, 520.0, 15.0)
    I = PL.quench_model(T, 1000, math.log(5e5), 0.4) * (1 + np.random.default_rng(3).normal(0, 0.005, len(T)))
    f = PL.fit_quench(T, I)
    assert f["status"] == "ok" and abs(f["Ea_meV"] - 400) < 3 * f["Ea_err_meV"] + 10 and np.isfinite(f["T50_err"])


# ------------------------------------------------------------------------------------------------------------------------------ colour, lines, files
def test_dominant_wavelength_of_narrow_lines():
    wl = np.arange(380.0, 781.0, 1.0)
    for lam in (470, 520, 610):
        c = PX.cie_xy(wl, np.exp(-0.5 * ((wl - lam) / 3) ** 2), False)
        d = PX.dominant_wavelength(c["x"], c["y"])
        assert abs(d["dominant_nm"] - lam) < 1.5 and d["purity_pct"] > 98
    c = PX.cie_xy(wl, np.exp(-0.5 * ((wl - 450) / 5) ** 2) + 1.5 * np.exp(-0.5 * ((wl - 650) / 5) ** 2), False)
    assert PX.dominant_wavelength(c["x"], c["y"])["complementary"]


def test_lamp_correction_recovers_the_response():
    wl = np.linspace(400, 900, 251)
    response = 0.4 + np.exp(-0.5 * ((wl - 650) / 120) ** 2)
    cert = 1e3 * (wl / 600) ** 4
    w, f = PX.correction_from_lamp(wl, cert * response, wl, cert, smooth_nm=0.1)
    corr = PX.apply_correction(wl, cert * response, w, f)
    assert np.allclose(corr / corr.mean(), cert / cert.mean(), rtol=0.02)
    assert abs(PX.correction_coverage(300, 900, w) - 5 / 6) < 0.01


def test_transition_library_and_coupled_pairs():
    assert ("Er³⁺", "²H₁₁/₂→⁴I₁₅/₂") in LN.find(525)
    assert LN.coupled_pair(525, 545)[0] == "Er³⁺" and LN.coupled_pair(545, 525)[0] == "Er³⁺" and LN.coupled_pair(980, 1020) is None


def test_fluoracle_and_ocean_headers_are_read():
    fl = ("Labels,S_20C,S_30C\nType,Emission Scan,Emission Scan\nFixed/Offset,980,980\nDwell Time,0.2,0.2\nRepeats,1,3\n"
          + "\n".join(f"{900 + i},{100 + i},{200 + i}" for i in range(20))).encode()
    sp = PL.read_pl(fl, "fl.csv")
    assert [s["name"] for s in sp] == ["20C", "30C"] and sp[0]["meta"]["dwell_s"] == 0.2 and sp[1]["meta"]["scans"] == 3 and sp[0]["meta"]["excitation_nm"] == 980
    oo = ("Integration Time (sec): 1.500000E-1\nScans to average: 4\n>>>>>Begin Spectral Data<<<<<\n" + "\n".join(f"{500 + i}\t{10 + i}" for i in range(20))).encode()
    s2 = PL.read_pl(oo, "oo.txt")[0]
    assert abs(s2["meta"]["dwell_s"] - 0.15) < 1e-9 and s2["meta"]["scans"] == 4


# ------------------------------------------------------------------------------------------------------------------------------ the page
def _page(monkeypatch, spectra_raw, name="multi.FS"):
    from streamlit.testing.v1 import AppTest
    warnings.filterwarnings("ignore")
    orig = PL.read_pl
    monkeypatch.setattr(PL, "read_pl", lambda src, nm="": orig(spectra_raw, name))
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600)
    at.session_state["_lv"] = {"global": "Advanced"}
    at.run()
    at.sidebar.radio[0].set_value("Photoluminescence").run()
    return at


def test_thermometry_page_shows_a_verdict_and_the_budget(monkeypatch):
    sp = _series(seed=8)
    head = "wavelength," + ",".join(f"{s['T']:.0f} K" for s in sp)
    body = "\n".join(f"{x:.3f}," + ",".join(f"{s['y'][i] + 300:.0f}" for s in sp) for i, x in enumerate(sp[0]["x"]))
    at = _page(monkeypatch, (head + "\n" + body).encode(), "series.csv")
    assert not at.exception, [e.value for e in at.exception][:2]
    tabs = [t.label for t in at.tabs]
    assert any("Thermometry" in t for t in tabs), tabs
    at.radio(key="pl_method").set_value("Integration windows (choose them on the graph)").run()
    assert not at.exception, [e.value for e in at.exception][:2]
    txt = " ".join(str(getattr(e, "value", "")) for e in list(at.success) + list(at.warning) + list(at.error))
    assert "Confidence" in txt, txt[:500]
    assert any(m.label == "ΔE" for m in at.metric)
    assert any("Uncertainty budget" in e.label for e in at.expander) and any("Quality report" in e.label for e in at.expander)


def test_decay_mode_runs_with_the_example():
    from streamlit.testing.v1 import AppTest
    warnings.filterwarnings("ignore")
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600)
    at.run()
    at.sidebar.radio[0].set_value("Photoluminescence").run()
    at.radio(key="pl_kind").set_value("Decay curves (lifetimes)").run()
    at.checkbox(key="dc_ex").check().run()
    assert not at.exception, [e.value for e in at.exception][:2]
    assert any("described well" in s.value for s in at.success) or len(at.warning)
    assert any("Eₐ" in m.label for m in at.metric)


@pytest.mark.parametrize("sep", [",", "\t", ";"])
def test_text_export_with_a_separator_at_the_end_of_every_line(sep):
    """v2.18.5: Fluoracle text exports end every line with the separator - this gave 'list index out of range'."""
    name = "Sample_20C"
    head = [f"Labels{sep}{name}{sep}", f"Type{sep}Emission Scan{sep}", f"Fixed/Offset{sep}980.00{sep}", f"Dwell Time{sep}0.10{sep}"]
    body = [f"{900 + i:.2f}{sep}{12000 + i}{sep}" for i in range(50)]
    sp = PL.read_pl(("\n".join(head + body)).encode("utf-8"), name + ".txt")
    assert len(sp) == 1 and sp[0]["name"] == name and len(sp[0]["wl"]) == 50
    assert PL.parse_temperature(sp[0]["T_text"])[0] == 20 and sp[0]["meta"]["excitation_nm"] == 980.0 and sp[0]["meta"]["dwell_s"] == 0.1


def test_rows_without_numbers_give_a_clear_message():
    with pytest.raises(ValueError, match="No spectrum"):
        PL.read_pl(b"a,b\nc,d\n", "x.txt")
