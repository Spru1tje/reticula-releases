"""v2.14: Zetasizer Excel import, per-figure font family, legend title size, DLS annotations, progress helper."""
import io
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "tests"))
from mfs import figctl as F, readers as R, style as S  # noqa: E402

CSV = os.path.join(HERE, "example_data", "Zeta-DLS", "Zetasizer_export_pD7.4_D2O.csv")


def _xlsx(preamble_rows=3, sheet="Results"):
    df = pd.read_csv(CSV)
    pre = pd.DataFrame([["Zetasizer export (test)"] + [None] * (df.shape[1] - 1)] * preamble_rows)
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as w:
        pre.to_excel(w, sheet_name=sheet, header=False, index=False, startrow=0)
        df.to_excel(w, sheet_name=sheet, index=False, startrow=preamble_rows)
    return buf.getvalue()


def test_zetasizer_excel_equals_csv():
    a = R.read_zetasizer(CSV, "x.csv")
    b = R.read_zetasizer(_xlsx(), "x.xlsx")
    pd.testing.assert_frame_equal(a["summary"].reset_index(drop=True), b["summary"].reset_index(drop=True), check_dtype=False)
    assert b["dist"] is not None and len(b["dist"]) == len(a["dist"])


def test_zetasizer_excel_without_results_table_explains():
    buf = io.BytesIO()
    pd.DataFrame({"a": [1, 2], "b": [3, 4]}).to_excel(buf, index=False)
    with pytest.raises(ValueError, match="No Zetasizer results table"):
        R.read_zetasizer(buf.getvalue(), "bad.xlsx")


def _fig_with_legend():
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], label="a"); ax.plot([0, 1], [1, 0], label="b")
    ax.legend(title="Samples")
    return fig, ax


def test_font_family_is_applied_to_everything():
    fam = next(f for f in S.available_fonts() if f not in ("auto", "DejaVu Sans"))
    fig, ax = _fig_with_legend(); ax.set_title("t"); ax.set_xlabel("x"); ax.text(0.5, 0.5, "note", transform=ax.transAxes)
    cfg = F.default_settings(); cfg["font"]["family"] = fam
    F.apply(fig, cfg)
    want = S._font_family(fam)
    got = {ax.xaxis.label.get_fontfamily()[0], ax.title.get_fontfamily()[0], ax.texts[0].get_fontfamily()[0],
           ax.get_legend().get_texts()[0].get_fontfamily()[0], ax.get_xticklabels()[0].get_fontfamily()[0]}
    assert got == {want}, got
    plt.close(fig)


def test_legend_title_size_and_family():
    fig, ax = _fig_with_legend()
    cfg = F.default_settings(); cfg["legend"]["title_size"] = 13.0; cfg["legend"]["family"] = "DejaVu Serif"
    F.apply(fig, cfg)
    lg = ax.get_legend()
    assert lg.get_title().get_fontsize() == 13.0
    assert lg.get_title().get_fontfamily()[0] == "DejaVu Serif" and lg.get_texts()[0].get_fontfamily()[0] == "DejaVu Serif"
    plt.close(fig)


def test_house_style_keeps_legend_title_size():
    cfg = F.default_settings(); cfg["legend"]["title_size"] = 11.0; cfg["legend"]["family"] = "Arial"
    h = F.house_style(cfg)
    assert h["legend"]["title_size"] == 11.0 and h["legend"]["family"] == "Arial"


def test_available_fonts_contains_the_classics():
    f = S.available_fonts()
    assert f[0] == "auto" and "Arial" in f and "Times New Roman" in f and len(f) == len(set(f))


@pytest.fixture
def at(tmp_path, monkeypatch):
    monkeypatch.setenv("MFS_STORE", str(tmp_path / "store.json"))
    from streamlit.testing.v1 import AppTest
    a = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=300); a.run()
    return a


def test_zeta_page_annotations_and_legend(at):
    at.button(key="nv_" + "Zeta / DLS").click().run()
    assert not at.exception, [e.value for e in at.exception][:2]
    for k in ("zeta_show_z", "zeta_show_pdi", "zeta_show_zp", "zeta_bar_leg"):
        at.checkbox(key=k).check()
    at.run()
    assert not at.exception, [e.value for e in at.exception][:2]


def test_progress_helper_learns_speed_and_returns_result():
    import numpy as np
    from mfs import progress as PG

    class _Box:
        def progress(self, *a, **k): pass
        def empty(self): pass
    p = PG.SegProgress.__new__(PG.SegProgress)
    p.n, p.label, p.box, p.done, p.rates, p._pending = 2, "x", _Box(), 0, {}, 0.0
    import time
    r = p.call(np.zeros((1000, 1000), np.uint8), "auto", lambda img, k: (time.sleep(0.6), k)[1], 7)
    assert r == 7 and p.done == 1 and "auto" in p.rates and p.rates["auto"] > 0.3
    # a fast (cache-like) call does not change the learned rate
    before = dict(p.rates)
    p.call(np.zeros((1000, 1000), np.uint8), "auto", lambda img: None)
    assert p.rates == before


# ---------------------------------------------------------------- volume / number distributions, quality notes, methods text
def _csv_with_vol_num():
    df = pd.read_csv(CSV)
    for k in range(1, 4):
        size_col = [c for c in df.columns if c.lower().startswith(f"sizes[{k}]")][0]
        int_col = [c for c in df.columns if c.lower().startswith(f"intensities[{k}]")][0]
        df[f"Volumes[{k}] (Percent)"] = df[int_col] * 0.5
        df[f"Numbers[{k}] (Percent)"] = df[int_col] * 0.25
    return df.to_csv(index=False)


def test_reader_volume_and_number_distributions():
    d = R.read_zetasizer(_csv_with_vol_num().encode("utf-8"), "x.csv")["dist"]
    assert {"intensity", "volume", "number"} <= set(d.columns)
    assert d.volume.notna().any() and d.number.notna().any()
    sub = d.dropna(subset=["volume"])
    assert ((sub.volume - 0.5 * sub.intensity).abs() < 1e-9).all()
    plain = R.read_zetasizer(CSV, "x.csv")["dist"]           # the plain export has no volume/number columns
    assert plain.volume.isna().all() and plain.number.isna().all()


def test_quality_notes():
    from mfs import methods as MT
    rows = [dict(sample="A", n=1, z_avg=100, pdi=0.8, zeta=-10, zeta_sd=1), dict(sample="B", n=3, z_avg=200, pdi=0.1, zeta=-45, zeta_sd=2),
            dict(sample="C", n=3, z_avg=300, pdi=0.5, zeta=40, zeta_sd=2)]
    notes = MT.dls_quality_notes(rows)
    text = " | ".join(t for _, t in notes)
    assert ("warn", ) == tuple({l for l, t in notes if t.startswith("A: PDI")})
    assert "A: one record only" in text and "A (-10 mV)" in text and "C: PDI 0.50" in text and "C (40" not in text
    assert "B:" not in text


def test_methods_text_uses_only_the_data():
    from mfs import methods as MT
    rows = [dict(sample="A", n=3, z_avg=120.4, pdi=0.12, zeta=-30.2, zeta_sd=1), dict(sample="B", n=3, z_avg=150.2, pdi=0.2, zeta=-22.0, zeta_sd=1)]
    t = MT.dls_methods(rows)
    assert "Z-average diameters were 120–150 nm" in t and "PDI 0.12–0.20" in t and "-30.2–-22.0 mV" in t
    assert "3 records per sample" in t and "[model]" in t and "[temperature]" in t
    assert MT.dls_methods([]) == ""


# ---------------------------------------------------------------- NMR: proton letters, figure layout
def test_proton_tags_are_literature_locants():
    """v2.17: a proton is named after the locant of the carbon it sits on (H-2, H-2/6), not by letters."""
    from mfs import nmrassign as A, nmrpredict as HP
    G = A.groups(HP.mol_from_any("CC(=O)Oc1ccccc1C(=O)O"))
    tags = [g["tag"] for g in G]
    assert all(t.startswith("-") and t[1:].replace("/", "").isdigit() for t in tags), tags
    assert len(set(tags)) == len(tags)
    assert A.proton_label(G[0]).startswith("H-") and not any(ch.isalpha() for ch in A.proton_label(G[0])[1:])


def _nmr_fig(**kw):
    import warnings
    warnings.filterwarnings("ignore")
    from mfs import nmrproc as NP, nmrpeaks as PK, nmrpredict as HP, nmrassign as A, nmrfig as NF
    fn = os.path.join(HERE, "example_data", "NMR_1H", "aspirin_300_XWIN35.zip")
    sp = NP.process(NP.load(open(fn, "rb").read(), "a.zip")[0], lb=None, p0=None, p1=None, reference="auto", use_topspin=True)
    ms = PK.analyse(sp, k=6.0, solvent="CDCl3"); PK.normalise(ms)
    mol = HP.mol_from_any("CC(=O)Oc1ccccc1C(=O)O"); G = A.groups(mol); A.assign(G, ms)
    fig, ax = plt.subplots(figsize=(3.25, 1.5)); NF.draw(ax, sp, ms, mol, G, {}, **kw)
    fig.canvas.draw()
    return fig, ax


def test_nmr_figure_texts_stay_inside_the_axes_and_off_the_axis():
    fig, ax = _nmr_fig()
    ab = ax.get_window_extent()
    for t in ax.texts:
        if not t.get_text() or t.get_transform() == ax.transAxes:
            continue
        bb = t.get_window_extent()
        assert bb.y1 <= ab.y1 + 1.0 and bb.y0 >= ab.y0 - 1.0, (t.get_text(), bb, ab)     # nothing clipped or through the x axis
    plt.close(fig)


def test_integral_numbers_can_be_hidden():
    f1, a1 = _nmr_fig(); f2, a2 = _nmr_fig(show_int_numbers=False)
    n1 = [t.get_text() for t in a1.texts if t.get_text().replace(".", "").isdigit()]
    n2 = [t.get_text() for t in a2.texts if t.get_text().replace(".", "").isdigit()]
    assert len(n1) >= 5 and n2 == []
    assert any(ln.get_color() == "#0072B2" for ln in a2.lines)            # the bars stay
    plt.close(f1); plt.close(f2)


# ---------------------------------------------------------------- NMR: satellites, integral uncertainty, assignment report
def test_13c_satellites_are_recognised_and_ordinary_peaks_are_not():
    from mfs import nmrpeaks as PK
    sf = 300.0
    main = PK.Multiplet(lo=2.95, hi=3.05, centre=3.0, integral=100.0, mult="s")
    sat1 = PK.Multiplet(lo=2.78, hi=2.80, centre=3.0 - 63 / sf, integral=0.55, mult="s")      # ±63 Hz = ¹J(CH 126 Hz)/2
    sat2 = PK.Multiplet(lo=3.20, hi=3.22, centre=3.0 + 63 / sf, integral=0.55, mult="s")
    other = PK.Multiplet(lo=1.0, hi=1.02, centre=1.01, integral=0.6, mult="s")               # an impurity, not symmetric
    ms = [main, sat1, other, sat2]
    assert PK.tag_satellites(ms, sf) == 2
    assert ms[1].note.startswith("13C satellite") and ms[3].note.startswith("13C satellite") and not ms[2].note and not ms[0].note
    lone = [PK.Multiplet(lo=2.95, hi=3.05, centre=3.0, integral=100.0, mult="s"),
            PK.Multiplet(lo=2.78, hi=2.80, centre=3.0 - 63 / sf, integral=0.55, mult="s")]       # one-sided: could be an impurity
    assert PK.tag_satellites(lone, sf) == 0


def test_integral_uncertainty_and_integer_flag():
    from mfs import nmrpeaks as PK
    from types import SimpleNamespace
    ppm = np.linspace(10, -1, 4000)
    rng = np.random.default_rng(1)
    y = rng.normal(0, 1.0, ppm.size)
    sp = SimpleNamespace(ppm=ppm, y=y)
    a = PK.Multiplet(lo=2.0, hi=2.5, centre=2.25, integral=300.0, nH=3.0, mult="s")
    b = PK.Multiplet(lo=7.0, hi=7.4, centre=7.2, integral=166.0, nH=1.66, mult="br s")          # 1.66 H: not an integer
    PK.add_integral_uncertainty([a, b], sp)
    assert a.nH_err is not None and 0 < a.nH_err < 0.2 and not a.flag
    assert "not an integer" in b.flag


def test_assignment_report_flags_doubtful_and_ambiguous():
    from mfs import nmrassign as A, nmrpeaks as PK, nmrpredict as HP
    sp, ms = _aspirin_for_report()
    mol = HP.mol_from_any("CC(=O)Oc1ccccc1C(=O)O"); G = A.groups(mol); asg = A.assign(G, ms)
    rep = A.assignment_report(G, ms, asg)
    assert set(rep) == {g["gid"] for g in G} and all(v[0] in ("ok", "check", "ambiguous", "doubtful") for v in rep.values())
    cooh = next(g for g in G if g["exch"])
    assert rep[cooh["gid"]][0] in ("doubtful", "check")                     # a 13 ppm prediction vs a signal at 8.2 ppm
    levels = [v[0] for v in rep.values()]
    assert "ok" in levels and ("ambiguous" in levels or "doubtful" in levels)


def _aspirin_for_report():
    import warnings
    warnings.filterwarnings("ignore")
    from mfs import nmrproc as NP, nmrpeaks as PK
    fn = os.path.join(HERE, "example_data", "NMR_1H", "aspirin_300_XWIN35.zip")
    sp = NP.process(NP.load(open(fn, "rb").read(), "a.zip")[0], lb=None, p0=None, p1=None, reference="auto", use_topspin=True)
    ms = PK.analyse(sp, k=6.0, solvent="CDCl3"); PK.normalise(ms); PK.add_integral_uncertainty(ms, sp)
    return sp, ms


# ---------------------------------------------------------------- particle sizer: parallel pool, edge refinement, shapes
def _synth_jobs(n=3):
    import particle_synth as PS
    jobs = []
    for seed in range(n):
        img, _ = PS.make("disc", n=8, size=40, hw=(400, 400), seed=seed)
        jobs.append(("adaptive", img, (True, 1.0, 4.0, 0.0, 0.25, 0.8, 2.0, True, 1536, None, 0.7, 10.0, "auto", 0.5)))
    return jobs


def test_segpool_parallel_equals_serial_and_caches():
    from mfs import segpool as SP, segjobs as SJ
    jobs = _synth_jobs(3)
    ser = [SJ.run(*j) for j in jobs]
    SP.clear_cache()
    seen = []
    last = dict(SP.prefetch(jobs, progress=lambda d, t, s: seen.append((d, t))))
    assert last["mode"] in ("parallel", "serial") and seen and seen[-1] == (3, 3)
    par = [SP.get(*j) for j in jobs]
    assert all(np.array_equal(a[0], b[0]) and a[1] == b[1] for a, b in zip(ser, par))
    assert all(SP.cached(*j) for j in jobs)
    assert SP.prefetch(jobs)["mode"] == "cache"                       # nothing left to do


def test_segpool_falls_back_to_serial_when_the_pool_breaks(monkeypatch):
    from mfs import segpool as SP
    SP.clear_cache(); SP.shutdown()

    def boom():
        raise RuntimeError("no processes allowed")
    monkeypatch.setattr(SP, "_get_pool", boom)
    jobs = _synth_jobs(2)
    last = dict(SP.prefetch(jobs))
    assert last["mode"] == "serial" and "no processes allowed" in last["error"]
    assert all(SP.cached(*j) for j in jobs)


def test_edge_refinement_removes_the_mask_bias_for_non_round_particles():
    import particle_synth as PS
    from mfs import particles as PA
    errs = {True: [], False: []}
    for sh in ("square", "ellipse"):
        for seed in range(2):
            img, truth = PS.make(sh, n=14, size=56, seed=seed)
            for er in (True, False):
                lab, f, t, info = PA.segment_adaptive(img, edge_refine=er)
                errs[er] += [(d / dt - 1) * 100 for _, dt, d in PS.evaluate(img, truth, lab, f, info)]
    refined, raw = np.array(errs[True]), np.array(errs[False])
    assert len(refined) >= 40
    assert abs(raw.mean()) > 6.0                                       # the raw mask is biased high
    assert abs(refined.mean()) < abs(raw.mean()) - 2.0                 # the calibrated edge is clearly closer


# ---------------------------------------------------------------- combined figure (several techniques in one figure)
def test_compose_letters_inside_caption_and_auto_columns():
    from mfs import composer as CP
    assert [CP.auto_columns(n) for n in (1, 2, 3, 4, 5, 6)] == [1, 2, 3, 2, 3, 3]
    assert [CP.panel_label(i, "A") for i in range(3)] == ["A", "B", "C"] and CP.panel_label(1, "(A)") == "(B)"
    fns = [lambda ax, k=k: ax.plot([0, 1], [0, k + 1]) for k in range(3)]
    fig, axs = CP.compose(fns, ncols=3, width=7.0, label_style="A", label_inside=True, caption="Fig. 2 (A) UV/Vis, (B) FTIR and (C) XRD of the same material.")
    assert len(axs) == 3
    letters = [t.get_text() for ax in axs for t in ax.texts]
    assert letters == ["A", "B", "C"]
    assert any("UV/Vis" in t.get_text() for t in fig.texts)                      # the caption is part of the figure
    fig.canvas.draw()
    cap = [t for t in fig.texts if "UV/Vis" in t.get_text()][0].get_window_extent()
    lowest_ax = min(ax.get_window_extent().y0 for ax in axs)
    assert cap.y1 <= lowest_ax + 1                                                # ... under the axes, not over them
    plt.close(fig)


def test_panel_settings_travel_with_the_panel():
    from mfs import composer as CP
    cfg = F.default_settings()
    cfg["panels"]["0"] = {"xlabel": "Wavelength (nm)", "xmin": 200.0, "xmax": 800.0}
    cfg["title"] = "UV/Vis"

    def draw(ax):
        ax.plot([100, 900], [0, 1]); ax.set_xlabel("original")
    fig, axs = CP.compose([lambda ax: (draw(ax), F.apply_panel(ax, cfg, None, 0))], ncols=1)
    ax = axs[0]
    assert ax.get_xlabel() == "Wavelength (nm)" and sorted(ax.get_xlim()) == [200.0, 800.0] and ax.get_title() == "UV/Vis"
    plt.close(fig)


def test_combined_figure_from_three_techniques_in_the_app(tmp_path, monkeypatch):
    monkeypatch.setenv("MFS_STORE", str(tmp_path / "store.json"))
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=400)
    at.run()
    for pg in ("PXRD", "FTIR / Raman", "UV/Vis"):
        at.button(key="nv_" + pg).click().run()
        assert not at.exception, [e.value for e in at.exception][:2]
        btn = [b for b in at.button if b.key and b.key.startswith("add_")][0]
        assert "Add to combined figure" in btn.label
        btn.click().run()
        assert not at.exception, [e.value for e in at.exception][:2]
    assert len(at.session_state["composer"]) == 3
    assert len(at.session_state["fig_registry"]) >= 3
    at.button(key="nv_" + "Figure composer").click().run()
    assert not at.exception, [e.value for e in at.exception][:2]
    fig = at.session_state["si_figs"]["combined_figure"]["fig"]
    assert len([a for a in fig.axes if a.axison]) == 3
    assert sorted(t.get_text() for a in fig.axes for t in a.texts if t.get_text() in ("A", "B", "C")) == ["A", "B", "C"]
    # the tray in the sidebar names the panels
    assert any("3 panels" in m.value for m in at.sidebar.markdown)


# ---------------------------------------------------------------- photoluminescence
def test_pl_temperature_parsing_and_units():
    from mfs import pl as PL
    cases = {"293.15 K": (293.15, "K"), "25C": (25.0, "C"), "50 °C": (50.0, "C"), "300K": (300.0, "K"), "T_300K": (300.0, "K"),
             "sample_25C_rep1": (25.0, "C"), "310": (310.0, None), "Er3+ sample": (None, None), "Wavelength (nm)": (None, None)}
    for t, want in cases.items():
        assert PL.parse_temperature(t) == want, (t, PL.parse_temperature(t))
    assert PL.to_kelvin([(20.0, None), (50.0, None)]) == [293.15, 323.15]                       # small plain numbers: degrees C
    assert PL.to_kelvin([(300.0, None), (400.0, None)]) == [300.0, 400.0]                        # large: K
    assert PL.to_kelvin([(25.0, "C"), (300.0, "K"), (None, None)]) == [298.15, 300.0, None]


def test_pl_reader_wide_table_two_column_and_excel(tmp_path):
    from mfs import pl as PL
    wide = "# comment\nWavelength (nm),293.15 K,303.15 K\n500,1,2\n501,2,3\n502,3,4\n503,4,5\n504,5,6\n505,6,7\n"
    sp = PL.read_pl(wide.encode(), "wide.csv")
    assert [s["name"] for s in sp] == ["293.15 K", "303.15 K"] and sp[0]["wl"][0] == 500 and sp[1]["I"][-1] == 7
    two = "\n".join(f"{500 + k}\t{k * 2.5}" for k in range(10)) + "\n"
    s2 = PL.read_pl(two.encode(), "sample_300K.txt")
    assert len(s2) == 1 and PL.parse_temperature(s2[0]["T_text"]) == (300.0, "K")
    import pandas as pd
    p = tmp_path / "pl.xlsx"
    pd.DataFrame({"Wavelength": np.arange(500, 510.0), "25 C": np.arange(10.0), "35 C": np.arange(10.0) * 2}).to_excel(p, index=False)
    s3 = PL.read_pl(p.read_bytes(), "pl.xlsx")
    assert [s["name"] for s in s3] == ["25 C", "35 C"] and s3[1]["I"][-1] == 18
    with pytest.raises(ValueError, match="No spectrum"):
        PL.read_pl(b"hello\nworld\n", "bad.txt")


def test_pl_fits_recover_known_parameters():
    from mfs import pl as PL
    T = np.linspace(80, 500, 20)
    I = 1000.0 / (1 + 2000.0 * np.exp(-0.25 / (PL.KB_EV * T)))
    q = PL.fit_quench(T, I)
    assert q["status"] == "ok" and abs(q["Ea_meV"] - 250) < 3 and abs(q["A"] / 2000 - 1) < 0.05
    assert abs(q["T50"] - 0.25 / (PL.KB_EV * np.log(2000))) < 3 and q["R2"] > 0.9999
    assert PL.fit_quench(T, np.full_like(T, 5.0))["status"].startswith("the intensity does not quench")
    Tt = np.linspace(293.15, 323.15, 7)
    D = 10.7 * np.exp(-(741.0 / PL.EV_TO_CM) / (PL.KB_EV * Tt))
    b = PL.fit_boltzmann(Tt, D)
    assert abs(b["dE_cm1"] - 741) < 1 and abs(b["alpha"] - 10.7) < 0.05
    sr = PL.relative_sensitivity(Tt, b["dE_eV"])
    num = 100 * np.abs(np.gradient(D, Tt) / D)                                    # numerical (1/D) dD/dT
    assert np.allclose(sr[1:-1], num[1:-1], rtol=0.01)
    assert np.allclose(PL.absolute_sensitivity(Tt, b), D * sr / 100, rtol=1e-6)
    st = PL.temperature_uncertainty(sr, 40000.0, 100000.0)
    assert np.allclose(st, np.sqrt(1 / 40000.0 + 1 / 100000.0) / (sr / 100))


def test_pl_peak_stats_axis_conversion_and_integrals():
    from mfs import pl as PL
    wl = np.linspace(400, 600, 801)
    y = 500.0 * np.exp(-0.5 * ((wl - 512.3) / 7.0) ** 2) + 3.0
    ps = PL.peak_stats(wl, y)
    assert abs(ps["pos"] - 512.3) < 0.1 and abs(ps["fwhm"] - 2.3548 * 7.0) < 0.2
    x, ye = PL.to_axis(wl, y, "eV")
    assert np.all(np.diff(x) > 0) and abs(np.trapezoid(ye, x) - np.trapezoid(y, wl)) / np.trapezoid(y, wl) < 5e-3     # Jacobian: area conserved
    assert PL.band_integral(wl, y, 400, 600, "linear") > 0.95 * 500 * 7 * np.sqrt(2 * np.pi)


def test_pl_page_runs_on_both_examples_and_gives_the_known_numbers(tmp_path, monkeypatch):
    monkeypatch.setenv("MFS_STORE", str(tmp_path / "store.json"))
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=400)
    at.run()
    at.button(key="nv_" + "Photoluminescence").click().run()
    assert not at.exception, [e.value for e in at.exception][:2]
    vals = {m.label: m.value for m in at.metric}
    assert abs(float(vals["ΔE"].split()[0]) - 741) < 25 and 1.1 < float(vals["Sᵣ at 293 K"].split()[0]) < 1.4
    figs = {k for k, r in at.session_state["si_figs"].items() if r["page"] == "Photoluminescence"}
    assert {"PL_spectra", "PL_temperature_dependence", "PL_thermometry"} <= figs
    at.radio(key="pl_example").set_value("Thermal quenching – phosphor, 80–500 K").run()
    assert not at.exception, [e.value for e in at.exception][:2]
    vals = {m.label: m.value for m in at.metric}
    assert abs(float(vals["Eₐ"].split()[0]) - 250) < 10 and abs(float(vals["T₅₀"].split()[0]) - 382) < 12


def test_rods_are_found_whole_and_squares_are_not_merged():
    import particle_synth as PS
    from mfs import particles as PA
    found = true = 0
    errs = []
    for seed in range(2):
        img, truth = PS.make("rod", n=14, size=56, seed=seed)
        lab, f, t, info = PA.segment_adaptive(img)
        assert info["_elongated"]["elongated"] and info["_pieces_merged"]               # the population is recognised as elongated
        ev = PS.evaluate(img, truth, lab, f, info)
        found += len(ev); true += len(truth); errs += [(d / dt - 1) * 100 for _, dt, d in ev]
    assert found / true >= 0.85, (found, true)                                          # was 44 % before the merge step
    assert abs(np.mean(errs)) < 13.0, np.mean(errs)
    img, truth = PS.make("square", n=14, size=56, seed=0)
    lab, f, t, info = PA.segment_adaptive(img)
    assert not info["_pieces_merged"]                                                   # squares / cubes keep their separate segments
    assert len(PS.evaluate(img, truth, lab, f, info)) == len(truth)


def test_zetasizer_dts_reads_names_and_intensity_distribution():
    """A real .dts (MFS_TEST_DTS=path) is read for the sample names and the intensity size distribution; xls files are not mistaken for .dts."""
    path = os.environ.get("MFS_TEST_DTS") or os.path.join(HERE, "tests", "data", "ZrO2-2-Aliquot1.dts")
    if not os.path.exists(path):
        pytest.skip("no .dts test file (private measurement data is not part of the repository)")
    z = R.read_zetasizer(open(path, "rb").read(), "x.dts")
    assert len(z["summary"]) >= 1 and z["summary"]["sample"].str.len().min() > 0
    sums = z["dist"].groupby("rec").intensity.sum()                       # every run (record) sums to 100 %
    assert np.allclose(sums, 100.0, atol=0.5)
    assert z["dist"].size_nm.min() == pytest.approx(0.4, rel=1e-3) and z["dist"].size_nm.max() == pytest.approx(10000.0, rel=1e-3)
    assert z["summary"]["zeta"].isna().all() and "ζ" in z["note"]
    za, pdi = z["summary"]["z_avg"], z["summary"]["pdi"]                    # confirmed against the Zetasizer software for ZrO2-aliquot2 (≈ 32 nm, PDI ≈ 0.52)
    if "ZrO2-2-Aliquot1" in path:
        assert za.between(30, 34).all() and pdi.between(0.48, 0.56).all(), (za.tolist(), pdi.tolist())
    else:
        assert za.isna().all() or (za > 0).all()                                 # other software versions: empty or positive, never invented
    # an .xls (also an OLE file) must not be taken for a .dts
    assert R._zeta_dts(b"\xd0\xcf\x11\xe0" + b"\x00" * 600) is None


def test_zeta_records_can_be_deselected_and_means_follow():
    z = R.read_zetasizer(CSV, "x.csv")
    assert "rec" in z["summary"] and "rec" in z["dist"] and z["summary"].rec.is_monotonic_increasing
    keep = [int(r) for r in z["summary"].rec][:-3]
    z2 = R.select_records(z, keep)
    assert len(z2["summary"]) == len(keep) and set(z2["dist"].rec) == set(keep)
    assert len(R.select_records(z, [])["summary"]) == 0
    # the mean of the remaining runs of a sample changes accordingly
    s = z["summary"].copy(); s["z_avg"] = [100.0, 110.0, 120.0] + list(s["z_avg"][3:])
    s["sample"] = ["A", "A", "A"] + list(s["sample"][3:])
    z3 = {"summary": s, "dist": z["dist"]}
    assert R.select_records(z3, s.rec)["summary"].query("sample == 'A'").z_avg.mean() == 110.0
    assert R.select_records(z3, [0, 2])["summary"].query("sample == 'A'").z_avg.mean() == 110.0 and \
        R.select_records(z3, [0, 1])["summary"].query("sample == 'A'").z_avg.mean() == 105.0


def test_dts_runs_are_split_into_sample_and_run_number():
    path = os.path.join(HERE, "tests", "data", "ZrO2-2-Aliquot1.dts")
    if not os.path.exists(path):
        pytest.skip("no .dts test file")
    z = R.read_zetasizer(open(path, "rb").read(), "x.dts")
    assert set(z["summary"]["sample"]) == {"ZrO2-aliquot2"} and z["summary"]["run"].tolist() == [1, 2, 3, 4, 5]
    assert z["summary"]["rec"].tolist() == [0, 1, 2, 3, 4]
    z2 = R.select_records(z, [0, 1, 4])
    assert z2["summary"]["run"].tolist() == [1, 2, 5] and set(z2["dist"]["rec"]) == {0, 1, 4}
