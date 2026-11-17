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
    at.sidebar.radio[0].set_value("Zeta / DLS").run()
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
def test_proton_tags_are_letters_from_the_most_downfield_group():
    from mfs import nmrassign as A, nmrpredict as HP
    G = A.groups(HP.mol_from_any("CC(=O)Oc1ccccc1C(=O)O"))
    tags = {g["tag"] for g in G}
    assert tags == set("abcdef")
    by_shift = sorted(G, key=lambda g: -g["shift"])
    assert [g["tag"] for g in by_shift] == list("abcdef")
    assert A._letters(0) == "a" and A._letters(25) == "z" and A._letters(26) == "aa" and A._letters(27) == "ab"
    assert A.proton_label(G[0]).startswith("H") and not any(ch.isdigit() for ch in A.proton_label(G[0]))


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
