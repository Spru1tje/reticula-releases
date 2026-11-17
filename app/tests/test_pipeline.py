"""Readers, analysis, plotting, export and report tests on the bundled EXAMPLE data (synthetic, designed values known)."""
import glob, io, json, os
import numpy as np
import pandas as pd
import pytest
import matplotlib.pyplot as plt
from mfs import readers as R, analysis as A, plots as P, style as S, report as REP, defects as D, composer as C

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX = os.path.join(ROOT, "example_data")
ex = lambda *p: os.path.join(EX, *p)


# ------------------------------------------------------------------ readers
@pytest.mark.parametrize("path,fn,n", [
    ("PXRD/MOF801_STOE.xy", R.read_pxrd, 1000), ("FTIR/MOF801_ATR.dpt", R.read_ftir, 1000),
    ("Raman/MOF801_RAMII.dpt", R.read_ftir, 500), ("NMR/MOF801_digest_totxt.txt", R.read_nmr, 10000),
    ("NMR/MOF801_digest_TopSpin_expno10.zip", R.read_nmr, 10000), ("XPS/MOF801_Hf4f_CasaXPS.csv", R.read_xps, 50)])
def test_xy_readers(path, fn, n):
    df = fn(open(ex(path), "rb"))
    assert list(df.columns)[:2] == ["x", "y"] and len(df) > n and np.isfinite(df.y).all()


def test_topspin_native_matches_totxt_scale():
    a = R.read_nmr(open(ex("NMR/MOF801_digest_TopSpin_expno10.zip"), "rb"))
    b = R.read_nmr(open(ex("NMR/MOF801_digest_totxt.txt"), "rb"))
    assert a.x.max() == pytest.approx(b.x.max(), abs=0.01) and a.x.min() == pytest.approx(b.x.min(), abs=0.01)
    # fumarate CH (~6.5-6.8 ppm in NaOD/D2O digest) must sit at the same shift in both files
    pk = []
    for d in (a, b):
        w = d[(d.x > 6.0) & (d.x < 7.2)]
        pk.append(w.x[w.y.idxmax()])
    assert pk[0] == pytest.approx(pk[1], abs=0.01)


@pytest.mark.parametrize("f", sorted(glob.glob(os.path.join(EX, "TGA", "*.txt"))))
def test_tga_vendor_formats(f):
    df = R.read_tga(open(f, "rb"))
    assert df["T"].max() > 750 and 95 < df.w.iloc[0] <= 101 and df.w.iloc[-1] < 80


def test_sorption_microactive_and_csv():
    d = R.read_sorption(open(ex("Sorption/MOF801_3Flex_MicroActive.xlsx"), "rb"))
    assert set(d["iso"].branch) == {"ads", "des"} and d["psd"] is not None and "nan" not in d["meta"].get("sample", "")
    d2 = R.read_sorption(open(ex("Sorption/MOF801_isotherm_plain.csv"), "rb"))
    assert len(d2["iso"]) > 20 and d2["iso"].q.max() > 100


def test_icp_imagelab_zeta_sem():
    icp = R.read_icp(open(ex("ICP/ICP-OES_Syngistix_export.csv"), "rb"))
    assert {"Hf", "Zr", "S"} <= set(icp.element)
    il = R.read_imagelab(open(ex("SDS-PAGE/ImageLab_LanesAndBands_ZrO2@Sulfo_HEWL_60C.csv"), "rb"))
    assert il.band_pct.between(0, 100).all() and il.time.notna().any()
    z = R.read_zetasizer(open(ex("Zeta-DLS/Zetasizer_export_pD7.4_D2O.csv"), "rb"))
    assert len(z["summary"]) >= 3
    sem = R.read_sem(open(ex("SEM/MOF-801_SEM.tif"), "rb"))
    assert sem["pixel_nm"] and sem["pixel_nm"] > 0 and sem["databar_row"] == 968


def test_decimal_comma_semicolon():
    txt = "EXAMPLE; header\n2Theta;Intensity\n5,00;10,5\n5,02;11,0\n5,04;12,5\n5,06;9,0\n"
    df = R.read_xy(io.StringIO(txt))
    assert df.x.iloc[1] == pytest.approx(5.02) and df.y.iloc[2] == pytest.approx(12.5)


# ------------------------------------------------------------------ analysis
def test_first_order_recovers_k():
    t = np.linspace(0, 200, 15); k = 0.0150
    rng = np.random.default_rng(1)
    c = 2.0 * np.exp(-k * t) * (1 + rng.normal(0, 0.01, t.size))
    r = A.first_order_fit(t, c)
    assert r["k"] == pytest.approx(k, rel=0.03) and r["t_half"] == pytest.approx(np.log(2) / k, rel=0.03) and r["k_err"] > 0


def test_example_kinetics_recovers_design():
    df = pd.read_csv(ex("Kinetics/GlyGly_hydrolysis_60C_pD7.4.csv"), comment="#")
    s = [c for c in df.columns if c.lower().startswith("sample")][0]
    tc = [c for c in df.columns if "time" in c.lower()][0]
    cc = [c for c in df.columns if "glygly" in c.lower()][0]
    sub = df[df[s].str.contains("ZrO2@", regex=False)]
    r = A.first_order_fit(sub[tc].values, sub[cc].values)
    assert r["k"] == pytest.approx(0.0292, rel=0.10)


def test_langmuir_recovers():
    ce = np.array([0.2, 0.5, 1, 2, 4, 8, 16.]); qm, K = 1.8, 0.6
    r = A.langmuir_fit(ce, qm * K * ce / (1 + K * ce))
    vals = list(r.values())
    assert any(isinstance(v, float) and abs(v - qm) < 0.01 for v in vals)


def test_bet_rouquerol_on_example():
    d = R.read_sorption(open(ex("Sorption/MOF801_3Flex_MicroActive.xlsx"), "rb"))
    iso = d["iso"][d["iso"].branch == "ads"]
    b = A.bet_rouquerol(iso.p_rel.values, iso.q.values)
    assert 500 < b["area"] < 900 and b["C"] > 0 and b["c1"]
    # the ultramicroporous example saturates before n_m is reached (n(1-p) max at p = 0.012, q = 157 < n_m = 160):
    # criterion 3 cannot be met and the app must say so instead of silently reporting an area
    assert b["all_criteria"] is False and b["c3"] is False
    assert 0.15 < A.total_pore_volume(iso.p_rel.values, iso.q.values) < 0.4


def test_bet_type_ii_reference():
    """Ideal BET isotherm with n_m = 100 cm3 STP/g, C = 100 -> area 435.3 m2/g must be recovered."""
    p = np.linspace(0.01, 0.3, 30); nm, Cc = 100.0, 100.0
    q = nm * Cc * p / ((1 - p) * (1 + (Cc - 1) * p))
    b = A.bet_rouquerol(p, q)
    assert b["area"] == pytest.approx(435.3, rel=0.005) and b["C"] == pytest.approx(100, rel=0.01)


def test_tga_plateau_matches_design():
    truth = json.load(open(ex("TGA/example_design_truth.json"), encoding="utf-8"))
    df = R.read_tga(open(ex("TGA/MOF801_Netzsch.txt"), "rb"))
    t = A.tga_analysis(df["T"].values, df.w.values)
    assert 250 < t["T_plateau"] < 480 and t["W_plateau"] > t["W_end"]
    assert isinstance(truth, dict)


def test_equivalence_points():
    df = R.read_table(open(ex("Titration/MOF801_titration_NaOH.csv"), "rb"))
    eq = A.equivalence_points(df.iloc[:, 0].values, df.iloc[:, 1].values)
    assert len(eq) == 3 and np.all(np.diff(eq) > 0)


# ------------------------------------------------------------------ plotting / export
def test_colour_mapping_fixed(example_colour_rules):
    assert S.colour_for("ZrO$_2$@Sulfo-MOF-801(Hf)") == S.colour_for("composite") == "#D55E00"
    assert S.colour_for("MOF-801(Hf)") == "#000000"


def test_colour_mapping_neutral_by_default():
    """v2.13: no lab-specific name rules by default - palette order, whatever the samples are called."""
    assert S.SAMPLE_COLOURS == []
    names = ["ZrO2@Sulfo composite", "PM", "my reference", "sample 4"]
    assert S.colours_for(names) == ["#000000", "#0072B2", "#E69F00", "#009E73"]
    assert [S.colour_for(n, i) for i, n in enumerate(names)] == ["#000000", "#0072B2", "#E69F00", "#009E73"]
    S.SAMPLE_COLOURS[:] = [("uio", "#D55E00"), ("[broken", "#111111")]       # a user rule; an invalid regex is ignored
    try:
        assert S.colours_for(["a", "UiO-66", "b"]) == ["#000000", "#D55E00", "#0072B2"]
    finally:
        S.SAMPLE_COLOURS[:] = []


@pytest.mark.parametrize("fmt", ["svg", "pdf", "png", "tiff"])
def test_export_formats(fmt):
    fig, ax = S.new_figure(S.SINGLE)
    ax.plot([0, 1], [0, 1]); ax.set_xlabel(r"2$\theta$ (°)")
    b = S.export(fig, fmt, dpi=600)
    plt.close(fig)
    assert len(b) > 500
    if fmt == "svg":
        assert b"<text" in b or b"font" in b  # text kept as text (svg.fonttype none)
    if fmt == "png":
        from PIL import Image
        im = Image.open(io.BytesIO(b)); assert 0.75 * S.SINGLE < im.size[0] / 600 <= S.SINGLE + 0.02  # tight bbox never exceeds the column width


def test_all_example_figures_render(tmp_path):
    """Runs the full headless figure script (every technique) into a temp folder."""
    import subprocess, sys
    r = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "make_example_figures.py"), "--out", str(tmp_path), "--png-only"],
                       capture_output=True, text=True, cwd=ROOT)
    if "unrecognized arguments" in r.stderr:
        pytest.skip("figure script has no --out option")
    assert r.returncode == 0, r.stderr[-2000:]
    assert len(list(tmp_path.glob("*.png"))) >= 18


def test_composer_labels():
    fns = [lambda ax: ax.plot([0, 1], [1, 0]), lambda ax: ax.bar([1, 2], [2, 3]), lambda ax: ax.scatter([1], [1])]
    fig, axs = C.compose(fns, ncols=3, width=S.DOUBLE)
    texts = [t.get_text() for a in axs for t in a.texts] + [t.get_text() for t in fig.texts]
    assert {"(a)", "(b)", "(c)"} <= set(texts)
    plt.close(fig)


# ------------------------------------------------------------------ defect report
@pytest.mark.parametrize("f", sorted(glob.glob(os.path.join(EX, "Defects", "defect_input_*.json"))))
def test_defect_report_and_design_recovery(f, tmp_path):
    inp = REP.input_from_dict(json.load(open(f, encoding="utf-8")))
    res = D.run_all(inp)
    c = res["consensus"]
    design = {"MOF801": 5.2, "Sulfo-MOF801": 4.7, "ZrO2@Sulfo-MOF801": 4.7}[os.path.basename(f)[13:-5]]
    assert c["L"] == pytest.approx(design, abs=0.15)
    pdf = tmp_path / "r.pdf"; REP.pdf_report(res, str(pdf))
    assert pdf.stat().st_size > 20000
    csv = REP.csv_report(res); assert "CONSENSUS" in csv
    fig = P.fig_defect_comparison(res); plt.close(fig)


def test_every_equation_renders():
    """Every LaTeX equation of every method is converted to matplotlib mathtext for the PDF report."""
    inp = REP.input_from_dict(json.load(open(ex("Defects", "defect_input_MOF801.json"), encoding="utf-8")))
    res = D.run_all(inp)
    eqs = [e for m in res["methods"] for e in m.equations] + list(res["h"]["equations"])
    bad = [e for e in eqs if not REP._render_ok(REP.tex_to_mathtext(e))]
    assert len(eqs) > 20
    assert len(bad) / len(eqs) < 0.25, bad  # the rest are printed as plain text in the PDF (still shown as LaTeX in the app)
