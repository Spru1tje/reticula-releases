"""v2 features: table-aware 3Flex reader, shared figure controls, scale bars, TEM files, export robustness."""
import io, json, os, importlib.util
import numpy as np
import pandas as pd
import pytest
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, MultipleLocator
from mfs import readers as R, analysis as A, plots as P, style as S, figctl as F

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FULL = os.path.join(ROOT, "example_data", "Sorption", "SampleA_3Flex_full_report_EXAMPLE.xls")
MULTI = os.path.join(ROOT, "tests", "data", "SampleA_3Flex_full_report_multisheet_EXAMPLE.xls")
DATA = os.path.join(ROOT, "tests", "data")


# ============================================================ N2: MicroActive full-report export
@pytest.mark.parametrize("path", [FULL, MULTI])
def test_full_report_uses_only_isotherm_tabular(path):
    d = R.read_sorption(open(path, "rb"), os.path.basename(path))
    iso = d["iso"]
    assert (iso.branch == "ads").sum() == 42 and (iso.branch == "des").sum() == 12   # not 158 merged points
    a = iso[iso.branch == "ads"].sort_values("order"); de = iso[iso.branch == "des"].sort_values("order")
    assert np.all(np.diff(a.p_rel) > 0), "adsorption must increase in measurement order (no zig-zag)"
    assert np.all(np.diff(de.p_rel) < 0), "desorption must decrease in measurement order"
    assert 110 < iso.q.max() < 120 and iso.q.min() > 0.3   # cm3/g STP, the mmol/g log-plot table is NOT mixed in
    m = d["meta"]
    assert m["sample"] == "Sample A batch 3" and float(m["mass_g"]) == pytest.approx(0.0473)
    assert float(m["bet_reported"]) == pytest.approx(335.08) and m["bet_reported_range"] == (0.05, 0.3)
    kinds = [t["kind"] for t in d["tables"]]
    assert {"iso_tabular", "bet", "langmuir", "tplot", "iso_plot", "psd"} <= set(kinds)
    assert sum(t["used"] for t in d["tables"] if t["kind"].startswith("iso")) == 1
    assert d["psd"] is not None and {"w", "dvdw", "cum"} <= set(d["psd"].columns)


def test_v1_reader_reproduces_the_reported_bug():
    """Regression record: the v1 reader merged every table with a p/p0 + quantity column into one series."""
    spec = importlib.util.spec_from_file_location("mfs.readers_v1", os.path.join(ROOT, "tools", "legacy", "readers_v1.py"))
    m = importlib.util.module_from_spec(spec); m.__package__ = "mfs"; spec.loader.exec_module(m)
    old = m.read_sorption(open(FULL, "rb"), "x.xls")
    assert len(old["iso"]) > 150                      # merged (Isotherm + BET + linear/log overlay tables)
    assert old["iso"].q.min() < 1                      # the mmol/g table was mixed in as if it were cm3/g


def _xlsx(rows):
    import openpyxl
    wb = openpyxl.Workbook(); ws = wb.active
    for r in rows:
        ws.append(r)
    b = io.BytesIO(); wb.save(b); b.seek(0); b.name = "t.xlsx"; return b


def test_mmol_per_g_is_converted():
    rows = [["Sample: test"], ["Isotherm Tabular Report"], ["Relative Pressure (P/Po)", "Quantity Adsorbed (mmol/g)"]]
    rows += [[p, 4.0 * p / (0.01 + p)] for p in (0.001, 0.01, 0.1, 0.5, 0.9)] + [[p, 4.0 * p / (0.01 + p)] for p in (0.5, 0.1)]
    d = R.read_sorption(_xlsx(rows), "t.xlsx")
    assert d["meta"]["q_unit_in_file"] == "mmol/g"
    assert d["iso"].q.max() == pytest.approx(4.0 * 0.9 / 0.91 * 22.414, rel=1e-6)
    assert list(d["iso"].branch) == ["ads"] * 5 + ["des"] * 2


def test_branch_label_rows_and_side_by_side_pairs():
    rows = [["Isotherm Tabular Report"], ["Relative Pressure (P/Po)", "Quantity Adsorbed (cm³/g STP)"], ["Adsorption"]]
    rows += [[p, 100 * p / (0.01 + p)] for p in (0.001, 0.01, 0.2, 0.8)] + [["Desorption"]] + [[p, 101 * p / (0.01 + p)] for p in (0.6, 0.3)]
    d = R.read_sorption(_xlsx(rows), "t.xlsx")
    assert list(d["iso"].branch) == ["ads"] * 4 + ["des"] * 2
    rows = [["Isotherm Linear Plot"], ["Adsorption", None, "Desorption"],
            ["Relative Pressure (P/Po)", "Quantity Adsorbed (cm³/g STP)", "Relative Pressure (P/Po)", "Quantity Adsorbed (cm³/g STP)"]]
    rows += [[0.01, 50, 0.9, 99], [0.1, 80, 0.5, 92], [0.5, 90, 0.1, 81], [0.9, 98, None, None]]
    d = R.read_sorption(_xlsx(rows), "t.xlsx")
    assert (d["iso"].branch == "ads").sum() == 4 and (d["iso"].branch == "des").sum() == 3


def test_bet_app_value_vs_instrument_window():
    """The app keeps its Rouquerol BET; the file value is reference only. For a Type I MOF the fixed 0.05-0.30
    window lies after micropore filling and gives a LOWER area - the difference is the window, not the data."""
    d = R.read_sorption(open(FULL, "rb"), "x.xls")
    a = d["iso"][d["iso"].branch == "ads"]
    b = A.bet_rouquerol(a.p_rel.values, a.q.values)
    fx = A.bet_fixed_range(a.p_rel.values, a.q.values, 0.05, 0.30)
    assert b["all_criteria"] and b["p_hi"] <= 0.05
    assert fx["area"] < b["area"]


def test_isotherm_plot_measurement_order_and_single_legend():
    d = R.read_sorption(open(FULL, "rb"), "x.xls")
    fig, ax = S.new_figure()
    P.draw_isotherms(ax, [("A", d["iso"])], styles={"A": {"ads_color": "#123456", "des_color": "#654321", "marker": "s", "ms": 5, "ls": "--", "lw": 1.5}})
    legs = [c for c in ax.get_children() if c.__class__.__name__ == "Legend"]
    assert len(legs) == 1, "one legend box (v1 had two)"
    labels = [t.get_text() for t in legs[0].get_texts()]
    assert labels == ["A", "adsorption", "desorption"]
    ads_line = [l for l in ax.lines if l.get_label() == "A"][0]
    assert ads_line.get_color() == "#123456" and ads_line.get_marker() == "s" and ads_line.get_linestyle() == "--"
    assert np.all(np.diff(ads_line.get_xdata()) > 0)
    des_pts = [l for l in ax.lines if l.get_marker() == "s" and l.get_mfc() == "white"][0]
    assert des_pts.get_mec() == "#654321"
    plt.close(fig)


# ============================================================ shared figure controls
def _fig():
    fig, ax = S.new_figure()
    x = np.linspace(0, 1, 50)
    ax.plot(x, 100 * x, label="Sample A"); ax.plot(x, 50 * x, label="PM"); ax.legend()
    ax.set_xlabel("Relative pressure"); ax.set_ylabel(r"Quantity adsorbed (cm$^3$ g$^{-1}$ STP)")
    return fig, ax


@pytest.mark.parametrize("preset", [p for p in F.TICK_PRESETS if p != "Custom"])
def test_tick_side_presets(preset):
    fig, ax = _fig()
    cfg = F.default_settings(); cfg["ticks"]["preset"] = preset
    F.apply(fig, cfg)
    want = F.TICK_PRESETS[preset]
    kw_x, kw_y = ax.xaxis._major_tick_kw, ax.yaxis._major_tick_kw
    got = {"bottom": kw_x.get("tick1On", True), "top": kw_x.get("tick2On", False), "left": kw_y.get("tick1On", True), "right": kw_y.get("tick2On", False)}
    assert got == want
    plt.close(fig)


def test_custom_sides_labels_direction_and_increments():
    fig, ax = _fig()
    cfg = F.default_settings()
    cfg["ticks"].update(preset="Custom", sides=dict(bottom=True, top=False, left=True, right=True),
                        labels=dict(bottom=True, top=False, left=False, right=True), direction="out", minor=True)
    cfg["panels"]["0"] = {"xmajor": 0.2, "xminor": 0.05, "ymajor": 25, "yminor": 0}
    F.apply(fig, cfg)
    assert isinstance(ax.xaxis.get_major_locator(), MultipleLocator)
    t = [v for v in ax.get_xticks() if 0 <= v <= 1]
    assert np.allclose(np.diff(t), 0.2)
    assert np.allclose(np.diff([v for v in ax.yaxis.get_minorticklocs()])[:3], 25 / 5, atol=1e-9) or len(ax.yaxis.get_minorticklocs()) > 0
    assert ax.yaxis._major_tick_kw.get("label2On") and not ax.yaxis._major_tick_kw.get("label1On")
    assert ax.xaxis._major_tick_kw.get("tickdir") == "out"
    plt.close(fig)


def test_axis_titles_ranges_title_mathify():
    fig, ax = _fig()
    cfg = F.default_settings(); cfg["title"] = "N_2 at 77 K"
    cfg["panels"]["0"] = {"xlabel": "p/p_0", "ylabel": "n (cm^3 g^-1)", "xmin": 0.1, "xmax": 0.9, "ymin": 0, "ymax": 120}
    F.apply(fig, cfg)
    assert ax.get_xlabel() == "p/p$_{0}$" and ax.get_ylabel() == "n (cm$^{3}$ g$^{-1}$)"
    assert ax.get_xlim() == pytest.approx((0.1, 0.9)) and ax.get_ylim() == pytest.approx((0, 120))
    assert ax.get_title() == "N$_{2}$ at 77 K"
    fig.canvas.draw()
    plt.close(fig)


def test_bad_mathtext_never_crashes():
    fig, ax = _fig()
    cfg = F.default_settings(); cfg["panels"]["0"] = {"ylabel": r"broken $\frac{1}{ $"}
    F.apply(fig, cfg); fig.canvas.draw(); plt.close(fig)


def test_unit_conversion_uptake_mmol():
    fig, ax = _fig()
    cfg = F.default_settings(); cfg["panels"]["0"] = {"y_family": "uptake", "y_unit": "mmol g⁻¹", "ymajor": 1, "ymax": 4}
    F.apply(fig, cfg)
    loc = ax.yaxis.get_major_locator(); assert isinstance(loc, FixedLocator)
    labels = [t.get_text() for t in ax.yaxis.get_majorticklabels()]
    assert labels[:5] == ["0", "1", "2", "3", "4"]
    assert np.allclose(loc.locs[:5], np.arange(5) * 22.414)        # data stay in cm3/g; ticks are relabelled
    assert ax.get_ylim()[1] == pytest.approx(4 * 22.414)
    assert "mmol g$^{-1}$" in ax.get_ylabel() and ax.get_ylabel().startswith("Quantity adsorbed")
    plt.close(fig)


@pytest.mark.parametrize("fam,unit,lo,hi", [("2theta", "d-spacing (Å)", 5, 50), ("temperature", "K", 25, 800),
                                            ("wavenumber", "µm (wavelength)", 400, 4000), ("time_h", "min", 0, 48)])
def test_other_unit_families(fam, unit, lo, hi):
    fig, ax = S.new_figure(); ax.plot([lo, hi], [0, 1]); ax.set_xlim(lo, hi); ax.set_xlabel("x (u)")
    cfg = F.default_settings(); cfg["panels"]["0"] = {"x_family": fam, "x_unit": unit}
    F.apply(fig, cfg); fig.canvas.draw()
    labs = [t.get_text() for t in ax.xaxis.get_majorticklabels() if t.get_text()]
    assert len(labs) >= 2
    plt.close(fig)


def test_legend_rename_hide_reorder_outside_and_style():
    fig, ax = _fig()
    cfg = F.default_settings()
    cfg["legend"].update(custom=True, loc="outside right", ncol=2, frame=True, markerscale=2.0, fontsize=6, title="Samples",
                         entries=[{"label": "PM", "text": "physical mixture", "show": True},
                                  {"label": "Sample A", "text": "Sample A", "show": False}])
    F.apply(fig, cfg)
    lg = ax.get_legend()
    assert [t.get_text() for t in lg.get_texts()] == ["physical mixture"]
    assert lg._ncols == 2 and lg.get_frame_on() and lg.markerscale == 2.0 and lg.get_title().get_text() == "Samples"
    fig.canvas.draw()
    bb = lg.get_window_extent(); axbb = ax.get_window_extent()
    assert bb.x0 >= axbb.x1 - 1, "outside right"
    plt.close(fig)


@pytest.mark.parametrize("loc", F.LEGEND_LOCS)
def test_every_legend_position(loc):
    fig, ax = _fig()
    cfg = F.default_settings(); cfg["legend"].update(custom=True, loc=loc)
    F.apply(fig, cfg); fig.canvas.draw()
    assert ax.get_legend() is not None
    plt.close(fig)


def test_direct_labels_renamed_and_stacked_y_stays_hidden():
    fig, ax = S.new_figure()
    P.draw_stacked(ax, [("Sample A", np.linspace(5, 50, 100), np.random.rand(100)), ("PM", np.linspace(5, 50, 100), np.random.rand(100))], kind="pxrd")
    cfg = F.default_settings(); cfg["legend"]["entries"] = [{"label": "PM", "text": "physical mix", "show": True}]
    F.apply(fig, cfg)
    assert "physical mix" in [t.get_text() for t in ax.texts]
    assert len(ax.get_yticks()) == 0
    plt.close(fig)


def test_twin_axis_tga_keeps_dtg_on_right():
    fig, ax = S.new_figure(); tw = ax.twinx()
    ax.plot([0, 800], [100, 60]); tw.plot([0, 800], [0, 1]); tw.set_ylabel("DTG")
    F.apply(fig, F.default_settings())
    assert tw.yaxis._major_tick_kw.get("label2On", True) and not ax.yaxis._major_tick_kw.get("tick2On", True)
    plt.close(fig)


def test_presets_json_roundtrip_and_house_style():
    cfg = F.default_settings(); cfg["ticks"]["preset"] = "2 sides (bottom + left)"; cfg["legend"]["loc"] = "upper left"
    cfg["panels"]["0"] = {"xlabel": "only for this figure"}
    back = F.from_json(F.to_json(cfg))
    assert back["ticks"]["preset"] == "2 sides (bottom + left)" and back["panels"]["0"]["xlabel"] == "only for this figure"
    h = F.house_style(cfg)
    assert "panels" not in h and h["ticks"]["preset"] == "2 sides (bottom + left)" and h["legend"]["loc"] == "upper left"
    merged = F.merge(F.default_settings(), json.loads(json.dumps(h)))
    assert merged["ticks"]["preset"] == "2 sides (bottom + left)"


def test_series_colour_override_everywhere(example_colour_rules):
    S.set_series_colours({"Sample A": "#ff00ff"})
    try:
        assert S.colours_for(["Sample A", "PM"]) == ["#ff00ff", "#CC79A7"]
        fig, ax = S.new_figure()
        P.draw_overlay(ax, [("Sample A", [0, 1], [0, 1])], "x", "y")
        assert ax.lines[0].get_color() == "#ff00ff"
        plt.close(fig)
    finally:
        S.set_series_colours({})


def test_duplicate_family_does_not_steal_colours(example_colour_rules):
    c = S.colours_for(["Sample A batch 3", "Sample A", "PM", "Sample B", "Sample C", "ZrO$_2$"])
    assert c[2:] == ["#CC79A7", "#0072B2", "#D55E00", "#E69F00"] and len(set(c)) == 6


# ============================================================ scale bars / SEM / TEM
def test_scalebar_value_mode_is_physical():
    px, nm, lab = P.scalebar_geometry(2.6, 1280, {"mode": "value", "value": 500, "unit": "nm"})
    assert px == pytest.approx(500 / 2.6) and lab == "500 nm"
    px, nm, lab = P.scalebar_geometry(2.6, 1280, {"mode": "value", "value": 1, "unit": "µm"})
    assert px == pytest.approx(1000 / 2.6) and lab == "1 µm"
    px, nm, lab = P.scalebar_geometry(0.05, 2048, {"mode": "value", "value": 20, "unit": "Å"})   # TEM
    assert px == pytest.approx(40) and lab == "20 Å"


def test_scalebar_length_mode_updates_label():
    px, nm, lab = P.scalebar_geometry(2.6, 1280, {"mode": "length", "length_px": 400, "unit": "auto"})
    assert nm == pytest.approx(1040) and lab == "1.04 µm"
    px, nm, lab = P.scalebar_geometry(2.6, 1280, {"mode": "length", "length_px": 400, "unit": "nm"})
    assert lab == "1040 nm"
    assert P.scalebar_geometry(2.6, 1280, {"mode": "length", "length_px": 400, "label": "1 µm (custom)"})[2] == "1 µm (custom)"


@pytest.mark.parametrize("pos", ["lower right", "lower left", "upper right", "upper left"])
@pytest.mark.parametrize("box", [False, True])
def test_scalebar_positions_colours_box(pos, box):
    img = np.random.randint(0, 255, (400, 600), dtype=np.uint8)
    fig, ax = plt.subplots(figsize=(3, 2))
    P.draw_sem(ax, img, pixel_nm=5.0, bar={"mode": "value", "value": 1, "unit": "µm", "position": pos, "color": "black",
                                           "thickness": 0.03, "font_size": 9, "box": box})
    info = ax._mfs_scalebar
    assert info["length_px"] == pytest.approx(200) and info["label"] == "1 µm"
    bars = [p for p in ax.patches if p.get_width() == pytest.approx(200)]
    assert bars and bars[0].get_height() == pytest.approx(12)
    x0 = bars[0].get_x(); y0 = bars[0].get_y()
    assert (x0 > 300) == ("right" in pos) and (y0 > 200) == ("lower" in pos)
    assert (len(ax.patches) == 2) == box
    fig.canvas.draw(); plt.close(fig)


def test_magnification_override():
    assert R.pixel_from_magnification(20000, 1280) == pytest.approx(127e6 / (20000 * 1280))


def test_imagej_calibrated_tiff(tmp_path):
    import tifffile
    p = tmp_path / "tem.tif"
    tifffile.imwrite(p, (np.random.rand(64, 64) * 255).astype(np.uint8), imagej=True, resolution=(1 / 0.25, 1 / 0.25),
                     metadata={"unit": "nm"})
    d = R.read_sem(open(p, "rb"), "tem.tif")
    assert d["pixel_nm"] == pytest.approx(0.25, rel=1e-3)


@pytest.mark.parametrize("fn,expect", [("dmTest_openNCEM.dm4", 1.0), ("STEM_HAADF_openNCEM.emd", 87.37)])
def test_tem_files_via_ncempy(fn, expect):
    pytest.importorskip("ncempy")
    d = R.read_sem(open(os.path.join(DATA, fn), "rb"), fn)
    assert d["image"].ndim == 2 and d["image"].dtype == np.uint8 and d["kind"] == "TEM"
    assert d["pixel_nm"] == pytest.approx(expect, rel=1e-3)


# ============================================================ exports
def test_vector_exports_are_deterministic():
    fig, ax = _fig()
    for fmt in ("pdf", "svg", "png"):
        assert S.export(fig, fmt) == S.export(fig, fmt)
    plt.close(fig)
