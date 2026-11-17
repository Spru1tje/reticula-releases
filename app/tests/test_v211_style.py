"""v2.11 global figure style: defaults identical to v2.10, journal presets, palette, persistence, applied on every page."""
import os, sys, json
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt


def test_default_global_style_equals_v210_rcparams():
    from mfs import style as S
    S.GLOBAL.clear(); S.GLOBAL.update(S.GLOBAL_DEFAULTS)
    S.apply_style(8.0)
    rc = plt.rcParams
    assert rc["lines.linewidth"] == 1.0 and rc["axes.linewidth"] == 0.75 and rc["xtick.direction"] == "in"
    assert rc["xtick.major.size"] == 3.0 and rc["legend.frameon"] is False and rc["font.size"] == 8.0
    assert rc["axes.prop_cycle"] == matplotlib.rcParamsDefault["axes.prop_cycle"]
    assert S.colours_for(["MOF-801", "Sulfo-MOF-801"]) == ["#000000", "#0072B2"]


def test_palette_and_lines_switch():
    from mfs import style as S
    try:
        S.set_global(palette="Tol bright", line_width=2.0, axes_width=1.5, tick_dir="out", minor=False)
        rc = plt.rcParams
        assert rc["lines.linewidth"] == 2.0 and rc["axes.linewidth"] == 1.5 and rc["xtick.major.size"] == 6.0
        assert rc["xtick.direction"] == "out" and rc["xtick.minor.visible"] is False
        assert S.colours_for(["MOF-801", "a", "b"]) == ["#4477AA", "#EE6677", "#228833"]
        S.SERIES_OVERRIDES["a"] = "#123456"
        assert S.colours_for(["MOF-801", "a", "b"]) == ["#4477AA", "#123456", "#EE6677"]
    finally:
        S.SERIES_OVERRIDES.pop("a", None); S.GLOBAL.clear(); S.GLOBAL.update(S.GLOBAL_DEFAULTS); S.apply_style()


def test_preset_applies_on_every_page_and_persists(tmp_path):
    from streamlit.testing.v1 import AppTest
    from mfs import style as S
    os.environ.update(MFS_NAV_GUARD="0", MFS_STORE=str(tmp_path / "s.json"), MFS_DB=str(tmp_path / "r.db"), MFS_EXPORT_DIR=str(tmp_path / "ex"))
    at = AppTest.from_file(os.path.join(ROOT, "app.py"), default_timeout=300); at.run()
    at.selectbox(key="gs_preset").set_value("Nature / Springer Nature").run()
    assert not at.exception, [e.value for e in at.exception][:1]
    assert at.number_input(key="gs_font").value == 7.0 and at.selectbox(key="gs_width").value.startswith("Nature single")
    at.selectbox(key="gs_palette").set_value("Tol muted").run()
    for pg in ("PXRD", "TGA", "FTIR / Raman"):
        at.sidebar.radio[0].set_value(pg).run()
        assert not at.exception, (pg, [e.value for e in at.exception][:1])
        assert S.GLOBAL["font_size"] == 7.0 and S.GLOBAL["palette"] == "Tol muted" and plt.rcParams["lines.linewidth"] == 0.8
    saved = json.load(open(tmp_path / "figure_style.json", encoding="utf-8"))
    assert saved["gs_font"] == 7.0 and saved["gs_palette"] == "Tol muted"
    b = AppTest.from_file(os.path.join(ROOT, "app.py"), default_timeout=300); b.run()      # next start: same style
    assert b.selectbox(key="gs_palette").value == "Tol muted" and b.number_input(key="gs_font").value == 7.0
    b.button(key="gs_reset").click().run()
    assert b.selectbox(key="gs_palette").value.startswith("Okabe–Ito (sample") and b.number_input(key="gs_font").value == 8.0
    S.GLOBAL.clear(); S.GLOBAL.update(S.GLOBAL_DEFAULTS); S.apply_style()
