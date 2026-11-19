"""v2.15: methods text + sample summary."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mfs import summary as M, ui


def test_methods_and_summary():
    skip = ("Home", "Figure composer", "Batch processing", "SI export bundle", "Sample database", "Export guide", "Reticon", "Help and info")
    assert all(p in M.TECHNIQUE_METHODS for p in ui.PAGES if p not in skip), "a technique has no methods text"
    t = M.methods_text(["PXRD", "TGA", "PXRD", "Home"], "2.15.0")
    assert t.count("**PXRD.**") == 1 and "Reticula 2.15.0" in t and "Home" not in t
    store = {"a": {"name": "MOF_1", "sections": {"tga": {"values": {"residue": 41.25, "x": {"y": 2}}}}},
             "b": {"name": "MOF 2", "sections": {"n2": {"values": {"bet": 1050.0}}}}}
    df = M.sample_summary(store)
    assert list(df["Sample"]) == ["MOF_1", "MOF 2"] and "tga: x.y" in df and df.loc[1, "n2: bet"] == 1050.0
    tex = M.summary_latex(df)
    assert "MOF\\_1" in tex and "\\toprule" in tex and tex.count("\\\\") == 3


def test_levels_standard_default_and_panel(monkeypatch):
    monkeypatch.setenv("MFS_DEFAULT_LEVEL", "Standard")
    import warnings
    warnings.filterwarnings("ignore")
    os.environ["MFS_NO_UPDATE_CHECK"] = "1"
    from streamlit.testing.v1 import AppTest
    HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600)
    at.run()
    at.sidebar.radio[0].set_value("TGA").run()
    assert not at.exception
    tabs = [t.label for t in at.tabs]
    assert "🖋️ Title & axes" in tabs and "📏 Axes & ticks" not in tabs and "💾 Preset" not in tabs
    assert not any("Sample → colour rules" in e.label for e in at.sidebar.expander)
    at.sidebar.selectbox(key="lv_mod|TGA").set_value("Advanced").run()      # per-module override
    assert "📏 Axes & ticks" in [t.label for t in at.tabs]
    assert any("Sample → colour rules" in e.label for e in at.sidebar.expander)


import pytest


@pytest.mark.parametrize("level", ["Standard", "Advanced", "Developer"])
def test_every_page_renders_at_every_level(level, monkeypatch):
    """The gated blocks run with a 'ghost' st below their level: no page may break at any level."""
    import warnings
    warnings.filterwarnings("ignore")
    from streamlit.testing.v1 import AppTest
    HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600)
    at.session_state["_lv"] = {"global": level}
    at.run()
    bad = {}
    for pg in ui.PAGES:
        at.sidebar.radio[0].set_value(pg).run()
        if at.exception:
            bad[pg] = str(at.exception[0].value)[:200]
    assert not bad, bad


def _overlaps(fig):
    from mfs import composer as C
    g = C._group_boxes(fig)
    bad = 0
    for a in g:
        for b in g:
            if a < b and set(range(a[0], a[1])) & set(range(b[0], b[1])) and (b[2] == a[3] or a[2] == b[3]):
                l, r = (a, b) if a[2] < b[2] else (b, a)
                bad += g[l][2] > g[r][0] + 1e-6
    return bad


def test_composer_colourbar_does_not_overlap_neighbour_and_ratios():
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mfs import composer as C

    def heat(ax):
        im = ax.imshow(np.random.default_rng(0).random((20, 20)) * 5e5, aspect="auto")
        ax.figure.colorbar(im, ax=ax).set_label("Intensity (counts)")
        ax.set_ylabel("Temperature (K)")

    def line(ax):
        ax.plot([0, 1], [0.3, 0.42]); ax.set_ylabel(r"$\Delta$ = I$_1$/I$_2$"); ax.set_xlabel("T")
    panels = [heat, line, heat, line, line, heat]
    fig0, _ = C.compose(panels, ncols=3, wspace=0.1, hspace=0.1, auto_space=False)
    assert _overlaps(fig0) > 0                                         # reproduces the problem
    fig, _ = C.compose(panels, ncols=3, wspace=0.1, hspace=0.1)
    assert _overlaps(fig) == 0
    fig2, axs = C.compose([line, line, line], ncols=3, width_ratios=[1, 2, 1], spans=None, label_size=12)
    w = [a.get_position().width for a in axs]
    assert w[1] > 1.7 * w[0]
    fig3, axs = C.compose([line, line, line], ncols=2, spans=[(1, 1), (1, 1), (1, 2)])
    assert axs[2].get_position().width > 1.8 * axs[0].get_position().width
    plt.close("all")
