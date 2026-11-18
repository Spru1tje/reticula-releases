"""v2.17 Schemes & Icons (phase 2/3): content from Reticula's own data and embedded figures."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest
from mfs import schemes as SC
from mfs.reticon import scene as RS

LIB_OK = os.path.exists(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "mfs", "icons_runtime", "specs.json"))


def test_svg_figure_is_embedded_with_its_viewbox_and_body():
    svg = '<?xml version="1.0"?><!-- c --><svg width="5mm" viewBox="0 0 10 20" xmlns="http://www.w3.org/2000/svg"><rect width="10" height="20"/></svg>'
    it = SC.svgart_from_svg(svg, "f1", 5, 6, 100)
    assert it["type"] == "svgart" and it["viewBox"] == "0 0 10 20" and it["body"].strip() == '<rect width="10" height="20"/>'
    assert (it["x"], it["y"], it["size"]) == (5.0, 6.0, 100.0)
    with pytest.raises(ValueError):
        SC.svgart_from_svg("not a figure", "x", 0, 0, 10)


def test_embedded_figure_renders_as_nested_svg_in_the_export():
    scene = {"version": 1, "width": 300, "height": 200, "background": "#FFFFFF", "theme": "default", "items": []}
    SC.add_svgart(scene, '<svg viewBox="0 0 10 10"><circle cx="5" cy="5" r="4"/></svg>', "fig", x=20, y=20, size=120)
    svg = RS.render_scene(scene, SC.library())
    assert 'viewBox="0 0 10 10"' in svg and "<circle" in svg and 'width="120"' in svg


def test_xy_parser_skips_headers_comments_and_reads_commas_or_spaces():
    x, y = SC.parse_xy(b"# comment\nangle,intensity\n10, 100\n20, 150\n30 , 120\n")
    assert list(x) == [10, 20, 30] and list(y) == [100, 150, 120]
    with pytest.raises(ValueError):
        SC.parse_xy(b"1 2\n")


@pytest.mark.skipif(not LIB_OK, reason="icon library not built")
def test_xy_plot_becomes_an_icon_that_embeds_in_a_scene():
    from mfs.reticon import generators as G
    x = np.linspace(10, 60, 40); y = np.sin(x / 5) ** 2 * 100
    ic = G.icon_from_xy(x, y, kind="line", name="t_plot", label="t")
    svg = SC.icon_svg(ic, prefix="t_plot_pref")
    item = SC.svgart_from_svg(svg, "p", 0, 0, 200)
    assert item["viewBox"] and "<" in item["body"]


@pytest.mark.skipif(not LIB_OK, reason="icon library not built")
def test_methods_schematic_is_a_scene_with_the_sample_and_its_techniques():
    scene = RS.methods_schematic("MOF-801", ["pxrd", "tga"])
    names = " ".join(str(it.get("text", "")) + str(it.get("label", "")) for it in scene["items"])
    assert "MOF-801" in names or any(it.get("type") == "icon" for it in scene["items"])
    assert any(it["type"] == "arrow" for it in scene["items"])
    assert scene["width"] > 0 and scene["height"] > 0


def test_page_shows_the_content_tabs_when_the_editor_is_open():
    from streamlit.testing.v1 import AppTest
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script = ("import streamlit as st\nimport sys\nsys.path.insert(0, %r)\nfrom mfs.pages.schemes import page_schemes\n"
              "st.session_state.setdefault('schemes', {'Demo': {'version': 1, 'width': 300, 'height': 200, 'background': '#FFFFFF', "
              "'theme': 'default', 'items': []}})\npage_schemes()\n" % root)
    at = AppTest.from_string(script, default_timeout=120).run()
    assert not at.exception
    at.toggle(key="optin|schemes").set_value(True).run()
    assert not at.exception
    labels = [t.label for t in at.tabs] if at.tabs else []
    assert any("Methods" in l for l in labels)


def test_a_scheme_panel_composes_into_a_combined_figure_next_to_a_graph():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mfs import composer as C
    scene = {"version": 1, "width": 200, "height": 120, "background": "#FFFFFF", "theme": "default",
             "items": [{"id": "t", "type": "text", "x": 10, "y": 60, "text": "A_{2}", "size": 18}]}
    png, _ = SC.export(scene, "png", width_mm=60, dpi=150)
    panel = SC.figure_panel(png)
    graph = lambda ax: (ax.plot([0, 1], [0, 1]), None)[1]
    fig, axes = C.compose([panel, graph], ncols=2)
    assert len(axes) >= 2
    plt.close(fig)


def test_a_bend_is_drawn_through_its_waypoints_in_the_export():
    scene = {"version": 1, "width": 400, "height": 200, "background": "#FFFFFF", "theme": "default",
             "items": [{"id": "p", "type": "text", "x": 20, "y": 40, "text": "a", "size": 14},
                       {"id": "q", "type": "text", "x": 300, "y": 40, "text": "b", "size": 14},
                       {"id": "ar", "type": "arrow", "from": {"x": 20, "y": 100}, "to": {"x": 300, "y": 100},
                        "route": "straight", "waypoints": [[100, 160], [200, 160]]}]}
    svg = RS.render_scene(scene, SC.library())
    assert 'd="M' in svg and "100.0 160.0" in svg and "200.0 160.0" in svg
