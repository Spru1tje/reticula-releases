"""v2.17 Reticon (phase 1): scenes, export and the page."""
import io, os, sys, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pytest
from mfs import schemes as SC, project as PRJ

SCENE = {"version": 1, "width": 400, "height": 300, "background": "#FFFFFF", "theme": "default",
         "items": [{"id": "f", "type": "icon", "icon": "round_bottom_flask", "shape": "base", "x": 40, "y": 40, "size": 120},
                   {"id": "t", "type": "text", "text": "A → B", "x": 200, "y": 150, "size": 18}]}


def test_templates_and_library_are_available():
    assert "toc_graphic" in SC.templates()
    assert SC.library() is not None


def test_scene_names_are_unique_and_scenes_travel_with_the_project():
    st_store = {}
    import streamlit as st
    st.session_state.clear()
    a = SC.add_scene("Fig 1", SCENE)
    b = SC.add_scene("Fig 1", SCENE)
    assert a == "Fig 1" and b == "Fig 1 (2)" and set(SC.scenes()) == {a, b}
    assert "schemes" in PRJ.GLOBAL_KEYS
    st.session_state.clear()


def test_export_formats_have_the_right_file_signature():
    svg, info = SC.export(SCENE, "svg")
    assert svg.lstrip().startswith(b"<?xml") or b"<svg" in svg[:400]
    pdf, _ = SC.export(SCENE, "pdf", width_mm=85)
    assert pdf.startswith(b"%PDF")
    png, info = SC.export(SCENE, "png", width_mm=85, dpi=300)
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    from PIL import Image
    im = Image.open(io.BytesIO(png))
    assert im.size[0] == round(85 / 25.4 * 300)
    assert info["vector_only"] in (True, False) and isinstance(info["warnings"], list)


def test_json_import_accepts_a_scheme_and_rejects_anything_else():
    assert SC.load_json_scene(json.dumps(SCENE).encode())["items"]
    with pytest.raises(ValueError):
        SC.load_json_scene(b'{"hello": 1}')


def test_page_renders_closed_and_open_without_errors():
    from streamlit.testing.v1 import AppTest
    script = ("import streamlit as st\nfrom mfs.pages.schemes import page_schemes\n"
              "import sys\nsys.path.insert(0, %r)\npage_schemes()\n" % os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    at = AppTest.from_string(script, default_timeout=120).run()
    assert not at.exception
    assert any("only when you open them" in i.value for i in at.info)
    at.toggle(key="optin|schemes").set_value(True).run()
    assert not at.exception
    assert any("Create a scheme" in c.value for c in at.caption)


def test_sub_and_superscript_markup_reaches_the_export():
    scene = {"version": 1, "width": 300, "height": 120, "background": "#FFFFFF", "theme": "default",
             "items": [{"id": "t", "type": "text", "x": 20, "y": 60, "text": "H_{2}O + CO^{2}", "size": 18}]}
    svg, _ = SC.export(scene, "svg")
    text = svg.decode("utf-8") if isinstance(svg, bytes) else svg
    assert 'baseline-shift="sub"' in text and 'baseline-shift="super"' in text
    assert "_{" not in text and "^{" not in text
