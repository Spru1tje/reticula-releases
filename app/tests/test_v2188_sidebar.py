"""v2.18.8: categorised sidebar. Every page is reachable through exactly one category button, has an icon file, and the
current page is kept in st.session_state['nav'] across clicks (Streamlit AppTest runs the real sidebar code)."""
import os
import sys

from streamlit.testing.v1 import AppTest

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, APP)

from mfs import sidebar_nav as SN  # noqa: E402
from mfs import ui as UI  # noqa: E402


def test_categories_cover_every_page_once():
    flat = [p for _, pages in SN.CATEGORIES for p in pages]
    assert sorted(flat) == sorted(UI.PAGES)
    assert len(flat) == len(set(flat))


def test_every_page_has_an_icon_file():
    for p in UI.PAGES:
        assert os.path.isfile(SN.icon_path(p)), p


_SCRIPT = f"""
import sys
sys.path.insert(0, {APP!r})
import streamlit as st
from mfs import pagestate as PS
from mfs import sidebar_nav as SN
from mfs import ui as UI
PS.begin()
page = SN.render(UI.PAGES)
st.write("PAGE=" + page)
"""


def test_default_page_is_home_and_click_switches_page():
    at = AppTest.from_string(_SCRIPT, default_timeout=60).run()
    assert not at.exception
    assert any(m.value == "PAGE=Home" for m in at.markdown) or any("PAGE=Home" in m.value for m in at.markdown)
    at.button(key="nv_TGA").click().run()
    assert not at.exception
    assert any("PAGE=TGA" in m.value for m in at.markdown)
    at.button(key="nv_DSC").click().run()
    assert any("PAGE=DSC" in m.value for m in at.markdown)
