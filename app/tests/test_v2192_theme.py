"""v2.19.2 sidebar and dark theme: every page has a visible (light) icon in the dark sidebar; the dark theme CSS covers the main area."""
import os

from mfs import sidebar_nav as SN
from mfs import ui as UI


def test_every_page_has_a_light_icon_for_the_dark_sidebar():
    for p in UI.PAGES:
        path = SN.icon_path(p)
        assert os.path.isfile(path), p
        assert "nav_icons_light" in path, p
        assert "#DCE7F5" in open(path, encoding="utf-8").read(), p


def test_dark_theme_css_covers_main_surfaces():
    css = SN.DARK_CSS
    for token in ('[data-testid="stMain"]', "#0F172A", "#1E293B", "stFileUploaderDropzone", "hm_intent_"):
        assert token in css, token
