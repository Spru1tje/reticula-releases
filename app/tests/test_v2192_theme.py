"""v2.19.2 sidebar and dark theme: every page has a visible (light) icon in the dark sidebar; the dark theme CSS covers the main area."""
import os

from mfs import sidebar_nav as SN
from mfs import ui as UI


def test_every_page_has_a_light_icon_for_the_dark_sidebar():
    # v2.19.4: one base icon per page, recoloured per row in CSS (light stroke on the dark sidebar)
    for p in UI.PAGES:
        path = SN.icon_path(p)
        assert os.path.isfile(path), p
        assert "nav_icons" in path, p
        css = SN.row_css(1, p, False)
        assert css, p
        assert "data:image/svg+xml;base64," in css, p
        uri = SN.icon_uri(path, SN.ICON_ON_DARK)
        assert uri in css, p


def test_dark_theme_css_covers_main_surfaces():
    css = SN.DARK_CSS
    for token in ('[data-testid="stMain"]', "#0F172A", "#1E293B", "stFileUploaderDropzone", "hm_intent_"):
        assert token in css, token
