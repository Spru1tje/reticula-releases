"""v2.19.2 (sidebar icon rows): the icon is recoloured per row and every category page has a readable base icon.

Known answers: a stroke of #1F2D3D becomes the given colour in the data URI; the selected row uses the navy colour."""
import base64
import os
import re

from mfs import sidebar_nav as SN


def test_icon_uri_recolours_the_stroke(tmp_path):
    f = tmp_path / "x.svg"
    f.write_text('<svg><path stroke="#1F2D3D" fill="none"/></svg>', encoding="utf-8")
    uri = SN.icon_uri(str(f), "#C9DCEB")
    svg = base64.b64decode(uri.split(",", 1)[1]).decode("utf-8")
    assert 'stroke="#C9DCEB"' in svg
    assert "#1F2D3D" not in svg


def test_row_css_uses_navy_on_selected_and_light_otherwise():
    page = "UV/Vis"
    def stroke(css):
        b64 = re.search(r"base64,([A-Za-z0-9+/=]+)", css).group(1)
        return re.search(r'stroke="([^"]+)"', base64.b64decode(b64).decode("utf-8")).group(1)
    assert stroke(SN.row_css(3, page, True)) == SN.ICON_ON_ACCENT
    assert stroke(SN.row_css(3, page, False)) == SN.ICON_ON_DARK


def test_every_category_page_has_a_base_icon_file():
    missing = []
    for _cat, pages in SN.CATEGORIES:
        for p in pages:
            path = SN.icon_path(p)
            if not (path and os.path.isfile(path)):
                missing.append(p)
    assert missing == []


def test_row_css_empty_for_page_without_icon():
    assert SN.row_css(1, "No such page", True) == ""
