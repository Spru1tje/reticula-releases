"""Reticon proof set: only the four palette tokens; no gradients, filters, text, images, masks or <style>.

Rules come from docs/icons/STYLE.md (palette, shadow, not-allowed list) and the 64-unit grid.
The allowed tokens are read from STYLE.md, so the guide is the single source (no copy in this test).
Scope: icon SVGs in assets/icons_proof only. The data palette (liquids, okabe_ito, pastel, vivid in
mfs/reticon/core.py) is deliberately not checked here; see STYLE.md "Data palette".
"""
import re
import xml.etree.ElementTree as ET
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PROOF = Path(__file__).resolve().parents[1] / "assets" / "icons_proof"
STYLE = REPO / "docs" / "icons" / "STYLE.md"
COLOUR_ATTRS = ("fill", "stroke", "stop-color")
FORBIDDEN = {"linearGradient", "radialGradient", "filter", "text", "image", "mask", "clipPath", "pattern", "style"}
HEX = re.compile(r"#[0-9A-Fa-f]{6}\b")


def _tokens_from_style() -> set[str]:
    """The four tokens are the bullet lines under '## Palette' that start with a name and a hex code."""
    text = STYLE.read_text(encoding="utf-8")
    section = text.split("## Palette", 1)[1].split("\n## ", 1)[0]
    found = set()
    for line in section.splitlines():
        m = re.match(r"\s*-\s+(ink|accent|muted|paper)\b.*?`(#[0-9A-Fa-f]{6})`", line)
        if m:
            found.add(m.group(2).upper())
    return found


ALLOWED = _tokens_from_style()


def _svgs():
    files = sorted(PROOF.glob("*.svg"))
    assert files, f"no proof SVGs in {PROOF}"
    return files


def _local(tag: str) -> str:
    return tag.split("}")[-1]


def test_style_guide_defines_exactly_four_tokens():
    assert ALLOWED == {"#1F2D3D", "#2E6F95", "#8FA3B5", "#FFFFFF"}, ALLOWED


def test_colours_are_only_palette_tokens():
    for f in _svgs():
        root = ET.parse(f).getroot()
        for el in root.iter():
            for attr in COLOUR_ATTRS:
                val = el.get(attr)
                if val in (None, "none"):
                    continue
                assert val.upper() in ALLOWED, f"{f.name}: {_local(el.tag)} {attr}={val}"
            style = el.get("style")
            if style:
                for hexcode in HEX.findall(style):
                    assert hexcode.upper() in ALLOWED, f"{f.name}: style colour {hexcode}"


def test_no_forbidden_elements():
    for f in _svgs():
        for el in ET.parse(f).getroot().iter():
            assert _local(el.tag) not in FORBIDDEN, f"{f.name}: <{_local(el.tag)}>"


def test_grid_is_64_units():
    for f in _svgs():
        assert ET.parse(f).getroot().get("viewBox") == "0 0 64 64", f.name
