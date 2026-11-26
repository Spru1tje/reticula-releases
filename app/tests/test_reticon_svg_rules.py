"""Reticon icon rules from docs/icons/STYLE.md, checked on every proof SVG.

1. One line weight: every stroke resolves to stroke-width 1.5 (own attribute or inherited).
2. Shadow: a shadow group is offset by (2, 2) and uses the muted colour at opacity 0.3.
3. One accent: at most one element outside the shadow uses the accent colour #2E6F95.
"""
import re
import xml.etree.ElementTree as ET
from pathlib import Path

PROOF = Path(__file__).resolve().parents[1] / "assets" / "icons_proof"
LINE_WIDTH = "1.5"
MUTED = "#8FA3B5"
ACCENT = "#2E6F95"
SHADOW_OFFSET = "translate(2 2)"


def _svgs():
    files = sorted(PROOF.glob("*.svg"))
    assert files, f"no proof SVGs in {PROOF}"
    return files


def _local(tag: str) -> str:
    return tag.split("}")[-1]


def _walk(el, inherited_width=None, in_shadow=False):
    """Yield (element, effective stroke-width, inside_shadow) for every element."""
    width = el.get("stroke-width", inherited_width)
    shadow = in_shadow or (_local(el.tag) == "g" and el.get("transform") == SHADOW_OFFSET)
    yield el, width, shadow
    for child in el:
        yield from _walk(child, width, shadow)


def test_stroke_width_is_one_and_a_half():
    for f in _svgs():
        for el, width, _ in _walk(ET.parse(f).getroot()):
            if el.get("stroke") not in (None, "none"):
                assert width == LINE_WIDTH, f"{f.name}: <{_local(el.tag)}> stroke-width={width}"


def test_shadow_group_offset_and_opacity():
    for f in _svgs():
        for el, _, _ in _walk(ET.parse(f).getroot()):
            # any offset group drawn in the muted shadow colour is a shadow and must follow the rule
            if _local(el.tag) == "g" and (el.get("transform") or "").startswith("translate(") and (
                el.get("stroke") == MUTED or el.get("fill") == MUTED
            ):
                assert el.get("transform") == SHADOW_OFFSET, f"{f.name}: shadow transform {el.get('transform')}"
                assert el.get("stroke") == MUTED, f"{f.name}: shadow stroke {el.get('stroke')}"
                fill_op = el.get("fill-opacity")
                stroke_op = el.get("stroke-opacity")
                assert "0.3" in (fill_op, stroke_op), f"{f.name}: shadow opacity fill={fill_op} stroke={stroke_op}"


def test_at_most_one_accent_element():
    for f in _svgs():
        accents = []
        for el, _, shadow in _walk(ET.parse(f).getroot()):
            if shadow or _local(el.tag) == "g":
                continue
            if ACCENT in (el.get("fill"), el.get("stroke")):
                accents.append(_local(el.tag))
        assert len(accents) <= 1, f"{f.name}: {len(accents)} accent elements ({accents})"


# Negative checks: each rule must fail on a broken SVG. The proof folder is swapped for a temp folder.
def _with_svg(monkeypatch, tmp_path, body):
    (tmp_path / "bad.svg").write_text(body, encoding="utf-8")
    monkeypatch.setattr(__import__(__name__), "PROOF", tmp_path)


def test_stroke_rule_detects_wrong_width(monkeypatch, tmp_path):
    _with_svg(monkeypatch, tmp_path,
              '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64">'
              '<path stroke="#1F2D3D" stroke-width="2" d="M0 0"/></svg>')
    try:
        test_stroke_width_is_one_and_a_half()
    except AssertionError as exc:
        assert "stroke-width=2" in str(exc)
    else:
        raise AssertionError("a wrong stroke width was not detected")


def test_shadow_rule_detects_wrong_offset(monkeypatch, tmp_path):
    _with_svg(monkeypatch, tmp_path,
              '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64">'
              '<g transform="translate(3 2)" stroke="#8FA3B5" stroke-opacity="0.3"><path d="M0 0"/></g></svg>')
    try:
        test_shadow_group_offset_and_opacity()
    except AssertionError as exc:
        assert "translate(3 2)" in str(exc)
    else:
        raise AssertionError("a wrong shadow offset was not detected")


def test_accent_rule_detects_two_accents(monkeypatch, tmp_path):
    _with_svg(monkeypatch, tmp_path,
              '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64">'
              '<path stroke="#2E6F95" d="M0 0"/><path stroke="#2E6F95" d="M1 1"/></svg>')
    try:
        test_at_most_one_accent_element()
    except AssertionError as exc:
        assert "2 accent elements" in str(exc)
    else:
        raise AssertionError("two accent elements were not detected")
