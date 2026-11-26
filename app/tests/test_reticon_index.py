"""docs/icons/INDEX.md must match what tools/build_icon_index.py generates from the proof SVGs."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import build_icon_index as idx  # noqa: E402


def test_index_matches_generated():
    assert idx.INDEX.read_bytes().decode("utf-8") == idx.build()  # byte-exact, LF


def test_index_row_count_equals_svg_count():
    """Every SVG in the proof set has exactly one row in INDEX.md, and no row lacks an SVG."""
    svgs = sorted(p.stem for p in idx.PROOF.glob("*.svg"))
    rows = [
        line for line in idx.INDEX.read_bytes().decode("utf-8").splitlines()
        if line.startswith("| ") and not line.startswith("| Icon ")
    ]
    names = [line.split("|")[1].strip() for line in rows]
    assert len(names) == len(set(names)), "duplicate index rows"
    assert sorted(names) == svgs, f"index rows differ from SVG files: rows={len(names)} svgs={len(svgs)}"
