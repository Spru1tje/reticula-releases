"""docs/icons/INDEX.md must match what tools/build_icon_index.py generates from the proof SVGs."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import build_icon_index as idx  # noqa: E402


def test_index_matches_generated():
    assert idx.INDEX.read_bytes().decode("utf-8") == idx.build()  # byte-exact, LF
