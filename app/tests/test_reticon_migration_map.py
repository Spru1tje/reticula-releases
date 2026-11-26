"""Migration map: every legacy 24-grid nav icon maps to an existing 64-grid icon, and the index says so.

The map lives in tools/build_icon_index.py (LEGACY_MAP). This test keeps it honest against the files.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import build_icon_index as idx  # noqa: E402


def test_every_legacy_icon_exists():
    for legacy in idx.LEGACY_MAP:
        assert (idx.NAV_LEGACY / f"{legacy}.svg").is_file(), f"legacy icon missing: {legacy}"


def test_every_target_exists_and_has_metadata():
    for new in idx.LEGACY_MAP.values():
        assert (idx.PROOF / f"{new}.svg").is_file(), f"target icon missing: {new}"
        assert new in idx.METADATA, f"no metadata for: {new}"


def test_index_lists_the_mapping():
    text = idx.INDEX.read_bytes().decode("utf-8")
    for legacy, new in idx.LEGACY_MAP.items():
        assert f"- {legacy} -> {new}" in text, f"index does not list {legacy} -> {new}"


def test_bogus_mapping_is_detected(monkeypatch):
    """A mapping to a missing icon must fail the target test (negative check, kept in the file)."""
    monkeypatch.setitem(idx.LEGACY_MAP, "bogus", "does_not_exist")
    try:
        test_every_target_exists_and_has_metadata()
    except AssertionError as exc:
        assert "does_not_exist" in str(exc)
    else:
        raise AssertionError("a missing target icon was not detected")
