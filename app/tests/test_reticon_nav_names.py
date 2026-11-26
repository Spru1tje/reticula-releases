"""Every navigation icon name used by the app must exist in assets/nav_icons.

The names come from ICON_FILE in mfs/sidebar_nav.py, read with ast (the app module is not imported).
A renamed or removed legacy icon fails here instead of breaking the sidebar at run time.
"""
import ast
import sys
from pathlib import Path

APP = Path(__file__).resolve().parents[1]
SIDEBAR = APP / "mfs" / "sidebar_nav.py"
NAV = APP / "assets" / "nav_icons"


def _icon_file_names(sidebar: Path) -> dict:
    tree = ast.parse(sidebar.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "ICON_FILE" for t in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError(f"ICON_FILE not found in {sidebar}")


def test_every_nav_icon_name_exists():
    mapping = _icon_file_names(SIDEBAR)
    assert mapping, "ICON_FILE is empty"
    missing = {page: name for page, name in mapping.items() if not (NAV / f"{name}.svg").is_file()}
    assert not missing, f"nav icon files missing for: {missing}"


def test_detects_a_missing_nav_icon(monkeypatch, tmp_path):
    """Runs the real test function against a fake sidebar that names a missing icon."""
    this = sys.modules[__name__]
    fake = tmp_path / "sidebar_nav.py"
    fake.write_text('ICON_FILE = {"Home": "home", "Nope": "does_not_exist"}\n', encoding="utf-8")
    monkeypatch.setattr(this, "SIDEBAR", fake)
    try:
        this.test_every_nav_icon_name_exists()
    except AssertionError as exc:
        assert "does_not_exist" in str(exc)
    else:
        raise AssertionError("a missing nav icon was not detected")
