"""The commit guard must block forbidden staged files and allow the rest (no git call in the test)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import check_staged as guard  # noqa: E402


def test_blocks_sample_store_and_test_data():
    assert guard.blocked(["app/sample_store.json"]) == ["app/sample_store.json"]
    assert guard.blocked(["app/tests/data/x.csv"]) == ["app/tests/data/x.csv"]


def test_blocks_core_and_nav_icons():
    assert guard.blocked(["app/mfs/reticon/core.py", "app/assets/nav_icons/home.svg"]) == [
        "app/assets/nav_icons/home.svg", "app/mfs/reticon/core.py"]


def test_allows_icon_work():
    names = ["app/tools/export_icons.py", "docs/icons/STYLE.md", "app/assets/icons_proof/contact_sheet.png"]
    assert guard.blocked(names) == []
