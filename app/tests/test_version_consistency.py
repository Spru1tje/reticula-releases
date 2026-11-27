"""The version in mfs/__init__.py matches the newest CHANGELOG heading and the release-notes file name (RELEASE_NOTES_v<version>.md)."""
import os
import re

from mfs import __version__

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(APP)


def _newest_changelog_version():
    with open(os.path.join(APP, "CHANGELOG.md"), encoding="utf-8") as fh:
        for line in fh:
            m = re.match(r"## (\d+\.\d+\.\d+)\b", line)
            if m:
                return m.group(1)
    raise AssertionError("no version heading found in CHANGELOG.md")


def test_version_matches_newest_changelog_heading():
    assert _newest_changelog_version() == __version__


def test_release_notes_file_exists_for_the_version():
    assert os.path.isfile(os.path.join(ROOT, f"RELEASE_NOTES_v{__version__}.md"))
