"""Presence check for numbers in the info documents (prose only; figure captions and formulas are not checked).

A number counts as sourced when the same number occurs somewhere in the code, the tests, or the notes (CHANGELOG, README,
docs). This is only a presence test: it does not check the unit, the quantity, or that a sentence matches its source.
Numbers 0 to 12 are skipped as structural.

The test counts the prose numbers that have no occurrence in the corpus. It stores only that count, as a ceiling
(MAX_UNSOURCED): it fails when the count rises above the ceiling. The numbers themselves are never stored or printed,
because some come from private test images."""
import glob
import os
import re

from mfs import helpdocs as HD

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(APP)
SELF = "test_helpdoc_numbers_sourced.py"

MAX_UNSOURCED = 3

NUM = re.compile(r"(?<![\w.])(\d+(?:[.,]\d+)?)(?!\d)")


def _corpus():
    files = []
    for pat in ["app/mfs/**/*.py", "app/tests/**/*.py", "app/tools/**/*.py", "app/app.py", "app/CHANGELOG.md",
                "app/README.md", "app/docs/*.md", "README.md", "RELEASE_NOTES*.md"]:
        files += glob.glob(os.path.join(ROOT, pat), recursive=True)
    out = []
    for f in files:
        base = os.path.basename(f)
        if base in ("helpdocs.py", SELF) or "icons_runtime" in f or "scheme_canvas" in f:
            continue
        with open(f, encoding="utf-8", errors="ignore") as fh:
            out.append(fh.read())
    return "\n".join(out)


def _prose_numbers_by_doc():
    for key, doc in HD.DOCS.items():
        for kind, val in doc["sections"]:
            if kind in ("p", "b", "n", "k", "h") and isinstance(val, str):
                for tok in NUM.findall(val):
                    if tok.isdigit() and int(tok) <= 12:
                        continue
                    yield key, tok


def test_prose_numbers_are_checked_against_the_corpus():
    corpus = _corpus()
    # distinct values per document: the same number repeated in one document counts once
    unknown = len({(key, tok) for key, tok in _prose_numbers_by_doc()
                   if not re.search(r"(?<![\d.])" + re.escape(tok) + r"(?!\d)", corpus)})
    assert unknown <= MAX_UNSOURCED, (
        f"prose numbers without an occurrence in the corpus rose to {unknown} (ceiling {MAX_UNSOURCED}); "
        "values are not shown: check the new number(s) in the info documents")
