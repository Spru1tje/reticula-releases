"""Developer tool: build the info PDFs of the Help and info page from mfs/helpdocs.py (the only source of the texts).

Uses mfs.helppdf.render_pdf, the same renderer as the Help page. Output: docs/info/<slug>.pdf (default), or the folder given with --out.
Keys in CONCEPT_KEYS are written with a CONCEPT_ prefix: their text is not yet final (see the module's review list).

Usage:
    python tools/build_info_pdf.py                       # the default modules (TEM and PXRD)
    python tools/build_info_pdf.py --all                 # every document in helpdocs.DOCS
    python tools/build_info_pdf.py --key "PXRD" --out /tmp/pdf
"""
from __future__ import annotations

import argparse
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # the app folder
sys.path.insert(0, ROOT)

DEFAULT_KEYS = ["TEM particle sizing (v2.21)", "PXRD"]
CONCEPT_KEYS = {"TEM particle sizing (v2.21)"}                           # waits for the TEM changes (tasks b and c)
DEFAULT_OUT = os.path.join(ROOT, "docs", "info")


def slug(key: str) -> str:
    s = re.sub(r"[^a-z0-9]+", "-", key.lower()).strip("-")
    return s or "doc"


def check_figures(doc: dict, key: str) -> list[str]:
    """Figure files that a document names but that are missing in mfs/assets/docs (returns the missing paths)."""
    missing = []
    for kind, val in doc["sections"]:
        if kind == "i":
            path = val["path"]
            if not os.path.isfile(path):
                missing.append(path)
    return missing


def build(keys: list[str], out_dir: str) -> int:
    from mfs import helpdocs as HD, helppdf

    os.makedirs(out_dir, exist_ok=True)
    status = 0
    for key in keys:
        doc = HD.DOCS.get(key)
        if doc is None:
            print(f"ERROR: no document '{key}' in helpdocs.DOCS", file=sys.stderr)
            status = 1
            continue
        missing = check_figures(doc, key)
        if missing:
            print(f"ERROR: {key}: figure file(s) missing: {', '.join(missing)}", file=sys.stderr)
            status = 1
            continue
        name = ("CONCEPT_" if key in CONCEPT_KEYS else "") + slug(key) + ".pdf"
        pdf = helppdf.render_pdf(doc)
        if not pdf.startswith(b"%PDF-"):
            print(f"ERROR: {key}: renderer did not return a PDF", file=sys.stderr)
            status = 1
            continue
        path = os.path.join(out_dir, name)
        with open(path, "wb") as f:
            f.write(pdf)
        print(f"{key} -> {path} ({len(pdf)} bytes)")
    return status


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--key", action="append", help="document key in helpdocs.DOCS (repeat for more)")
    ap.add_argument("--all", action="store_true", help="every document in helpdocs.DOCS")
    ap.add_argument("--out", default=DEFAULT_OUT, help="output folder (default docs/info)")
    a = ap.parse_args()
    if a.all:
        from mfs import helpdocs as HD
        keys = list(HD.DOCS)
    else:
        keys = a.key or DEFAULT_KEYS
    return build(keys, a.out)


if __name__ == "__main__":
    raise SystemExit(main())
