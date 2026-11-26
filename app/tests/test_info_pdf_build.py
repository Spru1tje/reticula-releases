"""Info-PDF build: tools/build_info_pdf.py renders PXRD and fails when a figure file is missing."""
import copy
import os

from mfs import helpdocs as HD
from tools import build_info_pdf as B


def test_pxrd_pdf_builds_into_a_folder(tmp_path):
    status = B.build(["PXRD"], str(tmp_path))
    assert status == 0
    pdf = tmp_path / "pxrd.pdf"
    assert pdf.is_file()
    assert pdf.read_bytes()[:5] == b"%PDF-"


def test_missing_figure_fails_the_build(tmp_path, monkeypatch):
    # PXRD has no figure; the TEM document has five, so it is the one to break
    doc = copy.deepcopy(HD.DOCS["TEM particle sizing (v2.21)"])
    idx = [i for i, (kind, _) in enumerate(doc["sections"]) if kind == "i"]
    assert idx, "TEM document should contain figures"
    kind, val = doc["sections"][idx[0]]
    doc["sections"][idx[0]] = (kind, {**val, "path": str(tmp_path / "no_such_figure.png")})
    monkeypatch.setitem(HD.DOCS, "PXRD-broken-figure", doc)
    assert B.check_figures(doc, "PXRD-broken-figure")
    assert B.build(["PXRD-broken-figure"], str(tmp_path)) == 1
    assert not (tmp_path / "pxrd-broken-figure.pdf").exists()


def test_unknown_key_fails_the_build(tmp_path):
    assert B.build(["no such document"], str(tmp_path)) == 1
