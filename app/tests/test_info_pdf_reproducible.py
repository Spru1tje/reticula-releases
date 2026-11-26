"""Info-PDF build is reproducible: two builds of the same document give identical bytes (SOURCE_DATE_EPOCH is fixed in the script)."""
import os

from tools import build_info_pdf as B


def test_two_builds_of_pxrd_are_byte_identical(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    assert B.build(["PXRD"], str(a)) == 0
    assert B.build(["PXRD"], str(b)) == 0
    pa, pb = a / "pxrd.pdf", b / "pxrd.pdf"
    assert pa.read_bytes() == pb.read_bytes()


def test_build_pins_the_pdf_timestamp():
    assert os.environ.get("SOURCE_DATE_EPOCH") == "1795651200"
