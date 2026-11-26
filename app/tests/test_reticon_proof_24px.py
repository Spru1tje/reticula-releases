"""Reticon proof set at 24 px: exported PNGs must match the stored 24 px goldens byte for byte.

What this checks: change detection only. The goldens are copies of the exports made with
tools/export_icons.py (headless Edge), so a wrong render becomes the new golden; this test does
NOT prove the PNG shows the SVG correctly (that is the visual contact-sheet review).
Byte equality depends on the renderer. The goldens were made with headless Microsoft Edge 158.0.4370.33
on Windows 11 (Microsoft Windows NT 10.0.26200, win32). On another Edge version or platform the
test can fail without the icon being wrong; then re-export, review visually, and update the goldens.
A deliberate change must update the golden in the same commit and say why.

Also checks the PNG header size (24 x 24), and that a changed PNG is detected.
"""
import shutil
import struct
import tempfile
from pathlib import Path

PNG_DIR = Path(__file__).resolve().parents[1] / "assets" / "icons_proof" / "png"
GOLDEN = Path(__file__).resolve().parent / "golden" / "reticon_proof_24"


def _png_size(path: Path) -> tuple[int, int]:
    head = path.read_bytes()[:24]
    assert head[:8] == b"\x89PNG\r\n\x1a\n", f"{path.name}: not a PNG"
    return struct.unpack(">II", head[16:24])


def _compare(png_dir: Path) -> None:
    golden = sorted(p.name for p in GOLDEN.glob("*.png"))
    assert golden, "no 24 px goldens"
    current = sorted(p.name for p in png_dir.glob("*_24.png"))
    assert current == golden, f"file set differs: current={current} golden={golden}"
    for name in golden:
        assert (png_dir / name).read_bytes() == (GOLDEN / name).read_bytes(), f"{name} differs from golden"


def test_24px_exports_match_goldens():
    _compare(PNG_DIR)


def test_24px_exports_are_24_by_24():
    for p in sorted(PNG_DIR.glob("*_24.png")):
        assert _png_size(p) == (24, 24), p.name


def test_changed_png_is_detected():
    tmp = Path(tempfile.mkdtemp())
    try:
        for p in GOLDEN.glob("*.png"):
            shutil.copy(p, tmp / p.name)
        target = tmp / "pxrd_pattern_24.png"
        target.write_bytes(target.read_bytes() + b"x")
        try:
            _compare(tmp)
        except AssertionError as exc:
            assert "pxrd_pattern_24.png differs" in str(exc)
        else:
            raise AssertionError("a changed PNG was not detected")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
