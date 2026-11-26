"""Reticon proof set at 24 px: exported PNGs must match the stored 24 px goldens pixel for pixel.

What this checks: change detection only. The goldens are copies of the exports made with
tools/export_icons.py (headless Edge), so a wrong render becomes the new golden; this test does
NOT prove the PNG shows the SVG correctly (that is the visual contact-sheet review).

Comparison is by decoded pixels (Pillow, RGBA) with a tolerance of 2 per channel, so a different
PNG encoder or a small renderer difference does not fail the test. The goldens were made with
headless Microsoft Edge 158.0.4370.33 on Windows 11 (Microsoft Windows NT 10.0.26200, win32). A
visible change in a render (more than the tolerance) fails the test; then re-export, review visually,
and update the goldens in the same commit with the reason.

Also checks the PNG size (24 x 24), and that a changed pixel is detected.
"""
import shutil
import tempfile
from pathlib import Path

from PIL import Image

PNG_DIR = Path(__file__).resolve().parents[1] / "assets" / "icons_proof" / "png"
GOLDEN = Path(__file__).resolve().parent / "golden" / "reticon_proof_24"
TOLERANCE = 2  # per channel, 0-255


def _rgba(path: Path) -> tuple[tuple[int, int], bytes]:
    with Image.open(path) as img:
        rgba = img.convert("RGBA")
        return rgba.size, rgba.tobytes()


def _same_pixels(a: Path, b: Path) -> bool:
    size_a, data_a = _rgba(a)
    size_b, data_b = _rgba(b)
    if size_a != size_b or len(data_a) != len(data_b):
        return False
    # RGBA bytes, so each byte is one channel value
    return all(abs(x - y) <= TOLERANCE for x, y in zip(data_a, data_b))


def _compare(png_dir: Path) -> None:
    golden = sorted(p.name for p in GOLDEN.glob("*.png"))
    assert golden, "no 24 px goldens"
    current = sorted(p.name for p in png_dir.glob("*_24.png"))
    assert current == golden, f"file set differs: current={current} golden={golden}"
    for name in golden:
        assert _same_pixels(png_dir / name, GOLDEN / name), (
            f"{name} differs from golden by more than {TOLERANCE} per channel. "
            "Check first whether the renderer changed: goldens were made with headless Edge "
            "158.0.4370.33 on Windows 11. If the Edge version or platform differs, re-export, "
            "review the PNG visually, and update the golden with a reason; if not, it is a real change."
        )


def test_24px_exports_match_goldens():
    _compare(PNG_DIR)


def test_24px_exports_are_24_by_24():
    for p in sorted(PNG_DIR.glob("*_24.png")):
        with Image.open(p) as img:
            assert img.size == (24, 24), p.name


def test_changed_pixel_is_detected():
    tmp = Path(tempfile.mkdtemp())
    try:
        for p in GOLDEN.glob("*.png"):
            shutil.copy(p, tmp / p.name)
        target = tmp / "pxrd_pattern_24.png"
        with Image.open(target) as img:
            rgba = img.convert("RGBA")
            rgba.putpixel((12, 12), (255, 0, 0, 255))
            rgba.save(target)
        try:
            _compare(tmp)
        except AssertionError as exc:
            assert "pxrd_pattern_24.png differs" in str(exc)
        else:
            raise AssertionError("a changed pixel was not detected")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
