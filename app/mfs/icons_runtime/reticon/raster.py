"""SVG -> PNG / PDF without a system program.

Uses `rsvg-convert` (librsvg) when it is on the PATH. Otherwise falls back to pure-Python packages:
  * PNG: `resvg-py` (resvg, a Rust renderer shipped as a wheel)
  * PDF: `svglib` + `reportlab` (vector output, no bitmaps)
Install the fallbacks with: pip install resvg-py svglib reportlab
"""
from __future__ import annotations

import io
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


def _rsvg() -> str | None:
    return shutil.which("rsvg-convert")


def has_png_backend() -> bool:
    if _rsvg():
        return True
    try:
        import resvg_py  # noqa: F401
        return True
    except ImportError:
        return False


def has_pdf_backend() -> bool:
    if _rsvg():
        return True
    try:
        import svglib  # noqa: F401
        import reportlab  # noqa: F401
        return True
    except ImportError:
        return False


def svg_to_png(svg: str, width: int, height: int | None = None, background: str | None = None) -> bytes:
    """PNG bytes of the SVG at the given pixel width (and height, if given)."""
    exe = _rsvg()
    if exe:
        args = [exe, "-w", str(width)] + (["-h", str(height)] if height else []) + (["-b", background] if background else [])
        return subprocess.run(args, input=svg.encode(), capture_output=True, check=True).stdout
    try:
        import resvg_py
    except ImportError as e:
        raise RuntimeError("PNG export needs rsvg-convert or the package resvg-py") from e
    # resvg rejects a root with physical sizes (e.g. width="43.2mm") as "invalid size". The viewBox already
    # carries the aspect ratio, so the physical size is dropped and the pixel size passed instead.
    return bytes(resvg_py.svg_to_bytes(svg_string=_strip_root_size(svg), width=width, height=height, background=background))


def _strip_root_size(svg: str) -> str:
    m = re.search(r"<svg\b[^>]*>", svg, re.S)
    if not m:
        return svg
    root = re.sub(r'\s(width|height)="[^"]*"', "", m.group(0), count=2)
    return svg[:m.start()] + root + svg[m.end():]


def svg_to_pdf(svg: str, path: str | Path) -> None:
    """Write a vector PDF of the SVG to `path`."""
    exe = _rsvg()
    if exe:
        subprocess.run([exe, "-f", "pdf", "-o", str(path)], input=svg.encode(), check=True)
        return
    try:
        from svglib.svglib import svg2rlg
        from reportlab.graphics import renderPDF
    except ImportError as e:
        raise RuntimeError("PDF export needs rsvg-convert or the packages svglib and reportlab") from e
    with tempfile.NamedTemporaryFile("w", suffix=".svg", delete=False, encoding="utf-8") as fh:
        fh.write(svg)
        tmp = fh.name
    try:
        drawing = svg2rlg(tmp)
        if drawing is None:
            raise RuntimeError("svglib could not read the scheme")
        renderPDF.drawToFile(drawing, str(path))
    finally:
        Path(tmp).unlink(missing_ok=True)


def svg_to_image(svg: str, px: int) -> "bytes":
    """Small helper for the QA and preview tools: PNG bytes, white background, square."""
    return svg_to_png(svg, px, px, "white")
