"""Export Reticon SVG sources to PNG at fixed sizes (24, 48, 96, 192 px).

Renders each SVG with headless Microsoft Edge (or Chrome) via a tiny HTML page,
so no extra Python package is needed. Output: <out>/<name>_<size>.png

Usage:
    python tools/export_icons.py --src assets/icons_proof --out assets/icons_proof/png
"""
from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys
import tempfile

SIZES = (24, 48, 96, 192)
BROWSERS = (
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
)


def find_browsers() -> list[str]:
    found = [path for path in BROWSERS if pathlib.Path(path).is_file()]
    if not found:
        raise SystemExit("No headless Edge or Chrome found; install one or edit BROWSERS.")
    return found


def render(browsers: list[str], svg: pathlib.Path, png: pathlib.Path, size: int) -> None:
    """Try each browser in order; a timeout or missing PNG falls back to the next."""
    errors: list[str] = []
    for browser in browsers:
        try:
            png.unlink(missing_ok=True)  # an old PNG must never count as a fresh render
            _render_once(browser, svg, png, size)
        except (OSError, subprocess.TimeoutExpired, subprocess.CalledProcessError) as exc:
            errors.append(f"{pathlib.Path(browser).name}: {type(exc).__name__}")
            continue
        if png.is_file() and png.stat().st_size > 0:
            return
        errors.append(f"{pathlib.Path(browser).name}: no PNG written")
    raise RuntimeError(f"{svg.name} {size}px failed: " + "; ".join(errors))


def _render_once(browser: str, svg: pathlib.Path, png: pathlib.Path, size: int) -> None:
    page = (
        '<html><head><meta charset="utf-8"></head>'
        '<body style="margin:0;background:#FFFFFF">'
        f'<img src="{svg.resolve().as_uri()}" width="{size}" height="{size}" '
        'style="display:block">'
        "</body></html>"
    )
    with tempfile.TemporaryDirectory() as tmp:
        page_path = pathlib.Path(tmp) / "page.html"
        page_path.write_text(page, encoding="utf-8")
        cmd = [
            browser,
            "--headless=new",
            "--disable-gpu",
            "--no-sandbox",
            "--hide-scrollbars",
            f"--window-size={size},{size}",
            f"--screenshot={png.resolve()}",
            page_path.resolve().as_uri(),
        ]
        subprocess.run(cmd, check=True, capture_output=True, timeout=60)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True, type=pathlib.Path)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    args = ap.parse_args()

    browsers = find_browsers()
    args.out.mkdir(parents=True, exist_ok=True)
    svgs = sorted(args.src.glob("*.svg"))
    if not svgs:
        print(f"no SVG files in {args.src}")
        return 1
    failures = 0
    for svg in svgs:
        for size in SIZES:
            png = args.out / f"{svg.stem}_{size}.png"
            try:
                render(browsers, svg, png, size)
            except RuntimeError as exc:
                failures += 1
                print(f"FAIL {exc}")
                continue
            print(f"ok {png.name}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
