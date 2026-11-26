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


def find_browser() -> str:
    for path in BROWSERS:
        if pathlib.Path(path).is_file():
            return path
    raise SystemExit("No headless Edge or Chrome found; install one or edit BROWSERS.")


def render(browser: str, svg: pathlib.Path, png: pathlib.Path, size: int) -> None:
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

    browser = find_browser()
    args.out.mkdir(parents=True, exist_ok=True)
    svgs = sorted(args.src.glob("*.svg"))
    if not svgs:
        print(f"no SVG files in {args.src}")
        return 1
    for svg in svgs:
        for size in SIZES:
            png = args.out / f"{svg.stem}_{size}.png"
            render(browser, svg, png, size)
            print(f"ok {png.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
