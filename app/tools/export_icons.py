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


CELL = 240  # px per row; a row shows one icon at each size
LABEL_W = 200


def build_contact_sheet(browsers: list[str], svgs: list[pathlib.Path], sheet_png: pathlib.Path) -> None:
    """Screenshot a review sheet: one row per icon, each shown at 24, 48, 96 and 192 px.

    Generated automatically after the PNG export, so the visual review is always current.
    Tries each browser in order, like render().
    """
    rows = []
    for svg in svgs:
        imgs = "".join(
            f'<img src="{svg.resolve().as_uri()}" width="{s}" height="{s}" style="display:block;margin-right:16px">'
            for s in SIZES
        )
        rows.append(
            '<div style="display:flex;align-items:flex-end;height:%dpx;border-bottom:1px solid #DDD;">'
            '<div style="width:%dpx;font:13px sans-serif;color:#1F2D3D">%s</div>%s</div>'
            % (CELL, LABEL_W, svg.stem, imgs)
        )
    page = (
        '<html><head><meta charset="utf-8"></head>'
        '<body style="margin:0;padding:12px;background:#FFFFFF">' + "".join(rows) + "</body></html>"
    )
    width = LABEL_W + sum(SIZES) + 16 * len(SIZES) + 24
    height = CELL * len(svgs) + 24
    errors: list[str] = []
    with tempfile.TemporaryDirectory() as tmp:
        page_path = pathlib.Path(tmp) / "sheet.html"
        page_path.write_text(page, encoding="utf-8")
        for browser in browsers:
            sheet_png.unlink(missing_ok=True)
            try:
                subprocess.run(
                    [
                        browser, "--headless=new", "--disable-gpu", "--no-sandbox", "--hide-scrollbars",
                        f"--window-size={width},{height}", f"--screenshot={sheet_png.resolve()}",
                        page_path.resolve().as_uri(),
                    ],
                    check=True, capture_output=True, timeout=120,
                )
            except (OSError, subprocess.TimeoutExpired, subprocess.CalledProcessError) as exc:
                errors.append(f"{pathlib.Path(browser).name}: {type(exc).__name__}")
                continue
            if sheet_png.is_file() and sheet_png.stat().st_size > 0:
                return
            errors.append(f"{pathlib.Path(browser).name}: no PNG written")
    raise RuntimeError("contact sheet failed: " + "; ".join(errors))


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
    sheet = args.out.parent / "contact_sheet.png"
    try:
        build_contact_sheet(browsers, svgs, sheet)
        print(f"ok {sheet.name}")
    except RuntimeError as exc:
        failures += 1
        print(f"FAIL {exc}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
