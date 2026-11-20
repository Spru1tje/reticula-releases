"""
Quality checks for the library. Run:  python -m reticon.qa  [--update-goldens]

  lint       grid bounds, stroke widths, 24 px legibility (ink coverage), missing tags
  lookalike  perceptual-hash pairs that are hard to tell apart (originals + CC0 sources)
  palette    colour-blind distinguishability of each theme's liquid colours (deuteranopia / protanopia)
  regress    compare every shape x style against goldens/dhash.json (catches silent renderer breakage)

Needs rsvg-convert + Pillow (same as the build).
"""
from __future__ import annotations
from .raster import svg_to_png

import io
import json
import math
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # only for the string annotation below
    from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
GOLDEN = ROOT / "tests" / "goldens" / "dhash.json"


def raster(svg: str, px: int) -> "Image.Image":
    from PIL import Image
    png = svg_to_png(svg, px, px, "white")
    return Image.open(io.BytesIO(png)).convert("L")


def dhash(img, n=8) -> int:
    im = img.resize((n + 1, n))
    p = list(im.get_flattened_data()) if hasattr(im, 'get_flattened_data') else list(im.getdata())
    bits = 0
    for r in range(n):
        for c in range(n):
            bits = (bits << 1) | (p[r * (n + 1) + c] > p[r * (n + 1) + c + 1])
    return bits


def ham(a: int, b: int) -> int:
    return bin(a ^ b).count("1")


# ----------------------------------------------------------------- lint
def lint(icons) -> list[dict]:
    from .core import render_svg, Style
    issues = []

    def one(ic):
        out = []
        wide = render_svg(ic, Style(mode="solid", color="#000000"), 96).replace('viewBox="0 0 64 64"', 'viewBox="-16 -16 96 96"')
        im = raster(wide, 96)
        px = im.load()
        outside = sum(1 for y in range(96) for x in range(96)
                      if px[x, y] < 128 and not (14 <= x <= 81 and 14 <= y <= 81))
        if outside > 4:
            out.append(("bounds", f"{outside} px of ink outside the 64 grid"))
        thin = [o.stroke for o in ic.ops if o.stroke and o.stroke < 1.0 and o.role != "x"]
        if thin:
            out.append(("thin", f"{len(thin)} strokes thinner than 1.0 (vanish at 24 px)"))
        img = raster(render_svg(ic, Style(mode="solid", color="#000000"), 24), 24)
        data = img.get_flattened_data() if hasattr(img, 'get_flattened_data') else img.getdata()
        ink = sum(1 for v in data if v < 128) / (24 * 24)
        if ink < 0.06:
            out.append(("faint", f"only {ink:.0%} ink at 24 px"))
        elif ink > 0.72:
            out.append(("blob", f"{ink:.0%} ink at 24 px: details merge"))
        if not ic.tags:
            out.append(("tags", "no search tags"))
        return [{"icon": ic.name, "check": c, "detail": d} for c, d in out]
    with ThreadPoolExecutor(8) as ex:
        for r in ex.map(one, icons):
            issues += r
    return issues


# ----------------------------------------------------------------- look-alikes
def lookalikes(entries: list[tuple[str, str]], threshold: int = 2) -> list[tuple[str, str, int]]:
    """entries: (name, svg). Returns pairs whose 32 px dHash differs by <= threshold bits."""
    with ThreadPoolExecutor(8) as ex:
        hashes = list(ex.map(lambda e: (e[0], dhash(raster(e[1], 32))), entries))
    pairs = []
    for i in range(len(hashes)):
        for j in range(i + 1, len(hashes)):
            d = ham(hashes[i][1], hashes[j][1])
            if d <= threshold:
                pairs.append((hashes[i][0], hashes[j][0], d))
    return sorted(pairs, key=lambda t: t[2])


# ----------------------------------------------------------------- palette
_M = {"deuteranopia": ((0.367, 0.861, -0.228), (0.280, 0.673, 0.047), (-0.012, 0.043, 0.969)),
      "protanopia": ((0.152, 1.053, -0.205), (0.115, 0.786, 0.099), (-0.004, -0.048, 1.052))}


def _lin(c):
    c /= 255
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def _lab(rgb):
    r, g, b = (_lin(v) for v in rgb)
    x = (0.4124 * r + 0.3576 * g + 0.1805 * b) / 0.9505
    y = 0.2126 * r + 0.7152 * g + 0.0722 * b
    z = (0.0193 * r + 0.1192 * g + 0.9505 * b) / 1.089
    f = lambda t: t ** (1 / 3) if t > 0.008856 else 7.787 * t + 16 / 116
    return 116 * f(y) - 16, 500 * (f(x) - f(y)), 200 * (f(y) - f(z))


def _sim(hexc, kind):
    rgb = [int(hexc[i:i + 2], 16) for i in (1, 3, 5)]
    m = _M[kind]
    return [max(0, min(255, sum(m[i][j] * rgb[j] for j in range(3)))) for i in range(3)]


def palette_report(themes) -> list[dict]:
    out = []
    for name, pal in themes.items():
        cols = list(pal.liquids)
        for kind in ("normal", "deuteranopia", "protanopia"):
            labs = [_lab([int(c[i:i + 2], 16) for i in (1, 3, 5)] if kind == "normal" else _sim(c, kind)) for c in cols]
            worst = min((math.dist(labs[i], labs[j]), cols[i], cols[j]) for i in range(len(cols)) for j in range(i + 1, len(cols))) \
                if len(cols) > 1 else (99, "", "")
            out.append({"theme": name, "vision": kind, "min_delta_e": round(worst[0], 1), "pair": worst[1:],
                        "ok": worst[0] >= 12})
    return out


# ----------------------------------------------------------------- regression
def load_specs(path) -> dict:
    """specs.json merged with specs_variants.json (all 6 shapes per icon)."""
    path = Path(path)
    specs = json.loads(path.read_text(encoding="utf-8"))
    vp = path.with_name("specs_variants.json")
    if vp.exists():
        for n, sh in json.loads(vp.read_text(encoding="utf-8")).items():
            if n in specs["icons"]:
                specs["icons"][n]["shapes"].update(sh)
    return specs


def regress(specs: dict, update: bool = False) -> list[str]:
    from .runtime import Library
    lib = Library(specs)
    jobs = [(n, sh, st) for n, v in specs["icons"].items() for sh in v["shapes"] for st in ("illustration", "solid")]
    with ThreadPoolExecutor(8) as ex:
        got = dict(zip([f"{n}|{sh}|{st}" for n, sh, st in jobs],
                       ex.map(lambda j: dhash(raster(lib.svg(j[0], j[1], j[2], metadata=False), 48), 12), jobs)))
    if update or not GOLDEN.exists():
        GOLDEN.parent.mkdir(parents=True, exist_ok=True)
        GOLDEN.write_text(json.dumps({k: format(v, "x") for k, v in got.items()}, indent=0), encoding="utf-8")
        return []
    gold = {k: int(v, 16) for k, v in json.loads(GOLDEN.read_text(encoding="utf-8")).items()}
    return [f"{k}: changed ({ham(v, gold[k])} bits)" for k, v in got.items() if k in gold and ham(v, gold[k]) > 6] + \
           [f"{k}: new (no golden)" for k in got if k not in gold]


def main():
    import argparse
    sys.path.insert(0, str(ROOT))
    from reticon import all_icons
    from reticon.core import THEMES
    ap = argparse.ArgumentParser()
    ap.add_argument("--update-goldens", action="store_true")
    ap.add_argument("--specs", default=str(ROOT / "dist" / "runtime" / "specs.json"))
    ap.add_argument("--out", default=str(ROOT / "dist" / "qa_report.json"))
    a = ap.parse_args()
    icons = all_icons()
    rep = {"lint": lint(icons), "palette": palette_report(THEMES)}
    specs = load_specs(a.specs)
    from .runtime import Library
    lib = Library(specs, Path(a.specs).parent)
    ent = [(n, lib.svg(n, metadata=False)) for n in specs["icons"]]
    ent += [(n, lib.svg(n)) for n in list(specs.get("external", {}))]
    rep["lookalike"] = lookalikes(ent)
    rep["regress"] = regress(specs, a.update_goldens)
    Path(a.out).write_text(json.dumps(rep, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"lint issues {len(rep['lint'])} | look-alike pairs {len(rep['lookalike'])} | "
          f"palette fails {sum(not p['ok'] for p in rep['palette'])} | regressions {len(rep['regress'])} -> {a.out}")


if __name__ == "__main__":
    main()
