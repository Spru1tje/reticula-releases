"""
Variants for external (CC0) SVG files: colour variants x shape variants,
done by rewriting the SVG text. Works on any SVG file you drop in.
"""
from __future__ import annotations
from .raster import svg_to_png

import colorsys
import io
import re

HEX = re.compile(r"#([0-9a-fA-F]{6}|[0-9a-fA-F]{3})\b")
RGB = re.compile(r"rgb\(\s*([\d.]+)(%?)\s*,\s*([\d.]+)%?\s*,\s*([\d.]+)%?\s*\)")
ROOT = re.compile(r"<svg\b[^>]*>", re.S)

ILLUSTRATION_COLOURS = ("original", "hue120", "hue240", "pastel", "greyscale")
GLYPH_COLOURS = {"ink": "#2F3B4A", "blue": "#3B82C4", "teal": "#1F9E8F", "coral": "#E06B5A", "purple": "#7C5CC4"}


def _hex(r, g, b):
    return "#{:02X}{:02X}{:02X}".format(*(max(0, min(255, round(x * 255))) for x in (r, g, b)))


def colour_fn(kind):
    def f(r, g, b):
        h, l, s = colorsys.rgb_to_hls(r, g, b)
        if kind == "greyscale":
            y = 0.299 * r + 0.587 * g + 0.114 * b
            return y, y, y
        if kind == "pastel":
            return colorsys.hls_to_rgb(h, l + (1 - l) * 0.35, s * 0.75)
        if kind.startswith("hue"):
            if s < 0.18 or l < 0.08 or l > 0.97:
                return r, g, b
            return colorsys.hls_to_rgb((h + int(kind[3:]) / 360) % 1, l, s)
        if kind.startswith("#"):
            return tuple(int(kind[i:i + 2], 16) / 255 for i in (1, 3, 5))
        return r, g, b
    return f


def recolour(text, fn):
    def hx(m):
        v = m.group(1)
        if len(v) == 3:
            v = "".join(c * 2 for c in v)
        return _hex(*fn(*(int(v[i:i + 2], 16) / 255 for i in (0, 2, 4))))

    def rg(m):
        pct = m.group(2) == "%"
        vals = [float(m.group(i)) for i in (1, 3, 4)]
        return _hex(*fn(*[v / 100 if pct else v / 255 for v in vals]))
    return RGB.sub(rg, HEX.sub(hx, text))


def colour_variants(svg: str, kind: str) -> dict[str, str]:
    if kind == "glyph":
        return {k: svg.replace("currentColor", c) for k, c in GLYPH_COLOURS.items()}
    return {k: svg if k == "original" else recolour(svg, colour_fn(k)) for k in ILLUSTRATION_COLOURS}


# ---------------------------------------------------------------- shapes

NS = re.compile(r'\sxmlns:([\w-]+)="([^"]+)"')
EDITOR_JUNK = re.compile(r"<sodipodi:namedview\b[^>]*?/>|<sodipodi:namedview\b.*?</sodipodi:namedview>", re.S)


def _split(svg):
    m = ROOT.search(svg)
    head = m.group(0)
    inner = EDITOR_JUNK.sub("", svg[m.end():svg.rfind("</svg>")])
    global _NS
    _NS = {k: v for k, v in NS.findall(head) if k != "xlink"}
    vb = re.search(r'viewBox="([^"]+)"', head)
    if vb:
        x, y, w, h = (float(v) for v in vb.group(1).replace(",", " ").split())
    else:
        def num(a):
            mm = re.search(rf'\b{a}="([\d.]+)', head)
            return float(mm.group(1)) if mm else 100.0
        x, y, w, h = 0.0, 0.0, num("width"), num("height")
    return inner, (x, y, w, h)


_NS: dict = {}


def _wrap(inner, vb, pad=0.0):
    x, y, w, h = vb
    px, py = w * pad, h * pad
    ns = "".join(f' xmlns:{k}="{v}"' for k, v in _NS.items())
    return (f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"{ns} '
            f'viewBox="{x - px:.2f} {y - py:.2f} {w + 2 * px:.2f} {h + 2 * py:.2f}" width="256" height="256">'
            f"{inner}</svg>")


REF = re.compile(r'(\bid="|url\(#|href="#)([^")]+)')


def _rename_ids(inner, suffix):
    ids = set(re.findall(r'\bid="([^"]+)"', inner))
    if not ids:
        return inner
    return REF.sub(lambda m: m.group(1) + m.group(2) + suffix if m.group(2) in ids else m.group(0), inner)


def mirrored(svg):
    inner, (x, y, w, h) = _split(svg)
    return _wrap(f'<g transform="translate({2 * x + w:.2f} 0) scale(-1 1)">{inner}</g>', (x, y, w, h))


def tilted(svg):
    inner, (x, y, w, h) = _split(svg)
    cx, cy = x + w / 2, y + h / 2
    return _wrap(f'<g transform="translate({cx:.2f} {cy:.2f}) rotate(-14) scale(0.88) '
                 f'translate({-cx:.2f} {-cy:.2f})">{inner}</g>', (x, y, w, h))


def group(svg):
    inner, (x, y, w, h) = _split(svg)
    parts = []
    for k, (ax, ay, s) in enumerate(((0.0, 0.03, 0.56), (0.42, 0.03, 0.56), (0.2, 0.28, 0.66))):
        body = inner if k == 2 else _rename_ids(inner, f"_g{k}")
        parts.append(f'<g transform="translate({x + ax * w:.2f} {y + ay * h:.2f}) scale({s}) '
                     f'translate({-x:.2f} {-y:.2f})">{body}</g>')
    return _wrap("".join(parts), (x, y, w, h))


def depth3d(svg):
    inner, (x, y, w, h) = _split(svg)
    grey = recolour(_rename_ids(inner, "_d"), colour_fn("#8592A3"))
    layers = "".join(f'<g transform="translate({w * d:.2f} {h * d:.2f})" opacity="0.55">{grey}</g>'
                     for d in (0.036, 0.024, 0.012))
    return _wrap(layers + inner, (x, y, w, h), pad=0.03)


SYM_CACHE: dict = {}


def is_symmetric(svg) -> bool:
    import hashlib
    key = hashlib.sha1(svg.encode()).hexdigest()
    if key in SYM_CACHE:
        return SYM_CACHE[key]
    SYM_CACHE[key] = r = _is_symmetric(svg)
    return r


def _is_symmetric(svg) -> bool:
    try:
        from PIL import Image, ImageChops, ImageOps
        png = svg_to_png(svg, 48, 48, "white")
        im = Image.open(io.BytesIO(png)).convert("L")
        diff = ImageChops.difference(im, ImageOps.mirror(im))
        return sum(diff.getdata()) / (im.width * im.height * 255) < 0.035
    except Exception:
        return False


def _badge_snippet(key: str, pid: str) -> str:
    """The Reticon state badge (64-grid) as plain coloured SVG, for embedding in external icons."""
    import copy as _c
    from .core import Icon, Op, render_svg, Style
    from . import variants as V
    ic = Icon("badge")
    ic.ops = [Op("m", "circle", {"cx": V.BX, "cy": V.BY, "r": 11.8, "_fill": "#FFFFFF"}),
              Op("m", "circle", {"cx": V.BX, "cy": V.BY, "r": 9.8, "_fill": V.BADGE_COL[key]})]
    for o in V._glyph(key):
        o = _c.deepcopy(o)
        o.role = "m"
        o.attrs["_fill"] = "#FFFFFF"
        ic.ops.append(o)
    out = render_svg(ic, Style(mode="color"), 64, id_prefix=pid, standalone=False)
    return out[out.index("</title>") + 8:out.rindex("</svg>")]


def _sparkles_snippet() -> str:
    from . import variants as V
    pts = [V._sparkle(53.5, 9.5, 5.2), V._sparkle(60, 19.5, 2.8), V._sparkle(8, 55, 3)]
    return "".join(f'<polygon points="{p.attrs["points"]}" fill="#F5A524"/>' for p in pts)


def highlighted(svg):
    inner, (x, y, w, h) = _split(svg)
    k = w / 64
    halo = f'<circle cx="{x + w / 2:.2f}" cy="{y + h / 2:.2f}" r="{w * 0.47:.2f}" fill="#F5A524" opacity="0.16"/>'
    body = (f'<g transform="translate({x + w * 0.06:.2f} {y + h * 0.06:.2f}) scale(0.88) translate({-x:.2f} {-y:.2f})">'
            f"{inner}</g>")
    spark = f'<g transform="translate({x:.2f} {y:.2f}) scale({k:.4f} {h / 64:.4f})">{_sparkles_snippet()}</g>'
    return _wrap(halo + body + spark, (x, y, w, h))


def badge(svg, key):
    inner, (x, y, w, h) = _split(svg)
    body = f'<g transform="translate({x:.2f} {y:.2f}) scale(0.84) translate({-x:.2f} {-y:.2f})">{inner}</g>'
    b = f'<g transform="translate({x:.2f} {y:.2f}) scale({w / 64:.4f} {h / 64:.4f})">{_badge_snippet(key, "xb-" + key)}</g>'
    return _wrap(body + b, (x, y, w, h))


def shape_variants(svg: str, allow_mirror: bool = True, category: str | None = None) -> dict[str, str]:
    """base + 5 variants, the same family size as the Reticon originals."""
    from . import variants as V
    from . import taxonomy as T
    out = {"base": svg}
    sym = is_symmetric(svg)
    if not allow_mirror:
        out["tilted"] = tilted(svg)
    elif sym:
        out["group"] = group(svg)
    else:
        out["mirrored"] = mirrored(svg)
    out["3d"] = depth3d(svg)
    out["highlighted"] = highlighted(svg)
    cat = T.unify(category) if category else ""
    for key in list(V.CAT_BADGES.get(cat, ("done", "timed"))) + ["group", "tilted"]:
        if len(out) > V.N_VARIANTS:
            break
        if key in out:
            continue
        out[key] = group(svg) if key == "group" else tilted(svg) if key == "tilted" else badge(svg, key)
    return out
