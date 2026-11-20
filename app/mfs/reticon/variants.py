"""
Shape variants: five real alternative drawings per icon (base + 5 = one "family").

Every variant has a stable id  <icon>--<key>  (e.g. round_bottom_flask--two-neck, beaker--crystals),
so Reticula can store it in a scene and the picker can open the whole family from any member.
The recipe depends on what the icon is:

  vessels (glassware with liquid)   -> "empty",     "suspension" (particles in liquid)
  data plots                        -> "alt-data",  "overlay"    (second series + fit)
  everything else                   -> "mirrored" (or "tilted" when mirroring is wrong),
                                       "3d"       (extruded, with depth)
Some icons get hand-written variants in OVERRIDES instead.
"""
from __future__ import annotations
from .raster import svg_to_png

import copy
import re

from .core import Icon, Op

NUM = re.compile(r"-?\d*\.?\d+(?:e-?\d+)?")
NO_MIRROR = {"sulfonate_group", "periodic_table", "heatmap"}   # letters / fixed layouts (+ every ghs_* / element_*)


# --------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------

def bbox(op: Op):
    a = op.attrs
    t = op.tag
    try:
        if t == "rect":
            x0, y0, x1, y1 = a["x"], a["y"], a["x"] + a["width"], a["y"] + a["height"]
        elif t == "circle":
            x0, y0, x1, y1 = a["cx"] - a["r"], a["cy"] - a["r"], a["cx"] + a["r"], a["cy"] + a["r"]
        elif t == "ellipse":
            m = max(a["rx"], a["ry"])
            x0, y0, x1, y1 = a["cx"] - m, a["cy"] - m, a["cx"] + m, a["cy"] + m
        elif t == "line":
            x0, x1 = sorted((a["x1"], a["x2"]))
            y0, y1 = sorted((a["y1"], a["y2"]))
        elif t in ("polygon", "polyline"):
            v = [float(n) for n in NUM.findall(a["points"])]
            xs, ys = v[0::2], v[1::2]
            x0, y0, x1, y1 = min(xs), min(ys), max(xs), max(ys)
        elif t == "path":
            d = a["d"]
            if re.search(r"[a-z]", d.replace("e-", "")):     # relative commands: skip
                return None
            toks = re.findall(r"[A-Za-z]|-?\d*\.?\d+", d)
            xs, ys, cmd, buf = [], [], "M", []
            cur = [0.0, 0.0]
            for tk in toks:
                if tk.isalpha():
                    cmd, buf = tk, []
                    continue
                buf.append(float(tk))
                if cmd == "H":
                    cur[0] = buf.pop(); xs.append(cur[0]); ys.append(cur[1])
                elif cmd == "V":
                    cur[1] = buf.pop(); xs.append(cur[0]); ys.append(cur[1])
                elif cmd == "A" and len(buf) == 7:
                    rx, ry = buf[0], buf[1]
                    for px, py in (cur, (buf[5], buf[6])):
                        xs += [px - rx, px + rx]; ys += [py - ry, py + ry]
                    big = buf[3] == 1
                    cur = [buf[5], buf[6]]; xs.append(cur[0]); ys.append(cur[1]); buf = []
                    if big:   # large arc: extend further
                        ys.append(cur[1] + 2 * ry)
                elif cmd in "MLCSQT" and len(buf) == 2:
                    cur = buf[:]; xs.append(cur[0]); ys.append(cur[1]); buf = []
            if not xs:
                return None
            x0, y0, x1, y1 = min(xs), min(ys), max(xs), max(ys)
        else:
            return None
    except (KeyError, ValueError, IndexError):
        return None
    if op.stroke:
        h = op.stroke / 2
        x0, y0, x1, y1 = x0 - h, y0 - h, x1 + h, y1 + h
    return x0, y0, x1, y1


def inside(b, outer, tol=2.5):
    return b and outer and b[0] >= outer[0] - tol and b[1] >= outer[1] - tol and \
        b[2] <= outer[2] + tol and b[3] <= outer[3] + tol


def union(bs):
    bs = [b for b in bs if b]
    if not bs:
        return None
    return min(b[0] for b in bs), min(b[1] for b in bs), max(b[2] for b in bs), max(b[3] for b in bs)


def with_transform(op: Op, tf: str) -> Op:
    o = copy.deepcopy(op)
    prev = o.attrs.get("transform")
    o.attrs["transform"] = f"{tf} {prev}" if prev else tf
    return o


def clone(icon: Icon, suffix: str, label: str) -> Icon:
    ic = Icon(f"{icon.name}", f"{icon.label} ({label})", icon.category, icon.tags)
    ic.ops = copy.deepcopy(icon.ops)
    ic.wrap = icon.wrap
    ic.anchors = copy.deepcopy(icon.anchors)
    ic.anchor_tf = icon.anchor_tf
    ic.defs = list(icon.defs)
    ic.ink = getattr(icon, "ink", False)
    return ic


# --------------------------------------------------------------------------
# classification
# --------------------------------------------------------------------------

def liquid_ops(icon: Icon):
    """s-role filled shapes that sit inside a glass vessel."""
    glass = [bbox(o) for o in icon.ops if o.glass]
    if not glass:
        return []
    res = []
    for i, o in enumerate(icon.ops):
        if o.role == "s" and o.stroke is None:
            b = bbox(o)
            if b and any(inside(b, g) for g in glass):
                res.append(i)
    return res


AXES_D = "M9 7 V55 H58"


def kind(icon: Icon) -> str:
    if icon.category in ("annotations", "symbols"):
        return "object"
    if icon.category == "data" or any(t in icon.name for t in ("pattern", "spectrum", "plot", "curve")):
        return "data"
    if liquid_ops(icon):
        return "vessel"
    return "object"


# --------------------------------------------------------------------------
# variant builders
# --------------------------------------------------------------------------

def v_empty(icon):
    ic = clone(icon, "empty", "empty")
    liq = set(liquid_ops(icon))
    glass = [bbox(o) for o in icon.ops if o.glass]
    keep = []
    for i, o in enumerate(ic.ops):
        if i in liq:
            continue
        if o.role == "a" and o.tag == "circle" and any(inside(bbox(o), g) for g in glass):
            continue                      # bubbles in the liquid
        keep.append(o)
    ic.ops = keep
    return ic


def v_suspension(icon):
    ic = clone(icon, "suspension", "with particles")
    liq = liquid_ops(icon)
    extra = []
    for i in liq:
        x0, y0, x1, y1 = bbox(icon.ops[i])
        w, h = x1 - x0, y1 - y0
        if w < 4 or h < 4:
            continue
        if icon.ops[i].tag == "path" and "A" in icon.ops[i].attrs.get("d", ""):
            y0 = y0 + h * 0.55          # arc bbox overestimates upward: use the lower band
            h = y1 - y0
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2 + h * 0.08
        r = max(0.9, min(w, h) * 0.075)
        pts = [(-0.22, -0.12), (0.18, -0.2), (0.02, 0.1), (-0.25, 0.2), (0.24, 0.16), (-0.05, -0.28), (0.3, -0.02)]
        n = 7 if w * h > 120 else 4
        for px, py in pts[:n]:
            extra.append(Op("a", "circle", {"cx": round(cx + px * w, 2), "cy": round(cy + py * h * 0.9, 2), "r": round(r, 2)}))
    # insert particles right after the last liquid op, so glass walls stay on top
    pos = max(liq) + 1
    ic.ops = ic.ops[:pos] + extra + ic.ops[pos:]
    return ic


def _trace_ops(icon):
    return [i for i, o in enumerate(icon.ops) if not (o.role == "m" and o.stroke)]


def v_altdata(icon):
    ic = clone(icon, "alt", "alternate data")
    for i in _trace_ops(icon):
        ic.ops[i] = with_transform(ic.ops[i], "translate(67 0) scale(-1 1)")
    return ic


def v_overlay(icon):
    ic = clone(icon, "overlay", "overlay series")
    extra = []
    for o in icon.ops:
        if o.stroke and o.role == "a":
            e = with_transform(o, "translate(3 7)")
            e.role = "s"
            e.attrs["stroke-dasharray"] = "3 2.5"
            e.stroke = max(1.8, o.stroke * 0.75)
            extra.append(e)
        elif o.stroke is None and o.role in ("a", "s") and o.tag == "rect" and o.attrs["height"] > 3:
            h = o.attrs["height"] * 0.62
            e = copy.deepcopy(o)
            e.role = "a" if o.role == "s" else "s"
            e.attrs = dict(o.attrs, x=o.attrs["x"] + o.attrs["width"] * 0.5, width=o.attrs["width"] * 0.62,
                           y=o.attrs["y"] + o.attrs["height"] - h, height=h)
            extra.append(("after", e))
    if not extra:  # dots only: add a fitted line through the plot
        extra.append(Op("s", "path", {"d": "M12 50 C24 40 36 26 56 14", "stroke-dasharray": "3 2.5"}, 2))
    after = [e for t, e in (x for x in extra if isinstance(x, tuple))]
    before = [e for e in extra if not isinstance(e, tuple)]
    first_trace = min(_trace_ops(icon)) if _trace_ops(icon) else len(icon.ops)
    ic.ops = ic.ops[:first_trace] + before + ic.ops[first_trace:] + after
    return ic


def v_mirror(icon):
    ic = clone(icon, "mirrored", "mirrored")
    ic.wrap = "translate(64 0) scale(-1 1)"
    return ic


def v_tilted(icon):
    ic = clone(icon, "tilted", "tilted")
    ic.wrap = "rotate(-14 32 32) translate(3.2 3.2) scale(0.9)"
    return ic


def v_group(icon):
    ic = clone(icon, "group", "set of three")
    ops = []
    for tf in ("translate(1 2) scale(0.56)", "translate(27 2) scale(0.56)", "translate(13 18) scale(0.66)"):
        ops += [with_transform(o, tf) for o in icon.ops]
    ic.ops = ops
    ic.anchor_tf = "translate(13 18) scale(0.66)"     # connect to the front copy
    return ic


SYM_CACHE: dict = {}


def is_symmetric(icon) -> bool:
    import hashlib
    from .core import render_svg, Style
    svg = render_svg(icon, Style(mode="solid", color="#000000"), 48)
    key = "r2:" + hashlib.sha1(svg.encode()).hexdigest()
    if key not in SYM_CACHE:
        SYM_CACHE[key] = _is_symmetric(icon)
    return SYM_CACHE[key]


def _is_symmetric(icon) -> bool:
    import io
    from PIL import Image, ImageChops, ImageOps
    from .core import render_svg, Style
    svg = render_svg(icon, Style(mode="solid", color="#000000"), 48)
    png = svg_to_png(svg, 48, 48, "white")
    im = Image.open(io.BytesIO(png)).convert("L")
    diff = ImageChops.difference(im, ImageOps.mirror(im))
    ink = sum(255 - v for v in im.getdata()) or 1
    return sum(diff.getdata()) / ink < 0.16          # relative to the ink, so thin icons (arrows) count too


def v_3d(icon):
    ic = clone(icon, "3d", "3D")
    depth = []
    for step in (2.4, 1.6, 0.8):
        for o in icon.ops:
            if o.role == "x":
                continue
            d = with_transform(o, f"translate({step} {step})")
            d.depth, d.glass = True, False
            depth.append(d)
    ic.ops = depth + ic.ops
    ic.wrap = "translate(-0.5 -0.5) scale(0.96)"
    return ic


# ==========================================================================
# v4.2: five variants per icon
# ==========================================================================
import functools
import json
import math

N_VARIANTS = 5                      # per icon, on top of the base drawing
SEP = "--"                          # variant id separator: <icon>--<key>

LABELS = {
    "base": "standard", "mirrored": "mirrored", "tilted": "tilted", "group": "set of three", "3d": "3D",
    "highlighted": "highlighted", "empty": "empty", "suspension": "with particles", "biphasic": "two phases",
    "crystals": "with crystals", "stirred": "stirred", "alt-data": "alternate data", "overlay": "overlay series",
    "stacked": "stacked / waterfall", "tile": "on a tile", "signal-word": "with signal word",
    "on-bottle": "on a bottle label", "annotated": "peak marked", "inset": "zoom inset",
    "heated": "heated", "cooled": "cooled", "irradiated": "irradiated (light)", "timed": "timed / ageing",
    "done": "done / OK", "failed": "failed / not OK", "electric": "under bias",
}


def variant_id(name: str, key: str) -> str:
    return name if key == "base" else f"{name}{SEP}{key}"


def split_id(vid: str) -> tuple[str, str]:
    """'beaker--crystals' -> ('beaker', 'crystals'); 'beaker' -> ('beaker', 'base')."""
    if SEP in vid:
        n, k = vid.split(SEP, 1)
        return n, k
    return vid, "base"


def _no_mirror(icon) -> bool:
    return icon.name in NO_MIRROR or icon.name.startswith(("ghs_", "element_")) or \
        any(o.tag == "text" for o in icon.ops)


# ----------------------------------------------------------------- accurate liquid box (pixels, cached)
_PXCACHE: dict = {}


def pixel_bbox(ops) -> tuple | None:
    """True ink bbox of some ops in 64-grid units (arc bboxes from the parser are too loose)."""
    from .core import Icon as _I, render_svg, Style
    key = json.dumps([(o.tag, o.attrs, o.stroke) for o in ops], sort_keys=True, default=str)
    if key in _PXCACHE:
        return _PXCACHE[key]
    import io
    from PIL import Image
    ic = _I("px")
    ic.ops = [Op("m", o.tag, {k: v for k, v in o.attrs.items() if k != "clip-path"}, o.stroke) for o in ops]
    svg = render_svg(ic, Style(mode="solid", color="#000000"), 256)
    png = svg_to_png(svg, 256, 256, "white")
    im = Image.open(io.BytesIO(png)).convert("L").point(lambda v: 255 if v < 128 else 0)
    b = im.getbbox()
    r = None if not b else tuple(round(v / 4, 2) for v in b)
    _PXCACHE[key] = r
    return r


def liquid_box(icon):
    idx = liquid_ops(icon)
    return pixel_bbox([icon.ops[i] for i in idx]) if idx else None


# ----------------------------------------------------------------- vessel variants
def v_biphasic(icon):
    idx = liquid_ops(icon)
    lb = liquid_box(icon)
    if not idx or not lb or lb[3] - lb[1] < 6:
        return None
    ic = clone(icon, "biphasic", LABELS["biphasic"])
    y_mid = lb[1] + (lb[3] - lb[1]) * 0.48
    ic.defs.append(f'<clipPath id="__PID__-phase"><rect x="-10" y="{y_mid:.2f}" width="84" height="80"/></clipPath>')
    lower = []
    for i in idx:
        o = copy.deepcopy(icon.ops[i])
        o.role = "a"
        o.attrs["clip-path"] = "url(#__PID__-phase)"
        lower.append(o)
    # little droplets of the lower phase still settling in the upper phase
    cx = (lb[0] + lb[2]) / 2
    w = lb[2] - lb[0]
    drops = [Op("a", "circle", {"cx": round(cx + dx * w, 2), "cy": round(y_mid - dy, 2), "r": r})
             for dx, dy, r in ((-0.18, 2.6, 0.9), (0.14, 4.2, 0.75))] if w > 10 else []
    pos = max(idx) + 1
    ic.ops = ic.ops[:pos] + lower + drops + ic.ops[pos:]
    return ic


def _crystal(cx, cy, s, rot):
    pts = []
    for k, (a, rr) in enumerate(((90, 1.0), (30, 0.62), (-30, 0.62), (-90, 1.0), (-150, 0.62), (150, 0.62))):
        x, y = cx + s * rr * math.cos(math.radians(a + rot)) * (0.75 if k % 3 else 1), \
            cy - s * rr * math.sin(math.radians(a + rot))
        pts.append(f"{x:.2f},{y:.2f}")
    return Op("a", "polygon", {"points": " ".join(pts)})


def _single_liquid(icon) -> bool:
    idx = liquid_ops(icon)
    bs = [pixel_bbox([icon.ops[i]]) for i in idx]
    bs = [b for b in bs if b]
    return all(not (a[2] < b[0] or b[2] < a[0]) for a in bs for b in bs)   # all liquid boxes overlap in x


def v_crystals(icon):
    if not _single_liquid(icon):
        return None
    idx = liquid_ops(icon)
    lb = liquid_box(icon)
    if not idx or not lb:
        return None
    w, h = lb[2] - lb[0], lb[3] - lb[1]
    if w < 7 or h < 6:
        return None
    ic = clone(icon, "crystals", LABELS["crystals"])
    s = max(1.8, min(4.2, w * 0.13, h * 0.32))
    cx, yb = (lb[0] + lb[2]) / 2, lb[3] - s * 0.9 - 0.8
    spots = [(-0.2, 0, 12), (0.04, -0.15, -18), (0.22, 0.05, 30), (-0.04, 0.12, 64)] if w > 14 else [(-0.12, 0, 10), (0.14, 0, -24)]
    crys = [_crystal(cx + dx * w, yb + dy * s, s * (0.85 if k % 2 else 1), r) for k, (dx, dy, r) in enumerate(spots)]
    pos = max(idx) + 1
    ic.ops = ic.ops[:pos] + crys + ic.ops[pos:]
    return ic


def v_stirred(icon):
    if not _single_liquid(icon):
        return None
    idx = liquid_ops(icon)
    lb = liquid_box(icon)
    if not idx or not lb:
        return None
    w, h = lb[2] - lb[0], lb[3] - lb[1]
    if w < 9 or h < 7:
        return None
    ic = clone(icon, "stirred", LABELS["stirred"])
    bw = min(12.0, w * 0.42)
    cx, yb = (lb[0] + lb[2]) / 2, lb[3] - 2.6
    bar = Op("x", "rect", {"x": round(cx - bw / 2, 2), "y": round(yb - 1.7, 2), "width": round(bw, 2), "height": 3.2, "rx": 1.6})
    # rotation hint: a curved white arrow in the liquid
    r = min(bw * 0.55, h * 0.32)
    ay = yb - 2.5 - r * 0.55
    arc = Op("x", "path", {"d": f"M{cx - r:.2f} {ay:.2f} A{r:.2f} {r * 0.42:.2f} 0 0 0 {cx + r:.2f} {ay:.2f}"}, 1.3)
    head = Op("x", "polygon", {"points": f"{cx + r + 1.6:.2f},{ay - 0.6:.2f} {cx + r - 1.4:.2f},{ay - 1.6:.2f} {cx + r - 0.2:.2f},{ay + 1.6:.2f}"})
    extra = [bar] + ([arc, head] if h > 11 else [])
    pos = max(idx) + 1
    ic.ops = ic.ops[:pos] + extra + ic.ops[pos:]
    return ic


# ----------------------------------------------------------------- data variants
def path_points(op) -> list[tuple[float, float]]:
    """End points of an absolute path / polyline / line (control points skipped)."""
    a = op.attrs
    if op.tag == "line":
        return [(a["x1"], a["y1"]), (a["x2"], a["y2"])]
    if op.tag in ("polyline", "polygon"):
        v = [float(n) for n in NUM.findall(a["points"])]
        return list(zip(v[0::2], v[1::2]))
    if op.tag == "circle":
        return [(a["cx"], a["cy"])]
    if op.tag != "path":
        return []
    d = a["d"]
    if re.search(r"[a-z]", d.replace("e-", "")):
        return []
    toks = re.findall(r"[A-Za-z]|-?\d*\.?\d+", d)
    need = {"M": 2, "L": 2, "T": 2, "C": 6, "S": 4, "Q": 4, "A": 7, "H": 1, "V": 1}
    pts, cmd, buf, cur = [], "M", [], [0.0, 0.0]
    for t in toks:
        if t.isalpha():
            cmd, buf = t.upper(), []
            continue
        buf.append(float(t))
        if cmd in need and len(buf) == need[cmd]:
            if cmd == "H":
                cur = [buf[0], cur[1]]
            elif cmd == "V":
                cur = [cur[0], buf[0]]
            else:
                cur = buf[-2:]
            pts.append(tuple(cur))
            buf = []
    return pts


def _traces(icon):
    return [o for i, o in enumerate(icon.ops) if i in set(_trace_ops(icon))]


def _scaled(o, tf, k):
    e = with_transform(o, tf)
    if e.stroke:
        e.stroke = round(e.stroke / k * 0.85, 2)
    return e


def v_stacked(icon):
    tr = _trace_ops(icon)
    if not tr:
        return None
    ic = clone(icon, "stacked", LABELS["stacked"])
    k = 0.72
    base_tf = f"translate(9 55) scale({k}) translate(-9 -55)"
    layers = []
    for n, (dx, dy) in enumerate(((10, -17), (5, -8.5), (0, 0))):
        for i in tr:
            o = _scaled(icon.ops[i], f"translate({dx} {dy}) {base_tf}", k)
            if n < 2 and o.role in ("a", "s"):
                o.role = "s" if n == 0 else ("a" if o.role == "s" else "s")
            if n == 0:
                o.attrs["opacity"] = 0.75
            layers.append(o)
    axes_ = [o for i, o in enumerate(icon.ops) if i not in set(tr)]
    ic.ops = axes_ + layers
    return ic


def _peak(icon):
    """Highest local maximum of the main trace (the stroked series with most points)."""
    cands = [o for o in _traces(icon) if o.role != "x" and o.tag in ("path", "polyline")]
    pts = max((path_points(o) for o in cands), key=len, default=[])
    if len(pts) < 3:
        return None
    best = None
    for i in range(1, len(pts) - 1):
        x, y = pts[i]
        if 13 <= x <= 55 and y <= pts[i - 1][1] and y <= pts[i + 1][1] and (best is None or y < best[1]):
            best = (x, y)
    return best


def v_annotated(icon):
    p = _peak(icon)
    if not p or p[1] > 46:
        return None
    ic = clone(icon, "annotated", LABELS["annotated"])
    x, y = p
    ic.ops.insert(min(_trace_ops(icon)), Op("s", "line", {"x1": x, "y1": round(y + 2, 2), "x2": x, "y2": 54,
                                                          "stroke-dasharray": "2.4 2"}, 1.5))
    ty = max(3.5, y - 4.2)
    ic.ops.append(Op("a", "polygon", {"points": f"{x - 3:.2f},{ty - 3.4:.2f} {x + 3:.2f},{ty - 3.4:.2f} {x:.2f},{ty + 0.6:.2f}"}))
    return ic


def v_inset(icon):
    p = _peak(icon)
    tr = _trace_ops(icon)
    if not p or not tr:
        return None
    ic = clone(icon, "inset", LABELS["inset"])
    px, py = p
    zw, zh = 12.0, 9.0
    zx, zy = min(max(px - zw / 2, 10), 56 - zw), min(max(py - 2, 6), 53 - zh)
    W, H = 22.0, 16.5
    pts = [q for o in _traces(icon) for q in path_points(o)]
    dense = pts + [((a[0] + b[0]) / 2, (a[1] + b[1]) / 2) for a, b in zip(pts, pts[1:])]

    def overlap(r):
        x0, y0 = r
        return sum(1 for qx, qy in dense if x0 - 1 <= qx <= x0 + W + 1 and y0 - 1 <= qy <= y0 + H + 1) + \
            (50 if x0 - 1 <= zx + zw and zx <= x0 + W + 1 and y0 - 1 <= zy + zh and zy <= y0 + H + 1 else 0)
    ix, iy = min(((12.5, 7.0), (35.0, 7.0), (12.5, 36.5), (35.0, 36.5), (23.0, 22.0)), key=overlap)
    s = W / zw
    ic.defs.append(f'<clipPath id="__PID__-inset"><rect x="{zx:.2f}" y="{zy:.2f}" width="{zw}" height="{zh}"/></clipPath>')
    tf = f"translate({ix - zx * s:.2f} {iy - zy * s:.2f}) scale({s:.3f})"
    zoomed = []
    for i in tr:
        o = _scaled(icon.ops[i], tf, s)
        if o.stroke:
            o.stroke = round(max(1.6, icon.ops[i].stroke * 0.8) / s, 2)
        o.attrs["clip-path"] = "url(#__PID__-inset)"   # clip lives in the element's own (pre-zoom) space
        zoomed.append(o)
    ic.ops += [Op("s", "rect", {"x": zx, "y": zy, "width": zw, "height": zh, "rx": 1, "stroke-dasharray": "2 1.6"}, 1.2),
               Op("s", "line", {"x1": zx + (zw if ix > zx else 0), "y1": zy, "x2": ix + (0 if ix > zx else W), "y2": iy + (H if iy > zy else 0),
                                "stroke-dasharray": "1.6 1.6"}, 1.0),
               Op("x", "rect", {"x": ix, "y": iy, "width": W, "height": H, "rx": 1.5})] + zoomed + \
              [Op("m", "rect", {"x": ix, "y": iy, "width": W, "height": H, "rx": 1.5}, 1.5)]
    return ic


# ----------------------------------------------------------------- generic object variants
def v_highlighted(icon):
    ic = clone(icon, "highlighted", LABELS["highlighted"])
    tf = "translate(3.84 3.84) scale(0.88)"
    ic.ops = [Op("a", "circle", {"cx": 32, "cy": 32, "r": 30, "opacity": 0.16})] + \
        [with_transform(o, tf) for o in icon.ops] + [_sparkle(53.5, 9.5, 5.2), _sparkle(60, 19.5, 2.8), _sparkle(8, 55, 3)]
    ic.anchor_tf = f"{tf} {icon.anchor_tf}" if icon.anchor_tf else tf
    return ic


def _sparkle(x, y, r):
    q = r * 0.26
    pts = [(x, y - r), (x + q, y - q), (x + r, y), (x + q, y + q), (x, y + r), (x - q, y + q), (x - r, y), (x - q, y - q)]
    return Op("a", "polygon", {"points": " ".join(f"{a:.2f},{b:.2f}" for a, b in pts)})


def v_tile(icon):
    """App-tile style: the drawing on a soft rounded square (good for symbols / legends)."""
    ic = clone(icon, "tile", "on a tile")
    tf = "translate(8.32 8.32) scale(0.74)"
    ic.ops = [Op("s", "rect", {"x": 3, "y": 3, "width": 58, "height": 58, "rx": 13, "opacity": 0.22}),
              Op("s", "rect", {"x": 3, "y": 3, "width": 58, "height": 58, "rx": 13, "opacity": 0.55}, 1.6)] + \
        [with_transform(o, tf) for o in icon.ops]
    ic.anchor_tf = f"{tf} {icon.anchor_tf}" if icon.anchor_tf else tf
    return ic


GHS_WORD = {"ghs_harmful": "WARNING", "ghs_environment": "WARNING", "ghs_gas": "WARNING"}   # others: DANGER (typical)


def v_ghs_word(icon):
    ic = clone(icon, "signal-word", "with signal word")
    tf = "translate(10.9 0.6) scale(0.66)"
    ic.ops = [with_transform(o, tf) for o in icon.ops]
    ic.ops.append(Op("m", "text", {"x": 32, "y": 58.5, "font-size": 10.5, "font-weight": "bold", "text-anchor": "middle",
                                   "font-family": "Helvetica, Arial, sans-serif", "_text": GHS_WORD.get(icon.name, "DANGER"),
                                   "_fill": "#111111"}))
    return ic


def v_ghs_bottle(icon):
    ic = clone(icon, "on-bottle", "on a bottle label")
    ic.ops = [Op("m", "rect", {"x": 24, "y": 3, "width": 16, "height": 8, "rx": 2}),
              Op("m", "rect", {"x": 26, "y": 10, "width": 12, "height": 7}),
              Op("m", "rect", {"x": 13, "y": 16, "width": 38, "height": 45, "rx": 7}, 3.2, True),
              Op("x", "rect", {"x": 17, "y": 25, "width": 30, "height": 30, "rx": 2.5})]
    tf = "translate(17.9 25.9) scale(0.44)"
    ic.ops += [with_transform(o, tf) for o in icon.ops]
    ic.anchors = []
    return ic


BADGE_COL = {"heated": "#F0643C", "cooled": "#3A8DDE", "irradiated": "#F2B705", "timed": "#64748B",
             "done": "#2FA36B", "failed": "#E5484D", "electric": "#7C5CDB"}
BX, BY = 51.5, 51.5


def _glyph(key, cx=BX, cy=BY):
    """White symbol on the badge (role x: white in colour styles, a cut-out in solid)."""
    o = []
    if key == "heated":
        o.append(Op("x", "path", {"d": f"M{cx:.2f} {cy - 6.2:.2f} C{cx + 4.6:.2f} {cy - 2.4:.2f} {cx + 5:.2f} {cy + 1.2:.2f} "
                                        f"{cx + 3.6:.2f} {cy + 3.8:.2f} A4.6 4.6 0 0 1 {cx - 3.8:.2f} {cy + 3.6:.2f} "
                                        f"C{cx - 4.9:.2f} {cy + 0.9:.2f} {cx - 3.4:.2f} {cy - 1.3:.2f} {cx - 2.2:.2f} {cy - 2.6:.2f} "
                                        f"C{cx - 1.9:.2f} {cy - 0.6:.2f} {cx - 1.1:.2f} {cy + 0.4:.2f} {cx - 0.2:.2f} {cy + 0.6:.2f} "
                                        f"C{cx - 0.9:.2f} {cy - 1.8:.2f} {cx - 0.6:.2f} {cy - 4:.2f} {cx:.2f} {cy - 6.2:.2f} Z"}))
    elif key == "cooled":
        for a in (90, 30, 150):
            dx, dy = 6 * math.cos(math.radians(a)), 6 * math.sin(math.radians(a))
            o.append(Op("x", "line", {"x1": round(cx - dx, 2), "y1": round(cy - dy, 2), "x2": round(cx + dx, 2), "y2": round(cy + dy, 2)}, 1.6))
        for a in range(0, 360, 60):
            tx, ty = cx + 4 * math.cos(math.radians(a + 90)), cy + 4 * math.sin(math.radians(a + 90))
            for s in (-1, 1):
                b = math.radians(a + 90 + s * 140)
                o.append(Op("x", "line", {"x1": round(tx, 2), "y1": round(ty, 2), "x2": round(tx + 1.9 * math.cos(b), 2),
                                          "y2": round(ty + 1.9 * math.sin(b), 2)}, 1.2))
    elif key == "irradiated":
        o.append(Op("x", "circle", {"cx": cx, "cy": cy, "r": 2.6}))
        for a in range(0, 360, 45):
            c, s_ = math.cos(math.radians(a)), math.sin(math.radians(a))
            o.append(Op("x", "line", {"x1": round(cx + 4.2 * c, 2), "y1": round(cy + 4.2 * s_, 2),
                                      "x2": round(cx + 6.4 * c, 2), "y2": round(cy + 6.4 * s_, 2)}, 1.4))
    elif key == "timed":
        o.append(Op("x", "path", {"d": f"M{cx - 5.4:.2f} {cy:.2f} A5.4 5.4 0 1 0 {cx + 5.4:.2f} {cy:.2f} A5.4 5.4 0 1 0 {cx - 5.4:.2f} {cy:.2f}"}, 1.5))
        o.append(Op("x", "path", {"d": f"M{cx:.2f} {cy - 3.4:.2f} V{cy:.2f} H{cx + 2.8:.2f}"}, 1.5))
    elif key == "done":
        o.append(Op("x", "path", {"d": f"M{cx - 4.6:.2f} {cy + 0.2:.2f} L{cx - 1.4:.2f} {cy + 3.4:.2f} L{cx + 4.8:.2f} {cy - 3.2:.2f}"}, 2.3))
    elif key == "failed":
        for s in (-1, 1):
            o.append(Op("x", "line", {"x1": round(cx - 3.8, 2), "y1": round(cy - 3.8 * s, 2), "x2": round(cx + 3.8, 2), "y2": round(cy + 3.8 * s, 2)}, 2.3))
    elif key == "electric":
        o.append(Op("x", "polygon", {"points": f"{cx + 1.2:.2f},{cy - 6.6:.2f} {cx - 4:.2f},{cy + 0.8:.2f} {cx - 0.4:.2f},{cy + 0.8:.2f} "
                                               f"{cx - 1.4:.2f},{cy + 6.6:.2f} {cx + 4:.2f},{cy - 0.9:.2f} {cx + 0.4:.2f},{cy - 0.9:.2f}"}))
    return o


def v_badge(icon, key):
    ic = clone(icon, key, LABELS[key])
    tf = "scale(0.84)"
    ic.ops = [with_transform(o, tf) for o in icon.ops]
    ic.ops.append(Op("x", "circle", {"cx": BX, "cy": BY, "r": 11.8}))
    ic.ops.append(Op("a", "circle", {"cx": BX, "cy": BY, "r": 9.8, "_fill": BADGE_COL[key]}))
    ic.ops += _glyph(key)
    ic.anchor_tf = f"{tf} {icon.anchor_tf}" if icon.anchor_tf else tf
    return ic


CAT_BADGES = {
    "glassware": ("heated", "cooled"), "synthesis": ("heated", "timed"), "characterisation": ("done", "timed"),
    "instruments": ("done", "timed"), "materials": ("irradiated", "heated"), "proteomics": ("timed", "done"),
    "biology": ("done", "failed"), "energy": ("electric", "heated"), "safety": ("done", "failed"),
    "infrastructure": ("done", "failed"), "computation": ("done", "timed"), "elements": ("irradiated", "electric"),
    "data": ("done", "timed"),
}
VESSEL_CATS = {"glassware", "synthesis", "energy", "proteomics", "biology", "materials"}
NOT_VESSEL = {"liebig_condenser", "glass_slide", "thermometer", "hourglass", "reference_electrode", "vacuum_desiccator",
              "cvd_reactor", "syringe_pump", "kugelrohr", "electrospinning", "tube_rack"}


def has_axes(icon) -> bool:
    return any(o.tag == "path" and o.role == "m" and o.stroke and o.attrs.get("d", "").startswith(AXES_D)
               for o in icon.ops)


def vessel_like(icon) -> bool:
    if icon.name in NOT_VESSEL or icon.category not in VESSEL_CATS:
        return False
    lb = liquid_box(icon)
    return bool(lb) and (lb[2] - lb[0]) * (lb[3] - lb[1]) >= 50
GROUPABLE = {"glassware", "materials", "biology", "proteomics", "symbols", "computation", "elements"}

BUILDERS = {
    "tile": v_tile, "signal-word": v_ghs_word, "on-bottle": v_ghs_bottle,
    "empty": v_empty, "suspension": v_suspension, "biphasic": v_biphasic, "crystals": v_crystals, "stirred": v_stirred,
    "alt-data": v_altdata, "overlay": v_overlay, "stacked": v_stacked, "annotated": v_annotated, "inset": v_inset,
    "mirrored": v_mirror, "tilted": v_tilted, "group": v_group, "3d": v_3d, "highlighted": v_highlighted,
    **{k: functools.partial(v_badge, key=k) for k in BADGE_COL},
}


def recipe(icon) -> list[str]:
    """Ordered candidate keys; the first N_VARIANTS that apply are used."""
    k = kind(icon)
    badges = list(CAT_BADGES.get(icon.category, ()))
    if icon.name.startswith("ghs_"):
        return ["signal-word", "on-bottle", "group", "3d", "highlighted"]
    if k == "vessel" and vessel_like(icon):
        return ["empty", "suspension", "biphasic", "crystals", "stirred", "highlighted", "3d", "heated"]
    if k == "data":
        if has_axes(icon):
            return ["alt-data", "overlay", "stacked", "annotated", "inset", "highlighted", "3d"]
        return ["alt-data", "overlay", "highlighted", "3d", "tile", "group"]
    if has_axes(icon):                    # a plot drawn inside an instrument icon (TGA, FTIR ...): never mirror
        return ["alt-data", "overlay", "3d", "highlighted"] + badges
    if k == "vessel" and icon.category not in VESSEL_CATS and icon.name not in NOT_VESSEL:
        return ["empty"] + (["mirrored"] if not _no_mirror(icon) and not is_symmetric(icon) else []) + \
            ["3d", "highlighted"] + badges
    if icon.category in ("annotations", "symbols"):
        return (["tilted"] if _no_mirror(icon) else ["mirrored"]) + ["group", "3d", "highlighted", "tile", "tilted"]
    keys = []
    if _no_mirror(icon):
        keys.append("tilted")
    elif not is_symmetric(icon):
        keys.append("mirrored")
    if icon.category in GROUPABLE or is_symmetric(icon):
        keys.append("group")
    keys += ["3d", "highlighted"] + badges + ["tile", "tilted", "group"]
    return keys


OVERRIDES: dict[str, tuple] = {}


def shape_variants(icon: Icon) -> list[tuple[str, Icon]]:
    """Return [(variant_key, Icon)]: the base drawing first, then exactly N_VARIANTS variants."""
    out = [("base", icon)]
    seen = {"base"}
    for key, fn in OVERRIDES.get(icon.name, ()):
        if len(out) > N_VARIANTS:
            break
        out.append((key, fn(icon)))
        seen.add(key)
    fill = ["highlighted", "3d", "tile", "group"] + ([] if icon.name.startswith("ghs_") else ["tilted"] + list(BADGE_COL))
    for key in recipe(icon) + fill:
        if len(out) > N_VARIANTS:
            break
        if key in seen:
            continue
        if key == "mirrored" and (_no_mirror(icon) or is_symmetric(icon)):
            continue
        v = BUILDERS[key](icon)
        if v is not None:
            out.append((key, v))
            seen.add(key)
    return out


def family(icon: Icon) -> list[dict]:
    """[{id, key, label}] for the picker: every member of the icon's variant family."""
    return [{"id": variant_id(icon.name, k), "key": k, "label": LABELS.get(k, k.replace("-", " "))}
            for k, _ in shape_variants(icon)]
