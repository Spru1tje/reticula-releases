"""
Reticon core: a tiny vector-icon DSL + SVG renderer.

Every icon is drawn on a 64x64 grid from a short list of primitives. Each
primitive has a *role*:

    m  main       (body of the instrument / glass outline)
    s  secondary  (liquid, sample, linker ...)
    a  accent     (beams, highlights, data)
    x  cut        (knock-out: removes everything drawn BEFORE it)

Because roles are abstract, the same icon renders as:

    solid    one colour, cuts are real transparency
    duotone  one colour, secondary at reduced opacity
    color    full palette (main / secondary / accent)

Cuts are order-aware: a cut only removes what was painted before it, so you
can do "halo" effects (cut a gap, then draw on top of it).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable

GRID = 64

# --------------------------------------------------------------------------
# Style
# --------------------------------------------------------------------------

@dataclass
class Palette:
    main: str = "#1E293B"
    secondary: str = "#14B8A6"
    accent: str = "#F59E0B"


@dataclass
class Style:
    mode: str = "solid"            # solid | duotone | color
    color: str = "#1E293B"         # used by solid / duotone
    palette: Palette = field(default_factory=Palette)
    duotone_opacity: float = 0.45
    knockout: str | None = None    # e.g. "#FFFFFF" -> paint cuts instead of masking
                                   # (use for renderers without <mask>, e.g. QtSvg)

    def fill_for(self, role: str) -> tuple[str, float]:
        if self.mode == "color":
            return {"m": self.palette.main, "s": self.palette.secondary,
                    "a": self.palette.accent}[role], 1.0
        if self.mode == "duotone" and role == "s":
            return self.color, self.duotone_opacity
        return self.color, 1.0


# --------------------------------------------------------------------------
# Icon builder
# --------------------------------------------------------------------------

def _f(v: float) -> str:
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


@dataclass
class Op:
    role: str
    tag: str
    attrs: dict
    stroke: float | None = None    # None -> filled shape, else stroke width
    glass: bool = False            # glass wall: illustration mode fills it translucent
    depth: bool = False            # extrusion copy (3D variant): painted as a darker shade


class Icon:
    def __init__(self, name: str, label: str = "", category: str = "misc",
                 tags: Iterable[str] = ()):
        self.name = name
        self.label = label or name.replace("_", " ").title()
        self.category = category
        self.tags = list(tags)
        self.ops: list[Op] = []
        self.wrap: str | None = None   # optional transform applied to the whole drawing
        self.anchors: list[dict] = []  # connection points: {"id","x","y","dir"} in 64-grid coords
        self.anchor_tf: str | None = None  # transform for anchors when it differs from wrap (group variant)
        self.defs: list[str] = []      # extra <defs> content (clip paths for fill level ...)

    def text(self, x, y, txt, size=14, weight="bold", anchor="middle", role="m", fill=None):
        a = {"x": x, "y": y, "font-size": size, "font-weight": weight, "text-anchor": anchor,
             "font-family": "Helvetica, Arial, sans-serif", "_text": txt}
        if fill:
            a["_fill"] = fill
        self.ops.append(Op(role, "text", a))
        return self

    def fixed(self, colour: str):
        """Give the last added shape a fixed colour (GHS red, emission colours …) in colour styles."""
        self.ops[-1].attrs["_fill"] = colour
        return self

    def anchor(self, aid: str, x: float, y: float, direction: str = "out"):
        """Connection point for arrows / tubing. direction: up | down | left | right | out."""
        self.anchors.append({"id": aid, "x": x, "y": y, "dir": direction})
        return self

    # ---- filled primitives ------------------------------------------------
    def rect(self, x, y, w, h, r=0, role="m"):
        a = {"x": x, "y": y, "width": w, "height": h}
        if r:
            a["rx"] = r
        self.ops.append(Op(role, "rect", a))
        return self

    def circle(self, cx, cy, r, role="m"):
        self.ops.append(Op(role, "circle", {"cx": cx, "cy": cy, "r": r}))
        return self

    def ellipse(self, cx, cy, rx, ry, role="m"):
        self.ops.append(Op(role, "ellipse", {"cx": cx, "cy": cy, "rx": rx, "ry": ry}))
        return self

    def path(self, d, role="m"):
        self.ops.append(Op(role, "path", {"d": d}))
        return self

    def poly(self, pts, role="m"):
        self.ops.append(Op(role, "polygon", {"points": " ".join(f"{_f(x)},{_f(y)}" for x, y in pts)}))
        return self

    # ---- stroked primitives -----------------------------------------------
    def line(self, x1, y1, x2, y2, w=3.5, role="m"):
        self.ops.append(Op(role, "line", {"x1": x1, "y1": y1, "x2": x2, "y2": y2}, w))
        return self

    def spath(self, d, w=3.5, role="m", glass=False):
        self.ops.append(Op(role, "path", {"d": d}, w, glass))
        return self

    def spoly(self, pts, w=3.5, role="m", glass=False):
        self.ops.append(Op(role, "polyline",
                           {"points": " ".join(f"{_f(x)},{_f(y)}" for x, y in pts)}, w, glass))
        return self

    def srect(self, x, y, w, h, r=0, sw=3.5, role="m", glass=False):
        a = {"x": x, "y": y, "width": w, "height": h}
        if r:
            a["rx"] = r
        self.ops.append(Op(role, "rect", a, sw, glass))
        return self

    def scircle(self, cx, cy, r, w=3.5, role="m", glass=False):
        self.ops.append(Op(role, "circle", {"cx": cx, "cy": cy, "r": r}, w, glass))
        return self

    # ---- convenience ------------------------------------------------------
    def cut(self):
        """Fluent helper: icon.cut().circle(...) == icon.circle(..., role='x')."""
        return _RoleProxy(self, "x")

    def sec(self):
        return _RoleProxy(self, "s")

    def acc(self):
        return _RoleProxy(self, "a")

    # ---- serialisation (JSON specs) ---------------------------------------
    def to_spec(self) -> dict:
        d = {"name": self.name, "label": self.label, "category": self.category,
             "tags": self.tags,
             "ops": [{"role": o.role, "tag": o.tag, "attrs": o.attrs,
                      **({"stroke": o.stroke} if o.stroke else {}),
                      **({"glass": True} if o.glass else {}),
                      **({"depth": True} if o.depth else {})} for o in self.ops]}
        if self.wrap:
            d["wrap"] = self.wrap
        if self.anchors:
            d["anchors"] = self.anchors
        if self.anchor_tf:
            d["anchor_tf"] = self.anchor_tf
        if getattr(self, "ink", False):
            d["ink"] = True
        if self.defs:
            d["defs"] = list(self.defs)
        return d

    @classmethod
    def from_spec(cls, spec: dict) -> "Icon":
        ic = cls(spec["name"], spec.get("label", ""), spec.get("category", "misc"),
                 spec.get("tags", []))
        for o in spec["ops"]:
            ic.ops.append(Op(o["role"], o["tag"], dict(o["attrs"]), o.get("stroke"), o.get("glass", False),
                             o.get("depth", False)))
        ic.wrap = spec.get("wrap")
        ic.anchors = [dict(a) for a in spec.get("anchors", [])]
        ic.anchor_tf = spec.get("anchor_tf")
        ic.ink = spec.get("ink", False)
        ic.defs = list(spec.get("defs", []))
        return ic


class _RoleProxy:
    def __init__(self, icon: Icon, role: str):
        self.icon, self.role = icon, role

    def __getattr__(self, item):
        fn = getattr(self.icon, item)

        def call(*a, **k):
            k["role"] = self.role
            return fn(*a, **k)
        return call


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------

def _cls(op: Op) -> str:
    return f'rt-{op.role}' + ("-glass" if op.glass else "") + ("-depth" if op.depth else "") + ("-line" if op.stroke else "")


def _fin(op: Op, s: str) -> str:
    """Text ops carry their content in attrs["_text"]: turn the self-closing tag into <text>…</text>."""
    if op.tag != "text" or not s.endswith("/>"):
        return s
    t = str(op.attrs.get("_text", "")).replace("&", "&amp;").replace("<", "&lt;")
    return s[:-2] + f">{t}</text>"


def _el(op: Op, colour: str, opacity: float = 1.0) -> str:
    own = op.attrs.get("opacity")                      # icon-level opacity multiplies the style opacity
    if own is not None:
        opacity = float(own) * opacity
    attrs = " ".join(f'{k}="{_f(v) if isinstance(v, (int, float)) else v}"'
                     for k, v in op.attrs.items() if not k.startswith("_") and k != "opacity") + f' class="{_cls(op)}"'
    op_attr = f' opacity="{_f(opacity)}"' if opacity < 1 else ""
    if op.stroke:
        return (f'<{op.tag} {attrs} fill="none" stroke="{colour}" stroke-width="{_f(op.stroke)}" '
                f'stroke-linecap="round" stroke-linejoin="round"{op_attr}/>')
    return f'<{op.tag} {attrs} fill="{colour}"{op_attr}/>'


def render_svg(icon: Icon, style: Style | None = None, size: int = 64,
               id_prefix: str | None = None, standalone: bool = True) -> str:
    style = style or Style()
    pid = id_prefix or icon.name
    defs: list[str] = []
    body: list[str] = []
    ops = icon.ops
    i, k = 0, 0
    while i < len(ops):
        op = ops[i]
        if op.role != "x":
            col, opa = style.fill_for(op.role)
            if style.mode == "color" and "_fill" in op.attrs:
                col = op.attrs["_fill"]
            body.append(_fin(op, _el(op, col, 0.4 if op.depth else opa)))
            i += 1
            continue
        # gather consecutive cuts
        cuts = []
        while i < len(ops) and ops[i].role == "x":
            cuts.append(ops[i])
            i += 1
        if style.knockout:
            body.extend(_fin(c, _el(c, style.knockout)) for c in cuts)
        else:
            mid = f"{pid}-k{k}"
            k += 1
            defs.append(
                f'<mask id="{mid}" maskUnits="userSpaceOnUse" x="0" y="0" width="{GRID}" height="{GRID}">'
                f'<rect x="0" y="0" width="{GRID}" height="{GRID}" fill="#fff"/>'
                + "".join(_fin(c, _el(c, "#000")) for c in cuts) + "</mask>")
            body = [f'<g mask="url(#{mid})">' + "".join(body) + "</g>"]

    head = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {GRID} {GRID}" '
            f'width="{size}" height="{size}">')
    title = f"<title>{icon.label}</title>"
    defs = list(icon.defs) + defs
    d = f"<defs>{''.join(defs)}</defs>" if defs else ""
    inner = "".join(body)
    if icon.wrap:
        inner = f'<g transform="{icon.wrap}">{inner}</g>'
    svg = (head + title + d + inner + "</svg>").replace("__PID__", pid)   # per-instance ids in icon defs
    if standalone:
        svg = '<?xml version="1.0" encoding="UTF-8"?>\n' + svg
    return svg


# --------------------------------------------------------------------------
# Geometry helpers for computed icons
# --------------------------------------------------------------------------

def polar(cx, cy, r, deg):
    a = math.radians(deg)
    return cx + r * math.cos(a), cy + r * math.sin(a)


def smooth_path(pts) -> str:
    """Catmull-Rom -> cubic bezier through points."""
    if len(pts) < 3:
        return "M" + " L".join(f"{_f(x)} {_f(y)}" for x, y in pts)
    d = [f"M{_f(pts[0][0])} {_f(pts[0][1])}"]
    for i in range(len(pts) - 1):
        p0 = pts[i - 1] if i else pts[i]
        p1, p2 = pts[i], pts[i + 1]
        p3 = pts[i + 2] if i + 2 < len(pts) else p2
        c1 = (p1[0] + (p2[0] - p0[0]) / 6, p1[1] + (p2[1] - p0[1]) / 6)
        c2 = (p2[0] - (p3[0] - p1[0]) / 6, p2[1] - (p3[1] - p1[1]) / 6)
        d.append(f"C{_f(c1[0])} {_f(c1[1])} {_f(c2[0])} {_f(c2[1])} {_f(p2[0])} {_f(p2[1])}")
    return " ".join(d)


def lorentz_trace(x0, x1, base, peaks, n=90, sign=-1):
    """Spectrum polyline: peaks = [(centre, height, width)] in grid units."""
    pts = []
    for i in range(n + 1):
        x = x0 + (x1 - x0) * i / n
        y = base
        for c, h, w in peaks:
            y += sign * h * (w * w) / ((x - c) ** 2 + w * w)
        pts.append((x, y))
    return pts


# --------------------------------------------------------------------------
# Illustration ("bio") mode: BioRender-like flat illustration
# --------------------------------------------------------------------------
import colorsys
import zlib


def shade(hex_color: str, dl: float = -0.30, ds: float = 0.05) -> str:
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))
    hh, l, s = colorsys.rgb_to_hls(r, g, b)
    l = min(1, max(0, l + dl))
    s = min(1, max(0, s + ds))
    r, g, b = colorsys.hls_to_rgb(hh, l, s)
    return "#{:02X}{:02X}{:02X}".format(round(r * 255), round(g * 255), round(b * 255))


@dataclass
class BioPalette:
    # 2.19.0 scientific palette (docs/icons/STYLE.md): restrained, low saturation, one steel-blue accent
    body: str = "#E9EDF1"         # instruments, hardware
    glass_fill: str = "#DCEAF3"
    glass_line: str = "#5F7A90"
    panel: str = "#F6F8FA"        # windows, displays, holes
    ink: str = "#1F2D3D"          # axes and text-like strokes
    accent: str = "#2E6F95"       # the one element the icon is about
    liquids: tuple = ("#8FB8CC", "#A7C4B5", "#B8B3D6", "#D9C49A", "#C9A7A7", "#9FB7A0")
    glass_opacity: float = 0.55
    outline: float = 1.2


def _attrs(op: Op) -> str:
    return " ".join(f'{k}="{_f(v) if isinstance(v, (int, float)) else v}"' for k, v in op.attrs.items()
                    if not k.startswith("_")) + f' class="{_cls(op)}"'


def render_bio(icon: Icon, pal: BioPalette | None = None, size: int = 64, id_prefix: str | None = None,
               standalone: bool = True, liquid: str | None = None, variant: str = "illustration") -> str:
    """variant: illustration | shaded | flat | line"""
    pal = pal or BioPalette()
    pid = id_prefix or icon.name
    liq = liquid or pal.liquids[zlib.crc32(icon.name.encode()) % len(pal.liquids)]
    data = getattr(icon, "ink", False) or icon.category in ("data", "annotations") or any(t in icon.name for t in ("pattern", "spectrum", "plot", "curve"))
    has_glass = any(o.glass for o in icon.ops)
    shaded, flat, line = variant == "shaded", variant == "flat", variant == "line"
    body = "#CDD6E0" if flat else pal.body
    col = {"m": pal.ink if data else body, "s": liq, "a": pal.accent}
    ow = 0 if flat else pal.outline
    grads: dict[str, str] = {}
    defs: list[str] = list(icon.defs)

    def fill_of(c):
        if not shaded:
            return c
        gid = grads.get(c)
        if not gid:
            gid = f"{pid}-g{len(grads)}"
            grads[c] = gid
            defs.append(f'<linearGradient id="{gid}" x1="0" y1="0" x2="0.35" y2="1">'
                        f'<stop offset="0" stop-color="{shade(c, 0.22, -0.04)}"/><stop offset="0.45" stop-color="{c}"/>'
                        f'<stop offset="1" stop-color="{shade(c, -0.20, 0.06)}"/></linearGradient>')
        return f"url(#{gid})"

    def glass_fill():
        if not shaded:
            return pal.glass_fill, pal.glass_opacity
        gid = f"{pid}-glass"
        if gid not in grads.values():
            grads["__glass"] = gid
            defs.append(f'<linearGradient id="{gid}" x1="0" y1="0" x2="1" y2="0">'
                        f'<stop offset="0" stop-color="{pal.glass_fill}" stop-opacity="0.95"/>'
                        f'<stop offset="0.22" stop-color="#FFFFFF" stop-opacity="0.9"/>'
                        f'<stop offset="0.5" stop-color="{pal.glass_fill}" stop-opacity="0.35"/>'
                        f'<stop offset="1" stop-color="{shade(pal.glass_fill, -0.08)}" stop-opacity="0.85"/>'
                        '</linearGradient>')
        return f"url(#{gid})", 1.0

    def filled(op, fill, line_c, opacity=1.0):
        o = f' fill-opacity="{_f(opacity)}"' if opacity < 1 else ""
        st = f' stroke="{line_c}" stroke-width="{_f(ow)}" stroke-linejoin="round"' if ow and op.tag != "text" else ""
        return f'<{op.tag} {_attrs(op)} fill="{fill}"{o}{st}/>'

    def stroke(op, c, w, extra=""):
        return (f'<{op.tag} {_attrs(op)} fill="none" stroke="{c}" stroke-width="{_f(w)}" '
                f'stroke-linecap="round" stroke-linejoin="round"{extra}/>')

    def paint_line(op):
        ink = pal.ink
        if op.role == "x":
            return f'<{op.tag} {_attrs(op)} fill="#FFFFFF" stroke="{ink}" stroke-width="1.6" stroke-linejoin="round"/>'
        c = pal.accent if op.role == "a" else ink
        if op.stroke is None and "_fill" in op.attrs:
            if op.tag == "text":
                return f'<{op.tag} {_attrs(op)} fill="{ink}"/>'
            return f'<{op.tag} {_attrs(op)} fill="{op.attrs["_fill"]}" stroke="{ink}" stroke-width="1.6" stroke-linejoin="round"/>'
        if op.stroke is None and op.tag == "text":
            return f'<{op.tag} {_attrs(op)} fill="{ink}"/>'
        if op.stroke is None:
            if op.role == "s":
                return (f'<{op.tag} {_attrs(op)} fill="{liq}" fill-opacity="0.35" stroke="{ink}" '
                        f'stroke-width="1.8" stroke-linejoin="round"/>')
            if op.role == "a":
                return f'<{op.tag} {_attrs(op)} fill="{pal.accent}" stroke="{ink}" stroke-width="1.6" stroke-linejoin="round"/>'
            return f'<{op.tag} {_attrs(op)} fill="#FFFFFF" stroke="{ink}" stroke-width="2" stroke-linejoin="round"/>'
        w = max(1.6, op.stroke * (0.7 if op.glass or op.role == "m" else 0.85))
        return stroke(op, c, w)

    def paint(op):
        if op.depth:
            dc = shade(pal.ink if line else (col[op.role] if op.role != "x" else body), -0.22 if not line else 0.25)
            if op.stroke:
                return stroke(op, dc, op.stroke + (0 if line else ow))
            return f'<{op.tag} {_attrs(op)} fill="{dc}" stroke="{dc}" stroke-width="{_f(max(ow, 1))}" stroke-linejoin="round"/>'
        if line:
            return paint_line(op)
        if op.role == "x":                      # filled cut -> panel (text cut may carry its own colour)
            if op.tag == "text":
                return f'<{op.tag} {_attrs(op)} fill="{op.attrs.get("_fill", pal.panel)}"/>'
            return filled(op, pal.panel, shade(body, -0.32))
        c = op.attrs.get("_fill") or col[op.role]
        if op.stroke is None:
            if op.tag == "text":
                return f'<{op.tag} {_attrs(op)} fill="{c}"/>'
            if data and op.role == "m":
                return f'<{op.tag} {_attrs(op)} fill="{c}"/>'
            return filled(op, fill_of(c), shade(c))
        if op.glass:
            gf, go = glass_fill()
            return (f'<{op.tag} {_attrs(op)} fill="{gf}" fill-opacity="{_f(go)}" stroke="none"/>'
                    + stroke(op, pal.glass_line, op.stroke * 0.8))
        if op.role == "m" and has_glass:
            return stroke(op, pal.glass_line, op.stroke * 0.8)
        if data or flat:
            return stroke(op, c, op.stroke)
        return stroke(op, shade(c), op.stroke + ow * 1.6) + stroke(op, c, op.stroke)

    out, k, i, ops = [], 0, 0, icon.ops
    while i < len(ops):
        op = ops[i]
        if not (op.role == "x" and op.stroke):
            out.append(_fin(op, paint(op)))
            i += 1
            continue
        cuts = []
        while i < len(ops) and ops[i].role == "x" and ops[i].stroke:
            cuts.append(ops[i])
            i += 1
        # painted in the panel colour like filled cuts: no <mask>, so PDF export stays 100 % vector
        # (cairo / rsvg rasterise masks at ~72 dpi)
        out += [stroke(c, pal.panel, c.stroke) for c in cuts]

    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="-1 -1 {GRID + 2} {GRID + 2}" '
           f'width="{size}" height="{size}"><title>{icon.label}</title>'
           + (f"<defs>{''.join(defs)}</defs>" if defs else "")
           + (f'<g transform="{icon.wrap}">' if icon.wrap else "") + "".join(out)
           + ("</g>" if icon.wrap else "") + "</svg>").replace("__PID__", pid)
    return ('<?xml version="1.0" encoding="UTF-8"?>\n' + svg) if standalone else svg


VARIANTS = ("illustration", "shaded", "flat", "line", "solid", "duotone")


THEMES: dict[str, BioPalette] = {
    "default": BioPalette(),
    "okabe_ito": BioPalette(liquids=("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7"),
                            accent="#D55E00"),
    "pastel": BioPalette(liquids=("#F7B6C8", "#A8E0D4", "#B4C9F5", "#F9D7A0", "#D2C0F0", "#C2E3A6"),
                         accent="#F5BE75", body="#EEF1F5"),
    "vivid": BioPalette(liquids=("#E8336B", "#00A88F", "#3D6FE0", "#F59E0B", "#8B5CF6", "#4CAF50"),
                        accent="#FF7A00", body="#DCE3EB"),
    "dark": BioPalette(body="#3B4656", glass_fill="#5B7A99", glass_line="#9CC3E6", panel="#1E2632",
                       ink="#D7DEE8", accent="#F7B955", glass_opacity=0.45),
    "mono_slate": BioPalette(liquids=("#7C8BA1",), accent="#4B5A70", body="#E5E9EF"),
}
SOLID_COLOURS = {"default": "#2F3B4A", "dark": "#E6EBF2"}


def render_variant(icon: Icon, variant: str, size: int = 64, id_prefix: str | None = None,
                   standalone: bool = True, liquid: str | None = None, theme: str = "default",
                   knockout: str | None = None) -> str:
    """knockout: paint cut-outs in this colour (the page background) instead of masking them.
    Use it for PDF / print export: masks are rasterised by cairo-based converters."""
    if variant in ("illustration", "shaded", "flat", "line"):
        return render_bio(icon, THEMES.get(theme, THEMES["default"]), size, id_prefix, standalone, liquid, variant)
    return render_svg(icon, Style(mode=variant, color=SOLID_COLOURS.get(theme, "#2F3B4A"), knockout=knockout), size,
                      id_prefix=id_prefix, standalone=standalone)
