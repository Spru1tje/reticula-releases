"""
Scenes: BioRender-style schemes as data (JSON), rendered to one vector SVG (PDF via rsvg-convert).

Scene format (see docs/SCENE_FORMAT.md for the full schema):
{
  "version": 1, "width": 900, "height": 420, "background": "#FFFFFF", "theme": "default", "font": "Helvetica, Arial, sans-serif",
  "items": [
    {"id": "flask", "type": "icon", "icon": "round_bottom_flask", "shape": "base", "style": "illustration",
     "x": 80, "y": 120, "size": 110, "rotate": 0, "flip": false, "liquid": null, "fill": null,
     "label": {"text": "ZrCl4 + linker", "pos": "below"}},
    {"id": "a1", "type": "arrow", "from": {"item": "flask", "anchor": "right"}, "to": {"item": "oven", "anchor": "left"},
     "route": "straight|elbow|curved", "label": "120 °C, 24 h", "head": "end|both|none", "dash": false},
    {"id": "t", "type": "text", "x": 450, "y": 30, "text": "Synthesis", "size": 18, "weight": "bold", "anchor": "middle"},
    {"id": "bg", "type": "box", "x": 20, "y": 60, "w": 400, "h": 300, "fill": "#F1F5F9", "radius": 14, "title": "Step 1"}
  ]
}
x, y = top-left of the icon box (size x size px). Anchors come from the icon (named or automatic) and are mapped
through rotate / flip, so arrows stay attached when icons move.
"""
from __future__ import annotations

import json
import math
import re
import shutil
import subprocess
from pathlib import Path

INK = "#435266"
VB = 66.0          # icon viewBox size (-1 .. 65)


def _esc(s: str) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _rich(s: str) -> str:
    """Text with _{sub} and ^{sup} (the same markup as the editor): escaped, then baseline-shifted runs."""
    out = _esc(s)
    out = re.sub(r"_\{([^}]*)\}", r'<tspan baseline-shift="sub" font-size="72%">\1</tspan>', out)
    return re.sub(r"\^\{([^}]*)\}", r'<tspan baseline-shift="super" font-size="72%">\1</tspan>', out)


class SceneRenderer:
    def __init__(self, lib, print_safe: bool = False):
        self.lib = lib
        self.print_safe = print_safe      # no <mask>: cut-outs painted in the colour behind the icon

    # ---- geometry
    def anchor_point(self, item: dict, anchor_id: str):
        anchors = {a["id"]: a for a in self.lib.anchors(item["icon"], item.get("shape", "base"))}
        a = anchors.get(anchor_id) or anchors["center"]
        s = item["size"] / VB
        lx, ly = (a["x"] + 1) * s, (a["y"] + 1) * s
        cx = cy = item["size"] / 2
        d = a["dir"]
        if item.get("flip"):
            lx = item["size"] - lx
            d = {"left": "right", "right": "left"}.get(d, d)
        r = math.radians(item.get("rotate", 0))
        dx, dy = lx - cx, ly - cy
        x = item["x"] + cx + dx * math.cos(r) - dy * math.sin(r)
        y = item["y"] + cy + dx * math.sin(r) + dy * math.cos(r)
        vec = {"up": (0, -1), "down": (0, 1), "left": (-1, 0), "right": (1, 0)}.get(d)
        if vec:
            vx, vy = vec
            vec = (vx * math.cos(r) - vy * math.sin(r), vx * math.sin(r) + vy * math.cos(r))
        return x, y, vec

    def _end(self, items, end):
        if "item" in end:
            it = items[end["item"]]
            if it["type"] == "icon":
                return self.anchor_point(it, end.get("anchor", "center"))
            if it["type"] == "box":
                return it["x"] + it["w"] / 2, it["y"] + it["h"] / 2, None
        return end["x"], end["y"], None

    # ---- items
    def icon(self, it, theme):
        svg = self.lib.svg(it["icon"], it.get("shape", "base"), it.get("style", "illustration"),
                           it.get("theme") or theme, it.get("liquid"), it.get("fill"), int(it["size"]),
                           metadata=False, id_prefix=f"s-{it['id']}",
                           print_safe=self._behind(it, self._scene) if self.print_safe else False)
        svg = re.sub(r"^<\?xml[^>]*>\s*", "", svg)
        svg = svg.replace("<svg ", f'<svg x="0" y="0" ', 1)
        s = it["size"]
        tf = [f"translate({it['x']:.2f} {it['y']:.2f})"]
        if it.get("rotate"):
            tf.append(f"rotate({it['rotate']} {s / 2} {s / 2})")
        if it.get("flip"):
            tf.append(f"translate({s} 0) scale(-1 1)")
        out = f'<g id="{_esc(it["id"])}" transform="{" ".join(tf)}">{svg}</g>'
        lab = it.get("label")
        if lab:
            lab = lab if isinstance(lab, dict) else {"text": lab}
            pos = lab.get("pos", "below")
            lx = it["x"] + s / 2
            ly = it["y"] + s + 16 if pos == "below" else it["y"] - 8
            out += self.text({"x": lx, "y": ly, "text": lab["text"], "size": lab.get("size", 13), "anchor": "middle",
                              "color": lab.get("color", INK), "weight": lab.get("weight", "normal")})
        return out

    def text(self, it, font=None):
        lines = str(it["text"]).split("\n")
        size = it.get("size", 13)
        f = font or self.font
        tsp = "".join(f'<tspan x="{it["x"]:.1f}" dy="{0 if i == 0 else size * 1.2:.1f}">{_rich(l)}</tspan>'
                      for i, l in enumerate(lines))
        return (f'<text x="{it["x"]:.1f}" y="{it["y"]:.1f}" font-family="{f}" font-size="{size}" '
                f'font-weight="{it.get("weight", "normal")}" fill="{it.get("color", INK)}" '
                f'text-anchor="{it.get("anchor", "start")}">{tsp}</text>')

    def svgart(self, it):
        """An embedded vector figure (a nested SVG with its own viewBox): the body is trusted only from mfs.schemes."""
        return (f'<svg x="{it["x"]:.2f}" y="{it["y"]:.2f}" width="{it["size"]:g}" height="{it["size"]:g}" viewBox="{it["viewBox"]}" '
                f'preserveAspectRatio="xMidYMid meet" overflow="visible">{it["body"]}</svg>')

    def box(self, it):
        out = (f'<rect x="{it["x"]}" y="{it["y"]}" width="{it["w"]}" height="{it["h"]}" rx="{it.get("radius", 12)}" '
               f'fill="{it.get("fill", "#F1F5F9")}" stroke="{it.get("stroke", "none")}" stroke-width="{it.get("stroke_width", 1.5)}"'
               + (' stroke-dasharray="6 4"' if it.get("dash") else "") + "/>")
        if it.get("title"):
            out += self.text({"x": it["x"] + 14, "y": it["y"] + 24, "text": it["title"], "size": it.get("title_size", 14),
                              "weight": "bold", "color": it.get("title_color", INK)})
        return out

    def arrow(self, it, items):
        x1, y1, v1 = self._end(items, it["from"])
        x2, y2, v2 = self._end(items, it["to"])
        gap = it.get("gap", 6)
        if v1:
            x1, y1 = x1 + v1[0] * gap, y1 + v1[1] * gap
        if v2:
            x2, y2 = x2 + v2[0] * gap, y2 + v2[1] * gap
        route = it.get("route", "straight")
        wps = [tuple(map(float, p)) for p in (it.get("waypoints") or [])]      # bends placed by hand or by the auto-router
        if route == "straight" and it.get("straighten", True) and not wps:
            if abs(y2 - y1) < 14 and abs(x2 - x1) > abs(y2 - y1):
                y1 = y2 = (y1 + y2) / 2
            elif abs(x2 - x1) < 14 and abs(y2 - y1) > abs(x2 - x1):
                x1 = x2 = (x1 + x2) / 2
        col, w = it.get("color", INK), it.get("width", 2.6)
        if route == "curved":
            k = max(40, math.dist((x1, y1), (x2, y2)) * 0.4)
            c1 = (x1 + (v1 or (1, 0))[0] * k, y1 + (v1 or (1, 0))[1] * k)
            c2 = (x2 + (v2 or (-1, 0))[0] * k, y2 + (v2 or (-1, 0))[1] * k)
            d = f"M{x1:.1f} {y1:.1f} C{c1[0]:.1f} {c1[1]:.1f} {c2[0]:.1f} {c2[1]:.1f} {x2:.1f} {y2:.1f}"
            tail_dir = (x2 - c2[0], y2 - c2[1])
            head_dir0 = (x1 - c1[0], y1 - c1[1])
            mid = (0.125 * x1 + 0.375 * c1[0] + 0.375 * c2[0] + 0.125 * x2, 0.125 * y1 + 0.375 * c1[1] + 0.375 * c2[1] + 0.125 * y2)
        elif route == "elbow":
            if v1 and abs(v1[1]) > abs(v1[0]):
                pts = [(x1, y1), (x1, (y1 + y2) / 2), (x2, (y1 + y2) / 2), (x2, y2)]
            else:
                pts = [(x1, y1), ((x1 + x2) / 2, y1), ((x1 + x2) / 2, y2), (x2, y2)]
            d = "M" + " L".join(f"{p[0]:.1f} {p[1]:.1f}" for p in pts)
            tail_dir = (pts[-1][0] - pts[-2][0], pts[-1][1] - pts[-2][1])
            head_dir0 = (pts[0][0] - pts[1][0], pts[0][1] - pts[1][1])
            mid = ((pts[1][0] + pts[2][0]) / 2, (pts[1][1] + pts[2][1]) / 2)
        elif wps:
            pts = [(x1, y1)] + wps + [(x2, y2)]
            d = "M" + " L".join(f"{p[0]:.1f} {p[1]:.1f}" for p in pts)
            tail_dir = (pts[-1][0] - pts[-2][0], pts[-1][1] - pts[-2][1])
            head_dir0 = (pts[0][0] - pts[1][0], pts[0][1] - pts[1][1])
            k = max(range(len(pts) - 1), key=lambda i: math.dist(pts[i], pts[i + 1]))
            mid = ((pts[k][0] + pts[k + 1][0]) / 2, (pts[k][1] + pts[k + 1][1]) / 2)
        else:
            d = f"M{x1:.1f} {y1:.1f} L{x2:.1f} {y2:.1f}"
            tail_dir, head_dir0 = (x2 - x1, y2 - y1), (x1 - x2, y1 - y2)
            mid = ((x1 + x2) / 2, (y1 + y2) / 2)
        hs = it.get("head_size", 11)
        out = (f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{w}" stroke-linecap="round" stroke-linejoin="round"'
               + (' stroke-dasharray="7 5"' if it.get("dash") else "") + "/>")

        def head(x, y, dv):
            a = math.atan2(dv[1], dv[0])
            p1 = (x - hs * math.cos(a - 0.42), y - hs * math.sin(a - 0.42))
            p2 = (x - hs * math.cos(a + 0.42), y - hs * math.sin(a + 0.42))
            return f'<polygon points="{x:.1f},{y:.1f} {p1[0]:.1f},{p1[1]:.1f} {p2[0]:.1f},{p2[1]:.1f}" fill="{col}"/>'
        hd = it.get("head", "end")
        if hd in ("end", "both"):
            out += head(x2, y2, tail_dir)
        if hd in ("start", "both"):
            out += head(x1, y1, head_dir0)
        if it.get("label"):
            horiz = abs(tail_dir[0]) >= abs(tail_dir[1])
            out += self.text({"x": mid[0] if horiz else mid[0] + 10, "y": mid[1] - 9 if horiz else mid[1] + 4,
                              "text": it["label"], "size": it.get("label_size", 12),
                              "anchor": "middle" if horiz else "start", "color": col})
        return out

    # ---- scene
    def _behind(self, it, scene) -> str:
        """Colour behind an icon: the top-most panel under its centre, else the page background."""
        cx, cy = it["x"] + it["size"] / 2, it["y"] + it["size"] / 2
        col = scene.get("background") or "#FFFFFF"
        for b in scene["items"]:
            if b["type"] == "box" and b["x"] <= cx <= b["x"] + b["w"] and b["y"] <= cy <= b["y"] + b["h"]:
                col = b.get("fill", "#F1F5F9")
        return col

    def render(self, scene: dict) -> str:
        self._scene = scene
        self.font = scene.get("font", "Helvetica, Arial, sans-serif")
        W, H = scene["width"], scene["height"]
        items = {it["id"]: it for it in scene["items"]}
        order = sorted(scene["items"], key=lambda it: (it.get("z", {"box": 0, "arrow": 2, "icon": 1, "svgart": 1, "text": 3}[it["type"]])))
        parts = []
        for it in order:
            t = it["type"]
            if t == "icon":
                parts.append(self.icon(it, scene.get("theme", "default")))
            elif t == "arrow":
                parts.append(self.arrow(it, items))
            elif t == "text":
                parts.append(self.text(it))
            elif t == "box":
                parts.append(self.box(it))
            elif t == "svgart":
                parts.append(self.svgart(it))
        bg = scene.get("background")
        return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">'
                + (f'<rect width="{W}" height="{H}" fill="{bg}"/>' if bg else "") + "".join(parts) + "</svg>")


def render_scene(scene: dict | str | Path, lib=None, print_safe: bool = False) -> str:
    """print_safe=True: no SVG masks (cut-outs painted in the colour behind each icon), so PDF export is 100 % vector."""
    from .runtime import Library
    if not isinstance(scene, dict):
        scene = json.loads(Path(scene).read_text(encoding="utf-8"))
    return SceneRenderer(lib or Library.load(), print_safe).render(scene)


def to_pdf(svg: str, path: str | Path) -> bool:
    from .raster import has_pdf_backend, svg_to_pdf
    if not has_pdf_backend():
        return False
    svg_to_pdf(svg, path)
    return True


MIN_ICON_MM = 11.0       # 64-grid icon: 1-unit lines stay >= 0.5 pt (0.176 mm) from 11.3 mm up; most lines are 2-3.5 units
JOURNAL_WIDTHS_MM = {"single column": 85, "1.5 column": 120, "double column": 178, "toc graphic": 84.7}


def export(scene: dict | str | Path, path: str | Path, dpi: int = 600, width_mm: float | None = None,
           lib=None) -> dict:
    """Publication export. Format from the suffix: .svg / .pdf / .eps (vector) or .png / .tif / .tiff (raster).

    width_mm: printed width (85 single column, 178 double; see JOURNAL_WIDTHS_MM). Default: the page width
    at 96 px per inch. Raster files get exactly width_mm/25.4*dpi pixels and the DPI written into the file,
    so journals read the right physical size. Vector files are resolution-free; dpi is ignored.
    Returns {"path", "format", "width_px", "height_px", "dpi", "width_mm", "vector_only"}."""
    from .runtime import Library
    lib = lib or Library.load()
    if not isinstance(scene, dict):
        scene = json.loads(Path(scene).read_text(encoding="utf-8"))
    path = Path(path)
    fmt = path.suffix.lower().lstrip(".")
    W, H = scene["width"], scene["height"]
    width_mm = width_mm or W / 96 * 25.4
    height_mm = width_mm * H / W
    svg = render_scene(scene, lib, print_safe=True)
    # physical size on the root element, so PDF / EPS pages have the right dimensions
    svg = re.sub(r'(<svg\b[^>]*?)\swidth="[^"]*"\s+height="[^"]*"',
                 lambda m: f'{m.group(1)} width="{width_mm:.3f}mm" height="{height_mm:.3f}mm"', svg, count=1)
    vector_only = all(lib.is_vector_only(it["icon"]) for it in scene["items"] if it["type"] == "icon")
    warnings = []
    mm_per_px = width_mm / W
    for it in scene["items"]:
        if it["type"] == "icon":
            mm = it["size"] * mm_per_px
            if mm < MIN_ICON_MM:      # thinnest icon lines (1 grid unit) would print below 0.5 pt
                warnings.append(f'{it["icon"]} prints {mm:.1f} mm wide: thin lines fall below 0.5 pt (keep icons >= {MIN_ICON_MM} mm)')
            if not lib.is_vector_only(it["icon"]):
                warnings.append(f'{it["icon"]} contains a bitmap or blur from its CC0 source: not fully vector')
        if it["type"] == "text" and it.get("size", 14) * mm_per_px / 0.3528 < 5:
            warnings.append(f'text "{it.get("text", "")[:20]}" prints below 5 pt')
    if fmt in ("png", "tif", "tiff") and dpi < 300:
        warnings.append(f"{dpi} dpi is below the usual journal minimum (300 for images, 600-1000 for line art)")
    info = {"path": str(path), "format": fmt, "dpi": dpi, "width_mm": round(width_mm, 2), "vector_only": vector_only,
            "warnings": warnings}
    if fmt == "svg":
        path.write_text(svg, encoding="utf-8")
        return info
    from .raster import svg_to_png, svg_to_pdf
    if fmt == "pdf":
        svg_to_pdf(svg, path)
        return info
    if fmt in ("eps", "ps"):
        raise RuntimeError("EPS / PS export needs rsvg-convert (librsvg); use PDF or SVG instead")
    wpx = round(width_mm / 25.4 * dpi)
    hpx = round(height_mm / 25.4 * dpi)
    png = svg_to_png(svg, wpx, hpx, scene.get("background") or "#FFFFFF")
    from PIL import Image
    import io
    im = Image.open(io.BytesIO(png))
    if fmt in ("tif", "tiff"):
        im.convert("RGB").save(path, dpi=(dpi, dpi), compression="tiff_lzw")
    else:
        im.save(path, dpi=(dpi, dpi), optimize=True)
    return info | {"width_px": wpx, "height_px": hpx}


# ----------------------------------------------------------------- automatic methods schematic
TECHNIQUE_ICONS = {
    "pxrd": ("xrd", "PXRD"), "n2": ("sorption_analyzer", "N₂ sorption"), "tga": ("tga_instrument", "TGA"),
    "dsc": ("dsc_instrument", "DSC"), "kinetics": ("conversion_time", "Kinetics"), "icp": ("icp_oes", "ICP-OES"),
    "xrf": ("xrf_spectrometer", "XRF"), "xps": ("xps", "XPS"), "nmr_1h": ("nmr", "¹H NMR"), "nmr": ("nmr", "NMR"),
    "pl": ("fluorimeter", "Photoluminescence"), "uvvis": ("uv_vis", "UV-Vis"), "ftir": ("ftir_atr", "FTIR"),
    "raman": ("raman_microscope", "Raman"), "sem": ("electron_microscope", "SEM"), "tem": ("electron_microscope", "TEM"),
    "zeta": ("dls_zeta", "Zeta / DLS"), "dls": ("dls_zeta", "DLS"), "defects": ("missing_linker_defect", "Defects"),
    "sdspage": ("sds_page", "SDS-PAGE"), "bands": ("sds_page", "SDS-PAGE"), "edx": ("edx_spectrum", "EDX"),
}


def methods_schematic(sample: str, sections: list[str], material_icon: str = "crystal_octahedron",
                      synthesis: list[tuple[str, str]] | None = None, width: int = 980) -> dict:
    """Scene: (optional synthesis steps →) the sample → every technique it went through.
    sections: Reticula sample-store section keys (pxrd, n2, tga, …). synthesis: [(icon, label), …]."""
    items, x = [], 30
    prev = None
    for k, (icn, lab) in enumerate(synthesis or []):
        iid = f"syn{k}"
        items.append({"id": iid, "type": "icon", "icon": icn, "x": x, "y": 150, "size": 96, "label": lab})
        if prev:
            items.append({"id": f"as{k}", "type": "arrow", "from": {"item": prev, "anchor": "right"},
                          "to": {"item": iid, "anchor": "left"}})
        prev, x = iid, x + 150
    items.append({"id": "sample", "type": "icon", "icon": material_icon, "x": x, "y": 140, "size": 116,
                  "label": {"text": sample, "weight": "bold"}})
    if prev:
        items.append({"id": "as_last", "type": "arrow", "from": {"item": prev, "anchor": "right"},
                      "to": {"item": "sample", "anchor": "left"}})
    techs, seen = [], set()
    for s in sections:
        key = re.sub(r"[^a-z0-9]", "", s.lower())
        for k, v in TECHNIQUE_ICONS.items():
            if key.startswith(k) and v[0] not in seen:
                techs.append(v)
                seen.add(v[0])
                break
    cols = 2 if len(techs) > 4 else 1
    rows = math.ceil(len(techs) / cols) or 1
    cell_h = 110
    H = max(400, 40 + rows * cell_h)
    x0 = x + 230
    for i, (icn, lab) in enumerate(techs):
        c, r = divmod(i, rows)
        iid = f"t{i}"
        y = (H - rows * cell_h) / 2 + r * cell_h + 6
        items.append({"id": iid, "type": "icon", "icon": icn, "x": x0 + c * 210, "y": y, "size": 74, "label": lab})
        items.append({"id": f"at{i}", "type": "arrow", "from": {"item": "sample", "anchor": "right"},
                      "to": {"item": iid, "anchor": "left"}, "route": "curved", "width": 2, "head_size": 9})
    for it in items:
        if it["type"] == "icon" and it["id"] in ("sample",) or it["id"].startswith("syn"):
            it["y"] = (H - it["size"]) / 2 - 10
    W = max(width, x0 + cols * 210 + 40)
    return {"version": 1, "width": W, "height": H, "background": "#FFFFFF", "theme": "default", "items": items}
