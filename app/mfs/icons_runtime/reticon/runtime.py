"""
Runtime API for Reticula (or any app): one object, no pre-rendered files needed.

    from reticon.runtime import Library
    lib = Library.load("specs.json")              # ~1 MB JSON with every original icon + all shapes
    svg = lib.svg("round_bottom_flask", shape="two-neck", style="shaded", theme="okabe_ito",
                  liquid="#E69F00", fill=0.6, size=96)
    lib.anchors("round_bottom_flask", "two-neck")  # connection points for arrows / tubing
    lib.search("kolf", lang="nl")                  # multilingual search

Rendering is pure Python string building (~0.1-0.5 ms per icon) and cached (LRU), so it can run inside a
Streamlit rerun without threads. External CC0 SVGs are recoloured / reshaped on demand the same way.
"""
from __future__ import annotations

import copy
import functools
import json
import math
import re
from pathlib import Path

from .core import Icon, render_variant, VARIANTS, THEMES, shade, _f
from .variants import SEP, LABELS

STYLES = VARIANTS
_TF = re.compile(r"(translate|scale|rotate)\(([^)]*)\)")


# ----------------------------------------------------------------- geometry
def _matrix(tf: str | None):
    """SVG transform list -> 2x3 affine (a, b, c, d, e, f)."""
    m = (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)

    def mul(p, q):
        a1, b1, c1, d1, e1, f1 = p
        a2, b2, c2, d2, e2, f2 = q
        return (a1 * a2 + c1 * b2, b1 * a2 + d1 * b2, a1 * c2 + c1 * d2, b1 * c2 + d1 * d2,
                a1 * e2 + c1 * f2 + e1, b1 * e2 + d1 * f2 + f1)
    for kind, args in _TF.findall(tf or ""):
        v = [float(x) for x in args.replace(",", " ").split()]
        if kind == "translate":
            q = (1, 0, 0, 1, v[0], v[1] if len(v) > 1 else 0)
        elif kind == "scale":
            sx = v[0]
            q = (sx, 0, 0, v[1] if len(v) > 1 else sx, 0, 0)
        else:
            a = math.radians(v[0])
            r = (math.cos(a), math.sin(a), -math.sin(a), math.cos(a), 0, 0)
            if len(v) == 3:
                q = mul(mul((1, 0, 0, 1, v[1], v[2]), r), (1, 0, 0, 1, -v[1], -v[2]))
            else:
                q = r
        m = mul(m, q)
    return m


def _apply(m, x, y):
    a, b, c, d, e, f = m
    return a * x + c * y + e, b * x + d * y + f


_DIR_VEC = {"up": (0, -1), "down": (0, 1), "left": (-1, 0), "right": (1, 0)}


def _dir_after(m, direction):
    if direction not in _DIR_VEC:
        return direction
    dx, dy = _DIR_VEC[direction]
    a, b, c, d, _, _ = m
    vx, vy = a * dx + c * dy, b * dx + d * dy
    if abs(vx) > abs(vy):
        return "right" if vx > 0 else "left"
    return "down" if vy > 0 else "up"


def auto_anchors(icon: Icon) -> list[dict]:
    from .variants import bbox, union
    bb = union([bbox(o) for o in icon.ops if o.role != "x" and not o.depth and "transform" not in o.attrs])
    if not bb:
        bb = (8, 8, 56, 56)
    x0, y0, x1, y1 = bb
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    return [{"id": "top", "x": cx, "y": y0, "dir": "up"}, {"id": "bottom", "x": cx, "y": y1, "dir": "down"},
            {"id": "left", "x": x0, "y": cy, "dir": "left"}, {"id": "right", "x": x1, "y": cy, "dir": "right"},
            {"id": "center", "x": cx, "y": cy, "dir": "out"}]


def resolved_anchors(icon: Icon) -> list[dict]:
    """Explicit anchors + automatic box anchors (explicit ids win), mapped through the variant transform."""
    base = {a["id"]: a for a in auto_anchors(icon)}
    for a in icon.anchors:
        base[a["id"]] = a
    m = _matrix(icon.anchor_tf or icon.wrap)
    out = []
    for a in base.values():
        x, y = _apply(m, a["x"], a["y"])
        out.append({"id": a["id"], "x": round(x, 2), "y": round(y, 2), "dir": _dir_after(m, a["dir"])})
    return out


def with_fill(icon: Icon, fill: float, pid: str) -> Icon:
    """Lower the liquid level to `fill` (0-1) of the drawn level by clipping the liquid shapes."""
    from .variants import liquid_ops, bbox, union
    idx = liquid_ops(icon)
    if not idx or fill >= 0.999:
        return icon
    ic = copy.copy(icon)
    ic.ops = copy.deepcopy(icon.ops)
    ic.defs = list(icon.defs)
    bb = union([bbox(icon.ops[i]) for i in idx])
    y_cut = bb[3] - max(0.0, min(1.0, fill)) * (bb[3] - bb[1])
    cid = f"{pid}-fill"
    ic.defs.append(f'<clipPath id="{cid}"><rect x="-10" y="{_f(y_cut)}" width="84" height="{_f(80 - y_cut)}"/></clipPath>')
    for i in idx:
        ic.ops[i].attrs["clip-path"] = f"url(#{cid})"
    # particles / bubbles inside the liquid follow the same clip
    for o in ic.ops:
        if o.role == "a" and o.tag == "circle" and o.stroke is None:
            b = bbox(o)
            if b and bb[0] - 1 <= b[0] and b[2] <= bb[2] + 1 and bb[1] - 1 <= b[1] and b[3] <= bb[3] + 1:
                o.attrs["clip-path"] = f"url(#{cid})"
    return ic


# ----------------------------------------------------------------- CSS variables (live recolouring in a browser)
_EL = re.compile(r"<(path|rect|circle|ellipse|line|polygon|polyline)\b([^>]*?)(/?)>")


def css_variables(svg: str, theme: str = "default", liquid: str | None = None) -> str:
    """Add style="fill:var(--rt-…, #hex)" so a browser can recolour instantly with CSS variables:
    --rt-body --rt-liquid --rt-accent --rt-glass --rt-glass-line --rt-panel --rt-ink.
    Outline shades follow via color-mix(). Other renderers ignore the vars and use the fallback colour."""
    pal = THEMES.get(theme, THEMES["default"])
    names = {pal.body: "body", pal.accent: "accent", pal.glass_fill: "glass", pal.glass_line: "glass-line",
             pal.panel: "panel", pal.ink: "ink"}
    if liquid:
        names[liquid] = "liquid"
    for c in pal.liquids:
        names.setdefault(c, "liquid")
    shades = {shade(c): n for c, n in names.items()}

    def var(col):
        cu = col.upper()
        for c, n in names.items():
            if c.upper() == cu:
                return f"var(--rt-{n}, {col})"
        for c, n in shades.items():
            if c.upper() == cu:
                return f"var(--rt-{n}-line, color-mix(in srgb, var(--rt-{n}, {col}) 70%, #1F2937))"
        return None

    def sub(m):
        tag, attrs, close = m.groups()
        st = []
        for prop in ("fill", "stroke"):
            mm = re.search(rf'\b{prop}="(#[0-9A-Fa-f]{{6}})"', attrs)
            if mm:
                v = var(mm.group(1))
                if v:
                    st.append(f"{prop}:{v}")
        if not st or "style=" in attrs:
            return m.group(0)
        return f'<{tag}{attrs} style="{";".join(st)}"{close}>'
    return _EL.sub(sub, svg)


def add_metadata(svg: str, meta: dict) -> str:
    desc = meta.get("desc") or f'{meta.get("label", "")} icon ({meta.get("category", "")})'
    md = json.dumps({k: meta[k] for k in ("id", "source", "license", "shape", "style", "theme") if k in meta})
    block = f"<desc>{_esc(desc)}</desc><metadata>{_esc(md)}</metadata>"
    return svg.replace("</title>", "</title>" + block, 1)


def _esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


# ----------------------------------------------------------------- library
class Library:
    def __init__(self, data: dict, root: Path | None = None):
        self.data = data
        self.root = root
        self.icons: dict = data["icons"]
        self.external: dict = data.get("external", {})
        self.aliases: dict = data.get("aliases", {})
        self._variants_loaded = all(len(v.get("shapes", {})) > 1 for v in list(self.icons.values())[:1])
        self.synonyms: dict = data.get("synonyms", {})
        self._icon_cache: dict = {}

    # ---- loading
    @classmethod
    def load(cls, path: str | Path | None = None) -> "Library":
        p = Path(path) if path else Path(__file__).resolve().parent.parent / "dist" / "specs.json"
        return cls(json.loads(p.read_text(encoding="utf-8")), p.parent)

    def load_variants(self) -> "Library":
        """Merge specs_variants.json (the 5 variants of every icon). Called automatically on first use of a
        non-base shape; call it up front (e.g. inside st.cache_resource) to pay the ~200 ms once at start-up."""
        if not self._variants_loaded:
            p = (self.root or Path(".")) / "specs_variants.json"
            if p.exists():
                for n, shapes in json.loads(p.read_text(encoding="utf-8")).items():
                    if n in self.icons:
                        self.icons[n]["shapes"].update(shapes)
            self._variants_loaded = True
        return self

    def _spec(self, n: str, shape: str) -> dict:
        sh = self.icons[n]["shapes"]
        if shape not in sh:
            self.load_variants()
        return sh[shape]

    # ---- lookup
    def resolve(self, name: str) -> str:
        name = name.split(":", 1)[-1] if name.startswith("reticon:") else name
        name = name.split(SEP, 1)[0]
        return self.aliases.get(name, name)

    @staticmethod
    def split(vid: str, shape: str = "base") -> tuple[str, str]:
        """'beaker--crystals' -> ('beaker', 'crystals'). An explicit shape argument wins over the id."""
        if SEP in vid:
            n, k = vid.split(SEP, 1)
            return n, (shape if shape != "base" else k)
        return vid, shape

    def family(self, name: str) -> list[dict]:
        """Every member of an icon's variant family: [{id, key, label}], base first.
        Use it to open the variant strip when the user clicks any member."""
        n = self.resolve(name)
        if n in self.icons:
            return list(self.icons[n].get("variants") or
                        [{"id": n if k == "base" else f"{n}{SEP}{k}", "key": k, "label": k} for k in self.icons[n]["shapes"]])
        e = self.external[n]
        return [{"id": n if k == "base" else f"{n}{SEP}{k}", "key": k, "label": LABELS.get(k, k)} for k in e["shapes"]]

    def info(self, name: str) -> dict:
        n = self.resolve(name)
        if n in self.icons:
            return {k: v for k, v in self.icons[n].items() if k != "shapes"} | {"name": n, "shapes": [m["key"] for m in self.family(n)]}
        if n in self.external:
            e = self.external[n]
            return e | {"name": n}
        raise KeyError(name)

    def names(self, category: str | None = None, source: str | None = None):
        out = [n for n, v in self.icons.items() if (not category or v["category"] == category)]
        if source != "reticon":
            out += [n for n, v in self.external.items() if (not category or v["category"] == category)]
        return out

    def icon(self, name: str, shape: str = "base") -> Icon:
        name, shape = self.split(name, shape)
        key = (self.resolve(name), shape)
        ic = self._icon_cache.get(key)
        if ic is None:
            meta = self.icons[key[0]]
            spec = dict(self._spec(key[0], shape))
            vlabel = spec.pop("label", shape)
            label = meta["label"] if shape == "base" else f'{meta["label"]} ({vlabel})'
            ic = Icon.from_spec({**spec, "name": key[0], "label": label, "category": meta["category"],
                                 "tags": meta.get("tags", [])})
            self._icon_cache[key] = ic
        return ic

    def anchors(self, name: str, shape: str = "base") -> list[dict]:
        name, shape = self.split(name, shape)
        n = self.resolve(name)
        if n in self.icons:
            return self._spec(n, shape).get("anchors_resolved") or resolved_anchors(self.icon(n, shape))
        return [{"id": k, "x": x, "y": y, "dir": d} for k, x, y, d in
                (("top", 32, 2, "up"), ("bottom", 32, 62, "down"), ("left", 2, 32, "left"),
                 ("right", 62, 32, "right"), ("center", 32, 32, "out"))]

    # ---- rendering
    def svg(self, name: str, shape: str = "base", style: str = "illustration", theme: str = "default",
            liquid: str | None = None, fill: float | None = None, size: int = 64, css_vars: bool = False,
            metadata: bool = True, id_prefix: str | None = None, print_safe: str | bool = False) -> str:
        """print_safe: True (white) or a background colour -> no <mask> elements, so PDF / EPS export stays
        fully vector. Only solid / duotone use masks; illustration styles never do."""
        name, shape = self.split(name, shape)          # "beaker--crystals" works as well as shape="crystals"
        ko = None if not print_safe else ("#FFFFFF" if print_safe is True else print_safe)
        return self._svg(self.resolve(name), shape, style, theme, liquid, None if fill is None else round(fill, 3),
                         size, css_vars, metadata, id_prefix, ko)

    def is_vector_only(self, name: str) -> bool:
        """False for the few CC0 source icons that embed a bitmap or a blur filter (they print at their own resolution)."""
        n = self.resolve(name.split(SEP, 1)[0])
        return self.external.get(n, {}).get("vector_only", True)

    @functools.lru_cache(maxsize=4096)
    def _svg(self, n, shape, style, theme, liquid, fill, size, css_vars, metadata, id_prefix, knockout=None):
        if n in self.external:
            return self._external_svg(n, shape, style, size)
        meta = self.icons[n]
        pid = id_prefix or f"{n}-{shape}-{style}"
        ic = self.icon(n, shape)
        if fill is not None:
            ic = with_fill(ic, fill, pid)
        s = render_variant(ic, style, size, id_prefix=pid, standalone=False, liquid=liquid, theme=theme,
                           knockout=knockout)
        if css_vars:
            s = css_variables(s, theme, liquid)
        if metadata:
            s = add_metadata(s, {"id": f"reticon:{n}", "label": meta["label"], "category": meta["category"],
                                 "source": "Reticon (original)", "license": "CC0-1.0", "shape": shape,
                                 "style": style, "theme": theme})
        return s

    def _external_svg(self, n, shape, colour, size):
        from . import external as X
        e = self.external[n]
        txt = (self.root / e["file"]).read_text(errors="ignore", encoding="utf-8")
        shapes = X.shape_variants(txt, e.get("allow_mirror", True), e.get("source_category"))
        if shape not in shapes:
            shape = next(k for k in shapes if k not in ("base", "3d")) if shape not in ("base", "3d") else shape
        cols = X.colour_variants(shapes[shape], e.get("kind", "illustration"))
        s = cols.get(colour) or cols[next(iter(cols))]
        return re.sub(r'width="[^"]*" height="[^"]*"', f'width="{size}" height="{size}"', s, count=1) \
            if 'width="' in s[:400] else s

    # ---- search
    def search(self, query: str, lang: str | None = None, limit: int = 40, category: str | None = None):
        q = [t for t in re.split(r"[\s,]+", query.lower().strip()) if t]
        if not q:
            return []
        expanded = []
        for t in q:
            alts = {t}
            for canon, tr in self.synonyms.items():
                words = {canon} | {w for ws in tr.values() for w in ws}
                if lang:
                    words |= set(tr.get(lang, []))
                if any(t == w or (len(t) > 3 and w.startswith(t)) for w in words):
                    alts |= words
            expanded.append(alts)
        scored = []

        def score(n, v, extra=""):
            hay_name = n.lower().replace("_", " ").replace("-", " ")
            hay = " ".join([hay_name, v.get("label", "").lower(), " ".join(v.get("tags", [])).lower(),
                            v.get("category", "").lower(), extra])
            s = 0
            for alts in expanded:
                best = 0
                for a in alts:
                    if a in hay_name.split():
                        best = max(best, 5)
                    elif a in hay:
                        best = max(best, 3 if a in v.get("label", "").lower() else 2)
                if not best:
                    return 0
                s += best
            return s

        for pool in (self.icons, self.external):
            for n, v in pool.items():
                if category and v["category"] != category:
                    continue
                s = score(n, v)
                if s:
                    scored.append((s + (1 if pool is self.icons else 0), n))
                    continue
                # "beaker crystals", "heated autoclave": match a variant of the family
                fam = v.get("variants") or [{"id": f"{n}{SEP}{k}", "key": k, "label": LABELS.get(k, k)}
                                            for k in (v["shapes"] if isinstance(v["shapes"], list) else v["shapes"].keys())]
                for m in fam:
                    if m["key"] == "base":
                        continue
                    s = score(n, v, f'{m["key"].replace("-", " ")} {m["label"].lower()}')
                    if s:
                        scored.append((s, m["id"]))
                        break
        return [n for _, n in sorted(scored, key=lambda t: (-t[0], t[1]))[:limit]]
