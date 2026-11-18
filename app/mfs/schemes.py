"""v2.17 Schemes & Icons: the icon library (reticon, vendored in mfs/reticon), the scheme canvas (a browser component, served from
mfs/components/scheme_canvas) and the scenes of a project. Scenes are plain JSON dicts; they travel with the project file.

The library is loaded from mfs/icons_runtime (built by the icon toolkit); nothing is downloaded and no external service is used."""
from __future__ import annotations

import copy
import io
import json
import re
import tempfile
from pathlib import Path

import streamlit as st

HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "icons_runtime" / "specs.json"
CANVAS_DIR = HERE / "components" / "scheme_canvas"
TEMPLATE_DIR = HERE / "schemes_templates"
SCENES_KEY = "schemes"          # project-global session key: {name: scene dict}

_COMP = None


@st.cache_resource(show_spinner=False)
def library():
    """The icon library (loaded once per process)."""
    from mfs.reticon.runtime import Library
    return Library.load(str(RUNTIME))


def canvas_component():
    """The browser canvas. The icon catalogue is served as a static file by the component, never sent through its arguments."""
    global _COMP
    if _COMP is None:
        import streamlit.components.v1 as components
        _COMP = components.declare_component("reticula_scheme_canvas", path=str(CANVAS_DIR))
    return _COMP


def canvas(scene: dict, key: str, height: int = 980):
    """Show the editor for one scene. Returns the edited scene (dict) when the user changed it, otherwise None."""
    value = canvas_component()(scene=scene, height=height, key=key, default=None)
    if value is None:
        return None
    try:
        return json.loads(value) if isinstance(value, str) else value
    except json.JSONDecodeError:
        return None


def templates() -> dict[str, dict]:
    """The built-in starting schemes, by name."""
    return {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in sorted(TEMPLATE_DIR.glob("*.json"))}


def empty_scene() -> dict:
    return {"version": 1, "width": 960, "height": 540, "background": "#FFFFFF", "theme": "default", "items": []}


def scenes() -> dict[str, dict]:
    """The scenes of the project (kept in the session and saved with the project file)."""
    return st.session_state.setdefault(SCENES_KEY, {})


def add_scene(name: str, scene: dict) -> str:
    """Store a scene under a unique name; returns the name used."""
    base = name.strip() or "scheme"
    store, n, out = scenes(), 2, base
    while out in store:
        out = f"{base} ({n})"
        n += 1
    store[out] = copy.deepcopy(scene)
    return out


def export(scene: dict, fmt: str, width_mm: float | None = None, dpi: int = 600) -> tuple[bytes, dict]:
    """Publication export of one scene (SVG, PDF, PNG or TIFF). Returns (file bytes, info incl. warnings).
    The vendored toolkit gets the library explicitly, so it does not look for its own dist folder."""
    from mfs.reticon import scene as RS
    lib = library()
    fmt = fmt.lower().lstrip(".")
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / f"scheme.{fmt}"
        info = RS.export(scene, path, dpi=dpi, width_mm=width_mm, lib=lib)
        return path.read_bytes(), info


def load_json_scene(data: bytes) -> dict:
    """A scene from an uploaded JSON file (the format of docs/SCENE_FORMAT.md)."""
    scene = json.loads(data.decode("utf-8"))
    if not isinstance(scene, dict) or "items" not in scene:
        raise ValueError("this file is not a scheme (no 'items')")
    return scene


# ---------------------------------------------------------------- content from the user's own data (embedded vector figures)
def svgart_from_svg(svg_text: str, item_id: str, x: float, y: float, size: float) -> dict:
    """A scene item that embeds an SVG figure (nested, with its own viewBox). The body is the SVG content without its root."""
    s = re.sub(r"<\?xml[^>]*\?>\s*", "", svg_text.strip())
    s = re.sub(r"<!--.*?-->", "", s, flags=re.S)
    m = re.search(r"<svg\b([^>]*)>", s, re.S)
    if not m:
        raise ValueError("not an SVG figure")
    vb = re.search(r'viewBox="([^"]+)"', m.group(1))
    viewbox = vb.group(1) if vb else f"0 0 {size} {size}"
    body = s[m.end(): s.rindex("</svg>")]
    return {"id": item_id, "type": "svgart", "x": float(x), "y": float(y), "size": float(size), "viewBox": viewbox, "body": body}


def add_svgart(scene: dict, svg_text: str, item_id: str, x: float = 40, y: float = 40, size: float = 240) -> dict:
    item = svgart_from_svg(svg_text, item_id, x, y, size)
    scene.setdefault("items", []).append(item)
    return item


def icon_svg(icon, prefix: str) -> str:
    """The SVG of an icon made by a generator (xy data, molecule, crystal structure)."""
    from mfs.reticon import core
    return core.render_svg(icon, core.Style(mode="solid", color="#1E293B"), 256, id_prefix=prefix, standalone=True)


def parse_xy(data: bytes):
    """Two numeric columns (x, y) from a text / CSV file. Header lines and comments are skipped."""
    import numpy as np
    rows = []
    for line in data.decode("utf-8", "replace").splitlines():
        t = line.strip()
        if not t or t.startswith("#"):
            continue
        try:
            rows.append([float(p) for p in re.split(r"[\s,;\t]+", t)[:2]])
        except ValueError:
            continue
    if len(rows) < 3:
        raise ValueError("need at least three rows of two numbers (x and y)")
    a = np.array(rows, float)
    return a[:, 0], a[:, 1]


def figure_panel(png: bytes):
    """A composer panel that draws a rendered scheme (raster, at the export dpi) into one panel of a combined figure."""
    import matplotlib.image as mpimg
    img = mpimg.imread(io.BytesIO(png))

    def draw(ax):
        ax.imshow(img)
        ax.set_axis_off()
    return draw
