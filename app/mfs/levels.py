"""v2.15: interface detail level - Standard (default; the bare essentials for a publishable figure), Advanced (more details) and
Developer (everything, incl. diagnostics). One global level in the sidebar, overridable per module.

Use in a page:   from mfs import levels as LV
                 if LV.at_least("Advanced"): ...            # show only from Advanced up
                 with LV.expander("Fit details", "Advanced"): ...   # an expander that exists only from that level
                 LV.opt_in("Rietveld refinement", key)       # heavy tool: hidden until the user clicks it open (any level)
A control without a level belongs to Advanced by convention, so Standard stays clean when something new is added."""
from __future__ import annotations

import contextlib
import os

import streamlit as st

LEVELS = ("Standard", "Advanced", "Developer")
GLOBAL_KEY = "ui_level"
USE_GLOBAL = "global"
HELP = {"Standard": "Only what a publishable figure needs: load data, pick samples, the main settings, export.",
        "Advanced": "Adds the details: fit options, corrections, ticks, presets, quality checks.",
        "Developer": "Everything, including diagnostics, raw tables and internal parameters."}


def _page() -> str:
    return str(st.session_state.get("nav") or st.session_state.get("_ps_page") or "")


def _store() -> dict:
    return st.session_state.setdefault("_lv", {"global": os.environ.get("MFS_DEFAULT_LEVEL", "Standard")})


def level_name(page: str | None = None) -> str:
    d = _store()
    ov = d.get(page or _page(), USE_GLOBAL)
    lv = d.get("global", "Standard") if ov == USE_GLOBAL else ov
    return lv if lv in LEVELS else "Standard"


def rank(name: str | None = None) -> int:
    return LEVELS.index(name or level_name())


def at_least(minimum: str = "Advanced") -> bool:
    return rank() >= LEVELS.index(minimum)


@contextlib.contextmanager
def expander(label: str, minimum: str = "Advanced", expanded: bool = False):
    """An expander that exists only from level `minimum` up (otherwise its content is not drawn; yields False)."""
    if at_least(minimum):
        with st.expander(label, expanded=expanded):
            yield True
    else:
        yield False


def opt_in(label: str, key: str, help: str | None = None) -> bool:
    """A heavy / optional tool stays hidden until the user clicks it open (a toggle that remembers its state). Returns True when open."""
    return bool(st.toggle(label, value=False, key=f"optin|{key}", help=help))


LEVEL_COLOUR = {"Standard": "#16A34A", "Advanced": "#2563EB", "Developer": "#9333EA"}
LEVEL_ICON = {"Standard": "🟢", "Advanced": "🔵", "Developer": "🟣"}


def sidebar(page: str | None = None):
    """The level switch at the top of the sidebar: a coloured box with the three levels as buttons (global) and an override for the current module
    underneath (kept in st.session_state['_lv'], a project-global key)."""
    d = _store()
    w = st.session_state.get("lv_global")                    # the button pressed in THIS run is already in the widget state: use it for the colours
    if w in LEVELS:
        d["global"] = w
    wm = st.session_state.get(f"lv_mod|{page or _page()}")
    if wm in (USE_GLOBAL,) + LEVELS and (page or _page()) and (page or _page()) != "Home":
        d[page or _page()] = wm
    page = page or _page()
    cur = level_name(page)
    col = LEVEL_COLOUR[cur]
    st.sidebar.markdown(
        "<style>"
        ".st-key-lv_box{background:linear-gradient(135deg,#17386B 0%,#0F2A52 100%);border:2px solid " + col + ";border-radius:14px;"
        "padding:.7rem .85rem .55rem;margin:.2rem 0 .8rem;box-shadow:0 2px 10px rgba(0,0,0,.35)}"
        ".st-key-lv_box p,.st-key-lv_box label,.st-key-lv_box span{color:#E2E8F0}"
        ".st-key-lv_box [data-testid='stSegmentedControl'] button{font-weight:700}"
        ".lv-pill{display:inline-block;padding:2px 10px;border-radius:999px;color:#fff;font-weight:700;font-size:.8rem;vertical-align:middle;background:" + col + "}"
        ".st-key-lv_box .lv-title{font-weight:800;color:#FFFFFF;font-size:1rem}"
        "</style>", unsafe_allow_html=True)
    with st.sidebar.container(border=False, key="lv_box"):
        st.markdown(f"<span class='lv-title'>🎚️ Detail level</span>&nbsp;&nbsp;<span class='lv-pill'>{cur}</span>", unsafe_allow_html=True)
        v = st.segmented_control("All modules", LEVELS, default=d.get("global", "Standard"), key="lv_global", label_visibility="collapsed",
                                 format_func=lambda x: f"{LEVEL_ICON[x]} {x}", help="Standard = bare essentials; Advanced = more details; Developer = all functions.")
        if v in LEVELS:
            d["global"] = v
        st.caption(HELP[level_name(page)])
        if page and page != "Home":
            opts = (USE_GLOBAL,) + LEVELS
            cur_o = d.get(page, USE_GLOBAL)
            with st.expander(f"Only for this module: {page}", expanded=cur_o != USE_GLOBAL):
                d[page] = st.selectbox(f"Level of {page}", opts, index=opts.index(cur_o if cur_o in opts else USE_GLOBAL), key=f"lv_mod|{page}",
                                       format_func=lambda x: "same as all modules" if x == USE_GLOBAL else f"{LEVEL_ICON[x]} {x}")


# ------------------------------------------------------------------------------------------------ gated blocks
class _Ghost:
    """Stands in for `st` while a block above the user's detail level is run: widgets are not drawn and return their DEFAULT value, so the
    code after the block gets the same variables as before. Nothing of the block is shown."""
    _PASS = ("column_config", "session_state", "cache_data", "cache_resource", "fragment", "secrets", "query_params", "stop", "set_page_config")

    def __init__(self, real):
        self._real = real

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def __iter__(self):
        return iter(())

    def __getitem__(self, i):
        return self

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        if name in self._PASS:
            return getattr(self._real, name)
        if name in ("checkbox", "toggle"):
            return lambda label=None, value=False, *a, **k: k.get("value", value)
        if name == "number_input":
            def f(label=None, min_value=None, max_value=None, value="min", step=None, *a, **k):
                value = k.get("value", value)
                mn = k.get("min_value", min_value)
                return (mn if mn is not None else 0) if isinstance(value, str) and value == "min" else value
            return f
        if name == "slider":
            def f(label=None, min_value=0, max_value=100, value=None, *a, **k):
                value = k.get("value", value)
                return k.get("min_value", min_value) if value is None else value
            return f
        if name in ("selectbox", "radio", "select_slider"):
            def f(label=None, options=(), index=0, *a, **k):
                options = list(k.get("options", options))
                index = k.get("index", index)
                v = k.get("value")
                if v is not None:
                    return v
                return options[index] if options and index is not None and index < len(options) else None
            return f
        if name == "multiselect":
            return lambda label=None, options=(), default=None, *a, **k: list(k.get("default", default) or [])
        if name in ("text_input", "text_area", "color_picker"):
            return lambda label=None, value="", *a, **k: k.get("value", value)
        if name == "data_editor":
            return lambda data=None, *a, **k: data
        if name in ("button", "download_button", "form_submit_button", "link_button"):
            return lambda *a, **k: False
        if name in ("file_uploader", "camera_input", "date_input", "time_input"):
            return lambda *a, **k: None
        if name == "columns":
            return lambda spec=2, *a, **k: [self] * (spec if isinstance(spec, int) else len(spec))
        if name == "tabs":
            return lambda labels, *a, **k: [self] * len(labels)
        if name in ("expander", "container", "empty", "status", "popover", "form", "sidebar", "spinner"):
            return lambda *a, **k: self
        return lambda *a, **k: self      # markdown, caption, dataframe, info, image, pyplot … : nothing is shown


@contextlib.contextmanager
def gate(minimum: str, g: dict):
    """`with LV.gate("Advanced", globals()):` around a block: shown from that level up; below it the block still runs (so the variables it
    defines exist with their default values) but draws nothing: `st` is a ghost, show()/dl() do nothing."""
    if at_least(minimum):
        yield True
        return
    saved = {k: g.get(k) for k in ("st", "show", "dl")}
    g["st"] = _Ghost(saved["st"] or st)
    for k in ("show", "dl"):
        if k in g:
            g[k] = lambda *a, **k_: None
    try:
        yield False
    finally:
        for k, v in saved.items():
            if v is not None:
                g[k] = v
