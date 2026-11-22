"""v2.11 global figure style: one set of fonts / sizes / line widths / palette / width for every page, with journal presets.
Values live in st.session_state['gs_*'] (saved in project files) and in <data dir>/figure_style.json (kept across starts)."""
from __future__ import annotations

import json
from pathlib import Path

import streamlit as st
from mfs import numentry as NE  # every slider takes a typed number too (2.19.6)

from . import style as S

DEFAULTS = {"gs_preset": "(custom)", "gs_width": list(S.WIDTHS)[0], "gs_aspect": 0.85, "gs_font": 8.0, "gs_family": "auto",
            "gs_lw": 1.0, "gs_alw": 0.75, "gs_ms": 3.5, "gs_palette": "Okabe–Ito (sample families)", "gs_tickdir": "in",
            "gs_minor": True, "gs_legframe": False, "gs_dpi": 600}
_MAP = {"gs_font": "font_size", "gs_family": "family", "gs_lw": "line_width", "gs_alw": "axes_width", "gs_ms": "marker_size",
        "gs_palette": "palette", "gs_tickdir": "tick_dir", "gs_minor": "minor", "gs_legframe": "legend_frame"}
_PATH = {"p": None}


def values(ss=None):
    ss = st.session_state if ss is None else ss
    return {k: ss.get(k, v) for k, v in DEFAULTS.items()}


def _save():
    p = _PATH["p"]
    if not p:
        return
    try:
        Path(p).parent.mkdir(parents=True, exist_ok=True)
        Path(p).write_text(json.dumps(values(), indent=1), encoding="utf-8")
    except Exception:
        pass


def _on_preset():
    ss = st.session_state
    pr = S.JOURNAL_PRESETS.get(ss.get("gs_preset"))
    if pr:
        ss["gs_width"] = pr["width"]; ss["gs_font"] = float(pr["font_size"]); ss["gs_family"] = pr.get("family", "auto")
        ss["gs_lw"] = float(pr["line_width"]); ss["gs_alw"] = float(pr["axes_width"]); ss["gs_ms"] = float(pr.get("marker_size", 3.5))
    _save()


def _reset():
    for k, v in DEFAULTS.items():
        st.session_state[k] = v
    _save()


def apply(vals: dict):
    S.set_global(**{_MAP[k]: vals[k] for k in _MAP if k in vals})


def sidebar(path=None):
    """Draw the '🎨 Figure style' box; returns (WIDTH in, ASPECT, FONT pt, DPI) and applies the style to matplotlib."""
    ss = st.session_state
    _PATH["p"] = path
    if "gs_font" not in ss:                                   # first run of the session: last used style from disk
        try:
            if path and Path(path).exists():
                d = json.loads(Path(path).read_text(encoding="utf-8"))
                for k, v in d.items():
                    if k in DEFAULTS and (k != "gs_width" or v in S.WIDTHS) and (k != "gs_palette" or v in S.PALETTES):
                        ss[k] = v
        except Exception:
            pass
    for k, v in DEFAULTS.items():
        ss.setdefault(k, v)
    box = st.sidebar.expander("🎨 Figure style – all pages (journal presets)", expanded=False)
    box.selectbox("Journal preset", ["(custom)"] + list(S.JOURNAL_PRESETS), key="gs_preset", on_change=_on_preset,
                  help="Sets width, font, font size and line widths from the journal's artwork guidelines. Every value stays editable.")
    box.selectbox("Width", list(S.WIDTHS), key="gs_width", on_change=_save)
    NE.slider("Height / width", 0.4, 1.6, key="gs_aspect", step=0.05, on_change=_save)
    c1, c2 = box.columns(2)
    c1.number_input("Font size (pt)", 5.0, 24.0, key="gs_font", step=0.5, on_change=_save)
    c2.selectbox("Font", S.FONT_FAMILIES, key="gs_family", on_change=_save,
                 help="'auto' = Arial / Helvetica / Liberation Sans, whichever is installed. Missing fonts fall back to a metric-compatible one.")
    c1, c2, c3 = box.columns(3)
    c1.number_input("Lines (pt)", 0.3, 5.0, key="gs_lw", step=0.1, on_change=_save)
    c2.number_input("Axes (pt)", 0.3, 3.0, key="gs_alw", step=0.05, on_change=_save)
    c3.number_input("Markers", 1.0, 12.0, key="gs_ms", step=0.5, on_change=_save)
    box.selectbox("Palette", list(S.PALETTES), key="gs_palette", on_change=_save,
                  help="Default: fixed colours for the known sample families, then Okabe–Ito (colour-blind safe). Colours set per sample "
                       "(sample list swatches) always win.")
    sw = S.PALETTES.get(ss["gs_palette"]) or S.palette()
    box.markdown("".join(f"<span style='display:inline-block;width:16px;height:12px;margin-right:2px;background:{c};"
                         f"border:1px solid #8884'></span>" for c in sw), unsafe_allow_html=True)
    c1, c2, c3 = box.columns(3)
    c1.selectbox("Ticks", ["in", "out", "inout"], key="gs_tickdir", on_change=_save)
    c2.checkbox("Minor ticks", key="gs_minor", on_change=_save)
    c3.checkbox("Legend frame", key="gs_legframe", on_change=_save)
    box.selectbox("Raster dpi", [300, 600, 1200], key="gs_dpi", on_change=_save)
    box.button("↺ Reset to the default style", key="gs_reset", on_click=_reset)
    vals = values()
    apply(vals)
    return S.WIDTHS[vals["gs_width"]], float(vals["gs_aspect"]), float(vals["gs_font"]), int(vals["gs_dpi"])
