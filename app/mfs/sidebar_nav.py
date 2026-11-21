"""v2.18.8: the technique navigation as collapsible categories (replaces the one long radio list).

- Each category is an expander; the category that holds the current page opens by itself.
- Every page is a button with a monochrome, scientific icon (assets/nav_icons/<name>.svg). The current page is highlighted.
- Navigation goes through pagestate.request(), so the unsaved-work guard and the 'nav' session key work as before.
- Single source of truth: st.session_state['nav']. The old radio widget with key 'nav' is gone; nothing else writes a widget under that key.
"""
from __future__ import annotations

import os

import streamlit as st

from mfs import pagestate as PS
from mfs import ui as UI

ICON_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "assets", "nav_icons")

# category label -> pages (every page of ui.NAV must appear exactly once; a test checks this)
CATEGORIES = [
    ("Start", ["Home"]),
    ("Structure and diffraction", ["PXRD", "Structure (CIF)"]),
    ("Porosity and sorption", ["N2 sorption (BET/PSD)"]),
    ("Spectroscopy and elements", ["FTIR / Raman", "UV/Vis", "Photoluminescence", "1H NMR (digest)", "XPS", "ICP-OES / XRF"]),
    ("Thermal", ["TGA", "DSC"]),
    ("Microscopy and particles", ["SEM / TEM / EDX", "Zeta / DLS"]),
    ("Kinetics and bio", ["Kinetics & adsorption", "SDS-PAGE (Image Lab)"]),
    ("Figures, data and tools", ["Defect calculator", "Figure composer", "Reticon", "Batch processing", "SI export bundle",
                                 "Sample database", "Export guide", "Help and info"]),
]

# page -> icon file (without .svg)
ICON_FILE = {"Home": "home", "PXRD": "pxrd", "Structure (CIF)": "structure_cif", "N2 sorption (BET/PSD)": "sorption",
             "FTIR / Raman": "ftir", "UV/Vis": "uvvis", "Photoluminescence": "photolum", "1H NMR (digest)": "nmr", "XPS": "xps",
             "ICP-OES / XRF": "icp", "TGA": "tga", "DSC": "dsc", "SEM / TEM / EDX": "sem", "Kinetics & adsorption": "kinetics",
             "SDS-PAGE (Image Lab)": "sds", "Zeta / DLS": "zeta", "Defect calculator": "defect", "Figure composer": "composer",
             "Reticon": "reticon", "Batch processing": "batch", "SI export bundle": "si_export", "Sample database": "sample_db",
             "Export guide": "export_guide", "Help and info": "help"}


def icon_path(page: str) -> str:
    name = ICON_FILE.get(page)
    return os.path.join(ICON_DIR, f"{name}.svg") if name else ""


def _label(page: str) -> str:
    return UI.DISPLAY.get(page, page)


def _pick(page: str):
    """on_click: go through the page-state guard (unsaved work asks first)."""
    PS.request(page)


def render(pages_in_nav: list[str]) -> str:
    """Draws the categories in the sidebar and returns the current page."""
    current = st.session_state.get("nav") or "Home"
    if current not in pages_in_nav:
        current = "Home"
    st.session_state["nav"] = current
    for cat, pages in CATEGORIES:
        has_current = current in pages
        with st.sidebar.expander(cat, expanded=has_current):
            for p in pages:
                c_icon, c_btn = st.columns([1, 6], vertical_alignment="center")
                path = icon_path(p)
                if path and os.path.isfile(path):
                    c_icon.image(path, width=20)
                c_btn.button(_label(p), key=f"nv_{p}", type="primary" if p == current else "tertiary",
                             on_click=_pick, args=(p,), width="stretch")
    return current
