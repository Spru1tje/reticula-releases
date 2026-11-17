"""Reticula look & feel: theme CSS, branded header, grouped technique navigation, home page with quick start + gallery.
Pure presentation - no analysis code here."""
from __future__ import annotations
import base64, functools, os
import streamlit as st

NAME = "Reticula"
TAGLINE = "Publication figures for framework materials"
TEAL, INDIGO = "#0F766E", "#312E81"

# (page, icon, group) - order defines the sidebar; group headers are drawn by CSS before the first page of a group
NAV = [("Home", "🏠", None),
       ("PXRD", "💠", "Structure"), ("Structure (CIF)", "💎", None),
       ("N2 sorption (BET/PSD)", "🫧", "Porosity"),
       ("FTIR / Raman", "〰️", "Spectroscopy"), ("UV/Vis", "🌈", None), ("1H NMR (digest)", "🧲", None), ("XPS", "⚛️", None), ("ICP-OES / XRF", "🧪", None),
       ("TGA", "🔥", "Thermal"), ("DSC", "🌡️", None),
       ("SEM / TEM / EDX", "🔬", "Microscopy"),
       ("Kinetics & adsorption", "⏱️", "Kinetics / Bio"), ("SDS-PAGE (Image Lab)", "🧬", None), ("Zeta / DLS", "⚡", None),
       ("Defect calculator", "🧮", "Tools"), ("Figure composer", "🧩", None), ("Batch processing", "🗂️", None),
       ("SI export bundle", "📦", None), ("Sample database", "🗃️", None), ("Export guide", "📘", None)]
PAGES = [p for p, _, _ in NAV]
ICON = {p: i for p, i, _ in NAV}

GALLERY = [("01_PXRD_stacked", "PXRD", "Stacked patterns + simulated CIF reference"),
           ("04_N2_isotherms_PSD", "N2 sorption (BET/PSD)", "N₂ isotherms, BET (Rouquerol) and PSD"),
           ("02_FTIR_stacked", "FTIR / Raman", "Offset FTIR / Raman spectra"),
           ("03_NMR_digest_stacked", "1H NMR (digest)", "¹H NMR: assignment and digested MOFs"),
           ("06_TGA_DTG", "TGA", "TGA + DTG, ash-normalised"),
           ("16_DSC", "DSC", "DSC: onset / peak / ΔH, Tg, exo arrow"),
           ("14_defects_MOF801", "Defect calculator", "Every defect method side by side"),
           ("10_SEM_panel", "SEM / TEM / EDX", "SEM/TEM panels with editable scale bars"),
           ("07_GlyGly_kinetics", "Kinetics & adsorption", "GlyGly hydrolysis kinetics"),
           ("13_XPS_Hf4f_Zr3d", "XPS", "Hf 4f / Zr 3d fits"),
           ("09_ICP_bars", "ICP-OES / XRF", "Metal / S content bars"),
           ("11_SDS-PAGE_bands", "SDS-PAGE (Image Lab)", "Protein cleavage from gel bands"),
           ("12_Zeta_DLS", "Zeta / DLS", "Zeta potential and DLS")]


@functools.lru_cache(None)
def _b64(path):
    try:
        with open(path, "rb") as fh:
            return base64.b64encode(fh.read()).decode()
    except OSError:
        return ""


def _group_css():
    rules = []
    for i, (_, _, g) in enumerate(NAV, start=1):
        if g:
            rules.append(f'section[data-testid="stSidebar"] div[role="radiogroup"] > div:nth-child({i})::before'
                         f'{{content:"{g.upper()}";}}')
    return "\n".join(rules)


CSS = """
<style>
:root { --rt-teal:%(teal)s; --rt-indigo:%(indigo)s; }
.block-container { padding-top: 1.2rem; max-width: 1350px; }
/* sidebar */
section[data-testid="stSidebar"] div[role="radiogroup"] { gap:0 !important; }
section[data-testid="stSidebar"] div[role="radiogroup"] > div { display:block; }
section[data-testid="stSidebar"] div[role="radiogroup"] > div > label { display:flex; width:100%%; padding:5px 10px; border-radius:8px; margin:1px 0; }
section[data-testid="stSidebar"] div[role="radiogroup"] > div > label:hover { background:rgba(94,234,212,.10); }
section[data-testid="stSidebar"] div[role="radiogroup"] > div[data-selected="true"] > label { background:rgba(94,234,212,.18); box-shadow: inset 3px 0 0 #5EEAD4; font-weight:600; }
section[data-testid="stSidebar"] div[role="radiogroup"] > div::before { display:block; font-size:.66rem; letter-spacing:.12em;
   font-weight:700; color:#7DD3C8; margin:12px 0 3px 8px; }
section[data-testid="stSidebar"] label[data-testid="stRadioOption"] > div > div:first-child { display:none; }  /* hide radio dots */
div[class*="st-key-exportbar_"] [data-testid="stDownloadButton"], div[class*="st-key-exportbar_"] [data-testid="stButton"] { width:100%%; }
.rt-side-brand { display:flex; align-items:center; gap:10px; margin:-8px 0 6px 0; }
.rt-side-brand img { width:42px; height:42px; border-radius:10px; }
.rt-side-brand .by { font-size:.72rem; color:#9FB4CF; margin-top:1px; }
.rt-side-brand .n { font-size:1.35rem; font-weight:800; letter-spacing:.02em; color:#fff; line-height:1.1; }
.rt-side-brand .t { font-size:.72rem; color:#9FB4CF; line-height:1.2; }
/* header */
.rt-header { display:flex; align-items:center; gap:16px; padding:14px 20px; margin:0 0 18px 0; border-radius:16px;
   background: linear-gradient(110deg, #0B1F3A 0%%, var(--rt-indigo) 55%%, var(--rt-teal) 100%%); color:#fff;
   box-shadow:0 6px 22px rgba(15,40,80,.18); }
.rt-header img { width:54px; height:54px; border-radius:12px; }
.rt-header .n { font-size:1.6rem; font-weight:800; line-height:1.05; }
.rt-header .t { font-size:.9rem; opacity:.85; }
.rt-header .chip { margin-left:auto; background:rgba(255,255,255,.14); border:1px solid rgba(255,255,255,.25); padding:6px 14px;
   border-radius:999px; font-size:.9rem; white-space:nowrap; }
/* cards */
[data-testid="stExpander"] details { border-radius:12px; border-color:#DCE3EE; background:#fff; }
[data-testid="stExpander"] summary { font-weight:600; }
div[class*="st-key-data_"] { background:#fff; border-radius:14px; }
.rt-sec { font-weight:700; color:var(--rt-indigo); font-size:.95rem; margin:2px 0 6px 0; }
/* export bar */
div[class*="st-key-exportbar_"] { background:linear-gradient(90deg,#ECFDF9 0%%,#EEF2FF 100%%); border:1px solid #99D5CC !important;
   border-radius:14px; padding:6px 10px 4px 10px; margin:6px 0 14px 0; }
div[class*="st-key-exportbar_"] button { width:100%%; font-weight:600; border-radius:10px; border:1px solid #0F766E55; background:#fff; }
div[class*="st-key-exportbar_"] [data-testid="stDownloadButton"] button:hover { background:var(--rt-teal); color:#fff; border-color:var(--rt-teal); }
div[class*="st-key-exportbar_"] [data-testid="stButton"] button { background:var(--rt-teal); color:#fff; }
/* home */
.rt-hero { padding:26px 30px; border-radius:18px; background:#fff; border:1px solid #DCE3EE; margin-bottom:16px; }
.rt-hero h2 { margin:0 0 6px 0; color:#0B1F3A; }
.rt-tiles { display:grid; grid-template-columns:repeat(4,1fr); gap:14px; margin:6px 0 22px 0; }
.rt-tile { background:#fff; border:1px solid #DCE3EE; border-radius:14px; padding:16px 16px 12px 16px; border-top:4px solid var(--rt-teal); }
.rt-tile .k { font-size:.75rem; font-weight:700; color:var(--rt-teal); letter-spacing:.1em; }
.rt-tile .h { font-weight:700; color:#0B1F3A; margin:4px 0; }
.rt-tile .b { font-size:.86rem; color:#475569; }
div[class*="st-key-gal_"] { background:#fff; border-radius:14px; }
div[class*="st-key-gal_"] img { border-radius:8px; }
</style>
""" % {"teal": TEAL, "indigo": INDIGO}


def inject(here):
    st.markdown(CSS.replace("</style>", _group_css() + "\n</style>"), unsafe_allow_html=True)


def sidebar_brand(here, version):
    b = _b64(os.path.join(here, "assets", "reticula_icon.png"))
    st.sidebar.markdown(f'<div class="rt-side-brand"><img src="data:image/png;base64,{b}"/><div><div class="n">{NAME}</div>'
                        f'<div class="t">{TAGLINE} · v{version}</div><div class="by">Built by Stijn Germonpré</div></div></div>',
                        unsafe_allow_html=True)


def about_box(version, data_dir=None):
    with st.sidebar.expander("ℹ️ About", expanded=False):
        st.markdown(f"**{NAME}** v{version}  \nBuilt by Stijn Germonpré")
        st.caption(TAGLINE + ". Runs locally; your data never leaves your computer (the optional update check only reads the public release list).")
        if data_dir:
            from mfs import update as UPD
            UPD.settings_box(data_dir, version)


def header(here, page):
    b = _b64(os.path.join(here, "assets", "reticula_icon.png"))
    st.markdown(f'<div class="rt-header"><img src="data:image/png;base64,{b}"/><div><div class="n">{NAME}</div>'
                f'<div class="t">{TAGLINE}</div></div><div class="chip">{ICON.get(page, "")} {DISPLAY.get(page, page)}</div></div>',
                unsafe_allow_html=True)


DISPLAY = {"1H NMR (digest)": "¹H NMR"}      # v2.13: page key kept (projects, saved state), label changed


def nav_label(p):
    return f"{ICON.get(p, '')}  {DISPLAY.get(p, p)}"


def _go(p):
    from mfs import pagestate as PS
    PS.request(p)          # goes through the unsaved-work guard


def home(here):
    st.markdown("""<div class="rt-hero"><h2>Journal-ready figures from raw instrument exports</h2>
ACS / Wiley-VCH style out of the box (Arial 7–9 pt, 3.25 in or 7 in wide, inward ticks, colour-blind-safe palette, units in parentheses),
plus a defect calculator for any MOF node (Zr₆, Hf₆, Ce₆, Cu paddlewheel, Zn₄O, Al-OH, M₃O, MOF-74, Ti₈, ZIF …) that runs every established method side by side. Everything stays on your computer.</div>""",
                unsafe_allow_html=True)
    st.markdown("""<div class="rt-tiles">
<div class="rt-tile"><div class="k">STEP 1</div><div class="h">Pick a technique</div><div class="b">Grouped on the left: structure, porosity, spectroscopy, thermal, microscopy, kinetics/bio and tools.</div></div>
<div class="rt-tile"><div class="k">STEP 2</div><div class="h">Load your data</div><div class="b">Drop the instrument export in <b>📁 Data</b>, or tick the bundled EXAMPLE data to try it.</div></div>
<div class="rt-tile"><div class="k">STEP 3</div><div class="h">Fine-tune</div><div class="b"><b>🎛️ Customise this figure</b>: title, axes &amp; ticks, legend, colours, presets.</div></div>
<div class="rt-tile"><div class="k">STEP 4</div><div class="h">Export</div><div class="b">The export bar under each figure: PNG/TIFF 600 dpi, SVG, PDF, CSV, ZIP or 💾 Save to folder.</div></div>
</div>""", unsafe_allow_html=True)
    st.markdown('<div class="rt-sec">Technique gallery · click to open (figures made from the synthetic EXAMPLE data)</div>',
                unsafe_allow_html=True)
    cols = st.columns(4)
    for i, (img, page, cap) in enumerate(GALLERY):
        with cols[i % 4].container(border=True, key=f"gal_{i}"):
            p = os.path.join(here, "assets", "gallery", img + ".jpg")
            if os.path.exists(p):
                st.image(p, width="stretch")
            st.caption(cap)
            st.button(f"{ICON.get(page, '')} Open {DISPLAY.get(page, page)}", key=f"gal_btn_{i}", on_click=_go, args=(page,), width="stretch")
