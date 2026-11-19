"""v2.16: the Home page - a guided start for any user: hero, 'what do you want to do?' cards, the five steps, one coloured tile per technique
(grouped, each with what you get out of it), the three detail levels and the gallery. Pure presentation."""
from __future__ import annotations

import os

import streamlit as st

from mfs import ui as UI

GROUP_COLOUR = {"Structure": "#6366F1", "Porosity": "#0EA5E9", "Spectroscopy": "#A855F7", "Thermal": "#F97316", "Microscopy": "#10B981",
                "Kinetics / Bio": "#EC4899", "Tools": "#64748B"}

# page -> (what you get)
WHAT = {"PXRD": "Stacked patterns, CIF reference, crystallite size, Rietveld refinement",
        "Structure (CIF)": "Simulated patterns and structure tables from a CIF",
        "N2 sorption (BET/PSD)": "BET (Rouquerol), pore size distribution, t-plot, isotherms",
        "FTIR / Raman": "Offset spectra, band assignment, multi-band fits",
        "UV/Vis": "Spectra, Tauc band gap, calibration curves",
        "Photoluminescence": "Emission, thermal quenching, luminescence thermometry",
        "1H NMR (digest)": "Assign spectra, digestion ratios, qNMR purity",
        "XPS": "Peak fits with doublets, atomic percent",
        "ICP-OES / XRF": "Compositions and metal ratios",
        "TGA": "Mass loss, DTG, residues for the defect calculator",
        "DSC": "Onset / peak / ΔH, glass transition, heat capacity",
        "SEM / TEM / EDX": "Particle sizing, scale bars, EDX maps",
        "Kinetics & adsorption": "Rate constants, TOF, isotherm fits",
        "SDS-PAGE (Image Lab)": "Gel band densitometry",
        "Zeta / DLS": "Size and ζ-potential, choose runs",
        "Defect calculator": "Defects of any MOF node, every method side by side",
        "Figure composer": "Several techniques in ONE figure (A, B, C …)",
        "Batch processing": "Many files, one preset, one ZIP",
        "SI export bundle": "Figures, data, methods text, cross-checks, dossiers",
        "Sample database": "Link samples and their data files",
        "Export guide": "How to export from each instrument",
        "Help and info": "An info PDF for every module: terms, tutorial, formulas, checks"}

INTENTS = [
    ("📈", "I have raw data and want a figure", "Pick the technique, drop the file, adjust, export a journal-ready figure.",
     [("PXRD", "💠 PXRD"), ("TGA", "🔥 TGA"), ("N2 sorption (BET/PSD)", "🫧 N₂"), ("FTIR / Raman", "〰️ FTIR")], "#6366F1"),
    ("🧩", "I want several techniques in one figure", "Collect graphs from different pages and arrange them as A, B, C.",
     [("Figure composer", "🧩 Combined figure"), ("Batch processing", "🗂️ Batch")], "#A855F7"),
    ("🔬", "I want to know what my material is", "Phases, structure and refinement; particle size; composition; defects.",
     [("PXRD", "💠 PXRD & Rietveld"), ("SEM / TEM / EDX", "🔬 Particles"), ("Defect calculator", "🧮 Defects")], "#10B981"),
    ("📦", "I need the numbers and the SI for my paper", "Methods text, sample summary with sources, cross-checks, one dossier per material.",
     [("SI export bundle", "📦 SI export"), ("Sample database", "🗃️ Samples")], "#F97316"),
]

HOME_CSS = """
<style>
[data-testid="stAppViewContainer"] { background: linear-gradient(180deg,#EEF2FF 0%,#F8FAFF 380px,#FFFFFF 900px); }
.hm-hero { position:relative; overflow:hidden; border-radius:22px; padding:34px 38px 30px; color:#F8FAFC;
  background: radial-gradient(900px 340px at 85% -10%, rgba(94,234,212,.35), transparent 60%), radial-gradient(700px 300px at 5% 110%, rgba(99,102,241,.55), transparent 60%),
              linear-gradient(135deg,#0B1F3A 0%,#16315F 55%,#1E3A8A 100%);
  box-shadow:0 18px 40px rgba(11,31,58,.28); margin-bottom:22px; }
.hm-hero h1 { margin:0 0 6px; font-size:2.35rem; font-weight:800; letter-spacing:-.02em; color:#fff; line-height:1.1; }
.hm-hero h1 span { background:linear-gradient(90deg,#5EEAD4,#A5B4FC); -webkit-background-clip:text; background-clip:text; color:transparent; }
.hm-hero p { margin:6px 0 16px; font-size:1.08rem; max-width:760px; color:#DCE7F7; line-height:1.5; }
.hm-pill { display:inline-block; margin:0 8px 6px 0; padding:5px 13px; border-radius:999px; font-size:.82rem; font-weight:600;
  background:rgba(255,255,255,.12); border:1px solid rgba(255,255,255,.28); color:#F1F5F9; }
.hm-h { font-size:1.15rem; font-weight:800; color:#0B1F3A; margin:26px 0 10px; display:flex; align-items:center; gap:10px; }
.hm-h small { font-weight:500; color:#64748B; font-size:.86rem; }
.hm-steps { display:grid; grid-template-columns:repeat(5,1fr); gap:12px; margin:6px 0 4px; }
.hm-step { position:relative; background:#fff; border-radius:16px; padding:16px 14px 14px; border:1px solid #E2E8F0; box-shadow:0 4px 14px rgba(15,23,42,.06); }
.hm-step .n { width:30px; height:30px; border-radius:50%; display:flex; align-items:center; justify-content:center; font-weight:800; color:#fff;
  background:linear-gradient(135deg,#6366F1,#0EA5E9); margin-bottom:8px; }
.hm-step .t { font-weight:700; color:#0B1F3A; margin-bottom:3px; }
.hm-step .d { font-size:.86rem; color:#475569; line-height:1.4; }
div[class*="st-key-hm_intent_"] { border-radius:18px; padding:16px 18px 12px; border:1px solid #E2E8F0; background:#fff; box-shadow:0 6px 18px rgba(15,23,42,.07); }
div[class*="st-key-hm_tile_"] { border-radius:16px; padding:12px 14px 8px; background:#fff; border:1px solid #E2E8F0; box-shadow:0 3px 10px rgba(15,23,42,.05); min-height:132px; }
div[class*="st-key-hm_tile_"]:hover { box-shadow:0 8px 22px rgba(15,23,42,.13); transform:translateY(-1px); transition:all .15s; }
div[class*="st-key-hm_lvl_"] { border-radius:16px; padding:14px 16px 10px; background:#fff; border:2px solid #E2E8F0; }
.hm-ic { font-size:1.5rem; line-height:1; }
.hm-tt { font-weight:800; color:#0B1F3A; font-size:1.02rem; margin:2px 0 2px; }
.hm-tx { color:#475569; font-size:.88rem; line-height:1.4; margin-bottom:6px; }
.hm-tag { display:inline-block; font-size:.68rem; font-weight:800; letter-spacing:.09em; text-transform:uppercase; padding:2px 8px; border-radius:6px; margin-bottom:4px; color:#fff; }
.hm-foot { margin:26px 0 6px; padding:16px 20px; border-radius:14px; background:#0B1F3A; color:#DCE7F7; font-size:.92rem; line-height:1.5; }
.hm-foot b { color:#5EEAD4; }
@media (max-width: 900px) { .hm-steps { grid-template-columns:1fr 1fr; } }
</style>
"""


def _tile_css(key: str, colour: str) -> str:
    return (f'<style>div.st-key-{key}{{border-left:5px solid {colour} !important;'
            f'background:linear-gradient(135deg,{colour}14 0%,#FFFFFF 55%) !important}}</style>')


def _project_box(version: str):
    """v2.18: continue a saved project (.reticula) straight from the home page, or save the current work to one."""
    from mfs import project as PRJ
    ss = st.session_state
    with st.container(border=True):
        st.markdown("**📂 Continue a project** · **💾 Save your work**")
        c1, c2 = st.columns(2, gap="large")
        with c1:
            up = st.file_uploader("Open a .reticula project", type=["reticula", "reticula", "zip"], key="home_proj_up",
                                  help="Replaces the current work with the saved project: every page, uploaded file and setting.")
            if up is not None and st.button(f"📂 Open {up.name}", key="home_proj_open", type="primary", width="stretch"):
                try:
                    PRJ.request_open(up.getvalue(), up.name)
                    st.rerun()
                except Exception as e:
                    st.error(f"{up.name}: {e}")
        with c2:
            nm = st.text_input("File name", ss.get("_proj_name", "my_project"), key="home_proj_name")
            frozen = {k: ss[k] for k in list(ss.keys())}      # the download runs outside the script thread
            cur = ss.get("nav")
            st.download_button("💾 Save as .reticula", data=lambda: PRJ.export_bytes(cur, version, frozen),
                               file_name=(nm or "project").strip().replace(" ", "_") + PRJ.EXT, mime="application/zip",
                               key="home_proj_save", width="stretch",
                               help="Every page's settings, uploaded files, sample store and figure style in one file.")


def render(here: str, version: str = ""):
    from mfs import levels as LV
    st.markdown(HOME_CSS, unsafe_allow_html=True)
    cur = LV.level_name("Home")
    st.markdown(f"""<div class="hm-hero">
<h1>Reticula<span> · journal-ready figures</span></h1>
<p>From the raw export of your instrument to a publication figure, the numbers for your Methods and a dossier per material. Powder diffraction, sorption,
spectroscopy, thermal analysis, microscopy, NMR and more – with the journal style, the statistics and the provenance built in.</p>
<span class="hm-pill">🔒 Runs on your computer</span><span class="hm-pill">🎨 ACS / Wiley-VCH style</span><span class="hm-pill">🧪 {len(UI.PAGES) - 1} tools</span>
<span class="hm-pill">🧾 Every figure records how it was made</span><span class="hm-pill">v{version}</span></div>""", unsafe_allow_html=True)

    _project_box(version)

    # ------------------------------------------------------------------ what do you want to do
    st.markdown('<div class="hm-h">👉 What do you want to do? <small>pick a starting point</small></div>', unsafe_allow_html=True)
    cols = st.columns(2)
    for i, (ic, title, text, buttons, colour) in enumerate(INTENTS):
        key = f"hm_intent_{i}"
        with cols[i % 2].container(key=key):
            st.markdown(_tile_css(key, colour) + f'<div class="hm-ic">{ic}</div><div class="hm-tt">{title}</div><div class="hm-tx">{text}</div>', unsafe_allow_html=True)
            bc = st.columns(len(buttons))
            for j, (page, label) in enumerate(buttons):
                bc[j].button(label, key=f"hm_i_{i}_{j}", on_click=UI._go, args=(page,), width="stretch")

    # ------------------------------------------------------------------ five steps
    st.markdown('<div class="hm-h">🧭 How it works <small>the same five steps on every page</small></div>', unsafe_allow_html=True)
    st.markdown("""<div class="hm-steps">
<div class="hm-step"><div class="n">1</div><div class="t">Pick a technique</div><div class="d">Grouped in the left menu: structure, porosity, spectroscopy, thermal, microscopy, kinetics / bio, tools.</div></div>
<div class="hm-step"><div class="n">2</div><div class="t">Load your data</div><div class="d">Drop the instrument export in the <b>📁 Data</b> box – or tick the bundled EXAMPLE data to try it.</div></div>
<div class="hm-step"><div class="n">3</div><div class="t">Fine-tune</div><div class="d"><b>🎛️ Customise this figure</b>: title, axes, legend, colours, presets. Everything has a sensible default.</div></div>
<div class="hm-step"><div class="n">4</div><div class="t">Export</div><div class="d">Under each figure: PNG / TIFF 600 dpi, SVG, PDF, CSV or a ZIP – with the settings stored inside.</div></div>
<div class="hm-step"><div class="n">5</div><div class="t">Combine &amp; report</div><div class="d"><b>🧩 Add to combined figure</b>, then <b>📦 SI export</b>: methods text, sample summary, checks, dossier.</div></div>
</div>""", unsafe_allow_html=True)

    # ------------------------------------------------------------------ the tools, grouped
    st.markdown('<div class="hm-h">🧰 All tools <small>click a tile to open it</small></div>', unsafe_allow_html=True)
    groups: dict = {}
    cur_g = None
    for page, icon, g in UI.NAV:
        if page == "Home":
            continue
        cur_g = g or cur_g
        groups.setdefault(cur_g, []).append((page, icon))
    n = 0
    for g, items in groups.items():
        colour = GROUP_COLOUR.get(g, "#64748B")
        st.markdown(f'<span class="hm-tag" style="background:{colour}">{g}</span>', unsafe_allow_html=True)
        tc = st.columns(4)
        for k, (page, icon) in enumerate(items):
            key = f"hm_tile_{n}"
            n += 1
            with tc[k % 4].container(key=key):
                st.markdown(_tile_css(key, colour) + f'<div class="hm-ic">{icon}</div><div class="hm-tt">{UI.DISPLAY.get(page, page)}</div><div class="hm-tx">{WHAT.get(page, "")}</div>',
                            unsafe_allow_html=True)
                st.button("Open →", key=f"hm_open_{n}", on_click=UI._go, args=(page,), width="stretch")

    # ------------------------------------------------------------------ the three levels
    st.markdown(f'<div class="hm-h">🎚️ How much do you want to see? <small>now: {cur} – switch it any time in the box at the top of the left menu</small></div>', unsafe_allow_html=True)
    lc = st.columns(3)
    info = {"Standard": ("🟢", "Standard", "Only what a publishable figure needs: load data, pick samples, the main settings, export. The default."),
            "Advanced": ("🔵", "Advanced", "Adds the details: fit options, corrections, ticks and presets, libraries, quality checks."),
            "Developer": ("🟣", "Developer", "Everything, including diagnostics, raw tables and internal parameters.")}
    for i, (lv, (ic, name, text)) in enumerate(info.items()):
        key = f"hm_lvl_{i}"
        col = LV.LEVEL_COLOUR[lv]
        with lc[i].container(key=key):
            st.markdown(f'<style>div.st-key-{key}{{border-color:{col if lv == cur else "#E2E8F0"} !important;'
                        f'background:{col + "12" if lv == cur else "#fff"} !important}}</style>'
                        f'<div class="hm-ic">{ic}</div><div class="hm-tt">{name}' + (" · current" if lv == cur else "") + f'</div><div class="hm-tx">{text}</div>', unsafe_allow_html=True)
            st.button(f"Use {name}", key=f"hm_setlvl_{i}", on_click=_set_level, args=(lv,), disabled=lv == cur, width="stretch")

    # ------------------------------------------------------------------ gallery
    st.markdown('<div class="hm-h">🖼️ Gallery <small>figures made from the synthetic EXAMPLE data – click to open the page</small></div>', unsafe_allow_html=True)
    cols = st.columns(4)
    for i, (img, page, cap) in enumerate(UI.GALLERY):
        with cols[i % 4].container(border=True, key=f"gal_{i}"):
            p = os.path.join(here, "assets", "gallery", img + ".jpg")
            if os.path.exists(p):
                st.image(p, width="stretch")
            st.caption(cap)
            st.button(f"{UI.ICON.get(page, '')} Open {UI.DISPLAY.get(page, page)}", key=f"gal_btn_{i}", on_click=UI._go, args=(page,), width="stretch")

    st.markdown("""<div class="hm-foot">💡 <b>Tips</b> · Everything stays on your computer. · Under every graph: <b>🧩 Add to combined figure</b>. ·
Something wrong? Open <b>🐞 Report a problem</b> in the left menu – it makes a small file without your data. · Your work is kept per page; <b>↶ Undo</b> is in the left menu.</div>""",
                unsafe_allow_html=True)


def _set_level(level: str):
    from mfs import levels as LV
    st.session_state.setdefault("_lv", {"global": "Standard"})["global"] = level
    st.session_state["lv_global"] = level
