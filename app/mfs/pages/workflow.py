"""Reticula page module: Figure composer, batch processing, SI bundle, sample database and export guide pages.

Moved verbatim out of app.py in v2.13. Names that the code takes from the app (helpers such as show(), dl(), files_input(),
the sidebar values WIDTH / ASPECT / FONT, the cached reader wrapper R …) are listed in INJECTED and filled in by
mfs.pages.bind(globals()) on every run of app.py, so the code behaves exactly as it did inside app.py."""
# ruff: noqa: F821
from __future__ import annotations
import io, os, re, json, glob, zipfile, warnings, datetime, hashlib
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

from mfs import analysis as A, plots as P, style as S, defects as D, composer as C, report as REP
from mfs import figctl as F, figctl_ui as FU, ui as UI, __version__
from mfs.chem import SPECIES, FRAMEWORKS, molar_mass

INJECTED = (
    'ASPECT',
    'DPI',
    'EX',
    'EXPORT_DIR',
    'FONT',
    'HERE',
    'STORE_PATH',
    'WIDTH',
    '_RENDER_LOCK',
    'dl',
    'open_with_default_app',
    'safe_name',
    'show',
)


def page_composer():
    st.header("🧩 Combined figure")
    st.caption("Bring the graphs of several techniques into ONE figure – for example (A) UV/Vis, (B) FTIR and (C) XRD of the same "
               "material – with panel letters, one font and one size, and a caption. Vector output (PDF / SVG) like every figure.")
    ss = st.session_state
    comp = ss.composer
    reg = ss.get("fig_registry", {})
    options = list(dict.fromkeys(list(reg) + list(comp)))
    if not options:
        st.info("**How it works**  \n1. Open a technique (PXRD, FTIR, UV/Vis, TGA …) and look at your graph.  \n"
                "2. Click **🧩 Add to combined figure** under it.  \n3. Do the same for the other techniques, then come back here "
                "(sidebar → *Open and arrange*): the panels appear below, ready to export.  \n\n"
                "Every graph you open in this session is remembered, so you can also just pick them here afterwards.")
        return

    def _lab(k):
        r = reg.get(k)
        return f"{r['page']} · {k}" if r else k
    sel = st.multiselect("1 · Panels – the order in which you pick them is the order in the figure", options,
                         default=[k for k in comp if k in options], format_func=_lab, key="comp_sel",
                         placeholder="Choose the graphs to combine …")
    if list(comp) != list(sel):                                   # keep the session list (and the sidebar tray) in step
        comp.clear()
        comp.update({k: (reg[k]["fn"] if k in reg else None) for k in sel})
    pages_seen = sorted({reg[k]["page"] for k in options if k in reg})
    if len(sel) == 0:
        st.info("Pick at least one graph above" + (f" (graphs from: {', '.join(pages_seen)})." if pages_seen else "."))
        return
    from mfs import templates as TP
    c1, c2, c3 = st.columns([1.4, 1.2, 1.2])
    lay = c1.radio("2 · Layout", ["Auto", "One row", "One column", "2 columns", "3 columns"], horizontal=True, key="comp_lay")
    ncols = {"Auto": C.auto_columns(len(sel)), "One row": len(sel), "One column": 1, "2 columns": 2, "3 columns": 3}[lay]
    ncols = max(1, min(ncols, len(sel)))
    lab = c2.selectbox("Panel letters", C.LABEL_STYLES, index=2, key="comp_lab", help="A, B, C … as in most papers")
    inside = c3.checkbox("Letters inside the panel", True, key="comp_inside")
    hide_long = c3.checkbox("Hide long legends (> 6 entries)", True, key="comp_hide_leg",
                            help="A legend with many samples takes over a small panel; the samples are then named by the colours / the caption.")
    with st.expander("⚙️ Size, journal template, font", expanded=False):
        t1, t2, t3 = st.columns([1.2, 1.2, 2.6])
        tpl = t1.selectbox("Journal template", ["custom"] + list(TP.TEMPLATES), key="comp_tpl")
        colw = t2.radio("Column", ["single", "1.5", "double"], index=2, horizontal=True, key="comp_col", disabled=tpl == "custom")
        if tpl != "custom":
            t3.caption(TP.TEMPLATES[tpl]["note"])
        c1, c2, c3 = st.columns(3)
        if tpl == "custom":
            width = S.WIDTHS[c1.selectbox("Width", list(S.WIDTHS), index=1)]
            fsz = c3.slider("Font (pt)", 5.0, 10.0, float(FONT), 0.5, key="comp_fs")
        else:
            T_ = TP.TEMPLATES[tpl]; width = TP.width(tpl, colw)
            fsz = c3.slider("Font (pt)", float(T_["min_font"]), 10.0, float(T_["font_size"]), 0.5, key=f"comp_fs_{tpl}")
            c1.metric("Width", f"{width:.2f} in / {width * 2.54:.1f} cm")
        asp = c2.slider("Panel height / width", 0.4, 1.5, 0.8, 0.05, key="comp_asp")
    letters = [C.panel_label(i, lab) for i in range(len(sel))]
    auto_cap = "Fig. X  " + ", ".join(f"{l} {k}" for l, k in zip(letters, sel)) + ". [material / sample]"
    cap_on = st.checkbox("3 · Caption under the figure", False, key="comp_cap_on")
    cap = st.text_area("Caption (edit freely)", auto_cap, height=80, key="comp_cap", label_visibility="collapsed" if not cap_on else "visible",
                       disabled=not cap_on)

    def _panel(k):
        r = reg.get(k) or {}
        fn, cfg, fams = r.get("fn") or comp.get(k), r.get("cfg"), r.get("fams")

        def draw(ax, fn=fn, cfg=cfg, fams=fams):
            fn(ax)
            if cfg:
                F.apply_panel(ax, cfg, fams, 0)                 # the settings made under the single figure travel with the panel
            lg = ax.get_legend()
            if hide_long and lg is not None and len(lg.get_texts()) > 6:
                lg.remove()
        return draw
    fig, axs = C.compose([_panel(k) for k in sel], ncols=ncols, width=width, panel_aspect=asp, label_style=lab, font_size=fsz,
                         label_inside=inside, caption=cap if cap_on else None)
    if tpl != "custom":
        for lvl, msg in TP.check(fig, tpl, colw):
            (st.success if lvl == "ok" else st.warning)(f"{tpl} check: {msg}")
    show(fig, f"combined_figure{'_' + tpl if tpl != 'custom' else ''}")
    if st.button("Clear the combined figure", key="comp_clear"):
        ss.composer = {}; ss.pop("comp_sel", None); st.rerun()


def page_guide():
    st.header("Per-instrument export guide")
    st.markdown(open(os.path.join(HERE, "docs", "EXPORT_GUIDE.md"), encoding="utf-8").read())


def page_batch():
    from mfs import batchproc as BP, templates as TP
    st.header("Batch processing")
    st.caption("All files of one folder (or an upload) → the standard figure of the technique with ONE preset → PNG/TIFF/PDF/SVG + CSV, "
               "zipped and optionally written to a folder. Preset = a figure-settings JSON saved under any figure (🎛️ → Presets) "
               "or the current house style.")
    c1, c2 = st.columns(2)
    kind = c1.selectbox("Technique", list(BP.KINDS), key="bp_kind")
    src = c2.radio("Files from", ["a folder on this computer", "upload", "the bundled EXAMPLE folder"], index=2, key="bp_src")
    exdir = {"PXRD": "PXRD", "FTIR": "FTIR", "Raman": "Raman", "XPS": "XPS", "TGA": "TGA", "DSC": "DSC", "N2 isotherm": "Sorption",
             "SEM/TEM particles": "SEM"}[kind]
    files = []
    if src.startswith("a folder"):
        folder = st.text_input("Folder (sub-folders included)", os.path.expanduser("~"), key="bp_folder")
        if os.path.isdir(folder):
            fl = BP.list_folder(folder, kind); files = [(p_, p_) for p_ in fl]
            st.caption(f"{len(fl)} {kind} file(s) found" + (": " + ", ".join(os.path.basename(p_) for p_ in fl[:8]) + (" …" if len(fl) > 8 else "") if fl else ""))
        else:
            st.warning("Folder not found.")
    elif src == "upload":
        up = st.file_uploader("Files", accept_multiple_files=True, key="bp_up")
        files = [(u, u.name) for u in up or []]
    else:
        fl = BP.list_folder(os.path.join(EX, exdir), kind); files = [(p_, p_) for p_ in fl]
        st.caption(f"{len(fl)} example file(s): " + ", ".join(os.path.basename(p_) for p_ in fl))
    popts = None
    if kind == "SEM/TEM particles":       # v2.11: the SEM page's adaptive detection + optional shape model
        with st.container(border=True, key="bp_particles"):
            q1, q2, q3, q4 = st.columns(4)
            pdark = q1.radio("Particles are", ["dark (TEM)", "bright (SEM)"], index=1 if src.startswith("the bundled") else 0,
                             key="bp_pdark", horizontal=True).startswith("dark")
            pprior = q2.toggle("🔷 Shape model", False, key="bp_prior",
                               help="Learn the outline of the clean isolated particles of each image and use it to split clumps.")
            pstr_ = q2.slider("Shape-model strength", 0.0, 1.0, 0.7, 0.05, key="bp_pstr", disabled=not pprior)
            pnm = q3.number_input("Fallback nm/px (files without metadata)", 0.0, 1e4, 0.0, format="%g", key="bp_pnm")
            plk = q3.slider("Sensitivity (× noise σ)", 2.0, 8.0, 4.0, 0.5, key="bp_plk")
            pagg = q4.checkbox("Exclude aggregates", True, key="bp_pagg")
            psph = q4.selectbox("Particle shape", ["auto", "spheres", "off"], key="bp_sphere",
                                format_func={"auto": "auto (round → circle fit)", "spheres": "spheres – circle fit", "off": "outlines only"}.get)
            pcm = q1.slider("Low-confidence threshold", 0.0, 1.0, 0.5, 0.05, key="bp_conf")
            popts = dict(dark=pdark, shape_prior=pprior, prior_strength=pstr_, nm_per_px=pnm or None, low_k=plk, exclude_aggregates=pagg,
                         sphere={"auto": "auto", "spheres": True, "off": False}[psph], conf_min=pcm)
            st.caption("One figure per file = outlines on the image (aggregates red dashed, splits dotted, low-confidence sphere fits orange); "
                       "'all files overlaid' = ECD histograms. CSV = every particle (ECD, Feret, kind, included).")
    k1, k2, k3, k4 = st.columns(4)
    pre = k1.radio("Preset", ["none (defaults)", "house style", "settings JSON"], key="bp_pre")
    cfg = None
    if pre == "house style":
        cfg = FU.house() or None
        if not cfg:
            k1.caption("No house style loaded – defaults are used.")
    elif pre == "settings JSON":
        pj = k1.file_uploader("Figure settings (JSON)", type=["json"], key="bp_json")
        if pj is not None:
            try:
                cfg = F.merge(F.default_settings(), json.loads(pj.getvalue().decode("utf-8")))
            except Exception as e:
                k1.error(f"Could not read the preset: {e}")
    layout = k2.radio("Layout", ["one figure per file", "all files overlaid", "stacked"], key="bp_lay")
    tpl = k3.selectbox("Journal template", ["sidebar width / font"] + list(TP.TEMPLATES), key="bp_tpl")
    colw = k3.radio("Column", ["single", "1.5", "double"], horizontal=True, key="bp_col") if tpl in TP.TEMPLATES else None
    fmts = k4.multiselect("Formats", ["png", "tiff", "pdf", "svg"], ["png", "pdf"], key="bp_fmt")
    dpi = k4.selectbox("dpi", [300, 600, 1200], index=1, key="bp_dpi")
    width = TP.width(tpl, colw) if tpl in TP.TEMPLATES else WIDTH
    fsz = TP.TEMPLATES[tpl]["font_size"] if tpl in TP.TEMPLATES else FONT
    out = st.text_input("Also write to folder (empty = ZIP download only)", "", key="bp_out",
                        placeholder=os.path.join(EXPORT_DIR, f"batch_{kind.split()[0]}"))
    if not files or not fmts:
        return
    if st.button(f"▶ Process {len(files)} file(s)", type="primary", key="bp_run"):
        pb = st.progress(0.0)
        z, log = BP.run(kind, files, cfg, layout, tuple(fmts), dpi, width, ASPECT, fsz, out or None, progress=lambda v: pb.progress(min(1.0, v)),
                        particle_opts=popts)
        st.session_state["bp_result"] = (z, log, kind)
    if st.session_state.get("bp_result"):
        z, log, k_ = st.session_state["bp_result"]
        lg = pd.DataFrame(log)
        ok = int((lg.status == "ok").sum()) if len(lg) else 0
        (st.success if ok == len(lg) else st.warning)(f"Batch done: {ok} / {len(lg)} file(s) read" + (f", written to {out}" if out else "") + ".")
        st.dataframe(lg, hide_index=True, width="stretch")
        dl(st, "⬇️ ZIP (figures + data + log + preset)", z, f"batch_{safe_name(k_)}.zip", "bp_zip")


def page_si():
    from mfs import si_bundle as SB
    st.header("SI export bundle")
    st.caption("Every figure shown while “Collect figures for the SI bundle” (sidebar) is on is kept here with its data and the "
               "settings used. One ZIP: figures (PDF vector, PNG at the sidebar dpi, SVG), data CSVs, SI_parameters.xlsx "
               "(figure list, analysis + figure settings, results per sample) and the defect SI table.")
    figs = st.session_state.get("si_figs", {})
    if not figs:
        st.info("No figures collected yet – open the technique pages (collection is on by default)."); return
    tab = pd.DataFrame([{"include": True, "figure": k, "page": v.get("page", ""), "captured": v.get("time", "")} for k, v in figs.items()])
    tab = st.data_editor(tab, hide_index=True, disabled=["figure", "page", "captured"], key="si_tab")
    sel = [r.figure for r in tab.itertuples() if r.include]
    inc = st.multiselect("Include", ["pdf", "png", "svg", "csv"], ["pdf", "png", "csv"], key="si_inc")
    c1, c2 = st.columns(2)
    if c2.button("🗑️ Clear the collection", key="si_clear"):
        st.session_state["si_figs"] = {}; st.rerun()
    if not sel:
        return

    def _bundle(sel=tuple(sel), inc=tuple(inc)):
        out = {}
        for k in sel:
            f = figs[k]; fig = f["fig"]
            d = {kk: f.get(kk) for kk in ("page", "time", "params", "cfg", "title", "n_rows")}
            with _RENDER_LOCK:
                for ext in ("pdf", "png", "svg"):
                    if ext in inc:
                        d[ext] = S.export(fig, ext, DPI)
            if "csv" in inc and f.get("data") is not None:
                d["csv"] = f["data"].to_csv(index=False).encode("utf-8-sig")
            out[k] = d
        return SB.build(out, st.session_state.get("sample_store", {}), __version__, inc)
    dl(c1, f"⬇️ SI bundle ({len(sel)} figures)", _bundle, "Reticula_SI_bundle.zip", "si_zip")
    figt, part, rest = SB.tables({k: figs[k] for k in sel}, st.session_state.get("sample_store", {}))
    st.markdown("**Parameters table (preview)**")
    st.dataframe(part, hide_index=True, width="stretch", height=260)


def _db_path():
    return os.environ.get("MFS_DB") or os.path.join(os.path.dirname(STORE_PATH), "reticula.db")


_SAMPLE_LIST_NAMES = ("Sample_overview.xlsx", "Master_sample_overview.xlsx")


_SAMPLE_LIST_EXAMPLE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "example_data", "SampleDB", "Sample_overview_EXAMPLE.xlsx")


def _sample_list_file(path):
    """An Excel sample list from what the user typed: the workbook itself, or a folder that contains one."""
    path = (path or "").strip().strip('"')
    if os.path.isdir(path):
        for nm in _SAMPLE_LIST_NAMES:
            if os.path.exists(os.path.join(path, nm)):
                return os.path.join(path, nm)
        return None
    return path if os.path.isfile(path) else None


def _pipeline_candidates():
    """The user's own sample list (data folder, Documents, Downloads, Desktop; also in a Synthesis_Pipeline sub-folder),
    else the bundled synthetic example."""
    h = os.path.expanduser("~")
    roots = [os.path.dirname(STORE_PATH)] + [os.path.join(h, d) for d in ("Documents", "Downloads", "Desktop", os.path.join("OneDrive", "Documents"))]
    c = []
    for r_ in roots:
        for sub in ("", "Synthesis_Pipeline"):
            f_ = _sample_list_file(os.path.join(r_, sub))
            if f_ and f_ not in c:
                c.append(f_)
    return c + [_SAMPLE_LIST_EXAMPLE]


def page_db():
    from mfs import sampledb as DB
    st.header("Sample database")
    path = _db_path()
    st.caption(f"Local SQLite database `{path}`: your sample IDs (protocol, family, parent, source), "
               "every result stored by the technique pages, and data files linked by the sample ID in their names "
               "(e.g. `EX-M801-F25-B1_PXRD.xy`).")
    con = DB.connect(path)
    cands = _pipeline_candidates()
    with st.expander("1 · Import sample IDs from an Excel sample list", expanded=not con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]):
        c1, c2 = st.columns([3, 1])
        root = c1.text_input("Sample list (.xlsx with a 'Sample ID' column) or the folder that contains it",
                             cands[0] if cands else "", key="db_root")
        if cands and cands[0] == _SAMPLE_LIST_EXAMPLE:
            c1.caption("This is the bundled EXAMPLE list (synthetic sample IDs). Point it to your own workbook: one row per sample with "
                       "a **Sample ID** column; optional columns *Name, Paper/Project, Protocol family, Parent / starting material, "
                       "Protocol basis, Main source (DOI), Protocol (docx), Protocol (PDF)*.")
        if c2.button("Import / update", key="db_imp"):
            try:
                f_ = _sample_list_file(root)
                if not f_:
                    raise FileNotFoundError("no Excel sample list found at that location")
                n = DB.import_master(con, f_); st.success(f"{n} sample IDs imported from {f_}.")
            except Exception as e:
                st.error(f"Import failed: {e}")
    nres = DB.sync_store(con, st.session_state.get("sample_store", {}))
    with st.expander("2 · Link data files"):
        c1, c2 = st.columns([3, 1])
        dfold = c1.text_input("Data folder (sub-folders included)", os.path.join(os.path.expanduser("~"), "Documents"), key="db_data")
        if c2.button("Scan and link", key="db_scan"):
            if os.path.isdir(dfold):
                n, un = DB.link_folder(con, dfold)
                st.success(f"{n} file(s) linked; {len(un)} without a known sample ID.")
            else:
                st.error("Folder not found.")
    ov = DB.overview(con)
    st.caption(f"{len(ov)} samples · {nres} result values synced from the sample store this run.")
    q = st.text_input("🔎 Search (ID, name, family, paper)", "", key="db_q")
    if q:
        m_ = ov.astype(str).apply(lambda c: c.str.contains(q, case=False, regex=False)).any(axis=1); ov = ov[m_]
    st.dataframe(ov, hide_index=True, width="stretch", height=300)
    if not len(ov):
        return
    sid = st.selectbox("Sample", list(ov["Sample ID"]), key="db_pick")
    s_ = DB.sample(con, sid)
    st.subheader(f"{sid} – {s_.get('name') or ''}")
    a, b = st.columns(2)
    a.markdown(f"**Paper:** {s_.get('paper') or '–'}  \n**Family:** {s_.get('family') or '–'}  \n**Parent:** {s_.get('parent') or '–'}  \n"
               f"**Children:** {', '.join(s_['children']) or '–'}  \n**Basis:** {s_.get('basis') or '–'}"
               + (f"  \n**Source:** [{s_['doi']}]({s_['doi']})" if s_.get("doi") and str(s_["doi"]).startswith("http") else ""))
    for lab, k_ in (("protocol PDF", "protocol_pdf"), ("protocol DOCX", "protocol_docx"), ("paper matrix", "matrix")):
        p_ = s_.get(k_)
        if p_ and os.path.exists(p_):
            dl(b, f"⬇️ {lab}: {os.path.basename(p_)}", open(p_, "rb").read(), os.path.basename(p_), f"db_{k_}_{sid}")
            if b.button(f"Open {lab}", key=f"db_open_{k_}_{sid}"):
                open_with_default_app(p_)
        elif p_:
            b.caption(f"{lab}: `{p_}` (not found on this computer)")
    if s_["results"]:
        st.markdown("**Results**")
        st.dataframe(pd.DataFrame(s_["results"]).astype({"value": str}), hide_index=True, width="stretch")
    if s_["files"]:
        st.markdown("**Data files**")
        st.dataframe(pd.DataFrame(s_["files"]), hide_index=True, width="stretch")
    nt = st.text_area("Notes", s_.get("notes") or "", key=f"db_notes_{sid}")
    if st.button("Save notes", key=f"db_savenotes_{sid}"):
        DB.set_notes(con, sid, nt); st.success("Saved.")
    with open(path, "rb") as fh:
        dl(st, "⬇️ Database file (reticula.db)", fh.read(), "reticula.db", "db_dl")
