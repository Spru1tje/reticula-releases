"""Reticula page module: TGA and DSC pages.

Moved verbatim out of app.py in v2.13. Names that the code takes from the app (helpers such as show(), dl(), files_input(),
the sidebar values WIDTH / ASPECT / FONT, the cached reader wrapper R …) are listed in INJECTED and filled in by
mfs.pages.bind(globals()) on every run of app.py, so the code behaves exactly as it did inside app.py."""
# ruff: noqa: F821
from __future__ import annotations
import io, os, re, json, glob, zipfile, warnings, datetime, hashlib
import numpy as np
import pandas as pd
import streamlit as st
from mfs import numentry as NE  # every slider takes a typed number too (2.19.6)
from mfs import levels as LV
import matplotlib.pyplot as plt

from mfs import analysis as A, plots as P, style as S, defects as D, composer as C, report as REP
from mfs import figctl as F, figctl_ui as FU, ui as UI, __version__
from mfs.chem import SPECIES, FRAMEWORKS, molar_mass

INJECTED = (
    'EX',
    'FONT',
    'R',
    'WIDTH',
    '_CachedReaders',
    '_F',
    '_res_lab',
    'files_input',
    'guess_framework',
    'guess_name',
    'newfig',
    'node_picker',
    'offset_controls',
    'sample_table',
    'set_even_offset',
    'show',
    'store_put',
)


def page_tga():
    st.header("Thermogravimetric analysis")
    st.caption("Netzsch Proteus ASCII, TA Universal Analysis/TRIOS text, Mettler STARe export, or any CSV with temperature + mass/weight columns.")
    files = files_input("TGA files", ["txt", "csv", "dat", "xlsx"], "TGA/*.txt", "tga")
    if not files:
        return
    samples = sample_table(files, "tga_tab")
    c1, c2, c3 = st.columns(3)
    dtg = c1.selectbox("DTG", ["panel", "twin", "none"])
    norm_end = c2.checkbox("Normalise to residue = 100 % (Shearer plot)", False)
    win = NE.slider_in(c3, "DTG smoothing window", 5, 101, 21, 2)
    res = []
    for nm, f in samples:
        try:
            t = R.read_tga(f, f.name); a = A.tga_analysis(t["T"], t["w"], win); res.append((nm, a))
        except Exception as e:
            st.error(f"{f.name}: {e}")
    # v2.17: the plateau point is automatic (the flattest stretch between the steps) or picked on the graph below
    ov = st.session_state.setdefault("tga_plat_ov", {})
    # EJIC 2020 (Lázaro et al., doi:10.1002/ejic.202000656): R_exp uses the mass at the START of the decomposition profile.
    rule = st.radio("Plateau rule", ["flattest stretch (automatic)", "start of decomposition (EJIC 2020)"], horizontal=True, key="tga_plat_rule",
                    help="Flattest stretch: the lowest-slope point of the plateau. Start of decomposition: where the mass first starts to fall "
                         "after the plateau, as in Lázaro et al., Eur. J. Inorg. Chem. 2020 (doi:10.1002/ejic.202000656).")
    for nm, a_ in res:
        info = a_.get("plateau_info") or {}
        if rule.startswith("start") and a_.get("T_onset") is not None:
            a_["T_plateau"] = a_["T_onset"]; a_["W_plateau"] = a_["W_onset"]
            a_["plateau_source"] = "start of decomposition (EJIC 2020)"
        elif rule.startswith("start"):
            a_["plateau_source"] = "start of decomposition – not found (no fall in window)"
        if nm in ov:
            a_["T_plateau"] = float(ov[nm])
            a_["W_plateau"] = float(np.interp(ov[nm], a_["T"], a_["w"]))
            a_["plateau_source"] = "picked on graph"
        elif not rule.startswith("start"):
            a_["plateau_source"] = "automatic (flat)" if info.get("flat") else "automatic – no clear plateau"
    with st.expander("🎯 Plateau point for the defect calculation – automatic, or pick it on the graph", expanded=False):
        st.caption("The point is taken from the flattest stretch between the solvent / dehydroxylation steps and the first linker-combustion step. "
                   "Click the curve to choose another point for that sample (the mass there is used for R_exp).")
        sel = st.selectbox("Sample", [nm for nm, _ in res], key="tga_plat_sample") if res else None
        if sel:
            a_sel = dict(res)[sel]
            import altair as alt
            df_p = pd.DataFrame({"T (°C)": a_sel["T"][::max(1, len(a_sel["T"]) // 900)], "mass (%)": a_sel["w"][::max(1, len(a_sel["T"]) // 900)]})
            pick = alt.selection_point(name="tga_pick", fields=["T (°C)"], nearest=True, on='click', clear="dblclick")
            ch = (alt.Chart(df_p).mark_line(color="#0F766E").encode(x=alt.X("T (°C):Q"), y=alt.Y("mass (%):Q", scale=alt.Scale(zero=False)),
                  tooltip=["T (°C)", "mass (%)"]).add_params(pick).properties(height=260))
            lo_w, hi_w = a_sel["plateau_info"].get("window", (None, None)) if a_sel.get("plateau_info") else (None, None)
            if lo_w:
                band = alt.Chart(pd.DataFrame({"a": [lo_w], "b": [hi_w]})).mark_rect(color="#14B8A6", opacity=0.12).encode(x="a:Q", x2="b:Q")
                ch = band + ch
            ev = st.altair_chart(ch, on_select="rerun", key=f"tga_pick_chart_{sel}", width="stretch")
            try:
                pts = ev.selection.get("tga_pick", {}).get("T (°C)", [])
                if pts:
                    t_pick = float(pts[0])
                    if abs(t_pick - float(ov.get(sel, -1e9))) > 0.5:
                        ov[sel] = t_pick
                        st.rerun()
            except Exception:
                pass
            c_a, c_b = st.columns([1, 3])
            if c_a.button("↺ Use the automatic point", key="tga_plat_reset"):
                ov.pop(sel, None); st.rerun()
            info_s = a_sel.get("plateau_info") or {}
            c_b.caption(f"Window {info_s.get('window', ('–', '–'))[0]}–{info_s.get('window', ('–', '–'))[1]} °C · "
                        + (f"local slope {info_s['slope_pct_per_C']:.4f} % per °C" if info_s.get("slope_pct_per_C") is not None else "")
                        + (f" · {info_s['reason']}" if info_s.get("reason") else "") + f" · now: {a_sel['plateau_source']}")
    tab = pd.DataFrame([{"sample": nm, "W(200 °C) %": round(a["W200"], 2) if a["W200"] else None, "plateau T (°C)": round(a["T_plateau"], 0) if a["T_plateau"] else None,
                         "plateau point": a.get("plateau_source", ""),
                         "W(plateau) %": round(a["W_plateau"], 2) if a["W_plateau"] else None, "W(end) %": round(a["W_end"], 2),
                         "plateau/end": round(a["W_plateau"] / a["W_end"], 4) if a["W_plateau"] else None,
                         "DTG peaks (°C)": ", ".join(f"{p:.0f}" for p in a["peaks"])} for nm, a in res])
    st.dataframe(tab, width="stretch")
    st.session_state["tga_table"] = tab
    for nm, a in res:
        store_put(nm, "tga", {"W_initial": a["W_initial"], "W200": a["W200"], "W_plateau": a["W_plateau"], "T_plateau": a["T_plateau"],
                              "W_end": a["W_end"], "DTG_peaks": list(a["peaks"])}, "TGA page")
    from mfs import tga as TG
    c1, c2, c3 = st.columns(3)
    lab_on = c1.checkbox("Automatic step labels (mass loss per DTG step)", True, key="tga_steps")
    min_loss = c2.number_input("min. step (%-points)", 0.1, 20.0, 1.0, 0.1, key="tga_minloss")
    lab_which = c3.multiselect("label which samples", [nm for nm, _ in res], default=[nm for nm, _ in res][:1], key="tga_labw")
    stp = {nm: TG.steps(a, min_loss) for nm, a in res}
    with st.expander("📉 Steps per sample", expanded=False):
        st.dataframe(pd.DataFrame([dict(sample=nm, step=i + 1, **{"T start (°C)": s_["T_start"], "DTG peak (°C)": s_["T_peak"],
                                                                   "T end (°C)": s_["T_end"], "loss (%-points)": s_["loss"],
                                                                   "loss (% of residue)": s_["loss"] * 100 / a["W_end"]})
                                   for nm, a in res for i, s_ in enumerate(stp[nm])]).style.format(precision=1), hide_index=True, width="stretch")
    with st.expander("⚖️ Plateaus normalised to the oxide residue (= 100 %) → linkers per node (any node / metal)", expanded=False):
        gfw, _gm = guess_framework(res[0][0] if res else "")
        c1, c2 = st.columns([1, 2])
        fw_ = c1.selectbox("Framework", list(FRAMEWORKS), index=list(FRAMEWORKS).index(gfw or "UiO-66"), key="tga_fw")
        nsel_ = node_picker(c2, fw_, "tga", metal_default=_gm)
        met = nsel_["metal"]
        st.caption(f"Ideal plateau: node with all ligand charge compensated by O²⁻ + L linkers → TGA residue {_res_lab(met, nsel_['def'])} "
                   "(M₆O₆(L)₆ → 6 MO₂ for Zr/Hf: Valenzano 2011; Shearer 2014). The defect calculator uses the same equations.")
        pr = []
        for nm, a in res:
            if not a["W_plateau"]:
                continue
            lk = FRAMEWORKS[fw_]["linker"]
            inp = D.DefectInput(sample=nm, framework=fw_, metal=met, species=[D.Species(lk, SPECIES[lk][0], SPECIES[lk][1], "linker", SPECIES[lk][3], None, False, True)],
                                tga_W_plateau=a["W_plateau"], tga_W_end=a["W_end"], node=nsel_["node"], node_custom=nsel_["node_custom"])
            try:
                r_ = D.method_a(inp)
                pr.append({"sample": nm, "plateau T (°C)": a["T_plateau"], "W(plateau) (% of residue)": 100 * a["W_plateau"] / a["W_end"],
                           "ideal (%)": r_.extra.get("W_plat_ideal_%") if r_.extra else None, "linkers / node L": r_.L, "missing x": r_.x})
            except Exception as e:
                pr.append({"sample": nm, "note": str(e)[:80]})
        st.dataframe(pd.DataFrame(pr).style.format(precision=2), hide_index=True, width="stretch")

    def _labels(ax):
        if lab_on:
            for nm, a in res:
                if nm in lab_which:
                    TG.label_steps(ax, a, stp[nm], normalise_end=norm_end)
    if dtg == "panel":
        fig = P.fig_tga(res, WIDTH, "panel", normalise_end=norm_end); _labels(fig.axes[0])
    else:
        fig, ax = newfig(0.8)
        tw = ax.twinx() if dtg == "twin" else None
        if tw is not None:
            ax.tick_params(right=False, which="both"); tw.tick_params(direction="in")
        P.draw_tga(ax, res, tw, normalise_end=norm_end); _labels(ax)
    def draw(ax, res=res):
        P.draw_tga(ax, res, normalise_end=norm_end); _labels(ax)
    show(fig, "TGA", pd.concat([pd.DataFrame({"sample": nm, "T": a["T"], "w": a["w"], "DTG": a["dtg"]}) for nm, a in res]), draw, "TGA",
         units={0: {"x": "temperature", "y": "mass_pct"}, 1: {"x": "temperature"}})
    st.subheader("TGA–MS overlay")
    c1, c2 = st.columns(2)
    msf = c1.file_uploader("Ion-current table (CSV/TXT: temperature or time + m/z columns)", type=["csv", "txt", "dat"], key="tga_ms")
    ms_ex = c2.checkbox("Use the TGA–MS EXAMPLE (MOF-801, synthetic m/z 18, 28, 44)", value=msf is None, key="tga_ms_ex")
    src = msf if msf is not None else (_F(os.path.join(EX, "TGA", "SampleA_QMS_ion_currents_EXAMPLE.csv")) if ms_ex else None)
    if src is not None:
        try:
            ms = TG.read_ms(src)
            chans = [c_ for c_ in ms.columns if c_.startswith("m/z")]
            c1, c2, c3, c4 = st.columns(4)
            tsel = c1.selectbox("TGA trace", [nm for nm, _ in res], key="tga_ms_s")
            chs = c2.multiselect("m/z channels", chans, default=chans[:4], key="tga_ms_ch")
            logy = c3.checkbox("log ion current", False, key="tga_ms_log")
            if "T" not in ms:
                T0 = c4.number_input("start T (°C)", -100.0, 1000.0, 30.0, key="tga_ms_T0"); beta = c4.number_input("rate (K/min)", 0.1, 100.0, 10.0, key="tga_ms_b")
                ms["T"] = T0 + beta * ms["t"]
            a = dict(res)[tsel]

            def draw_ms(ax, a=a, ms=ms, chs=chs):
                ax.plot(a["T"], a["w"] / a["W_end"] * 100 if norm_end else a["w"], color="k", lw=1.0)
                ax.set_ylabel("Weight (%)"); ax.set_xlabel("Temperature (°C)")
                tw = ax.twinx(); TG.draw_ms(tw, ms, chs, log=logy, norm=not logy)
                ax.set_xlim(max(a["T"].min(), ms["T"].min()), min(a["T"].max(), ms["T"].max()))
            fig, ax = newfig(0.8); draw_ms(ax)
            show(fig, "TGA_MS", ms, draw_ms, "TGA–MS")
        except Exception as e:
            st.warning(f"MS file: {e}")


def page_dsc():
    from mfs import dsc as DSC
    DSC = _CachedReaders(DSC)
    st.header("Differential scanning calorimetry")
    st.caption("Mettler-Toledo STARe ‘Curve Values’ export, TA Instruments TRIOS / Universal Analysis text, Netzsch Proteus ASCII, "
               "PerkinElmer Pyris text, or any CSV/XLSX with temperature + heat-flow columns (decimal commas ok). Heating/cooling "
               "segments are split automatically. Peaks: onset / peak / endset and ΔH by integration over time; Tg by ISO 11357-2.")
    files = files_input("DSC files", ["txt", "csv", "dat", "xlsx", "xls"], "DSC/MOF-PEG_TA-TRIOS_EXAMPLE.txt", "dsc")
    if st.checkbox("➕ Use the full-cycle example instead: REAL heating (SG_DSC_1.1, STARe) + SIMULATED cooling (separate export)",
                   False, key="dsc_cycle_ex", help="Two files of one sample that are combined into Heat 1 → Cool 1. The cooling "
                   "curve and its crystallisation exotherm are synthetic (labelled EXAMPLE)."):
        files = [_F(p_) for p_ in sorted(glob.glob(os.path.join(EX, "DSC", "full_cycle", "SG_DSC_*ing_Mettler_*.txt")))]
    if st.checkbox("➕ Use the heat–cool–heat EXAMPLE (REAL heating 1; SIMULATED cooling and heating 2) – cycle-to-cycle comparison",
                   False, key="dsc_hch_ex"):
        files = [_F(os.path.join(EX, "DSC", "full_cycle", "SG_DSC_heat-cool-heat_EXAMPLE.csv"))]
    if not files:
        return
    samples = sample_table(files, "dsc_tab")
    data = []
    for nm, f in samples:
        try:
            data.append((nm, f, DSC.read_dsc(f, f.name)))
        except Exception as e:
            st.error(f"{f.name}: could not read ({type(e).__name__}: {e})")
    if not data:
        return
    def _base(nm):   # 'X_heating' / 'X cooling' / 'X_H1' -> 'X' : files of one sample share a cycle group
        return DSC.cycle_group_name(nm)
    with st.expander("⚖️ Files – instrument, sample mass, sign, cycle", expanded=True):
        info = pd.DataFrame([{"sample": nm, "instrument / format": d["meta"]["vendor"],
                              "mass (mg)": d["meta"].get("mass_mg"), "curves": len(d["curves"]),
                              "file time": d["meta"].get("timestamp").strftime("%Y-%m-%d %H:%M") if d["meta"].get("timestamp") else "",
                              "cycle group": _base(nm), "flip sign": False} for nm, f, d in data])
        info = st.data_editor(info, key="dsc_info", hide_index=True, disabled=["sample", "instrument / format", "curves", "file time"],
                              column_config={"flip sign": st.column_config.CheckboxColumn(help="Tick if exothermic events point the wrong way for this file"),
                                             "cycle group": st.column_config.TextColumn(help="Files with the same group are one sample; with 'Combine' they are merged into one heating/cooling cycle")})
        groups = list(dict.fromkeys(info["cycle group"].astype(str)))
        multi = any((info["cycle group"].astype(str) == g).sum() > 1 for g in groups)
        k1, k2 = st.columns(2)
        combine = k1.checkbox("Combine files of one sample into one cycle (e.g. separate heating and cooling exports)", multi, key="dsc_combine_" + hashlib.md5("|".join(info["sample"].astype(str)).encode()).hexdigest()[:6],
                              help="Merges the files of each cycle group in order, makes the time axis continuous and numbers the segments Heat 1, Cool 1, Heat 2 …")
        order = k2.radio("Order", ["by file time stamp", "as in the sample table"], horizontal=True, key="dsc_order",
                         help="Time stamps come from the file header; if a file has none, the table order is used.") if combine else None
    curves = []
    packs = []
    for (nm, f, d), (_, row) in zip(data, info.iterrows()):
        mass = float(row["mass (mg)"]) if pd.notna(row["mass (mg)"]) and row["mass (mg)"] else None
        cs = []
        for c in d["curves"]:
            c = dict(c)
            if row["flip sign"]:
                c["q"] = -c["q"]
            c["mass"] = mass; c["file"] = nm
            cs.append(c)
        packs.append((str(row["cycle group"]) if combine else nm, {"curves": cs, "meta": d["meta"]}))
    for g in dict.fromkeys(p_[0] for p_ in packs):
        ds = [p_[1] for p_ in packs if p_[0] == g]
        segs = DSC.combine_cycle(ds, "time" if order == "by file time stamp" else "user") if combine else DSC.label_segments(
            [c for d_ in ds for c in d_["curves"]])
        for c in segs:
            c["sample"] = g; c["uid"] = f"{g}|{c['seg']}"
            curves.append(c)
    one = len({c["sample"] for c in curves}) == 1
    for c in curves:
        c["label"] = c["seg"] if one else f"{c['sample']} · {c['seg']}"
    has_cycle = any(c["kind"] == "cooling" for c in curves) and any(c["kind"] == "heating" for c in curves)
    with LV.gate("Advanced", globals()):
        with st.expander("🔁 Segments (detected from the sign of dT/dt)", expanded=False):
            st.dataframe(pd.DataFrame([{"sample": c["sample"], "segment": c["seg"], "file": c.get("file", ""), "kind": c["kind"],
                                        "T start (°C)": round(float(c["T"][0]), 2), "T end (°C)": round(float(c["T"][-1]), 2),
                                        "rate (K/min)": round(DSC.heating_rate(c), 3) if DSC.heating_rate(c) is not None else None,
                                        "t start (min)": round(float(c["t"][0]), 2) if c.get("t") is not None else None,
                                        "points": len(c["T"])} for c in curves]), hide_index=True, width="stretch")
            if not has_cycle:
                st.caption("Only one direction found. For the full cycle upload the cooling export too (same cycle group) – or "
                           "export all segments from the instrument software (STARe: all curves of the method).")
    with st.expander("⚙️ Plot and evaluation settings", expanded=True):
        c1, c2, c3, c4 = st.columns(4)
        labels = [c["label"] for c in curves]
        default = [l for l, c in zip(labels, curves) if c["kind"] != "isothermal"]
        import hashlib as _hl   # key tied to the loaded curves: new files -> fresh selection (no stale/empty choice)
        sel = c1.multiselect("Segments / cycles", labels, default=default,
                             key="dsc_sel_" + _hl.md5("|".join(labels).encode()).hexdigest()[:8])
        exo_up = c1.radio("Convention", ["exo up", "exo down"], horizontal=True, key="dsc_exo") == "exo up"
        VIEWS = ["Curves", "Full cycle vs temperature", "Full cycle vs time", "Stacked: heating above cooling"]
        view = c1.selectbox("View", VIEWS, index=1 if has_cycle else 0, key="dsc_view_" + ("c" if has_cycle else "s"),
                            help="Full cycle vs T: heating and cooling as one closed loop with direction arrows. Vs time: with the "
                                 "temperature programme on a twin axis and shaded segments. Stacked: heating panel above the cooling panel.")
        yunit = c2.selectbox("Heat flow", ["W/g (= mW/mg)", "mW/mg", "mW"], key="dsc_unit",
                             help="W/g needs normalised data or a sample mass; mW multiplies normalised data by the mass.")
        xax = c2.radio("x axis", ["temperature", "time"], horizontal=True, key="dsc_x")
        off = offset_controls(c3, "dsc_off", 0.0, 0.0, 5.0, 0.05, "Vertical offset", help="Stacks the curves (y ticks are hidden when > 0)")
        bl = c3.selectbox("Baseline", ["linear", "sigmoidal"], key="dsc_bl")
        do_pk = c4.checkbox("Detect peaks (onset / peak / endset / ΔH)", True, key="dsc_pk")
        prom = c4.number_input("Min. peak prominence (0 = auto)", 0.0, 100.0, 0.0, format="%g", key="dsc_prom")
        ann = c4.checkbox("Annotate peaks on the figure", True, key="dsc_ann")
        showbl = c4.checkbox("Show integration baselines", True, key="dsc_showbl")
        tg_on = c1.checkbox("Glass transition Tg", False, key="dsc_tg")
        if tg_on:
            Tall = np.concatenate([c["T"] for c in curves])
            lo, hi = float(np.nanmin(Tall)), float(np.nanmax(Tall))
            heat = [c for c in curves if c["kind"] == "heating"] or curves
            try:
                h0 = heat[0]
                sug = DSC.auto_tg_window(h0, DSC.evaluate_all(dict(h0, id=0), "linear", h0["mass"]))
            except Exception:
                sug = None
            sug = sug or (lo + 0.15 * (hi - lo), lo + 0.45 * (hi - lo))
            tgw = NE.slider_in(c2, "Tg window (°C)", lo, hi, (round(sug[0], 1), round(sug[1], 1)), key="dsc_tgw",
                            help="Auto-suggested around the steepest step outside the peaks; tangents use the outer quarters.")
            tg_cool = c2.checkbox("Tg also on cooling segments", False, key="dsc_tgcool")
    styles = DSC.default_styles(list(dict.fromkeys(c["sample"] for c in curves)))
    if view != "Curves":
        with LV.gate("Advanced", globals()):
            with st.expander("🎨 Segment styles – colour and line style for heating / cooling, per sample", expanded=False):
                inv = {v: k for k, v in DSC.LINESTYLES.items()}
                sty = pd.DataFrame([{"sample": smp, "segment kind": kd, "colour": col, "line style": inv.get(ls, "solid")}
                                    for (smp, kd), (col, ls) in styles.items() if any(c["sample"] == smp and c["kind"] == kd for c in curves)])
                sty = st.data_editor(sty, key="dsc_sty_" + str(abs(hash(tuple(sty["sample"]))) % 10**6), hide_index=True,
                                     disabled=["sample", "segment kind"],
                                     column_config={"colour": st.column_config.TextColumn(help="Hex colour, e.g. #C0392B"),
                                                    "line style": st.column_config.SelectboxColumn(options=list(DSC.LINESTYLES))})
                for _, r in sty.iterrows():
                    col = str(r["colour"]).strip()
                    if re.fullmatch(r"#?[0-9a-fA-F]{6}", col):
                        col = col if col.startswith("#") else "#" + col
                        styles[(r["sample"], r["segment kind"])] = (col, DSC.LINESTYLES.get(r["line style"], "-"))
                q1, q2, q3 = st.columns(3)
                arrows = q1.checkbox("Direction arrows on the curves", True, key="dsc_arrows")
                connect = q2.checkbox("Close the loop (join consecutive segments)", True, key="dsc_connect")
                prog = q3.radio("Temperature trace (vs time)", ["Ts (sample)", "Tr (programme)"], horizontal=True, key="dsc_prog")
    else:
        arrows, connect, prog = True, True, "Ts (sample)"
    use = [c for c in curves if c["label"] in sel]
    if not use:
        st.info("Select at least one segment."); return
    items, peaks, tgs, rows, warn = [], {}, {}, [], []
    for c in use:
        spec = DSC.to_specific(c, c["mass"])
        if yunit == "mW":
            y = c["q"] if c["unit"] in ("mW", "?") else (spec * c["mass"] if (spec is not None and c["mass"]) else None)
            if y is None:
                warn.append(f"{c['label']}: no sample mass – shown in W/g"); y = spec
        else:
            y = spec if spec is not None else c["q"]
            if spec is None:
                warn.append(f"{c['label']}: heat flow in {c['unit']} and no sample mass – enter it in the table to normalise")
        cc = dict(c); cc["id"] = c["uid"]
        if do_pk:
            try:
                prs = DSC.evaluate_all(cc, bl, c["mass"], prominence=prom or None)
            except Exception as e:
                prs = []; warn.append(f"{c['label']}: peak evaluation failed ({e})")
            if yunit == "mW" and c["mass"] and spec is not None:
                for r in prs:
                    r["baseline"] = r["baseline"] * (c["mass"] if c["unit"] == "W/g" else 1)
            peaks[cc["id"]] = prs
            for r in prs:
                rows.append({"curve": c["label"], "segment": c["seg"], "file": c.get("file", ""), "event": r["kind"], "onset (°C)": r["T_onset"], "peak (°C)": r["T_peak"],
                             "endset (°C)": r["T_endset"], "limits (°C)": f"{r['T_start']:.1f}–{r['T_end']:.1f}",
                             "ΔH (J/g)": r["dH"], "area (mJ)": r["area_mJ"], "baseline": bl})
        if tg_on and (c["kind"] == "heating" or tg_cool):
            g = DSC.glass_transition(cc, *tgw, mass_mg=c["mass"])
            if g:
                tgs[cc["id"]] = [g]
                rows.append({"curve": c["label"], "segment": c["seg"], "file": c.get("file", ""), "event": "Tg", "onset (°C)": g["Tg_onset"], "peak (°C)": g["Tg_mid"],
                             "endset (°C)": g["Tg_endset"], "limits (°C)": f"{g['window'][0]:.1f}–{g['window'][1]:.1f}",
                             "ΔH (J/g)": None, "area (mJ)": None, "baseline": f"Tg midpoint; inflection {g['Tg_inflection']:.1f} °C"
                             + (f"; Δcp {g['dcp_J_gK']:.3f} J g⁻¹ K⁻¹" if g["dcp_J_gK"] else "")})
        items.append((c["label"], cc, y))
    for w in warn:
        st.warning(w)
    if rows:
        tab = pd.DataFrame(rows)
        st.dataframe(tab.round(2), width="stretch", hide_index=True)
        st.session_state["dsc_table"] = tab
    for nm, f, d in data:
        rep = d["meta"].get("reported")
        if rep:
            mine = [r for r in rows if r.get("file") == nm and r["event"] != "Tg"]
            best = min(mine, key=lambda r: abs((r["peak (°C)"] or 0) - (rep.get("peak_C") or 0))) if mine else None
            st.info(f"**{nm}** – values reported by {d['meta']['vendor']}: "
                    + ", ".join(f"{k.replace('_C', ' (°C)').replace('_Jg', ' (J/g)').replace('_mJ', ' (mJ)').replace('_', ' ')} {v:g}" for k, v in rep.items())
                    + (f".  Reticula: onset {best['onset (°C)']:.2f} °C, peak {best['peak (°C)']:.2f} °C, ΔH {best['ΔH (J/g)']:.2f} J/g." if best and best["ΔH (J/g)"] is not None else ""))
    ylab = {"W/g (= mW/mg)": "Heat flow (W g$^{-1}$)", "mW/mg": "Heat flow (mW mg$^{-1}$)", "mW": "Heat flow (mW)"}[yunit]
    pk_ = peaks if do_pk else None
    if any(c.get("seg") for _, c, _ in items) and do_pk:
        summ = DSC.cycle_summary([c for _, c, _ in items], peaks)
        if len(summ) and summ[["Tm peak (°C)", "Tc peak (°C)"]].notna().any().any():
            st.markdown("**Cycle summary – melting (heating) vs crystallisation (cooling), supercooling ΔT = Tm − Tc**")
            st.dataframe(summ.round(2), width="stretch", hide_index=True)
            st.session_state["dsc_cycle_summary"] = summ
        cmp_ = DSC.cycle_compare([c for _, c, _ in items], peaks)
        kinds_ = pd.Series([r_["kind"] for r_ in cmp_]).value_counts() if cmp_ else pd.Series(dtype=int)
        if len(kinds_) and kinds_.max() >= 2:
            with st.expander("🔄 Cycle-to-cycle comparison (Heat 1 vs Heat 2 …, Cool 1 vs Cool 2 …)", expanded=True):
                st.dataframe(pd.DataFrame(cmp_).round(3), width="stretch", hide_index=True)
                ck = st.radio("overlay", [k_ for k_ in ("heating", "cooling") if kinds_.get(k_, 0) >= 2], horizontal=True, key="dsc_cmp_k")
                its_ = [it for it in items if it[1].get("kind") == ck]
                pks_ = peaks

                def draw_cmp(ax, its_=its_, pks_=pks_):
                    DSC.draw_dsc(ax, its_, "T", exo_up, 0.0, peaks=pks_ if do_pk else None, annotate=ann, show_baseline=False, ylabel=ylab,
                                 linestyles=["-", "--", ":", "-."] * 3)
                    ax.legend(frameon=False)
                figc, axc = newfig(0.75); draw_cmp(axc)
                show(figc, f"DSC_cycle_compare_{ck}", None, draw_cmp, f"DSC {ck} cycles")
    if view == "Full cycle vs time" and any(c.get("t") is None for _, c, _ in items):
        st.info("Some segments have no time column – showing the full cycle vs temperature instead."); view = "Full cycle vs temperature"
    progk = "Tr" if prog.startswith("Tr") else "Ts"
    if view == "Stacked: heating above cooling":
        S.apply_style(FONT)
        fig, axs = plt.subplots(2, 1, figsize=(WIDTH, WIDTH * 1.25), sharex=True)
        hi_ = [it for it in items if it[1]["kind"] == "heating"]; co_ = [it for it in items if it[1]["kind"] == "cooling"]
        DSC.draw_cycle_stacked(axs, hi_, co_, exo_up, styles, pk_, ann, showbl, ylab, arrows)
        for a_, its in zip(axs, (hi_, co_)):
            if len(its) > 1:
                a_.legend(frameon=False, fontsize=FONT - 1)
        fig.tight_layout()
        draw = None
    else:
        set_even_offset("dsc_off", [(it[1], it[2]) if not isinstance(it[1], dict) else (it[1].get("T"), it[2]) for it in items])
        fig, ax = newfig(0.75)
        if view == "Full cycle vs temperature":
            DSC.draw_cycle_T(ax, items, exo_up, styles, pk_, ann, showbl, ylab, arrows, connect)
        elif view == "Full cycle vs time":
            tw = DSC.draw_cycle_time(ax, items, exo_up, styles, pk_, ann, showbl, ylab, program=progk)
            tw.tick_params(which="both", direction="in")
        else:
            DSC.draw_dsc(ax, items, "T" if xax == "temperature" else "t", exo_up, off, peaks=pk_,
                         tgs=tgs, annotate=ann, show_baseline=showbl, ylabel=ylab)
        if len(items) > 1:
            ax.legend(frameon=False, loc="best")
        def draw(ax, items=items, peaks=peaks, tgs=tgs, view=view):
            if view == "Full cycle vs temperature":
                DSC.draw_cycle_T(ax, items, exo_up, styles, peaks if do_pk else None, ann, showbl, ylab, arrows, connect)
            elif view == "Full cycle vs time":
                DSC.draw_cycle_time(ax, items, exo_up, styles, peaks if do_pk else None, ann, showbl, ylab, program=progk)
            else:
                DSC.draw_dsc(ax, items, "T" if xax == "temperature" else "t", exo_up, off, peaks=peaks if do_pk else None, tgs=tgs,
                             annotate=ann, show_baseline=showbl, ylabel=ylab)
    exp = pd.concat([pd.DataFrame({"curve": lab, "sample": c.get("sample", ""), "segment": c.get("seg", ""), "kind": c.get("kind", ""),
                                   "T (°C)": c["T"], "Tr (°C)": c["Tr"] if c.get("Tr") is not None else np.nan,
                                   "t (min)": c["t"] if c.get("t") is not None else np.nan,
                                   "heat flow (exo +)": y}) for lab, c, y in items])
    stem = {"Curves": "DSC", "Full cycle vs temperature": "DSC_cycle_T", "Full cycle vs time": "DSC_cycle_time"}.get(view, "DSC_cycle_stacked")
    show(fig, stem, exp, draw, "DSC", units={0: {"x": "temperature" if view != "Full cycle vs time" else "time"}},
         series=[lab for lab, _, _ in items])
    sta = [c for c in use if c.get("w") is not None]
    if sta:
        st.subheader("TGA–DSC overlay (STA / SDT)")
        c = sta[0] if len(sta) == 1 else next(cc for cc in sta if cc["label"] == st.selectbox("Curve", [s_["label"] for s_ in sta], key="dsc_sta_sel"))
        y = DSC.to_specific(c, c["mass"]); y = y if y is not None else c["q"]
        fig2, ax2 = newfig(0.75)
        tw = DSC.draw_sta(ax2, c, y, exo_up)
        ax2.tick_params(right=False, which="both")
        show(fig2, "TGA_DSC_overlay", pd.DataFrame({"T (°C)": c["T"], "heat flow (exo +)": y, "mass (%)": c["w"]}), None, "TGA–DSC",
             units={0: {"x": "temperature"}})
    _dsc_cp()


def _dsc_cp():
    from mfs import dsc as DSC
    DSC = _CachedReaders(DSC)
    st.subheader("Specific heat capacity – sapphire method (ASTM E1269)")
    with LV.gate("Advanced", globals()):
        with st.expander("🌡️ cp from blank + sapphire + sample runs (same pans, rate and purge)", expanded=False):
            st.caption("cp(T) = cp,sapphire(T) · (q_sample − q_blank)/(q_sapphire − q_blank) · m_sapphire/m_sample. Sapphire values: "
                       "ASTM E1269 Table 1 (150–600 K), Ginnings & Furukawa above. Use the heating segment of each run (isothermal holds removed).")
            ex_ = st.checkbox("Use the Cp EXAMPLE (synthetic: blank, 25.30 mg sapphire, 9.80 mg MOF-801; truth cp = 0.62 + 0.0016 (T − 25))",
                              True, key="cp_ex")
            c1, c2, c3 = st.columns(3)
            fb = c1.file_uploader("blank (empty pans)", key="cp_b"); fs_ = c2.file_uploader("sapphire", key="cp_s")
            fx = c3.file_uploader("sample(s)", accept_multiple_files=True, key="cp_x")
            D_ = os.path.join(EX, "DSC", "Cp_sapphire")
            fb = fb or (_F(os.path.join(D_, "Cp_blank_EXAMPLE.csv")) if ex_ else None)
            fs_ = fs_ or (_F(os.path.join(D_, "Cp_sapphire_EXAMPLE.csv")) if ex_ else None)
            fx = fx or ([_F(os.path.join(D_, "Cp_sample_MOF-801_EXAMPLE.csv"))] if ex_ else [])
            if not (fb and fs_ and fx):
                st.caption("Load the three runs."); return
            try:
                db, ds = DSC.read_dsc(fb, fb.name), DSC.read_dsc(fs_, fs_.name)
                heat = lambda d: max([c_ for c_ in d["curves"] if (DSC.heating_rate(c_) or 1) > 0], key=lambda c_: len(c_["T"]))
                msap = c2.number_input("sapphire mass (mg)", 0.1, 200.0, float(ds["meta"].get("mass_mg") or 25.0), key="cp_ms_" + fs_.name)
                res_ = []
                for f_ in fx:
                    dx = DSC.read_dsc(f_, f_.name)
                    mx = c3.number_input(f"mass {f_.name} (mg)", 0.01, 200.0, float(dx["meta"].get("mass_mg") or 10.0), key="cp_mx_" + f_.name)
                    r_ = DSC.cp_from_sapphire(heat(dx), heat(ds), heat(db), mx, msap)
                    res_.append((guess_name(f_.name), r_))
                show_ref = st.checkbox("show the sapphire reference curve", False, key="cp_ref")

                def draw_cp(ax, res_=res_):
                    DSC.draw_cp(ax, res_, show_ref=show_ref)
                fig, ax = newfig(0.75); draw_cp(ax)
                exp_ = pd.concat([pd.DataFrame({"sample": nm, "T (°C)": r_["T"], "cp (J/g/K)": r_["Cp"], "cp direct (J/g/K)": r_.get("Cp_direct"),
                                                "cp sapphire (J/g/K)": r_["Cp_sapphire"]}) for nm, r_ in res_])
                show(fig, "DSC_cp", exp_, draw_cp, "cp")
                Tt = [25, 50, 100, 150, 200, 250]
                st.dataframe(pd.DataFrame([{"sample": nm, **{f"cp {T_} °C": float(np.interp(T_, r_["T"], r_["Cp"])) for T_ in Tt if r_["T"].min() <= T_ <= r_["T"].max()},
                                            "sapphire check (direct/table)": float(np.median(r_["Cp_sapphire_check"] / r_["Cp_sapphire"])) if r_.get("Cp_sapphire_check") is not None else None}
                                           for nm, r_ in res_]).round(4), hide_index=True, width="stretch")
                for nm, r_ in res_:
                    store_put(nm, "dsc", {"cp_25C": float(np.interp(25, r_["T"], r_["Cp"])) if r_["T"].min() <= 25 else None,
                                          "cp_100C": float(np.interp(100, r_["T"], r_["Cp"])) if r_["T"].min() <= 100 <= r_["T"].max() else None}, "DSC cp")
            except Exception as e:
                st.warning(f"cp evaluation failed: {e}")
