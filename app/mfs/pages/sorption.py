"""Reticula page module: N2 sorption (isotherm, BET/Rouquerol, t-plot, PSD) page.

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
    'FONT',
    'R',
    'WIDTH',
    'dl',
    'files_input',
    'newfig',
    'safe_name',
    'sample_table',
    'show',
    'store_put',
)


MARKERS = {"circle ●": "o", "square ■": "s", "triangle ▲": "^", "down-triangle ▼": "v", "diamond ◆": "D", "pentagon ⬟": "p",
           "hexagon ⬢": "h", "star ★": "*", "plus +": "P", "cross ×": "X", "none": "none"}


LINESTYLES = {"solid —": "-", "dashed – –": "--", "dotted ···": ":", "dash-dot –·": "-.", "none": "none"}


def page_sorption():
    st.header("N₂ physisorption at 77 K – Micromeritics 3Flex / TriStar / ASAP")
    st.caption("MicroActive 'Export' (.xls via xlrd / .xlsx) of the report, or CSV/TXT with p/p₀ and quantity adsorbed. "
               "v2 splits the export into its individual reports and uses only the Isotherm Tabular Report "
               "(the BET, Langmuir, t-plot and plot-overlay tables are no longer merged into the isotherm). TriStar/ASAP "
               "'Full Report Set' exports with all reports side by side on one sheet are split into their report columns.")
    files = files_input("Isotherm files", ["xls", "xlsx", "csv", "txt"], "Sorption/*_3Flex_*.xls*", "sorp")
    if not files:
        return
    samples = sample_table(files, "sorp_tab")
    data = {}
    for nm, f in samples:
        try:
            data[nm] = R.read_sorption(f, f.name)
        except Exception as e:
            st.error(f"{f.name}: {e}")
    if not data:
        return
    with st.expander("🔎 What was read from each file (tables found, which one is used, units, sample name)"):
        for nm, d in data.items():
            m = d["meta"]
            st.markdown(f"**{nm}** – sample in file: `{m.get('sample', '–')}` · mass {m.get('mass_g', '–')} g · "
                        f"{int((d['iso'].branch == 'ads').sum())} adsorption + {int((d['iso'].branch == 'des').sum())} desorption points · "
                        f"quantity unit in file: {m.get('q_unit_in_file')} (plotted as cm³ g⁻¹ STP) · isotherm from: {m.get('isotherm_table', '–')}")
            if d.get("tables"):
                st.dataframe(pd.DataFrame(d["tables"]), width="stretch", hide_index=True)
    # ---------------------------------------------------------------- per-sample style
    with st.expander("🎨 Series style per sample (adsorption / desorption colour, marker, line)", expanded=False):
        styles = {}
        defaults = S.colours_for(list(data))
        for nm, dc in zip(data, defaults):
            h = hashlib.md5(nm.encode()).hexdigest()[:10]
            c = st.columns([2, 1, 1, 2, 2, 1, 2, 1])
            c[0].markdown(f"**{re.sub(r'[$_{{}}]|mathrm', '', nm)}**")
            ca = c[1].color_picker("ads", dc, key=f"iso_ca_{h}")
            cd = c[2].color_picker("des", ca, key=f"iso_cd_{h}")
            mk = c[3].selectbox("marker", list(MARKERS), key=f"iso_mk_{h}")
            mkd = c[4].selectbox("des marker", ["same"] + list(MARKERS), key=f"iso_mkd_{h}")
            ms = c[5].number_input("size", 0.5, 12.0, 3.0, 0.5, key=f"iso_ms_{h}")
            ls = c[6].selectbox("line", list(LINESTYLES), key=f"iso_ls_{h}")
            lw = c[7].number_input("width", 0.0, 4.0, 0.8, 0.1, key=f"iso_lw_{h}")
            styles[nm] = {"ads_color": ca, "des_color": cd, "marker": MARKERS[mk], "des_marker": None if mkd == "same" else MARKERS[mkd],
                          "ms": ms, "ls": LINESTYLES[ls], "lw": lw}
        st.caption("Adsorption = filled symbols, desorption = open symbols (MOF convention). Lines join the points in measurement order.")
    c1, c2, c3, c4 = st.columns(4)
    logx = c1.checkbox("log p/p₀ axis (micropore filling)", False)
    both = c1.checkbox("also show linear + log side by side", False)
    show_des = c2.checkbox("show desorption", True)
    psd_kind = c3.selectbox("PSD representation", ["dvdw", "dvdlogw", "cum"], format_func={"dvdw": "dV/dw", "dvdlogw": "dV/dlog(w)", "cum": "cumulative"}.get)
    wmax = c4.number_input("PSD max width (nm)", 0.5, 50.0, 2.0)
    with st.expander("✨ Figure options (v2.13 redesign)", expanded=False):
        o1, o2, o3, o4 = st.columns(4)
        hyst = o1.checkbox("shade hysteresis loop", False, key="sorp_hyst")
        dlab = o1.checkbox("sample names at curve ends (no legend box)", False, key="sorp_dlab")
        thin = o2.checkbox("thin out markers (≈18 per branch)", False, key="sorp_thin")
        inset = False  # user (t102u): no graphs inside graphs
        cum = o3.checkbox("cumulative pore volume on a twin axis (PSD)", False, key="sorp_cum")
        bet_all = o3.checkbox("BET / Rouquerol: overlay all samples", False, key="sorp_bet_all")
        bet_box = o4.checkbox("BET annotation box", True, key="sorp_bet_box")
        st.caption("ACS styling with Okabe–Ito colours: adsorption filled, desorption open, lines in measurement order. "
                   "Labels and legends can be dragged in the figure editor below each plot.")
    # ---------------------------------------------------------------- BET
    rows = []
    notes = []                                   # v2.14: warnings / infos about the fits: shown in the collapsed details, not up front
    rq_ok = {}
    st.subheader("BET surface area")
    for nm, d in data.items():
        a = d["iso"][d["iso"].branch == "ads"]
        bet = A.bet_rouquerol(a.p_rel.values, a.q.values)
        tp = A.t_plot(a.p_rel.values, a.q.values)
        vt = A.total_pore_volume(a.p_rel.values, a.q.values, 0.95)
        m = d["meta"]
        rr = m.get("bet_reported_range")
        fixed = A.bet_fixed_range(a.p_rel.values, a.q.values, *(rr or (0.05, 0.30)))
        for w in A.isotherm_quality(a.p_rel.values, a.q.values):
            notes.append(("warn", f"{nm}: {w}"))
        if m.get("bet_C_reported"):
            try:
                if float(m["bet_C_reported"]) <= 0:
                    notes.append(("info", f"{nm}: the file's BET fit has C = {float(m['bet_C_reported']):.1f} (negative). A negative C has no physical "
                                  "meaning (ISO 9277 requires C > 0), so the reported area is not a valid BET area – the app's own value is used."))
            except ValueError:
                pass
        if bet is None:
            st.error(f"{nm}: fewer than 3 adsorption points below p/p₀ 0.35 – no BET area can be computed from this file.")
        if bet:
            rows.append({"sample": nm, "S_BET (m2/g) – this app": round(bet["area"], 0), "C": round(bet["C"], 0),
                         "p/p0 range": f"{bet['p_lo']:.2e}–{bet['p_hi']:.2e}", "points": bet["n"], "R2": round(bet["R2"], 5),
                         "area range over valid windows (m2/g)": (lambda sp_: (f"{sp_['min']:.0f}–{sp_['max']:.0f} ({sp_['n_valid'] if sp_['strict'] else sp_['n_windows']} windows"
                                                                                + ("" if sp_["strict"] else ", none fully valid") + ")") if sp_ else None)(A.bet_range_spread(bet)),
                         "all 4 Rouquerol criteria": bet["all_criteria"],
                         "reported in file (m2/g)": float(m["bet_reported"]) if m.get("bet_reported") else None,
                         "C in file": float(m["bet_C_reported"]) if m.get("bet_C_reported") else None,
                         "file's p/p0 range": f"{rr[0]:.3g}–{rr[1]:.3g}" if rr else None,
                         "same data, file's range (m2/g)": round(fixed["area"], 0) if fixed else None,
                         "V_total @0.95 (cm3/g)": round(vt, 3), "V_micro t-plot (cm3/g)": round(tp["V_micro"], 3) if tp else None,
                         "S_ext t-plot (m2/g)": round(tp["S_ext"], 0) if tp else None})
            d["bet"] = bet
            store_put(nm, "n2", {"S_BET": bet["area"], "C": bet["C"], "p_lo": bet["p_lo"], "p_hi": bet["p_hi"], "points": bet["n"],
                                 "all_criteria": bet["all_criteria"], "V_total": vt, "V_micro": tp["V_micro"] if tp else None,
                                 "S_ext": tp["S_ext"] if tp else None}, "N2 sorption page")
            rq_ok[nm] = bool(bet["all_criteria"]) and not bet.get("few_points")
            if bet.get("few_points"):
                notes.append(("warn", f"{nm}: only {bet['n']} adsorption points lie in the Rouquerol range (p/p₀ {bet['p_lo']:.3g}–{bet['p_hi']:.3g}); "
                              "the BET area is computed from these 3 points (the app normally requires ≥ 4). Treat it as an estimate and "
                              "measure more points below p/p₀ 0.05."))
            if not bet["all_criteria"]:
                miss = [n for n, k in (("1 (n(1-p) increasing)", "c1"), ("3 (n_m inside the fitted range)", "c3"), ("4 (1/(√C+1) ≈ p at n_m, ±20 %)", "c4")) if not bet[k]]
                notes.append(("warn", f"{nm}: no p/p0 range satisfies every Rouquerol criterion – failed: {', '.join(miss)}. "
                              "Common for ultramicroporous MOFs: report a BET-equivalent area with the p/p0 range and show the Rouquerol plot in the SI."))
    if rows:
        tab = pd.DataFrame(rows)
        main = pd.DataFrame({"sample": tab["sample"], "S_BET (m² g⁻¹)": tab["S_BET (m2/g) – this app"],
                             "p/p₀ range": tab["p/p0 range"], "V_total @0.95 (cm³ g⁻¹)": tab["V_total @0.95 (cm3/g)"],
                             "V_micro t-plot (cm³ g⁻¹)": tab["V_micro t-plot (cm3/g)"],
                             "Rouquerol": ["✅ consistent" if rq_ok.get(n) else "⚠️ see details" for n in tab["sample"]]})
        st.dataframe(main, width="stretch", hide_index=True)
        n_warn = sum(1 for lv, _ in notes if lv == "warn")
        with st.expander("🔎 Rouquerol consistency, comparison with the file's own BET value, notes"
                         + (f" ({n_warn} note{'s' if n_warn != 1 else ''})" if n_warn else ""), expanded=False):
            for lv, txt in notes:
                (st.warning if lv == "warn" else st.info)(txt)
            crit = []
            for nm, d in data.items():
                b_ = d.get("bet")
                if not b_:
                    continue
                ok = lambda v: "✅" if v else "❌"
                crit.append({"sample": nm, "p/p₀ range (auto)": f"{b_['p_lo']:.2e}–{b_['p_hi']:.2e}", "points": b_["n"],
                             "1 · n(1−p/p₀) increases": ok(b_["c1"]), "2 · C > 0": ok(b_["C"] > 0) + f" C = {b_['C']:.0f}",
                             "3 · n_m inside the range": ok(b_["c3"]) + (f" p(n_m) = {b_['p_nm']:.2e}" if b_["c3"] else ""),
                             "4 · 1/(√C+1) ≈ p(n_m) ±20 %": ok(b_["c4"]) + f" {b_['p_th']:.2e}",
                             "S_BET (m²/g)": round(b_["area"], 0)})
            if crit:
                st.markdown("**Rouquerol criteria panel** – the BET range is chosen automatically as the range with the highest R² among "
                            "all ranges that satisfy every criterion (ISO 9277:2022, Rouquerol et al. 2007)")
                st.dataframe(pd.DataFrame(crit), width="stretch", hide_index=True)
            st.markdown("**All fit values** (C, R², area range over the valid windows, the file's own value and the same data in the file's window)")
            st.dataframe(tab, width="stretch", hide_index=True)
            st.caption("The value **reported in the file** is shown for reference only. Instrument software fits a fixed window "
                       "(typically p/p₀ 0.05–0.30); for microporous MOFs that window lies after micropore filling and usually gives a lower, "
                       "Rouquerol-inconsistent area. This app's value uses the Rouquerol-consistent range (ISO 9277:2022). The column "
                       "'same data, file's range' shows that the difference is the pressure window, not the data.")
        dl(st, "CSV BET/t-plot table", tab.to_csv(index=False).encode("utf-8-sig"), "BET_table.csv", "bet_tab_csv")
    elif notes:
        with st.expander("🔎 Notes", expanded=False):
            for lv, txt in notes:
                (st.warning if lv == "warn" else st.info)(txt)
    sbet_leg = c2.checkbox("S_BET in the legend", False, key="sorp_sbet_leg")
    lab_ = {nm: (f"{nm} ({d['bet']['area']:.0f} m$^2$ g$^{{-1}}$)" if sbet_leg and d.get("bet") else nm) for nm, d in data.items()}
    styles = {lab_[k]: v for k, v in styles.items()}
    isos = [(lab_[nm], d["iso"] if show_des else d["iso"][d["iso"].branch == "ads"]) for nm, d in data.items()]
    psds = [(nm, d["psd"]) for nm, d in data.items() if d["psd"] is not None]

    def draw_iso(ax, isos=isos, logx=logx, styles=styles):
        P.draw_isotherms(ax, isos, logx=logx, styles=styles, hysteresis=hyst, direct_labels=dlab, thin=thin,
                         inset_psd=psds if (inset and psds and not logx) else None, inset_kind="dvdw" if psd_kind == "cum" else psd_kind,
                         inset_xlim=(0.3, wmax))

    def draw_psd(ax, psds=psds):
        P.draw_psd(ax, psds, kind=psd_kind, xlim=(0.3, wmax), legend=len(psds) > 1, cumulative=cum and psd_kind != "cum")
    csv = pd.concat([df.assign(sample=n) for n, df in isos])
    if both:
        fig, axs = C.compose([lambda ax: draw_iso(ax, logx=False), lambda ax: P.draw_isotherms(ax, isos, logx=True, styles=styles, legend=False, hysteresis=hyst, thin=thin)],
                             ncols=2, width=max(WIDTH, 6.0), panel_aspect=0.85, wspace=0.4)
        show(fig, "N2_isotherms_lin_log", csv, None, units={"y": "uptake"})
    else:
        fig, ax = newfig(); draw_iso(ax)
        show(fig, "N2_isotherms", csv, draw_iso, "N2 isotherm", units={"y": "uptake"})
    no_psd = [nm for nm, d in data.items() if d["psd"] is None]
    if no_psd:
        st.info("No pore-size-distribution table in: " + ", ".join(no_psd) + ". These files contain at most BJH/D-H summary values; "
                "for a PSD panel export the DFT, HK or BJH *tabular* report as well.")
    if psds:
        st.subheader("Pore-size distribution (separate panel)")
        fig, ax = newfig(); draw_psd(ax)
        show(fig, "DFT_PSD", pd.concat([df.assign(sample=n) for n, df in psds]), draw_psd, "PSD", units={"x": "pore_nm"})
    with st.expander("📈 BET / Rouquerol plots and t-plot (micropore volume) – click to open", expanded=False):
        sel = st.selectbox("BET plots for", [n for n in data if "bet" in data[n]])
        if sel:
            bet = data[sel]["bet"]
            withbet = [n for n in data if data[n].get("bet")]
            ov = [(n, data[n]["bet"], c) for n, c in zip(withbet, S.colours_for(withbet))] if (bet_all and len(withbet) > 1) else None
            fig, axs = C.compose([lambda ax: P.draw_rouquerol(ax, bet, name=sel, overlay=ov),
                                  lambda ax: P.draw_bet(ax, bet, sel, overlay=ov, box=bet_box)],
                                 ncols=2, width=max(WIDTH, 5.0), panel_aspect=0.8, wspace=0.5)
            show(fig, "BET_all" if ov else f"BET_{re.sub(r'[^A-Za-z0-9]+', '_', sel)}")
            st.subheader("t-plot (Harkins–Jura) – micropore volume")
            trng = st.slider("t range for the linear fit (Å)", 3.0, 10.0, (3.5, 5.0), 0.1, key="sorp_trng")
            a_ = data[sel]["iso"][data[sel]["iso"].branch == "ads"]
            tp = A.t_plot(a_.p_rel.values, a_.q.values, trng)
            if tp:
                st.write(f"**{sel}**: V_micro = **{tp['V_micro']:.3f} cm³ g⁻¹** (intercept {tp['intercept']:.1f} cm³ STP g⁻¹), "
                         f"S_ext = {tp['S_ext']:.0f} m² g⁻¹ (t {trng[0]}–{trng[1]} Å). A negative intercept means no micropores in this t window.")

                def draw_t(ax, tp=tp, trng=trng):
                    P.draw_tplot(ax, tp, trng, name=sel)
                fig, ax = newfig(); draw_t(ax)
                show(fig, f"tplot_{safe_name(sel)}", pd.DataFrame({"t (A)": tp["t"], "q (cm3/g STP)": tp["q"]}), draw_t, "t-plot")
            else:
                st.info("Fewer than 3 points in this t range.")
    st.subheader("Pore-size distribution computed from the isotherm")
    from mfs import psd as PSDm
    METH = list(PSDm.KERNELS) + ["Horváth–Kawazoe slit (carbon parameters)", "Saito–Foley cylinder (oxide parameters)",
                                 "BJH, desorption branch (mesopores)", "BJH, adsorption branch (mesopores)"]
    c1, c2, c3 = st.columns([3, 1, 1])
    meth = c1.selectbox("Method", METH, key="sorp_psdm",
                        help="The DFT kernel is the openly licensed pyGAPS N2/77 K carbon-slit kernel (MIT) – no free NLDFT/QSDFT "
                             "kernel exists for oxides or cylindrical pores, so for those Saito–Foley (micropores) and BJH (mesopores) are "
                             "offered. Commercial NLDFT/QSDFT PSDs: import the instrument's DFT tabular report (panel above).")
    pkind = c2.selectbox("Show", ["dV/dw", "cumulative"], key="sorp_psdk")
    wmax2 = c3.number_input("max width (nm)", 0.5, 100.0, 5.0 if "BJH" not in meth else 50.0, key="sorp_psdw_" + ("b" if "BJH" in meth else "m"))
    comp = []
    for nm, d in data.items():
        iso = d["iso"]; a_ = iso[iso.branch == "ads"]; de_ = iso[iso.branch == "des"]
        try:
            if meth in PSDm.KERNELS:
                r_ = PSDm.dft_fit(a_.p_rel.values, a_.q.values, meth)
            elif meth.startswith("Horv"):
                r_ = PSDm.hk_psd(a_.p_rel.values, a_.q.values, "slit", "Carbon (HK)")
            elif meth.startswith("Saito"):
                r_ = PSDm.hk_psd(a_.p_rel.values, a_.q.values, "cylinder", "AlSi oxide ion (SF)")
            else:
                br = de_ if ("desorption" in meth and len(de_) > 5) else a_
                r_ = PSDm.bjh_psd(br.p_rel.values, br.q.values, "desorption" if br is de_ else "adsorption")
            comp.append((lab_[nm], r_))
        except Exception as e:
            st.warning(f"{nm}: PSD failed ({e})")
    if comp:
        st.caption("Method: **" + comp[0][1]["method"] + "**" + (f" · fit residual " + ", ".join(f"{n}: {100*r['residual_rel']:.1f} %" for n, r in comp)
                   if "residual_rel" in comp[0][1] else "") + ". Shown on the figure legend title.")
        cols_ = S.colours_for([n for n, _ in comp])

        def draw_cpsd(ax, comp=comp, pkind=pkind, wmax2=wmax2, cols_=cols_):
            for (n, r_), col in zip(comp, cols_):
                yv = r_["dVdw"] if pkind == "dV/dw" else r_["V_cum"]
                ax.plot(r_["w"], yv, "-", color=col, lw=1.0, label=n)
            ax.set_xlim(0, wmax2); ax.set_ylim(0, None)
            ax.set_xlabel("Pore width (nm)")
            ax.set_ylabel("d$V$/d$w$ (cm$^3$ g$^{-1}$ nm$^{-1}$)" if pkind == "dV/dw" else "Cumulative pore volume (cm$^3$ g$^{-1}$)")
            ax.legend(frameon=False, title=comp[0][1]["method"].split(" (")[0], title_fontsize=FONT - 2, fontsize=FONT - 2)
        fig, ax = newfig(); draw_cpsd(ax)
        show(fig, "PSD_computed", pd.concat([pd.DataFrame({"sample": n, "w (nm)": r_["w"], "dV/dw": r_["dVdw"], "V_cum": r_["V_cum"],
                                                            "method": r_["method"]}) for n, r_ in comp]), draw_cpsd, "PSD (computed)", units={"x": "pore_nm"})
        if meth in PSDm.KERNELS:
            with st.expander("Kernel fit of the isotherm (check the fit quality)"):
                fig, ax = newfig()
                for (n, r_), col in zip(comp, cols_):
                    ax.semilogx(r_["p"], r_["q"], "o", ms=2.5, mfc="white", color=col, label=n)
                    ax.semilogx(r_["p"], r_["fit"], "-", color=col, lw=0.9)
                ax.set_xlabel("$p/p_0$"); ax.set_ylabel("Quantity adsorbed (cm$^3$ g$^{-1}$ STP)"); ax.legend(frameon=False, fontsize=FONT - 2)
                show(fig, "PSD_kernel_fit")
