"""Reticula page module: SDS-PAGE, Zeta / DLS and XPS pages.

Moved verbatim out of app.py in v2.13. Names that the code takes from the app (helpers such as show(), dl(), files_input(),
the sidebar values WIDTH / ASPECT / FONT, the cached reader wrapper R …) are listed in INJECTED and filled in by
mfs.pages.bind(globals()) on every run of app.py, so the code behaves exactly as it did inside app.py."""
# ruff: noqa: F821
from __future__ import annotations
import os, re, glob
import numpy as np
import pandas as pd
import streamlit as st
from mfs import numentry as NE  # every slider takes a typed number too (2.19.6)
from mfs import levels as LV
import matplotlib.pyplot as plt

from mfs import plots as P, style as S

INJECTED = (
    'EX',
    'R',
    '_F',
    '_apply_sample_list',
    'dl',
    'files_input',
    'newfig',
    'safe_name',
    'sample_names_list',
    'sample_table',
    'show',
    'store_put',
)


def page_sds():
    st.header("SDS-PAGE – Bio-Rad Gel Doc EZ / Image Lab")
    st.caption("Image Lab: Analysis Table → 'Lane and Bands' → Export → CSV/Excel. Add a time column (or map lane → time below).")
    files = files_input("Image Lab export", ["csv", "txt", "xlsx", "xls"], "SDS-PAGE/*.csv", "sds", multiple=False)
    if not files:
        return
    b = R.read_imagelab(files[0], files[0].name)
    if "time" not in b or b["time"].isna().all():
        lanes = sorted(b.lane.unique())
        lt = st.data_editor(pd.DataFrame({"lane": lanes, "time": [float(i) for i in range(len(lanes))]}), key="lt")
        b = b.merge(lt, on="lane")
    renorm = st.checkbox("Re-normalise band % to 100 % per lane", True)
    if renorm:
        b["band_pct"] = 100 * b.band_pct / b.groupby("lane").band_pct.transform("sum")
    st.dataframe(b, width="stretch")
    tunit = st.selectbox("time unit", ["h", "min", "d"])
    def draw(ax, b=b):
        P.draw_bands(ax, b, time_label=f"Time ({tunit})")
    fig, ax = newfig(); draw(ax)
    show(fig, "SDS-PAGE_bands", b, draw, "SDS-PAGE", units={"x": "time_h"})
    st.caption("Tip: fit the intact-protein band with the Kinetics tab (first-order) to report k_obs for HEWL.")


def page_zeta():
    st.header("Zeta potential / DLS – Malvern Zetasizer")
    st.caption("Zetasizer Ultra / Nano: export the records with Sample Name, Z-Average, PDI, Zeta Potential, Zeta Deviation and the size distribution "
               "(Sizes / Intensities) as **CSV or Excel (.xlsx / .xls)** and drop the file here. "
               "ℹ️ `.dts` files (Zetasizer software) are read for the sample names, the size distribution, Z-average / PDI and – for zeta records – "
               "the ζ-potential with its deviation, conductivity, count rate, the zeta peaks and the zeta distribution (software 8.x).")
    with st.expander("⚡ Zeta potential from the electrophoretic mobility (Henry / Smoluchowski)"):
        from mfs import zeta_conv as _ZC
        st.caption("v2.18.12. Enter the mobility from a laser-Doppler or phase-analysis light-scattering instrument. "
                   "Smoluchowski needs a large particle; the Henry form covers small ones when you give the radius.")
        _c1, _c2, _c3, _c4 = st.columns(4)
        _mob = _c1.number_input("mobility (µm·cm/(V·s))", value=-2.0, step=0.1, format="%.3f", key="zm_mob")
        _ion = _c2.number_input("ionic strength (mM)", value=1.0, min_value=0.001, step=0.1, format="%.3f", key="zm_ion")
        _rad = _c3.number_input("particle radius (nm, 0 = unknown)", value=0.0, min_value=0.0, step=1.0, key="zm_rad")
        _tmp = _c4.number_input("temperature (°C)", value=25.0, step=0.5, key="zm_tmp")
        _zr = _ZC.zeta_from_mobility(_mob, _ion, radius_nm=_rad or None, temp_C=_tmp, model="auto")
        if _zr.get("status") == "ok":
            st.markdown(f"**ζ = {_zr['zeta_mV']:.2f} mV** (model: {_zr['model']}; Smoluchowski value {_zr['zeta_smoluchowski_mV']:.2f} mV; "
                        f"Debye length {_zr['debye_length_nm']:.2f} nm)" if _zr.get("debye_length_nm") is not None else
                        f"**ζ = {_zr['zeta_mV']:.2f} mV** (model: {_zr['model']})")
            st.caption(_zr["note"])
        else:
            st.info(_zr.get("status"))
    files = files_input("Zetasizer export (CSV / Excel / .dts)", ["csv", "txt", "xlsx", "xls", "dts"], "Zeta-DLS/*.csv", "zeta", multiple=False)
    if not files:
        return
    z = dict(R.read_zetasizer(files[0], files[0].name))
    if z.get("note"):
        st.info(z["note"])
    if "rec" in z["summary"].columns and len(z["summary"]) > 1:        # v2.14: every run (record) can be switched off; means and graphs follow
        _s0 = z["summary"]
        _tab = pd.DataFrame({"use": True, "sample": _s0["sample"].astype(str), "run": _s0["run"] if "run" in _s0 else _s0["rec"] + 1,
                             "Z-avg (nm)": _s0["z_avg"].round(1), "PDI": _s0["pdi"].round(3), "ζ (mV)": _s0["zeta"].round(1)})
        with st.expander(f"🔁 Runs / records ({len(_tab)}) – untick a run to leave it out of the averages and graphs", expanded=len(_tab) <= 12):
            _ed = st.data_editor(_tab, key="zeta_runs|" + files[0].name, hide_index=True, width="stretch", height=min(360, 38 + 35 * len(_tab)),
                                 disabled=["sample", "run", "Z-avg (nm)", "PDI", "ζ (mV)"])
            _keep = _s0["rec"][_ed["use"].to_numpy()].tolist()
            st.caption(f"{len(_keep)} of {len(_tab)} runs used.")
        if not _keep:
            st.info("No run selected."); return
        z = R.select_records(z, _keep)
    _sel = sample_names_list(list(pd.unique(z["summary"]["sample"].astype(str))), "zeta_samples")
    z["summary"] = _apply_sample_list(z["summary"], "sample", None, _sel)
    if z.get("dist") is not None:
        z["dist"] = _apply_sample_list(z["dist"], "sample", None, _sel)
    _agg = dict(zeta=("zeta", "mean"), zeta_sd=("zeta", "std"), z_avg=("z_avg", "mean"), pdi=("pdi", "mean"), n_runs=("sample", "size"))
    for _c, _n_ in (("cond_mS_cm", "conductivity_mS_cm"), ("count_kcps", "count_rate_kcps"), ("zeta_sd", "zeta_deviation_mV")):   # v2.18.2: from .dts zeta records
        if _c in z["summary"].columns and z["summary"][_c].notna().any() and _c != "zeta_sd":
            _agg[_n_] = (_c, "mean")
    if "zeta_sd" in z["summary"].columns and z["summary"]["zeta_sd"].notna().any():
        _agg["zeta_deviation_mV"] = ("zeta_sd", "mean")        # width (SD) of the zeta distribution of one record, averaged; zeta_sd above = spread between records
    zs = z["summary"].groupby("sample", sort=False).agg(**_agg).reset_index()
    st.dataframe(zs, width="stretch")
    for _r in zs.itertuples():                                   # v2.16: the cross-technique checks compare DLS with the image sizes
        store_put(str(_r.sample), "dls", {"z_avg_nm": None if pd.isna(_r.z_avg) else float(_r.z_avg), "pdi": None if pd.isna(_r.pdi) else float(_r.pdi),
                                          "zeta_mV": None if pd.isna(_r.zeta) else float(_r.zeta), "n_runs": int(_r.n_runs)}, "Zeta / DLS page")
    zs["label"] = zs["sample"].str.replace("ZrO2", "ZrO$_2$")
    zs_sd = z["summary"].groupby("sample", sort=False).agg(z_avg_sd=("z_avg", "std"), pdi_sd=("pdi", "std")).reset_index()
    zs = zs.merge(zs_sd, on="sample", how="left")
    from mfs import methods as MT
    _n = z["summary"].groupby("sample", sort=False).size().rename("n").reset_index()
    _q = zs.merge(_n, on="sample", how="left")
    for lvl_, txt_ in MT.dls_quality_notes(_q.to_dict("records")):
        (st.warning if lvl_ == "warn" else st.caption)(("⚠️ " if lvl_ == "warn" else "ℹ️ ") + txt_)
    st.markdown("**Show in the figures** (each text can be dragged anywhere with *🏷️ Move / edit labels on the figure* under the figure)")
    c1, c2, c3, c4 = st.columns(4)
    show_z = c1.checkbox("Z-average", False, key="zeta_show_z", help="Z-average hydrodynamic diameter (d.nm) per sample")
    show_p = c2.checkbox("PDI", False, key="zeta_show_pdi", help="Polydispersity index per sample")
    show_zp = c3.checkbox("ζ-potential", False, key="zeta_show_zp", help="Zeta potential (mV) per sample in the DLS figure")
    with_sd = c4.checkbox("± SD (records)", True, key="zeta_show_sd", help="Mean ± standard deviation over the records of a sample (needs ≥ 2)")
    c1, c2 = st.columns(2)
    show_leg_z = c1.checkbox("Legend in the zeta-potential figure", False, key="zeta_bar_leg")
    _annot = show_z or show_p or show_zp
    # with the coloured texts on, the legend is redundant (each line names its sample): off by default then, tick it to bring it back
    show_leg_d = c2.checkbox("Legend in the DLS figure", not _annot, key="zeta_dls_leg_" + ("a" if _annot else "n"))

    def _fmt(v, sd, nd, unit=""):
        if not np.isfinite(v):
            return "–"
        s = f"{v:.{nd}f}" + (f" ± {sd:.{nd}f}" if with_sd and np.isfinite(sd) else "")
        return s + (f" {unit}" if unit else "")

    def draw_zeta(ax, zs=zs):
        cols = S.colours_for(zs["label"]); xs = np.arange(len(zs))
        ax.bar(xs, zs.zeta, 0.65, yerr=zs.zeta_sd, color=cols, edgecolor="k", lw=0.6, error_kw=dict(lw=0.7, capsize=2.5))
        ax.axhline(0, color="k", lw=0.6); ax.set_xticks(xs, zs.label, rotation=30, ha="right")
        ax.tick_params(axis="x", which="minor", bottom=False, top=False); ax.set_ylabel(r"$\zeta$-potential (mV)")
        if show_leg_z:
            from matplotlib.patches import Patch
            ax.legend([Patch(facecolor=c, edgecolor="k", lw=0.6) for c in cols], list(zs.label), loc="best")
    fig, ax = newfig(); draw_zeta(ax)
    show(fig, "zeta", zs, draw_zeta, "Zeta potential")
    if z.get("dist") is not None:
        kinds = [(k_, lab_) for k_, lab_ in (("intensity", "Intensity"), ("volume", "Volume"), ("number", "Number"))
                 if k_ in z["dist"].columns and z["dist"][k_].notna().any()]
        kind = st.radio("Size distribution by", [k_ for k_, _ in kinds], horizontal=True, key="zeta_dist_kind",
                        format_func=dict(kinds).get) if len(kinds) > 1 else "intensity"
        dd = z["dist"].groupby(["sample", "size_nm"], sort=False)[kind].mean().reset_index().rename(columns={kind: "value"})
        tr = [(s_.replace("ZrO2", "ZrO$_2$"), g_.size_nm.values, g_.value.values) for s_, g_ in dd.groupby("sample", sort=False)]
        info = {r_.label: r_ for r_ in zs.itertuples()}

        def draw_dls(ax, tr=tr):
            P.draw_overlay(ax, tr, "Hydrodynamic diameter (nm)", f"{dict(kinds).get(kind, 'Intensity')} (%)", logx=True, legend=show_leg_d); ax.set_ylim(bottom=0)
            if show_z or show_p or show_zp:           # v2.14: per-sample Z-average / PDI / ζ texts (movable like every label)
                cols_ = S.colours_for([t[0] for t in tr]); k = 0
                for c_, (nm_, _x, _y) in zip(cols_, tr):
                    r_ = info.get(nm_)
                    if r_ is None:
                        continue
                    parts = []
                    if show_z:
                        parts.append("$Z_{avg}$ = " + _fmt(r_.z_avg, r_.z_avg_sd, 0, "nm"))
                    if show_p:
                        parts.append("PDI = " + _fmt(r_.pdi, r_.pdi_sd, 2))
                    if show_zp:
                        parts.append(r"$\zeta$ = " + _fmt(r_.zeta, r_.zeta_sd, 1, "mV"))
                    ax.text(0.98, 0.97 - 0.075 * k, f"{nm_}: " + ", ".join(parts), transform=ax.transAxes, ha="right", va="top", color=c_,
                            fontsize=plt.rcParams["font.size"] - 2.5)
                    k += 1
                lo_, hi_ = ax.get_ylim()                    # head-room so the texts do not sit on the peaks
                ax.set_ylim(0, hi_ * (1.0 + 0.11 * k + 0.04))
        fig, ax = newfig(); draw_dls(ax)
        show(fig, "DLS", dd, draw_dls, "DLS")
    if z.get("zdist") is not None and len(z["zdist"]):                   # v2.18.2: zeta potential distribution (counts vs apparent zeta potential)
        zd = z["zdist"].groupby(["sample", "zeta_mV"], sort=False)["counts"].mean().reset_index()
        trz = [(s_.replace("ZrO2", "ZrO$_2$"), g_.zeta_mV.values, g_.counts.values) for s_, g_ in zd.groupby("sample", sort=False)]

        def draw_zd(ax, trz=trz):
            P.draw_overlay(ax, trz, "Apparent zeta potential (mV)", "Total counts", legend=True); ax.set_ylim(bottom=0)
        fig, ax = newfig(); draw_zd(ax)
        show(fig, "zeta_distribution", zd, draw_zd, "Zeta potential distribution")
        st.caption("The distribution is shown for the records whose stored block reproduces the zeta potential of the software; other records have none.")
    with st.expander("📝 Methods / SI text (edit the [brackets])", expanded=False):
        _txt = MT.dls_methods(_q.to_dict("records"), with_sd=with_sd)
        st.caption("Generated from this data only; fill in instrument, temperature, medium and concentration yourself.")
        _txt = st.text_area("Text", _txt, height=170, key="zeta_methods_txt")
        dl(st, "⬇️ Save as .txt", _txt.encode("utf-8"), "DLS_zeta_methods.txt", "zeta_methods_dl")


def page_xps():
    from mfs import xps as XP
    st.header("XPS")
    st.caption("CasaXPS / Avantage exported region (BE, counts) as CSV/TXT; simple VAMAS (.vms) regions. Shirley or Tougaard "
               "background, doublet fitting (Zr 3d, Hf 4f, S 2p with fixed spin–orbit splitting and area ratio), atomic-% with RSFs.")
    multi = st.checkbox("➕ Use the multi-region EXAMPLE instead (ZrO₂@Sample B: Zr 3d, Hf 4f, S 2p, O 1s, C 1s – synthetic)", False, key="xps_multi")
    files = files_input("Spectra", ["csv", "txt", "vms", "dat"], "XPS/*Hf4f*.csv", "xps")
    if multi:
        files = [_F(p_) for p_ in sorted(glob.glob(os.path.join(EX, "XPS", "SampleC_*.csv")))]
    if not files:
        return
    samples = sample_table(files, "xps_tab" + ("_m" if multi else ""))
    spec = []
    for nm, f in samples:
        try:
            d = R.read_xps(f, f.name)
        except Exception as e:
            st.error(f"{f.name}: {e}"); continue
        reg = XP.guess_region(d.x.values)
        grp = re.sub(r"\s*\(\d+\)$", "", nm)
        # v2.18.5: '\b' failed for names like 'Sample_S2p_CasaXPS' ('_' is a word character), so every region became its own sample
        grp = re.sub(r"[ _-]*(Zr ?3d|Hf ?4[fd]|S ?2p|O ?1s|C ?1s|N ?1s|F ?1s|Ti ?2p)(?![A-Za-z0-9])", "", grp, flags=re.I).strip(" _-") or nm
        spec.append({"name": nm, "group": grp, "region": reg, "x": d.x.values, "y": d.y.values})
    st.dataframe(pd.DataFrame([{"spectrum": d_["name"], "sample": d_["group"], "region (from BE range)": d_["region"],
                                "BE range (eV)": f"{d_['x'].min():.1f}–{d_['x'].max():.1f}"} for d_ in spec]), hide_index=True, width="stretch")
    regions = sorted({d_["region"] for d_ in spec if d_["region"]})
    c1, c2, c3 = st.columns(3)
    show_reg = c1.selectbox("Overview: region", regions or ["(all)"], key="xps_ovreg")
    bgk = c2.radio("Background", ["Shirley", "Tougaard", "linear", "none"], horizontal=True, key="xps_bg")
    region = c3.text_input("Region label", show_reg if regions else "Hf 4f", key="xps_lab_" + str(show_reg))
    tr = [(d_["name"], d_["x"], d_["y"]) for d_ in spec if not regions or d_["region"] == show_reg]
    lo_, hi_ = (min(t[1].min() for t in tr), max(t[1].max() for t in tr)) if tr else (13.0, 23.0)
    lo, hi = NE.slider("Binding energy (eV)", 0.0, 1200.0, (float(np.floor(lo_)), float(np.ceil(hi_))), 0.5, key="xps_rng_" + str(show_reg))

    def draw(ax, tr=tr):
        P.draw_xps(ax, tr, (lo, hi), bgk == "Shirley"); ax.text(0.03, 0.96, region, transform=ax.transAxes, va="top")
    fig, ax = newfig(); draw(ax)
    show(fig, "XPS", None, draw, "XPS")

    st.subheader("Peak fitting")
    names = [d_["name"] for d_ in spec]
    c1, c2, c3 = st.columns(3)
    pick = c1.selectbox("Spectrum", names, key="xps_fit_pick", format_func=lambda n_: re.sub(r"\$_(\d)\$", lambda m_: chr(0x2080 + int(m_.group(1))), n_))
    d_ = spec[names.index(pick)]
    regs = list(XP.DOUBLETS) + list(XP.SINGLETS) + ["other (singlets)"]
    rsel = c2.selectbox("Region / doublet constraint", regs, index=regs.index(d_["region"]) if d_["region"] in regs else len(regs) - 1,
                        key="xps_fit_reg_" + pick)
    eta = c3.slider("Lorentzian fraction (GL mix)", 0.0, 1.0, 0.3, 0.05, key="xps_eta")
    dd = XP.DOUBLETS.get(rsel)
    if dd:
        st.caption(f"{rsel}: {dd['major']} + {dd['minor']} at +{dd['split']} eV, area {dd['minor']}/{dd['major']} = {dd['ratio']:.3f}, same FWHM.")
    x, y = d_["x"], d_["y"]
    bg0 = XP.shirley(x, y)
    xpk = float(x[np.argmax(y - bg0)])
    defaults = {"Zr 3d": [("Zr–O (MOF / ZrO₂)", xpk)], "Hf 4f": [("Hf–O", xpk)], "S 2p": [("sulfonate / sulfate", xpk)],
                "O 1s": [("O lattice (M–O)", 530.2), ("O carboxylate / OH", 531.8)], "C 1s": [("C–C / C–H", 284.8), ("O–C=O", 288.6)]}
    comp0 = pd.DataFrame([{"name": n_, "BE guess (eV)": round(c_, 2), "doublet": bool(dd), "fix BE": False}
                          for n_, c_ in defaults.get(rsel, [("component 1", xpk)])])
    comps = st.data_editor(comp0, num_rows="dynamic", key="xps_comp_" + pick + rsel, hide_index=True)
    try:
        fr = XP.fit(x, y, [{"name": r_["name"], "center": float(r_["BE guess (eV)"]), "region": rsel if (r_["doublet"] and dd) else None,
                            "fix_center": bool(r_["fix BE"])} for _, r_ in comps.dropna(subset=["BE guess (eV)"]).iterrows()],
                    background=bgk.lower(), eta=eta)
        tab = pd.DataFrame(fr["rows"]); st.dataframe(tab.style.format(precision=2), hide_index=True, width="stretch")
        st.caption(f"Background: {bgk}{' (universal cross-section B·T/(C+T²)², C = 1643 eV², B scaled to the high-BE end)' if bgk == 'Tougaard' else ''}. "
                   f"Residual RMS = {100 * fr['rms_rel']:.1f} % of the peak height.")
        st.session_state.setdefault("xps_areas", {}).setdefault(d_["group"], {})[rsel] = float(tab["area (total)"].sum())

        def draw_f(ax, fr=fr):
            XP.draw_fit(ax, fr); ax.legend(frameon=False, fontsize=6); ax.text(0.03, 0.96, f"{pick}\n{rsel}", transform=ax.transAxes, va="top", fontsize=7)
        fig, ax = newfig(0.85); draw_f(ax)
        show(fig, f"XPS_fit_{safe_name(pick)}", pd.DataFrame({"BE": fr["x"], "counts": fr["y"], "background": fr["bg"], "envelope": fr["fit"]}),
             draw_f, f"XPS fit {pick}")
    except Exception as e:
        st.warning(f"Fit failed: {e}")

    with LV.gate("Advanced", globals()):
        with st.expander("🗂️ Batch: apply this fit model to ALL spectra of the same region", expanded=False):
            same = [s_ for s_ in spec if s_["region"] == d_["region"] or d_["region"] is None]
            st.caption(f"{len(same)} spectrum/spectra of region **{d_['region'] or 'other'}**; the components above (names, BE guesses, doublet, fixed BE) "
                       "are fitted to each of them with the background and Lorentzian fraction chosen above.")
            if st.button("▶ Fit all", key="xps_batch_go", type="primary"):
                comp_def = [{"name": r_["name"], "center": float(r_["BE guess (eV)"]), "region": rsel if (r_["doublet"] and dd) else None,
                             "fix_center": bool(r_["fix BE"])} for _, r_ in comps.dropna(subset=["BE guess (eV)"]).iterrows()]
                out, bad = [], []
                for s_ in same:
                    try:
                        o_ = XP.fit(s_["x"], s_["y"], comp_def, background=bgk.lower(), eta=eta)
                        for r_ in o_["rows"]:
                            out.append({"sample": s_["group"], "spectrum": s_["name"], **{k: v for k, v in r_.items() if k != "region"},
                                        "RMS (% of peak)": 100 * o_["rms_rel"]})
                    except Exception as e:
                        bad.append(f"{s_['name']}: {e}")
                st.session_state["xps_batch_res"] = (pd.DataFrame(out), bad)
            got = st.session_state.get("xps_batch_res")
            if got:
                tb, bad = got
                for b_ in bad:
                    st.warning(b_)
                if len(tb):
                    st.dataframe(tb.style.format(precision=2), hide_index=True, width="stretch")
                    dl(st, "⬇️ Batch fit table (CSV)", lambda: tb.to_csv(index=False).encode("utf-8-sig"), "XPS_batch_fit.csv", "xps_batch_csv")
                    pivot = tb.pivot_table(index="spectrum", columns="component", values="% of envelope", aggfunc="sum", sort=False)

                    def draw_b(ax, pivot=pivot):
                        pivot.plot(kind="bar", stacked=True, ax=ax, width=0.7, edgecolor="none")
                        ax.set_ylabel("Share of the envelope (%)"); ax.set_xlabel("")
                        ax.legend(frameon=False, fontsize=6)
                        plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
                    fig, ax = newfig(0.8); draw_b(ax)
                    show(fig, "XPS_batch_components", pivot.reset_index(), draw_b, "XPS components per sample")

    st.subheader("Atomic % (RSF-corrected)")
    st.caption("Areas above the chosen background, per region, per sample (whole doublets). Default RSFs: Scofield cross-sections "
               "(Al Kα, C 1s = 1) without transmission/IMFP correction – replace them with your instrument library's RSFs.")
    rsf = st.data_editor(pd.DataFrame({"region": list(XP.RSF_SCOFIELD), "RSF": list(XP.RSF_SCOFIELD.values())}), key="xps_rsf", hide_index=True)
    rsfd = dict(zip(rsf.region, rsf.RSF))
    areas = {}
    for s_ in spec:
        if not s_["region"]:
            continue
        b_ = XP.shirley(s_["x"], s_["y"]) if bgk == "Shirley" else (XP.tougaard(s_["x"], s_["y"])[0] if bgk == "Tougaard" else
                                                                      np.interp(s_["x"], [s_["x"][0], s_["x"][-1]], [s_["y"][:3].mean(), s_["y"][-3:].mean()]))
        o = np.argsort(s_["x"]); A_ = float(np.trapezoid(np.clip((s_["y"] - b_)[o], 0, None), s_["x"][o]))
        areas.setdefault(s_["group"], {})[s_["region"]] = A_
    if areas:
        at = []
        for g_, ar_ in areas.items():
            pc = XP.atomic_percent(ar_, rsfd)
            row = {"sample": g_}; row.update({f"{k} at.%": v for k, v in pc.items()})
            M_ = {k: pc.get(k, 0) for k in ("Zr 3d", "Hf 4f")}
            if sum(M_.values()):
                row["Hf/(Zr+Hf)"] = M_["Hf 4f"] / sum(M_.values())
                if "S 2p" in pc:
                    row["S/(Zr+Hf)"] = pc["S 2p"] / sum(M_.values())
            at.append(row)
            store_put(g_, "xps", {k: v for k, v in row.items() if k != "sample"}, "XPS page")
        att = pd.DataFrame(at); st.dataframe(att.style.format(precision=2), hide_index=True, width="stretch")
        dl(st, "CSV atomic %", att.to_csv(index=False).encode("utf-8-sig"), "XPS_atomic_percent.csv", "xps_at_csv")
