"""Reticula page module: Kinetics and adsorption page.

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
    'EX',
    'FONT',
    'R',
    '_F',
    '_apply_sample_list',
    'files_input',
    'newfig',
    'show',
    'store_put',
)


def page_kinetics():
    st.header("Hydrolysis kinetics (¹H NMR) and adsorption isotherms")
    st.caption("Table with columns sample, time, concentration (or integrals vs TMSP-d4). Replicates are averaged (mean ± SD) and all points are fitted.")
    files = files_input("Kinetics table (CSV/XLSX)", ["csv", "txt", "xlsx"], "Kinetics/GlyGly_hydrolysis_60C*.csv", "kin", multiple=False)
    if st.checkbox("➕ Use the temperature-series EXAMPLE instead (40–70 °C, for Arrhenius / Eyring)", False, key="kin_tser"):
        files = [_F(os.path.join(EX, "Kinetics", "GlyGly_hydrolysis_Tseries_EXAMPLE.csv"))]
    if st.session_state.get("_uv_kin_csv") and st.checkbox("➕ Use the table sent from UV/Vis (A at one λ vs time)", True, key="kin_from_uv"):
        _b = io.BytesIO(st.session_state["_uv_kin_csv"]); _b.name = "uvvis_kinetics.csv"; files = [_b]
        st.caption(f"Time unit from UV/Vis: {st.session_state.get('_uv_kin_unit', 'min')} – set the same unit below.")
    from mfs import kinetics as K
    if files:
        df = R.read_csv_table(files[0], (r"sample", r"time"))
        st.dataframe(df.head(), width="stretch")
        cols = list(df.columns)
        c1, c2, c3, c4 = st.columns(4)
        sc = c1.selectbox("sample column", cols, index=cols.index(R.find_col(df, "sample") or cols[0]))
        tc = c2.selectbox("time column", cols, index=cols.index(R.find_col(df, "time") or cols[1]))
        df = _apply_sample_list(df, sc, "kin_samples")
        mode = c3.selectbox("y is", ["concentration", "integral vs internal standard"])
        yc = c4.selectbox("y column", cols, index=len(cols) - 1)
        if mode != "concentration":
            c1, c2, c3 = st.columns(3)
            stdc = c1.selectbox("standard integral column", cols, index=len(cols) - 1)
            nh = c2.number_input("protons in peptide signal", 1, 20, 4)
            cstd = c3.number_input("[TMSP-d4] (mM)", 0.01, 100.0, 1.0)
            df["_c"] = A.nmr_concentration(df[yc], nh, df[stdc], 9, cstd); yc = "_c"
        c1, c2, c3, c4 = st.columns(4)
        tunit = c1.selectbox("time unit", ["h", "min", "s", "d"])
        offset = c2.checkbox("fit with plateau C∞", False)
        norm = c3.checkbox("plot C/C₀", True)
        errbar = c4.radio("error bars on k", ["95 % CI", "1 SE"], horizontal=True, key="kin_eb")
        Tcol = R.find_col(df, r"^t\s*\(°?c\)|temp|^t_c$"); pDcol = R.find_col(df, r"^pd\b|^ph\b|pd$")
        gcols = [sc] + [c_ for c_ in (Tcol, pDcol) if c_ is not None and df[c_].nunique() > 1]
        sets, rows, fits = [], [], []
        for key_, g_ in df.groupby(gcols, sort=False):
            key_ = key_ if isinstance(key_, tuple) else (key_,)
            lab = str(key_[0]) + "".join(f", {v:g} °C" if c_ == Tcol else f", pD {v:g}" for c_, v in zip(gcols[1:], key_[1:]))
            try:
                fit = K.first_order_ci(g_[tc].values, g_[yc].values, offset=offset)
            except Exception as e:
                st.warning(f"{lab}: fit failed ({e})"); fit = None
            agg = g_.groupby(tc)[yc].agg(["mean", "std"]).reset_index()
            sets.append((lab, agg[tc].values, agg["mean"].values, agg["std"].fillna(0).values, fit))
            if fit:
                Tv = float(g_[Tcol].iloc[0]) if Tcol is not None else None; pDv = float(g_[pDcol].iloc[0]) if pDcol is not None else None
                fits.append({"label": lab, "sample": str(key_[0]), "T": Tv, "pD": pDv, "fit": fit})
                rows.append({"sample": lab, f"k_obs ({tunit}-1)": fit["k"], "± 95 % CI": fit["k_ci"], "SE": fit["k_se"],
                             f"t1/2 ({tunit})": fit["t_half"], "± t1/2 (95 %)": fit["t_half_ci"], "C0": fit["C0"], "R2": fit["R2"], "n": fit["n"]})
                store_put(str(key_[0]), "kinetics", {"k": fit["k"], "k_ci95": fit["k_ci"], "unit": f"{tunit}-1", "T_C": Tv, "pD": pDv,
                                                      "R2": fit["R2"], "label": lab}, "Kinetics page")
        tab = pd.DataFrame(rows); st.dataframe(tab.style.format(precision=4), width="stretch")
        st.caption("95 % CI = t(n − p, 0.975) × standard error from the fit covariance (replicates fitted point by point).")
        def draw(ax, sets=sets):
            P.draw_kinetics(ax, sets, ylabel=r"$C/C_0$" if norm else "[peptide] (mM)", time_label=f"Time ({tunit})", normalise=norm,
                            annotate_k=False, k_unit=f"{tunit}$^{{-1}}$")
        fig, ax = newfig(); draw(ax)
        show(fig, "kinetics", tab, draw, "Kinetics", units={"x": "time_h"})
        def drawk(ax, sets=sets):
            ok = [s for s in sets if s[4]]
            P.draw_k_bars(ax, [s[0] for s in ok], [s[4]["k"] for s in ok], [s[4]["k_ci" if errbar.startswith("95") else "k_se"] for s in ok],
                          ylabel=rf"$k_\mathrm{{obs}}$ (10$^{{-3}}$ {tunit}$^{{-1}}$)")
        fig, ax = newfig(); drawk(ax)
        show(fig, "k_obs_bars", None, drawk, "k_obs bars")
        with st.expander("⏱️ Initial rates (linear fit up to a chosen conversion)"):
            frac = st.slider("use points up to conversion", 0.02, 0.5, 0.10, 0.01, key="kin_ir_frac")
            ir = []
            for key_, g_ in df.groupby(gcols, sort=False):
                try:
                    r_ = K.initial_rate(g_[tc].values, g_[yc].values, frac)
                    ir.append({"sample": str(key_ if not isinstance(key_, tuple) else ", ".join(map(str, key_))), f"r0 (mM {tunit}-1)": r_["r0"],
                               "± 95 % CI": r_["r0_ci"], "points": r_["n"], f"up to t ({tunit})": r_["t_max"], "R2": r_["R2"]})
                except Exception as e:
                    st.caption(f"{key_}: {e}")
            st.dataframe(pd.DataFrame(ir).style.format(precision=5), width="stretch", hide_index=True)
        with st.expander("⚗️ TOF per open site (open sites per node from the defect calculator)"):
            from mfs import samplestore as SS
            dstore = {k_: v_ for k_, v_ in SS.ids(st.session_state).items() if SS.get(st.session_state, v_, "defects")}
            q1, q2, q3, q4 = st.columns(4)
            src = q1.selectbox("Catalyst (defect result)", ["(type values)"] + list(dstore), format_func=lambda k_: dstore.get(k_, k_), key="tof_src")
            dv = SS.get(st.session_state, dstore[src], "defects") if src in dstore else {}
            osn = q2.number_input("open sites per node", 0.0, 12.0, float(dv.get("OS") or 1.0), key="tof_os_" + str(src))
            Mn = q3.number_input("M per node (g/mol)", 100.0, 5000.0, float(dv.get("M") or 1650.0), key="tof_M_" + str(src))
            mcat = q4.number_input("catalyst (mg)", 0.01, 1000.0, 2.0, key="tof_m")
            q1, q2 = st.columns(2)
            c0 = q1.number_input("substrate C0 (mM)", 0.001, 1000.0, 2.0, key="tof_c0")
            vol = q2.number_input("volume (mL)", 0.01, 100.0, 0.6, key="tof_v")
            per_s = {"h": 3600, "min": 60, "s": 1, "d": 86400}[tunit]
            tr = []
            for f_ in fits:
                r_ = K.tof(f_["fit"]["k"] / per_s, c0, vol, mcat, osn, Mn, f_["fit"]["k_ci"] / per_s)
                tr.append({"sample": f_["label"], "open sites (µmol)": r_["n_open_sites_umol"], "TOF (h-1)": r_["TOF_h"],
                           "± 95 % (h-1)": r_["TOF_err_s"] * 3600 if r_["TOF_err_s"] else None})
            st.dataframe(pd.DataFrame(tr).style.format(precision=4), width="stretch", hide_index=True)
            st.caption("TOF = k_obs·C0·V / (m_cat / M_node × open sites per node) – initial turnover per open site. "
                       "Pull OS and M from the defect calculator (sample store) or type them.")
        with st.expander("ln(C/C₀) linearisation"):
            fig, ax = newfig()
            cols_ = S.colours_for([s[0] for s in sets])
            for c_, (nm, t, cc, e, fit) in zip(cols_, sets):
                ax.plot(t, np.log(np.asarray(cc) / (fit["C0"] if fit else cc[0])), "o", color=c_, ms=3, label=nm)
                if fit: ax.plot(t, -fit["k"] * np.asarray(t), "-", color=c_, lw=0.8)
            ax.set_xlabel(f"Time ({tunit})"); ax.set_ylabel(r"ln($C/C_0$)"); ax.legend(fontsize=FONT - 2)
            show(fig, "ln_kinetics")
    st.subheader("Adsorption isotherm (Langmuir / Freundlich)")
    files = files_input("Adsorption table (sample, Ce, qe)", ["csv", "txt", "xlsx"], "Kinetics/GlyGly_adsorption_*.csv", "ads", multiple=False)
    if files:
        df = R.read_csv_table(files[0], (r"sample", r"ce"))
        cols = list(df.columns)
        c1, c2, c3 = st.columns(3)
        sc = c1.selectbox("sample", cols, key="as"); cc = c2.selectbox("Ce", cols, index=1, key="ac"); qc = c3.selectbox("qe", cols, index=2, key="aq")
        sets, rows = [], []
        for s_, g_ in df.groupby(sc, sort=False):
            fit = A.langmuir_fit(g_[cc], g_[qc]); fr = A.freundlich_fit(g_[cc].values, g_[qc].values)
            sets.append((str(s_), g_[cc].values, g_[qc].values, fit))
            rows.append({"sample": s_, "q_max": fit["q_max"], "± q_max": fit["q_max_err"], "K_L": fit["K"], "± K_L": fit["K_err"], "R2 Langmuir": fit["R2"],
                         "K_F": fr["K_F"], "n_F": fr["n"]})
        st.dataframe(pd.DataFrame(rows).style.format(precision=4), width="stretch")
        def drawl(ax, sets=sets):
            P.draw_langmuir(ax, sets)
        fig, ax = newfig(); drawl(ax)
        show(fig, "langmuir", pd.DataFrame(rows), drawl, "Langmuir", units={"x": "conc_mM"})
    st.subheader("Arrhenius and Eyring analysis (several temperatures)")
    fT = [f_ for f_ in (locals().get("fits") or []) if f_["T"] is not None]
    if len({f_["T"] for f_ in fT}) >= 2:
        st.caption("Rates from the fitted kinetics table above (temperature column found). Edit if needed.")
        ey0 = pd.DataFrame({"T_C": [f_["T"] for f_ in fT], "k": [f_["fit"]["k"] for f_ in fT], "k_err (95 %)": [f_["fit"]["k_ci"] for f_ in fT]})
    else:
        ey0 = pd.read_csv(os.path.join(EX, "Kinetics", "Eyring_composite.csv")).rename(columns={"k_h-1": "k"})
        ey0["k_err (95 %)"] = np.nan
    ey = st.data_editor(ey0, num_rows="dynamic", key="ey2_" + str(len(ey0)))
    if len(ey) >= 2:
        unit = st.selectbox("k unit", ["h-1", "min-1", "s-1"], key="ar_unit")
        per = {"h-1": 3600, "min-1": 60, "s-1": 1}[unit]
        kerr = ey["k_err (95 %)"].values if ey["k_err (95 %)"].notna().all() else None
        ar = K.arrhenius(ey["T_C"].values, ey["k"].values / per, None if kerr is None else kerr / per)
        eyr = K.eyring(ey["T_C"].values, ey["k"].values / per, None if kerr is None else kerr / per)
        st.markdown(f"**Arrhenius**: Ea = **{ar['Ea_kJ']:.1f} ± {ar['Ea_ci']:.1f} kJ mol⁻¹** (95 % CI), ln(A/s⁻¹) = {ar['lnA']:.1f} ± {ar['lnA_ci']:.1f}, R² {ar['R2']:.4f}  \n"
                    f"**Eyring**: ΔH‡ = **{eyr['dH_kJ']:.1f} ± {eyr['dH_ci']:.1f} kJ mol⁻¹**, ΔS‡ = **{eyr['dS_J']:.0f} ± {eyr['dS_ci']:.0f} J mol⁻¹ K⁻¹**, "
                    f"ΔG‡(298 K) = {eyr['dG298_kJ']:.1f} kJ mol⁻¹")

        def draw_ar(ax, r=ar, lab=r"ln($k$/s$^{-1}$)"):
            ax.plot(1000 * r["X"], r["Y"], "o", color="#0072B2", ms=4)
            xx = np.linspace(r["X"].min(), r["X"].max(), 10); ax.plot(1000 * xx, r["slope"] * xx + r["intercept"], "-", color="#0072B2", lw=0.9)
            ax.set_xlabel(r"1000/$T$ (K$^{-1}$)"); ax.set_ylabel(lab)

        def draw_ey(ax, r=eyr):
            ax.plot(1000 * r["X"], r["Y"], "o", color="#D55E00", ms=4)
            xx = np.linspace(r["X"].min(), r["X"].max(), 10); ax.plot(1000 * xx, r["slope"] * xx + r["intercept"], "-", color="#D55E00", lw=0.9)
            ax.set_xlabel(r"1000/$T$ (K$^{-1}$)"); ax.set_ylabel(r"ln($k$/$T$)")
        fig, ax = newfig(); draw_ar(ax)
        show(fig, "arrhenius", ey, draw_ar, "Arrhenius")
    if len(ey) >= 2:
        fig, ax = newfig(); draw_ey(ax)
        show(fig, "eyring", ey, draw_ey, "Eyring")
    st.subheader("pD (pH)–rate profile")
    fpd = [f_ for f_ in (locals().get("fits") or []) if f_["pD"] is not None]
    if len({f_["pD"] for f_ in fpd}) >= 3:
        pdt0 = pd.DataFrame({"pD": [f_["pD"] for f_ in fpd], "k": [f_["fit"]["k"] for f_ in fpd], "± k": [f_["fit"]["k_ci"] for f_ in fpd]})
    else:
        pdt0 = pd.read_csv(os.path.join(EX, "Kinetics", "pD_rate_profile_EXAMPLE.csv"), comment="#")
        pdt0.columns = ["pD", "k", "± k"]
        st.caption("EXAMPLE values (synthetic) – replace with your own k_obs vs pD (or add a pD column to the kinetics table).")
    pdt = st.data_editor(pdt0, num_rows="dynamic", key="pdt_" + str(len(pdt0)))
    c1, c2 = st.columns(2)
    pmodel = c1.radio("Model", ["sigmoidal (one pKa)", "bell-shaped (two pKa)"], horizontal=True, key="pd_model")
    if len(pdt.dropna()) >= 3:
        try:
            pdt_ = pdt.dropna(subset=["pD", "k"])
            fr = K.fit_pd_profile(pdt_["pD"].values, pdt_["k"].values, "sigmoidal" if pmodel.startswith("sig") else "bell",
                                  pdt_["± k"].values if pdt_["± k"].notna().all() and (pdt_["± k"] > 0).all() else None)
            st.markdown(" · ".join(f"**{k_} = {v_:.3g} ± {fr['ci'][k_]:.2g}**" for k_, v_ in fr["params"].items()) + f" (95 % CI), R² {fr['R2']:.4f}. "
                        "pKa in D₂O; pKa(H₂O) ≈ pKa(D₂O) − 0.4–0.5 for carboxylic acids / Zr–OH₂.")

            def draw_pd(ax, pdt_=pdt_, fr=fr):
                ax.errorbar(pdt_["pD"], pdt_["k"], yerr=pdt_["± k"] if pdt_["± k"].notna().all() else None, fmt="o", color="k", ms=4, capsize=2, lw=0.8)
                xx = np.linspace(pdt_["pD"].min() - 0.3, pdt_["pD"].max() + 0.3, 200); ax.plot(xx, fr["curve"](xx), "-", color="#D55E00", lw=1.0)
                for k_, v_ in fr["params"].items():
                    if k_.startswith("pKa"):
                        ax.axvline(v_, color="#D55E00", lw=0.6, ls="--")
                ax.set_xlabel("pD"); ax.set_ylabel(r"$k_\mathrm{obs}$")
            fig, ax = newfig(); draw_pd(ax)
            show(fig, "pD_rate_profile", pdt_, draw_pd, "pD–rate profile")
        except Exception as e:
            st.warning(f"Fit failed: {e}")
