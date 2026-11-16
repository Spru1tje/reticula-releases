"""Reticula page module: ICP-OES / XRF page.

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
    'R',
    'WIDTH',
    'dl',
    'files_input',
    'newfig',
    'show',
    'store_put',
)


def page_icp():
    st.header("ICP-OES (PerkinElmer Optima 8300 / Syngistix) and XRF")
    kind = st.radio("Data", ["ICP-OES (Syngistix CSV)", "XRF / any table (sample, element, value, sd)"], horizontal=True, key="icp_kind")
    if kind.startswith("ICP"):
        files = files_input("Syngistix export", ["csv", "txt", "xlsx"], "ICP/*.csv", "icp", multiple=False)
        if not files:
            return
        df = R.read_icp(files[0], files[0].name)
        c1, c2, c3 = st.columns(3)
        mass = c1.number_input("Digested sample mass (mg)", 0.01, 1000.0, 10.0)
        vol = c2.number_input("Final volume (mL)", 0.1, 1000.0, 50.0)
        dil = c3.number_input("Extra dilution factor", 1.0, 1e4, 1.0)
        df["wt%"] = df.conc * dil * vol / 1000 / mass * 100
        df["sd wt%"] = df.sd * dil * vol / 1000 / mass * 100
        val, err, ylabel = "wt%", "sd wt%", "wt%"
    else:
        files = files_input("Table", ["csv", "txt", "xlsx"], "XRF/*.csv", "xrf", multiple=False)
        if not files:
            return
        df = R.read_csv_table(files[0], (r"sample", r"element"))
        cols = list(df.columns)
        c1, c2, c3, c4 = st.columns(4)
        sc = c1.selectbox("sample", cols); ec = c2.selectbox("element", cols, index=1); val = c3.selectbox("value", cols, index=2); err = c4.selectbox("error", cols, index=min(3, len(cols) - 1))
        df = df.rename(columns={sc: "sample", ec: "element"}); ylabel = val
    st.dataframe(df, width="stretch")
    els = sorted(df.element.dropna().unique())
    chosen = st.multiselect("Elements", els, default=els)
    layout = st.radio("Layout", ["one panel per element", "grouped bars"], horizontal=True)
    # molar ratios
    piv = df.pivot_table(index="sample", columns="element", values=val, aggfunc="mean")
    from mfs.chem import AW
    mol = piv.apply(lambda c: c / AW.get(c.name, np.nan))
    if {"Hf", "Zr"} <= set(mol.columns):
        piv["Hf/(Hf+Zr) mol"] = mol.Hf / (mol.Hf + mol.Zr)
    if "S" in mol.columns and ({"Hf"} <= set(mol.columns)):
        piv["S per M6 node"] = 6 * mol.S / (mol.Hf + (mol.Zr.fillna(0) if "Zr" in mol.columns else 0))
    st.dataframe(piv, width="stretch")
    ratio_vals = _icp_ratios(df, val, err, "icp" if kind.startswith("ICP") else "xrf")
    _icp_generic_ratios(df, val, err, "icp" if kind.startswith("ICP") else "xrf")
    for smp, r_ in piv.iterrows():
        vals = {f"{el}_wt" if kind.startswith("ICP") else f"{el}": (None if pd.isna(v) else float(v)) for el, v in r_.items()}
        if not kind.startswith("ICP") and {"S", "Hf"} <= set(mol.columns) and pd.notna(mol.loc[smp, "S"]):
            M_ = (mol.loc[smp, "Hf"] if pd.notna(mol.loc[smp, "Hf"]) else 0) + (mol.loc[smp, "Zr"] if "Zr" in mol.columns and pd.notna(mol.loc[smp, "Zr"]) else 0)
            vals["xrf_S_per_M"] = float(mol.loc[smp, "S"] / M_) if M_ else None
        vals.update(ratio_vals.get(str(smp), {}))
        store_put(str(smp), "icp" if kind.startswith("ICP") else "xrf", vals, "ICP page" if kind.startswith("ICP") else "XRF page")
    sub = df[df.element.isin(chosen)].copy()
    sub["sample"] = sub["sample"].str.replace("ZrO2", "ZrO$_2$")
    if layout.startswith("one"):
        fns = [(lambda ax, el=el: P.draw_grouped_bars(ax, sub[sub.element == el], val, err, ylabel=f"{el} ({ylabel})", legend=False)) for el in chosen]
        fig, axs = C.compose(fns, ncols=min(3, len(fns)), width=max(WIDTH, 2.4 * min(3, len(fns))), panel_aspect=0.95)
    else:
        fig, ax = newfig(); P.draw_grouped_bars(ax, sub, val, err, ylabel=ylabel, colors_by_group=False)
    show(fig, "elemental", piv.reset_index())


def _icp_generic_ratios(df, val, err, sec):
    """v2.13: molar ratios for ANY elements: choose the framework metal(s) as reference, every other element is given
    per reference atom (or per `n` reference atoms, e.g. 6 for an M6 node), σ propagated."""
    from mfs import ratios as RT
    from mfs.chem import AW
    els = [e for e in sorted({str(x).split()[0] for x in df["element"].dropna().unique()}) if e in AW]
    with st.expander("Molar ratios for any elements (choose the reference metal)", expanded=False):
        if len(els) < 2:
            st.info("Needs at least two elements with known atomic weights."); return
        c1, c2 = st.columns([3, 1])
        refs = c1.multiselect("Reference element(s) (the framework metal, e.g. Zr, or Zr + Hf)", els, default=els[:1], key=f"{sec}_gr_ref")
        per = c2.number_input("per n reference atoms", 1, 100, 1, key=f"{sec}_gr_per", help="6 gives 'per M6 node' for a hexanuclear cluster.")
        if not refs:
            return
        d = df.copy(); d[val] = pd.to_numeric(d[val], errors="coerce"); d[err] = pd.to_numeric(d[err], errors="coerce").fillna(0)
        rows = []
        for smp, g in d.dropna(subset=[val]).groupby("sample", sort=False):
            wt, sd = {}, {}
            for el, ge in g.groupby("element"):
                el = str(el).split()[0]
                if el not in AW:
                    continue
                wt[el] = float(ge[val].mean())
                sd[el] = float(max(np.sqrt((ge[err] ** 2).mean()), ge[val].std(ddof=1) if len(ge) > 1 else 0))
            r = RT.molar_ratios(wt, sd, refs, per=float(per))
            row = {"sample": smp}
            for k, (v_, s_) in r.items():
                row[k], row[f"σ {k}"] = v_, s_
            rows.append(row)
        if rows:
            tab = pd.DataFrame(rows)
            st.dataframe(tab.style.format(precision=3), hide_index=True, width="stretch")
            dl(st, "CSV molar ratios ± σ", tab.to_csv(index=False).encode("utf-8-sig"), f"{sec.upper()}_molar_ratios.csv", f"{sec}_gr_csv")


def _icp_ratios(df, val, err, sec):
    """Zr:Hf:S molar ratios with propagated uncertainties; Hf fraction and S per M6 node go to the defect calculator."""
    from mfs import ratios as RT
    st.subheader("Zr : Hf : S molar ratios ± σ")
    st.caption("σ propagated (Gaussian, uncorrelated) from the per-element SD; with replicate rows the larger of the replicate SD "
               "and the pooled instrument SD is used. Hf/(Zr+Hf) and S per M₆ node are stored for the defect calculator.")
    d = df.copy(); d[val] = pd.to_numeric(d[val], errors="coerce"); d[err] = pd.to_numeric(d[err], errors="coerce").fillna(0)
    rows, out = [], {}
    for smp, g in d.dropna(subset=[val]).groupby("sample", sort=False):
        wt, sd = {}, {}
        for el, ge in g.groupby("element"):
            el = str(el).split()[0]
            if el not in RT.M:
                continue
            wt[el] = float(ge[val].mean())
            sd[el] = float(max(np.sqrt((ge[err] ** 2).mean()), ge[val].std(ddof=1) if len(ge) > 1 else 0))
        r = RT.zr_hf_s(wt, sd)
        if not r:
            continue
        row = {"sample": smp}
        for k in ("Zr:Hf", "Hf/(Zr+Hf)", "S/(Zr+Hf)", "S per M6 node"):
            if k in r:
                row[k], row[f"σ {k}"] = r[k]
        if "Zr:Hf:S (per M6)" in r:
            row["Zr:Hf:S (per M6)"] = " : ".join(f"{v:.2f}" for v in r["Zr:Hf:S (per M6)"])
        rows.append(row)
        sv = {"Hf_frac": row.get("Hf/(Zr+Hf)"), "Hf_frac_sd": row.get("σ Hf/(Zr+Hf)"), "S_per_M6": row.get("S per M6 node"),
              "S_per_M6_sd": row.get("σ S per M6 node")}
        if "S/(Zr+Hf)" in row:
            sv.update(xrf_S_per_M=row["S/(Zr+Hf)"], xrf_S_per_M_sd=row["σ S/(Zr+Hf)"])
        out[str(smp)] = {k: (None if v is None or (isinstance(v, float) and not np.isfinite(v)) else float(v)) for k, v in sv.items()}
    if not rows:
        st.info("Needs at least two of Zr, Hf, S."); return out
    tab = pd.DataFrame(rows)
    st.dataframe(tab.style.format(precision=3), hide_index=True, width="stretch")
    dl(st, "CSV ratios ± σ", tab.to_csv(index=False).encode("utf-8-sig"), f"{sec.upper()}_ratios_sigma.csv", f"{sec}_ratio_csv")
    if "S per M6 node" in tab:
        t2 = tab.dropna(subset=["S per M6 node"])

        def draw_r(ax, t2=t2):
            xs = np.arange(len(t2)); cols = S.colours_for(list(t2["sample"]))
            ax.bar(xs, t2["S per M6 node"], 0.65, yerr=t2["σ S per M6 node"], color=cols, edgecolor="k", lw=0.6, error_kw=dict(lw=0.7, capsize=2.5))
            ax.set_xticks(xs, [str(x).replace("ZrO2", "ZrO$_2$") for x in t2["sample"]], rotation=30, ha="right")
            ax.tick_params(axis="x", which="minor", bottom=False, top=False); ax.set_ylabel(r"S per M$_6$ node")
        if len(t2):
            fig, ax = newfig(); draw_r(ax)
            show(fig, f"{sec.upper()}_S_per_node", t2, draw_r, "S per M6 node")
    return out
