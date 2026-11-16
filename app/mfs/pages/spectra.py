"""Reticula page module: PXRD and FTIR / Raman pages.

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
    'REF',
    'files_input',
    'newfig',
    'proc_controls',
    'sample_table',
    'show',
    'store_put',
)


def page_pxrd():
    st.header("Powder X-ray diffraction")
    st.caption("STOE WinXPOW ASCII export (.xy / .xye / .dat), any 2θ–intensity text. Binary STOE .raw only if the optional `xylib` is installed.")
    files = files_input("PXRD files", ["xy", "xye", "dat", "txt", "csv", "raw", "asc", "rasx", "ras"], "PXRD/*.xy", "pxrd")
    if not files:
        return
    samples = sample_table(files, "pxrd_tab")
    opts, off = proc_controls("pxrd", "ALS", "max")
    c1, c2, c3 = st.columns(3)
    lo, hi = c1.slider("2θ range (°)", 2.0, 90.0, (4.0, 50.0), 0.5)
    ref_mode = c2.selectbox("Reference pattern", ["Simulate from bundled MOF-801 CIF (COD 4121459)", "Simulate from bundled UiO-66 CIF (COD 4512074)",
                                                  "Upload CIF", "Upload .xy reference", "none"])
    subst = c3.selectbox("Metal in simulation", ["Hf (replace Zr)", "Zr (as in CIF)"])
    fwhm = c3.number_input("Simulated FWHM (°)", 0.03, 1.0, 0.12, 0.01)
    sticks_on = c2.checkbox("Show reflection sticks", True)
    ref, sticks = None, None
    try:
        if ref_mode.startswith("Simulate") or ref_mode == "Upload CIF":
            if ref_mode == "Upload CIF":
                cf = st.file_uploader("CIF", type=["cif"], key="cif")
                cif = cf.read() if cf else None
            else:
                cif = os.path.join(REF, "MOF-801_Zr_COD4121459.cif" if "801" in ref_mode else "UiO-66_Zr_COD4512074.cif")
            if cif:
                sim, sticks, formula = A.simulate_pxrd(cif, (lo, hi), fwhm=fwhm, substitute={"Zr": "Hf"} if subst.startswith("Hf") else None)
                name = st.text_input("Reference label", "Simulated " + ("MOF-801" if "801" in ref_mode else "UiO-66" if "UiO" in ref_mode else formula) + ("(Hf)" if subst.startswith("Hf") else "(Zr)"))
                ref = (name, sim.x.values, sim.y.values)
        elif ref_mode == "Upload .xy reference":
            rf = st.file_uploader("Reference .xy", key="refxy")
            if rf:
                d = R.read_xy(rf); ref = (st.text_input("Reference label", "Reference"), d.x.values, d.y.values)
    except Exception as e:
        st.warning(f"Reference simulation failed: {e}")
    traces, rows = [], []
    for nm, f in samples:
        try:
            d = R.read_pxrd(f, f.name)
        except Exception as e:
            st.error(f"{f.name}: {e}"); continue
        rg = d.attrs.get("rigaku")
        extra = []
        if rg:                                   # Rigaku SmartLab .rasx / .ras: metadata → store, flag non-2θ/θ scans
            from mfs import rigaku as RG
            allsc = d.attrs.get("rigaku_all") or [rg]
            m_ = rg["meta"]
            st.caption(f"{f.name}: Rigaku SmartLab – {len(allsc)} scan(s); target {m_.get('target', '?')}, Kα1 "
                       f"{m_.get('wavelength_ka1', '?')} Å, axis {rg['axis']}, {m_.get('start', '?')}–{m_.get('end', '?')}°, step "
                       f"{m_.get('step', '?')}°, speed {m_.get('speed', '?')} {m_.get('speed_unit', '')}; intensity = counts × attenuator "
                       f"(max factor {float(np.max(rg['att'])):.4g}).")
            for s_ in allsc:
                if not s_["is_2theta_theta"]:
                    st.warning(f"{f.name} · {s_['name']}: scan axis is **{s_['axis']}**, not 2θ/θ – "
                               + ("not plotted as a powder pattern." if s_ is not rg else "plotted anyway (no 2θ/θ scan in the file), check before use."))
            extra = [s_ for s_ in allsc if s_ is not rg and s_["is_2theta_theta"]]
            store_put(nm, "pxrd", RG.store_values(rg), "PXRD (Rigaku)")
        y = A.process(d.x.values, d.y.values, **opts)
        traces.append((nm, d.x.values, y))
        rows.append(pd.DataFrame({"sample": nm, "two_theta": d.x, "I_raw": d.y, "I_processed": y}))
        for k_, s_ in enumerate(extra, 2):
            ye = A.process(s_["x"], s_["y"], **opts)
            traces.append((f"{nm} · scan {k_}", s_["x"], ye))
            rows.append(pd.DataFrame({"sample": f"{nm} · scan {k_}", "two_theta": s_["x"], "I_raw": s_["y"], "I_processed": ye}))
    st.session_state["_even|pxrdoff"] = 1.05      # traces are min-max scaled to 1 in the window -> even = 1.05
    def draw(ax, traces=traces, ref=ref, sticks=sticks):
        P.draw_stacked(ax, traces, "pxrd", offset=off, xlim=(lo, hi), reference=ref, sticks=sticks if sticks_on else None)
    fig, ax = newfig(); draw(ax)
    show(fig, "PXRD", pd.concat(rows) if rows else None, draw, "PXRD", units={"x": "2theta"})
    from mfs import xrd as X
    with st.expander("📏 Peak width / Scherrer crystallite size (with instrumental correction)"):
        if traces:
            c1, c2, c3, c4 = st.columns(4)
            pks = c1.text_input("Peaks (°2θ, comma separated)", "8.6, 9.9", key="sch_pk")
            Kf = c2.number_input("Shape factor K", 0.5, 1.5, 0.9, 0.01, key="sch_K")
            shape = c3.selectbox("Width correction", ["gaussian", "lorentzian", "voigt"], key="sch_shape",
                                 help="gaussian: β² = B² − b² · lorentzian: β = B − b · voigt: β = B − b²/B")
            st.caption("Defaults: MOF-801 (111) 8.6° and (200) 9.9° (Cu Kα1); UiO-66: 7.4° and 8.5°.")
            inst_src = c4.radio("Instrumental FWHM", ["constant", "from a standard (LaB₆ / Si) pattern"], key="sch_src")
            pk_list = [float(v) for v in pks.replace(";", ",").split(",") if v.strip()]
            inst = lambda tt: 0.05
            if inst_src == "constant":
                b0 = st.number_input("Instrumental FWHM (°)", 0.0, 1.0, 0.05, 0.005, key="sch_b")
                inst = lambda tt, b0=b0: b0
            else:
                sf = st.file_uploader("Standard pattern (.xy)", key="sch_std")
                spk = st.text_input("Standard peaks (°2θ)", "21.36, 30.38, 37.44, 43.51" , key="sch_spk",
                                    help="LaB₆ SRM 660 (Cu Kα1): 21.36 (100), 30.38 (110), 37.44 (111), 43.51 (200) …")
                if sf:
                    try:
                        dstd = R.read_xy(sf)
                        inst, pos, fw = X.instrument_fwhm_from_standard(dstd.x.values, dstd.y.values, [float(v) for v in spk.split(",") if v.strip()])
                        st.caption("Standard fitted: " + ", ".join(f"{p_:.2f}° → {f_:.3f}°" for p_, f_ in zip(pos, fw)) + " (linear FWHM(2θ))")
                    except Exception as e:
                        st.warning(f"Standard fit failed: {e}")
                else:
                    st.caption("No standard loaded – 0.05° used.")
            srows = []
            for nm, x, y in traces:
                for g in pk_list:
                    try:
                        r = X.fit_peak(x, y, g)
                        sz = X.scherrer_size(r["center"], r["fwhm"], inst(r["center"]), K=Kf, shape=shape)
                        srows.append({"sample": nm, "2θ (°)": r["center"], "FWHM obs (°)": r["fwhm"], "± (°)": r["fwhm_err"],
                                      "FWHM inst (°)": inst(r["center"]), "β sample (°)": sz["beta_sample"], "D (nm)": sz["size_nm"], "note": sz["note"]})
                    except Exception as e:
                        srows.append({"sample": nm, "2θ (°)": g, "note": f"fit failed ({e})"})
            st.dataframe(pd.DataFrame(srows).style.format(precision=3), width="stretch", hide_index=True)
            st.caption("pseudo-Voigt fit with linear background; D = Kλ/(β cos θ), λ = 1.540598 Å (Cu Kα1). Sizes > ~100 nm are not reliable.")
    with st.expander("💎 Crystallinity index vs a reference sample"):
        if len(traces) >= 2:
            names = [t[0] for t in traces]
            c1, c2 = st.columns(2)
            rname = c1.selectbox("Reference (100 % crystalline)", names, key="ci_ref")
            win = c2.text_input("Reflection windows (°2θ, lo-hi; …)", "8.2-9.0; 9.5-10.3; 14.6-15.4", key="ci_win")
            ws = [tuple(float(v) for v in w.split("-")) for w in win.split(";") if "-" in w]
            rt = traces[names.index(rname)]
            raw = {r_["sample"].iloc[0]: r_ for r_ in rows}
            use_raw = st.checkbox("use the raw intensities (recommended: same mass / holder / time); untick to use the processed curves", True, key="ci_raw")
            ci = []
            for nm, x, y in traces:
                yy = raw[nm].I_raw.values if use_raw else y; yr = raw[rname].I_raw.values if use_raw else rt[2]
                r_ = X.crystallinity_index(x, yy, rt[1], yr, ws)
                ci.append({"sample": nm, "CI (%)": r_["CI_pct"], "Σ area sample": r_["A_sample"], "Σ area reference": r_["A_ref"]})
            st.dataframe(pd.DataFrame(ci).style.format(precision=1), width="stretch", hide_index=True)
        else:
            st.caption("Load at least two patterns (sample(s) + a reference).")
    with st.expander("↔️ Peak-shift tracking (Δ2θ, d-spacing, cubic lattice parameter)"):
        if traces:
            c1, c2, c3 = st.columns(3)
            tp = c1.text_input("Peaks (°2θ)", "8.6, 9.9", key="sh_pk")
            hk = c2.text_input("hkl (optional, same order)", "111, 200", key="sh_hkl")
            rsel = c3.selectbox("Reference trace", [t[0] for t in traces], key="sh_ref")
            pl = [float(v) for v in tp.split(",") if v.strip()]
            hl = [tuple(int(ch) for ch in h.strip()) if h.strip().isdigit() and len(h.strip()) == 3 else None for h in hk.split(",")]
            tr = X.track_shifts(traces, pl, ref_index=[t[0] for t in traces].index(rsel), hkl=hl)
            tdf = pd.DataFrame(tr)
            st.dataframe(tdf.style.format(precision=4), width="stretch", hide_index=True)

            def draw_sh(ax, tdf=tdf):
                for i_, (g, gg) in enumerate(tdf.groupby("peak (°2θ, guess)", sort=False)):
                    ax.errorbar(range(len(gg)), 1000 * gg["Δ2θ vs ref (°)"], yerr=1000 * gg["± (°)"], fmt="o-", ms=4, lw=0.8, capsize=2, label=f"{g:g}°")
                ax.axhline(0, color="0.6", lw=0.5)
                ax.set_xticks(range(len(traces))); ax.set_xticklabels([t[0] for t in traces], rotation=30, ha="right")
                ax.set_ylabel(r"$\Delta 2\theta$ vs reference (m°)"); ax.legend(frameon=False)
            fig, ax = newfig(); draw_sh(ax)
            show(fig, "PXRD_peak_shifts", tdf, draw_sh, "Peak shifts")


def page_ftir():
    st.header("FTIR (ATR) / FT-Raman")
    st.caption("Bruker OPUS binary (.0, .1, …) via brukeropusreader, or OPUS *Save As → Data Point Table* (.dpt) / CSV.")
    mode = st.radio("Spectrum type", ["ATR-FTIR (transmittance)", "ATR-FTIR (absorbance)", "FT-Raman (1064 nm)"], horizontal=True)
    files = files_input("Spectra", None, "FTIR/*.dpt" if "FTIR" in mode else "Raman/*.dpt", "ftir")
    if not files:
        return
    samples = sample_table(files, "ftir_tab")
    opts, off = proc_controls("ftir", "none" if "FTIR" in mode else "ALS", "min-max")
    c1, c2 = st.columns(2)
    lo, hi = c1.slider("Range (cm⁻¹)", 100.0, 4000.0, (400.0, 2000.0) if "FTIR" in mode else (150.0, 1800.0), 10.0)
    conv = c2.checkbox("Convert transmittance → absorbance", False) if "transmittance" in mode else False
    hl = c2.text_input("Highlight bands (lo-hi:label; …)", "1700-1740:C=O; 1190-1230:S=O; 1030-1055:" if "FTIR" in mode else "")
    highlight = []
    for part in [p for p in hl.split(";") if p.strip()]:
        try:
            rng_, lab = part.split(":") if ":" in part else (part, "")
            a, b = [float(v) for v in rng_.split("-")]; highlight.append((a, b, lab.strip()))
        except Exception:
            pass
    traces, rows = [], []
    for nm, f in samples:
        try:
            d = R.read_ftir(f, f.name)
        except Exception as e:
            st.error(f"{f.name}: {e}"); continue
        y = A.transmittance_to_absorbance(d.y.values) if conv else d.y.values
        y = A.process(d.x.values, y, **opts)
        traces.append((nm, d.x.values, y)); rows.append(pd.DataFrame({"sample": nm, "wavenumber": d.x, "y_processed": y}))
    kind = "raman" if "Raman" in mode else ("ftir_A" if ("absorbance" in mode or conv) else "ftir_T")
    from mfs import spectra as SP
    is_T = kind == "ftir_T"
    with st.expander("🧪 Reference band, peak labels, zoomed inset, difference spectra", expanded=False):
        c1, c2, c3, c4 = st.columns(4)
        nb = c1.checkbox("Normalise to a reference band", False, key="ft_nb")
        dlo, dhi, dlab = (1540.0, 1620.0, "νas(COO⁻)") if kind != "raman" else (1630.0, 1690.0, "ν(C=C)")
        blo = c2.number_input("band from (cm⁻¹)", 50.0, 4000.0, dlo, key="ft_blo")
        bhi = c3.number_input("band to (cm⁻¹)", 50.0, 4000.0, dhi, key="ft_bhi")
        bmode = c4.radio("by", ["height", "area"], horizontal=True, key="ft_bm")
        c1, c2, c3, c4 = st.columns(4)
        pl = c1.checkbox("Label peaks", False, key="ft_pl")
        prom = c2.slider("min. prominence (fraction of range)", 0.01, 0.5, 0.08, 0.01, key="ft_prom")
        npk = c3.number_input("max labels per trace", 1, 40, 8, key="ft_npk")
        npk0 = {t[0]: len(SP.pick_peaks(t[1], t[2], prom, 40, (lo, hi), minima=is_T)) for t in traces}
        pl_which = c4.multiselect("traces", [t[0] for t in traces], default=[max(npk0, key=npk0.get)] if npk0 else [], key="ft_plw")
        c1, c2, c3 = st.columns(3)
        ins = c1.checkbox("Zoomed inset", False, key="ft_ins")
        zr = c2.text_input("inset range (cm⁻¹, lo-hi)", "1000-1250" if kind != "raman" else "200-700", key="ft_zr")
        zloc = c3.selectbox("inset position", ["upper left", "upper center", "upper right", "center left", "lower left"], key="ft_zloc")
        st.caption(f"Reference band suggestion: {dlab}. Normalisation keeps the relative intensities between traces (one common scale).")
    if nb:
        nt, weak = [], []
        for nm, x, y in traces:
            m_ = (x >= lo) & (x <= hi)
            try:
                v_ = SP.band_value(x, -y if is_T else y, blo, bhi, bmode if bmode == "height" else "height")
                if v_ < 0.05 * np.nanmax(np.abs(np.ptp(y[m_]))):
                    raise ValueError("band absent (< 5 % of the trace range)")
                nt.append((nm, x, SP.normalise_to_band(x, y, blo, bhi, bmode, minimum=is_T)))
            except Exception as e:
                weak.append(f"{nm} ({e})")
                r0 = np.ptp(y[m_]) or 1.0          # fallback: own range = 1 so the trace stays readable
                nt.append((nm, x, (y - np.nanmax(y[m_])) / r0 + 1 if is_T else (y - np.nanmin(y[m_])) / r0))
        if weak:
            st.caption("Not band-normalised (scaled to their own range instead): " + "; ".join(weak))
        traces = nt
    peaks_ = {nm: SP.pick_peaks(x, y, prom, int(npk), (lo, hi), minima=is_T) for nm, x, y in traces if pl and nm in pl_which}
    try:
        zlo, zhi = sorted(float(v) for v in zr.split("-")[:2])
    except Exception:
        zlo, zhi = lo, hi

    def draw(ax, traces=traces):
        st.session_state["_even|ftiroff"] = 1.05
        P.draw_stacked(ax, traces, kind, offset=off, xlim=(lo, hi), highlight=highlight, invert=(kind != "raman"),
                       scale="common" if nb else "each")
        info = {e[0]: e for e in getattr(ax, "_mfs_stack", [])}
        for nm, x, y in traces:
            if nm in peaks_ and nm in info:
                _, y0, r_, base, _k = info[nm]
                yl = (np.asarray(y) - y0) / r_ + base
                if is_T:     # transmittance: write the labels just above the trace's top, at the dip positions
                    m_ = (np.asarray(x) >= lo) & (np.asarray(x) <= hi); yl = np.full_like(yl, np.nanmax(yl[m_]) + 0.01)
                SP.label_peaks(ax, x, yl, peaks_[nm])
        if ins and zhi > zlo:
            SP.add_inset(ax, lambda ia: P.draw_stacked(ia, traces, kind, offset=off, xlim=(zlo, zhi), invert=(kind != "raman"),
                                                      scale="common" if nb else "each"), (zhi, zlo) if kind != "raman" else (zlo, zhi), loc=zloc)
    fig, ax = newfig(1.1); draw(ax)
    show(fig, "FTIR" if kind != "raman" else "Raman", pd.concat(rows) if rows else None, draw, "FTIR" if kind != "raman" else "Raman",
         units={"x": "wavenumber"} if kind != "raman" else None)
    if peaks_:
        st.dataframe(pd.DataFrame([{"sample": k, "peaks (cm⁻¹)": ", ".join(f"{v:.0f}" for v in vv)} for k, vv in peaks_.items()]),
                     hide_index=True, width="stretch")
    if len(traces) >= 2:
        with st.expander("➖ Difference spectra (sample − f × reference)", expanded=False):
            names = [t[0] for t in traces]
            c1, c2, c3 = st.columns(3)
            refn = c1.selectbox("reference", names, index=0, key="ft_dref")
            sams = c2.multiselect("samples", [n for n in names if n != refn], default=[n for n in names if n != refn][:2], key="ft_dsam")
            fmode = c3.radio("factor", ["1", "least squares in a window", "manual"], key="ft_dfm")
            fw = c3.text_input("window (cm⁻¹)", "1350-1450", key="ft_dfw") if fmode.startswith("least") else None
            fman = c3.number_input("f", 0.0, 10.0, 1.0, 0.05, key="ft_dfman") if fmode == "manual" else 1.0
            rt = traces[names.index(refn)]
            diffs, drows = [], []
            for sn in sams:
                stt = traces[names.index(sn)]
                f_ = fman
                if fw:
                    try:
                        a_, b_ = sorted(float(v) for v in fw.split("-")[:2]); f_ = SP.auto_factor((rt[1], rt[2]), (stt[1], stt[2]), a_, b_)
                    except Exception as e:
                        st.caption(f"factor: {e}")
                g_, d_ = SP.difference((rt[1], rt[2]), (stt[1], stt[2]), f_)
                diffs.append((f"{sn} − {f_:.2f}×{refn}", g_, d_)); drows.append(pd.DataFrame({"difference": diffs[-1][0], "x": g_, "y": d_}))

            def draw_d(ax, diffs=diffs):
                for k, (nm, g_, d_) in enumerate(diffs):
                    sc = np.nanmax(np.abs(d_[(g_ >= lo) & (g_ <= hi)])) or 1
                    ax.plot(g_, d_ / sc + 2.2 * k, lw=0.9, color=P.colours_for([n for n, _, _ in diffs])[k])
                    ax.axhline(2.2 * k, color="0.7", lw=0.5, ls=":")
                    ax.text(0.02, (2.2 * k + 0.6), nm, transform=ax.get_yaxis_transform(), fontsize=6)
                ax.set_xlim((hi, lo) if kind != "raman" else (lo, hi)); ax.set_yticks([])
                ax.set_xlabel(P.LABELS[kind][0]); ax.set_ylabel("Difference (norm.)")
            if diffs:
                fig, ax = newfig(0.8); draw_d(ax)
                show(fig, ("FTIR" if kind != "raman" else "Raman") + "_difference", pd.concat(drows), draw_d, "Difference spectra")
