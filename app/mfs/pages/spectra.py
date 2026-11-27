"""Reticula page module: PXRD and FTIR / Raman pages.

Moved verbatim out of app.py in v2.13. Names that the code takes from the app (helpers such as show(), dl(), files_input(),
the sidebar values WIDTH / ASPECT / FONT, the cached reader wrapper R …) are listed in INJECTED and filled in by
mfs.pages.bind(globals()) on every run of app.py, so the code behaves exactly as it did inside app.py."""
# ruff: noqa: F821
from __future__ import annotations
import io, os, re
import numpy as np
import pandas as pd
import streamlit as st
from mfs import numentry as NE  # every slider takes a typed number too (2.19.6)
from mfs import levels as LV

from mfs import analysis as A, plots as P, style as S
from mfs import __version__

INJECTED = (
    'R',
    'REF',
    'files_input',
    'newfig',
    'proc_controls',
    'sample_table',
    'show',
    'store_put',
    'dl',
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
    lo, hi = NE.slider_in(c1, "2θ range (°)", 2.0, 90.0, (4.0, 50.0), 0.5)
    ref_mode = c2.selectbox("Reference pattern", ["Simulate from bundled MOF-801 CIF (COD 4121459)", "Simulate from bundled UiO-66 CIF (COD 4512074)",
                                                  "Upload CIF", "Upload .xy reference", "none"], index=4)     # v2.18.4: no reference unless you choose one
    subst = c3.selectbox("Metal in simulation", ["As in the CIF (no substitution)", "Hf (replace Zr)"])           # v2.18.4: Hf is an option, not the default
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
                name = st.text_input("Reference label", "Simulated " + ("MOF-801" if "801" in ref_mode else "UiO-66" if "UiO" in ref_mode else formula) + ("(Hf)" if subst.startswith("Hf") else ""))
                ref = (name, sim.x.values, sim.y.values)
        elif ref_mode == "Upload .xy reference":
            rf = st.file_uploader("Reference .xy", key="refxy")
            if rf:
                d = R.read_xy(rf); ref = (st.text_input("Reference label", "Reference"), d.x.values, d.y.values)
    except Exception as e:
        st.warning(f"Reference simulation failed: {e}")
    traces, rows, raws = [], [], {}
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
        raws[nm] = (d.x.values, d.y.values)
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
    if LV.opt_in("📐 Crystallite size and microstrain (Scherrer · Williamson–Hall · Halder–Wagner)", "pxrd_size",
                 help="Off until you switch it on: fits the peaks of every sample and gives size and strain with uncertainties."):
        _crystal_size_ui(traces, lo, hi)
    if LV.opt_in("🧮 Rietveld / Le Bail refinement (cell, phases, weight fractions)", "pxrd_rv",
                 help="Off until you switch it on: choose the phases (CIF) and press Refine – everything else is automatic."):
        _rietveld_ui(raws, lo, hi, ref_mode, subst)
    if LV.opt_in("🧬 Solid solution: Vegard shift of the reference peaks (cubic end members)", "pxrd_vegard",
                 help="Off until you switch it on: choose two cubic end members and the composition x; the reference peaks move by Vegard's law."):
        _vegard_ui(raws, lo, hi)
    with LV.gate("Developer", globals()):
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
                st.caption("pseudo-Voigt fit with linear background; D = Kλ/(β cos θ), λ = 1.540598 Å (Cu Kα1). Sizes > ~100 nm: treat the size as a lower limit.")
    with LV.gate("Advanced", globals()):
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
    with LV.gate("Advanced", globals()):
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


def _vegard_ui(raws, lo, hi):
    """v2.24.0: solid solution of two end members (Vegard, linear per axis). Each end member comes from one of four sources:
    the library, cell numbers, a pattern of the end member (indexing + Le Bail cell fit), or an own CIF. Both end members must be of
    the same system (cubic, tetragonal or hexagonal). The reference peaks move with x. Source: Vegard (1921), Z. Phys. 5, 17–26,
    doi:10.1007/BF01349680. No bowing in this panel."""
    import math
    from mfs import rietveld as RV
    from mfs import phaseid as PID
    from mfs import vegard as VG
    from mfs import readers as RD
    if not raws:
        st.info("Load a pattern first.")
        return
    lib = [e for e in PID.library() if not e.get("exclude_matching")]
    by_label = {f"{e['name']} ({e['source']})": e for e in lib}
    names = list(by_label)
    box = st.container(border=True, key="vg_box")
    c1, c2, c3 = box.columns([2, 2, 3])
    nm = c1.selectbox("Pattern", list(raws), key="vg_sample")
    wl = RV.RADIATION["Cu Kα1"][0][0]
    dmin = wl / (2.0 * math.sin(math.radians(hi / 2.0)))
    SRC = ["Library", "Cell numbers", "Pattern (indexing + Le Bail cell fit)", "Own CIF"]

    def _end(which, col, default_index):
        src = col.selectbox(f"Source of end member {which}", SRC, index=default_index, key=f"vg_src_{which}")
        if src == "Library":
            if not names:
                col.warning("No library entries.")
                return None, None, None
            # default: the first cubic library entry (the default of cell numbers is cubic too), so the panel opens with a valid pair
            cubic_idx = next((i for i, n in enumerate(names) if VG.system_of(tuple(by_label[n]["cell"])) == "cubic"), 0)
            pick = col.selectbox(f"End member {which}", names, index=cubic_idx, key=f"vg_lib_{which}")
            try:
                ph = RV.Phase.from_cif(by_label[pick]["path"]).build(dmin)
            except Exception as exc:
                col.error(f"Could not read {pick}: {type(exc).__name__}: {exc}")
                return None, None, None
            return tuple(float(v) for v in ph.cell0), [tuple(int(v) for v in h) for h in ph.hkl], pick
        if src == "Cell numbers":
            sysn = col.selectbox(f"System {which}", ["cubic", "tetragonal", "hexagonal"], key=f"vg_sys_{which}")
            a = col.number_input(f"a {which} (Å)", min_value=1.0, max_value=80.0, value=5.4113, step=0.0001, format="%.4f", key=f"vg_a_{which}")
            c = a
            if sysn != "cubic":
                c = col.number_input(f"c {which} (Å)", min_value=1.0, max_value=80.0, value=float(a), step=0.0001, format="%.4f", key=f"vg_c_{which}")
            cell = VG.symmetric_cell(sysn, a, c)
            ph = RV.Phase.from_cell(cell, f"{which} (numbers)").build(dmin)
            return cell, [tuple(int(v) for v in h) for h in ph.hkl], f"{which}: {sysn} a={a:.4f}" + ("" if sysn == "cubic" else f" c={c:.4f}") + " (cell-only)"
        if src.startswith("Pattern"):
            up = col.file_uploader(f"Pattern of end member {which}", type=["xy", "xye", "dat", "txt", "csv", "asc", "ras", "raw"], key=f"vg_pat_{which}")
            if up is None:
                return None, None, None
            try:
                df = RD.read_pxrd(up.getvalue(), up.name)
                xe = df.iloc[:, 0].to_numpy(float); ye = df.iloc[:, 1].to_numpy(float)
            except Exception as exc:
                col.error(f"Could not read the pattern: {type(exc).__name__}: {exc}")
                return None, None, None
            res = VG.cell_from_pattern(xe, ye, lo, hi, wl)
            if res["status"] != "ok":
                col.warning(f"Cell not found for {which}: {res.get('reason', res['status'])}")
                return None, None, None
            cell = res["cell"]
            col.caption(f"{which}: {res['system']} a={cell[0]:.4f} Å" + ("" if res["system"] == "cubic" else f" c={cell[2]:.4f} Å")
                        + f" (Le Bail Rwp {res['rwp']:.2f} %; cell-only reflections)")
            ph = RV.Phase.from_cell(cell, f"{which} (pattern)").build(dmin)
            return cell, [tuple(int(v) for v in h) for h in ph.hkl], f"{which}: from pattern {up.name}"
        up = col.file_uploader(f"CIF of end member {which}", type=["cif"], key=f"vg_up_{which}")
        if up is None:
            return None, None, None
        try:
            ph = RV.Phase.from_cif(up.getvalue(), name=up.name).build(dmin)
        except Exception as exc:
            col.error(f"Could not read the CIF: {type(exc).__name__}: {exc}")
            return None, None, None
        return tuple(float(v) for v in ph.cell0), [tuple(int(v) for v in h) for h in ph.hkl], f"uploaded: {up.name}"

    cell_A, hkl_A, la = _end("A", c2, 0)
    cell_B, hkl_B, lb = _end("B", c2, 1)
    x = c3.slider("x (fraction of A)", 0.0, 1.0, 0.5, 0.01, key="vg_x")
    if cell_A is None or cell_B is None:
        box.info("Choose two end members.")
        return
    sA, sB = VG.system_of(cell_A), VG.system_of(cell_B)
    if not sA or sA != sB:
        box.warning(f"The end members must be of the same system (cubic, tetragonal or hexagonal): A is {sA or 'other'}, B is {sB or 'other'}.")
        return
    hkl = hkl_A
    cell_x = VG.cell_from_composition(x, cell_A, cell_B)
    th_x = [VG.two_theta_for_cell(cell_x, h, wl) for h in hkl]
    th_A = [VG.two_theta_for_cell(cell_A, h, wl) for h in hkl]
    th_B = [VG.two_theta_for_cell(cell_B, h, wl) for h in hkl]
    xo, yo = raws[nm]
    xo, yo = np.asarray(xo, float), np.asarray(yo, float)
    m = (xo >= lo) & (xo <= hi)
    box.caption(f"System: {sA}. A: {la}. B: {lb}. a(x = {x:.2f}) = {cell_x[0]:.4f} Å. Dashed: end members. Solid: predicted at x. Cu Kα1 assumed.")

    def draw(ax, xo=xo, yo=yo, th_x=th_x, th_A=th_A, th_B=th_B):
        ax.plot(xo[m], yo[m], lw=0.8, color="0.2", label="measured")
        for series, style in ((th_A, dict(color="C0", ls="--", lw=0.7)), (th_B, dict(color="C1", ls="--", lw=0.7)),
                              (th_x, dict(color="C3", lw=0.9))):
            for t in series:
                if math.isfinite(t) and lo <= t <= hi:
                    ax.axvline(t, **style)
        ax.set_xlim(lo, hi)
        ax.set_xlabel(r"2$\theta$ (°)"); ax.set_ylabel("intensity")

    fig, ax = newfig(aspect=0.75)
    draw(ax)
    show(fig, f"PXRD_Vegard_{nm}", pd.DataFrame({"two_theta": xo[m], "observed": yo[m]}), draw, "Solid solution (Vegard shift)",
         units={"x": "2theta"})
    with box.expander("Estimate x from observed peaks (several reflections; a spread is shown)", expanded=False):
        labels_r = [f"({h[0]} {h[1]} {h[2]})" for h in hkl]
        picked = st.multiselect("Reflections (of end member A)", labels_r, default=labels_r[:min(3, len(labels_r))], key="vg_refl")
        peaks = []
        for lab in picked:
            i = labels_r.index(lab)
            obs = st.number_input(f"Observed 2θ of {lab} (°)", min_value=0.0, max_value=180.0, value=float(lo), step=0.01,
                                  format="%.3f", key=f"vg_obs_{lab}")
            peaks.append((float(obs), hkl[i]))
        if st.button("Estimate x", key="vg_est"):
            r = VG.estimate_x_grid(peaks, wl, cell_A, cell_B)
            if r["status"] in ("ok", "end member"):
                st.success(f"x(A) ≈ {r['x']:.3f} (rms residual {r['rms_deg']:.4f}°, {r['n']} reflections, status: {r['status']})")
            else:
                st.warning(f"Not estimated: {r.get('reason', r['status'])}" + (f" (rms {r['rms_deg']:.3f}°)" if "rms_deg" in r else ""))


def _crystal_size_ui(traces, lo, hi):
    """v2.15: crystallite size + microstrain of every sample. Opt-in (the caller shows it only after a click)."""
    from mfs import crystsize as CS, xrd as X
    if not traces:
        st.info("Load a pattern first.")
        return
    box = st.container(border=True, key="cs_box")
    c1, c2 = box.columns([3, 2])
    inst_mode = c1.radio("Instrument broadening", ["Typical laboratory diffractometer (0.05°)", "I measured a standard (LaB₆ / Si)", "I know the value"],
                         key="cs_inst", help="The peak width of the instrument itself is subtracted. A standard measured on YOUR instrument is the correct way.")
    wl_name = c2.selectbox("Radiation", ["Cu Kα1 (1.5406 Å)", "Co Kα1 (1.7889 Å)", "Mo Kα1 (0.7093 Å)", "Cr Kα1 (2.2897 Å)", "Other"], key="cs_wl")
    wl = {"Cu": 1.5406, "Co": 1.7889, "Mo": 0.7093, "Cr": 2.2897}.get(wl_name[:2])
    if wl is None:
        wl = c2.number_input("Wavelength (Å)", 0.2, 5.0, 1.5406, 0.0001, format="%.4f", key="cs_wl_o")
    inst = None
    if inst_mode.startswith("I measured"):
        sf = box.file_uploader("Standard pattern (.xy, e.g. LaB₆ SRM 660)", key="cs_std")
        if sf:
            try:
                dstd = R.read_xy(sf)
                gp = CS.auto_peaks(dstd.x.values, dstd.y.values, n_max=8)
                inst, pos, fw = X.instrument_fwhm_from_standard(dstd.x.values, dstd.y.values, gp)
                box.caption("Standard: " + ", ".join(f"{p_:.2f}° → {f_:.3f}°" for p_, f_ in zip(pos, fw)))
            except Exception as e:
                box.warning(f"Standard could not be fitted: {e}")
        else:
            box.caption("Upload the pattern of the standard – until then 0.05° is used.")
    elif inst_mode.startswith("I know"):
        b0 = box.number_input("Instrumental FWHM (°2θ)", 0.0, 1.0, 0.05, 0.005, format="%.3f", key="cs_b0")
        inst = lambda tt, b0=b0: b0
    K, shape, manual = 0.9, "voigt", None
    if LV.at_least("Advanced"):
        a1, a2, a3 = box.columns(3)
        K = a1.number_input("Shape factor K", 0.5, 1.5, 0.9, 0.01, key="cs_K", help="0.9 for spherical crystallites; 0.94 for cubic.")
        shape = a2.selectbox("Width correction", ["voigt", "gaussian", "lorentzian"], key="cs_shape",
                             help="voigt (default): β = √((βo−βi)·√(βo²−βi²)) · gaussian: β² = βo²−βi² · lorentzian: β = βo−βi")
        txt = a3.text_input("Peaks to use (°2θ; empty = automatic)", "", key="cs_peaks")
        try:
            manual = [float(v) for v in txt.replace(";", ",").split(",") if v.strip()] or None
        except ValueError:
            box.warning("Peak list not understood – automatic detection is used.")
    results, rows, notes = {}, [], []
    for nm, x, y in traces:
        r = CS.analyse(x, y, wl, inst, K, shape, manual, lo=lo, hi=hi)
        results[nm] = r
        sc, wh, hw = r["scherrer"], r["wh"], r["hw"]
        rows.append({"sample": nm, "peaks": sc["n"] if sc else 0,
                     "D Scherrer (nm)": sc["D_mean"] if sc else np.nan, "± SD": sc["D_sd"] if sc else np.nan,
                     "D W–H (nm)": wh["D_nm"] if wh else np.nan, "± err": wh["D_err_nm"] if wh else np.nan,
                     "strain ε (%)": wh["strain"] * 100 if wh else np.nan, "± ε err (%)": wh["strain_err"] * 100 if wh else np.nan,
                     "D H–W (nm)": hw["D_nm"] if hw else np.nan, "R² (W–H)": wh["r2"] if wh else np.nan})
        notes += [f"{nm}: {n}" for n in r["notes"]]
        if sc:
            store_put(nm, "pxrd", {"crystallite_size_scherrer_nm": sc["D_mean"], "crystallite_size_scherrer_sd_nm": sc["D_sd"],
                                   "crystallite_size_WH_nm": wh["D_nm"] if wh else None, "microstrain_WH_percent": wh["strain"] * 100 if wh else None},
                      "PXRD crystallite size")
    tab = pd.DataFrame(rows)
    box.dataframe(tab.style.format(precision=2, na_rep="–"), hide_index=True, width="stretch")
    for n in notes:
        box.caption("⚠️ " + n)
    if not inst:
        box.caption("ℹ️ The instrumental width was assumed (0.05°) – measure a standard for a result you can defend in a paper.")
    if any(r["wh"] for r in results.values()):
        def draw(ax, results=results):
            CS.draw_wh(ax, results, {nm: S.SERIES_COLOURS.get(nm) if hasattr(S, "SERIES_COLOURS") else None for nm in results})
        fig, ax = newfig(); draw(ax)
        show(fig, "PXRD_Williamson-Hall", tab, draw, "Williamson–Hall plot")
    if LV.at_least("Advanced"):
        with st.expander("🔬 Fitted peaks (position, width, size of every peak)"):
            pk = pd.concat([pd.DataFrame(r["peaks"]).assign(sample=nm) for nm, r in results.items() if r["peaks"]]) if results else pd.DataFrame()
            if len(pk):
                st.dataframe(pk[[c for c in ["sample", "two_theta", "fwhm_obs", "fwhm_err", "fwhm_inst", "beta", "D_scherrer_nm", "D_err_nm", "ok", "note"] if c in pk]]
                             .style.format(precision=3, na_rep="–"), hide_index=True, width="stretch")
    first = next((r for r in results.values() if r["scherrer"]), None)
    if first:
        st.caption("📝 Methods: Crystallite sizes were estimated from the broadening of the reflections (pseudo-Voigt fits, instrumental broadening "
                   f"{'from a measured standard' if inst_mode.startswith('I measured') and inst else 'assumed / given'} removed) with the Scherrer equation (K = {K:g}) and by "
                   "Williamson–Hall analysis of ≥ 3 reflections; uncertainties are from the weighted regression.")


def _add_cell_phase(sol):
    """Callback: the indexed cell becomes an 'unknown phase' for Le Bail."""
    ss = st.session_state
    lst = ss.setdefault("rv_cellph", [])
    lst.append((f"unknown {len(lst) + 1}", tuple(sol["cell"])))


def _use_ident(labels):
    """Callback of 'Use these phases': the identified phases become the selection of the refinement."""
    st.session_state["rv_phases"] = list(labels)


def _rv_progress(bar):
    """v2.17: progress bar of a refinement with the elapsed time and an estimate of the time left (from the fraction done)."""
    import time as _tm
    t0 = _tm.time()

    def cb(f, text):
        f = min(max(float(f), 0.0), 1.0)
        el = _tm.time() - t0
        left = el * (1 - f) / f if f > 0.05 else None
        tail = f" · {el:.0f} s elapsed" + (f" · about {left:.0f} s left" if left is not None and f < 0.99 else "")
        bar.progress(f, text=text + tail)
    return cb


def _std_profile(up):
    """v2.18.2: the instrument profile fitted on an uploaded Si standard pattern (cached per file content).
    Returns ({U, V, W, X, Y, fixed: True} for Refinement, the full result of the standard fit)."""
    import hashlib
    from mfs import rietveld as RV
    from mfs import readers as R
    raw = up.getvalue()
    key = hashlib.md5(raw).hexdigest()
    cache = st.session_state.setdefault("rv_std_cache", {})
    if key not in cache:
        df = R.read_pxrd(io.BytesIO(raw), up.name)
        cache[key] = RV.instrument_profile_from_standard(df["x"].values, df["y"].values)
    s_ = cache[key]
    return {k: s_[k] for k in ("U", "V", "W", "X", "Y")} | {"fixed": True}, s_


def _solve_ui(box, raws, nm):
    """v2.18.4: guided structure solution of an unknown compound from this powder pattern (mfs.solve): cell, space group, charge flipping, atoms, refinement."""
    import os as _os
    from mfs import solve as SV
    from mfs import rietveld as RV
    from mfs import refout as RO
    ss = st.session_state
    with box.expander("🧪 Solve the crystal structure of an unknown compound (from this powder pattern)", expanded=False):
        st.caption("The route of the field, automated: **data check → unit cell (indexing + Le Bail on your pattern) → space group (systematic absences) → structure solution "
                   "(charge flipping, many random starts in parallel) → atoms (assigned from your formula by the R factor) → Rietveld refinement**. Every decision is shown and can be "
                   "overridden. It works best with a high-quality pattern to 2θ ≈ 80° (Cu Kα) or more; the result is a **starting model** to continue with, not a finished structure.")
        c1, c2 = st.columns([3, 1])
        formula = c1.text_input("Chemical formula (one formula unit), e.g. Zr6O4(OH)4(C8H4O4)6, or the elements, e.g. Sn O", "", key="sv_formula",
                                help="The atoms are named from this: the electron-density peaks cannot tell the elements apart by themselves. Hydrogen is ignored (it is invisible to X-rays).")
        zin = c2.text_input("Formula units per cell (Z)", "", key="sv_z", help="Empty = estimated from the cell volume; the best fit of the intensities decides between the plausible values.")
        c3, c4, c5 = st.columns([2, 1, 1])
        si_up = c3.file_uploader("Silicon (or other standard) pattern measured on the same instrument – optional", key="sv_si",
                                 help="Gives the zero shift of the instrument and the peak profile. A zero error of 0.05° spoils indexing; with a reference the cell is more reliable.")
        sgn = c4.text_input("Space group No. (empty = determine)", "", key="sv_sg")
        mono = c5.checkbox("Also try monoclinic cells (slower)", False, key="sv_mono")
        n_tr, n_it, n_thr = 12, 300, max(1, (_os.cpu_count() or 2))
        if LV.at_least("Advanced"):
            a1, a2, a3 = st.columns(3)
            n_tr = NE.slider_in(a1, "Random starts per space group", 4, 48, 12, key="sv_trials")
            n_it = NE.slider_in(a2, "Charge-flipping iterations", 100, 1500, 300, 50, key="sv_iter")
            n_thr = NE.slider_in(a3, "Parallel threads", 1, max(2, (_os.cpu_count() or 2)), max(1, (_os.cpu_count() or 2)), key="sv_threads")
        go = st.button("▶ Solve the structure", key="sv_go", type="primary", width="stretch", help="Typically one to several minutes.")
        store = ss.setdefault("sv_res", {})
        if go:
            try:
                comp = SV.parse_formula(formula)
            except Exception as e:
                comp = {}
                st.error(f"The formula was not understood: {e}")
            if not comp:
                st.warning("Give the chemical formula or the elements so that the peaks can be named.")
            else:
                x, y = raws[nm]
                rad_ = ss.get("rv_rad", "Cu Kα1")
                inst, zero = None, 0.0
                if si_up is not None:
                    try:
                        from mfs import readers as R2
                        d_si = R2.read_pxrd(io.BytesIO(si_up.getvalue()), si_up.name)
                        std2 = RV.instrument_profile_from_standard(d_si["x"].values, d_si["y"].values, rad_, offset="zero")        # separate standard sample: instrument zero, not its displacement
                        inst = {k_: std2[k_] for k_ in ("U", "V", "W", "X", "Y")} | {"fixed": True}
                        zero = float(std2["zero"])
                        st.caption(f"Reference pattern: zero shift {zero:+.4f}° applied, profile U, V, W, X, Y fixed; its Rwp is {100 * std2['rwp']:.1f} % and a = {std2['cell_a']:.4f} Å (SRM 640: 5.4312 Å).")
                    except Exception as e:
                        st.warning(f"The reference pattern could not be used: {e}")
                bar = st.progress(0.0, text="Starting …")

                def cb(f, t=""):
                    bar.progress(min(1.0, max(0.0, float(f))), text=str(t)[:110])
                systems = ("cubic", "tetragonal", "hexagonal", "orthorhombic") + (("monoclinic",) if mono else ())
                try:
                    out = SV.solve(x, y, rad_, comp, int(zin) if zin.strip().isdigit() else None, None, int(sgn) if sgn.strip().isdigit() else None, inst, zero, systems,
                                   None, n_tr, n_it, 2, 6, 0.10, n_thr, cb, True)
                    out["_formula"] = formula
                    out["_zero"] = zero
                    store[nm] = out
                except Exception as e:
                    st.error(f"The structure solution stopped: {e}")
                bar.empty()
        out = store.get(nm)
        if out is None:
            st.caption("Give the formula and press **Solve the structure**. The decision log appears here.")
            return
        icon = {"ok": "✅", "warning": "⚠️", "failed": "❌"}
        for s_ in out["steps"]:
            with st.expander(f"{icon.get(s_['status'], 'ℹ️')} {s_['title']}", expanded=s_["key"] in ("cf", "atoms", "rietveld") or s_["status"] != "ok"):
                st.write(s_["text"])
                tb = s_.get("data", {}).get("table")
                if tb:
                    st.dataframe(pd.DataFrame(tb), hide_index=True, width="stretch")
        if not out.get("ok"):
            st.warning("No structure was obtained. Read the steps above: the first step that is marked ⚠️ or ❌ tells what to change (more peaks / higher angle range, a zero-shift reference, the formula, "
                       "the space group, or more random starts in the Advanced settings).")
            return
        res_r = out.get("refinement")
        m = st.columns(4)
        m[0].metric("Space group", out["space_group"]["label"], f"No. {out['space_group']['number']}", delta_color="off")
        m[1].metric("Cell volume", f"{RV.cell_volume(tuple(out['cell'])):.1f} Å³", f"{out['system']}", delta_color="off")
        m[2].metric("R(F) of the amplitudes", f"{100 * out['R_F']:.1f} %", help="Agreement of the model with the Le Bail intensities (a model of the heavy atoms from powder data: about 10–30 %).")
        if res_r:
            m[3].metric("Rwp after refinement", f"{100 * res_r['stats']['Rwp']:.1f} %", "good" if res_r["stats"]["Rwp"] < 0.15 else "check", delta_color="off" if res_r["stats"]["Rwp"] < 0.15 else "inverse")
        for w_ in out.get("warnings", [])[-2:]:
            st.caption("⚠️ " + w_)
        ph_o = (res_r or {}).get("_phase_objs") or []
        if res_r:
            fig, ax = newfig(aspect=0.75)
            RV.draw(ax, res_r)
            show(fig, f"PXRD_solved_{nm}", pd.DataFrame({"two_theta": res_r["x"], "observed": res_r["y"], "calculated": res_r["calc"], "background": res_r["bg"], "difference": res_r["resid"]}),
                 lambda ax_, r_=res_r: RV.draw(ax_, r_), "Rietveld fit of the solved structure", units={"x": "2theta"})
            for _i, _p in enumerate(res_r["phases"]):
                if _p.get("structure"):
                    from mfs import struct_ui as SU
                    SU.render(res_r, _i, nm, {"show": show, "newfig": newfig, "dl": dl}, __version__)
                else:
                    RO.render(res_r, _i, ph_o[_i] if _i < len(ph_o) else None, nm, {"dl": dl}, __version__)
        _fn = re.sub(r"[^A-Za-z0-9_.-]+", "_", nm)
        dl(st, "⬇️ Solved structure (CIF, starting model)", out["cif"].encode("utf-8"), f"{_fn}_solved.cif", "sv_cif",
           help="Symmetry written by the program; light atoms may be missing and hydrogen is not placed – refine with restraints before use.")
        extra = {"solution_report.txt": SV.report_text(out, nm, __version__), "solved_structure_before_refinement.cif": out["cif"]}
        if res_r:
            dl(st, "⬇️ Everything to continue with (ZIP: CIFs, pattern files, reflection lists, decision report)", lambda: RO.bundle(res_r, ph_o, nm, __version__, extra_files=extra),
               f"{_fn}_structure_solution.zip", "sv_zip")
        else:
            dl(st, "⬇️ Decision report (text)", SV.report_text(out, nm, __version__).encode("utf-8"), f"{_fn}_solution_report.txt", "sv_rep")


def _rietveld_ui(raws, lo, hi, ref_mode="", subst=""):
    """v2.15: Rietveld / Le Bail refinement - pick phases, press Refine. Opt-in (shown only after a click)."""
    import math
    from mfs import rietveld as RV
    if not raws:
        st.info("Load a pattern first.")
        return
    box = st.container(border=True, key="rv_box")
    ss = st.session_state
    c1, c2, c3 = box.columns([2, 3, 2])
    nm = c1.selectbox("Pattern", list(raws), key="rv_sample")
    from mfs import phaseid as PID
    _lib = PID.library()
    bundled = {f"{e['name']} ({e['source']})": e["path"] for e in _lib}
    dflt = [k for k in bundled if ("801" in k and "801" in ref_mode) or ("UiO" in k and "UiO" in ref_mode)]
    pick = c2.multiselect("Phases in the sample", list(bundled), default=dflt, key="rv_phases",
                          help="Choose every crystalline phase that is present. A missing phase shows up as extra peaks in the difference curve.")
    with box.expander("🔎 Don't know the phases? Identify them automatically", expanded=False):
        st.caption(f"Compares the peaks of the pattern with {len(_lib)} reference structures (MOFs, oxides, standards, minerals; COD) – with a small search over the cell size and "
                   "the zero shift – and proposes the phases. The refinement afterwards decides.")
        if st.button("🔎 Identify the phases", key="rv_ident"):
            _x, _y = raws[nm]
            with st.spinner("Searching the library …"):
                ss["rv_ident_res"] = (nm, PID.match(_x, _y, RV.RADIATION[ss.get("rv_rad", "Cu Kα1")][0][0], progress=None))
        _ir = ss.get("rv_ident_res")
        if _ir and _ir[0] == nm:
            if not _ir[1]:
                st.warning("No phase of the library matches this pattern well – the sample is probably something else (upload its CIF below), or a mixture with many unknown peaks.")
            else:
                st.dataframe(pd.DataFrame([{"phase": r_["name"], "confidence": r_["confidence"], "score": round(r_["score"], 2),
                                             "reflections found": f"{r_['n_matched']}/{r_['n_expected']}", "cell vs library (%)": round(r_["cell_change_pct"], 2),
                                             "pattern explained (%)": round(100 * r_.get("explained", 0), 0)} for r_ in _ir[1]]), hide_index=True, width="stretch")
                st.button("➕ Use these phases", key="rv_use_ident", type="primary", on_click=_use_ident, args=([f"{r_['name']} ({r_['source']})" for r_ in _ir[1]],))
                st.caption("The cell change shows how far the measured cell is from the library cell (solid solution, another metal, temperature …).")
    with box.expander("🧭 Unknown phase? Find its unit cell (indexing)", expanded=False):
        from mfs import indexing as IX
        st.caption("For a phase that is not in the library and has no CIF: the unit cell is searched from the peak positions (cubic, tetragonal, hexagonal, orthorhombic; monoclinic on request). "
                   "The quality is the de Wolff M20 (> 10 reliable); the zero shift of the pattern is refined together with the cell. "
                   "The cell is then fitted with Le Bail (no structure needed) – a good cell gives a good fit.")
        _mono = st.checkbox("Also try monoclinic cells (slower)", False, key="rv_index_mono",
                            help="Monoclinic cells have four free parameters and many more candidates; the search takes longer and is only reliable with a good pattern.")
        if st.button("🧭 Index the pattern", key="rv_index"):
            _x, _y = raws[nm]
            _lam = RV.RADIATION[ss.get("rv_rad", "Cu Kα1")][0][0]
            _pk = PID.exp_peaks(_x, _y, 0.015)[0]
            _pb = st.progress(0.0, text="Searching unit cells …")

            def _cb(f, t=""):
                _pb.progress(min(1.0, max(0.0, float(f))), text=f"Searching unit cells – {t}" if t else "Searching unit cells …")
            _syst = ("cubic", "tetragonal", "hexagonal", "orthorhombic") + (("monoclinic",) if _mono else ())
            ss["rv_index_res"] = (nm, IX.index(_pk, _lam, systems=_syst, progress=_cb), len(_pk))
            _pb.empty()
        _xr = ss.get("rv_index_res")
        if _xr and _xr[0] == nm:
            if not _xr[1]:
                st.warning(f"No cell found from {_xr[2]} peaks (triclinic cells, impurities or too few peaks; tick 'Also try monoclinic cells' if you have not). Dedicated programs (DICVOL, McMaille) go further.")
            else:
                st.dataframe(pd.DataFrame([{"#": k + 1, "system": s_["system"], "a (Å)": s_["cell"][0], "b (Å)": s_["cell"][1], "c (Å)": s_["cell"][2], "γ (°)": s_["cell"][5],
                                             "volume (Å³)": s_["volume"], "M20": s_["M20"], "unindexed": s_["n_unindexed"], "quality": s_["quality"]} for k, s_ in enumerate(_xr[1])]),
                             hide_index=True, width="stretch")
                _sel = st.selectbox("Cell to use", range(len(_xr[1])), format_func=lambda k: f"#{k + 1}: {_xr[1][k]['system']} · V = {_xr[1][k]['volume']} Å³ · M20 = {_xr[1][k]['M20']}", key="rv_index_pick")
                st.button("➕ Add as unknown phase (Le Bail)", key="rv_index_add", on_click=_add_cell_phase, args=(_xr[1][_sel],))
        if ss.get("rv_cellph"):
            st.markdown("**Unknown phases (cell only):** " + ", ".join(f"{n_} ({c_[0]:.3f}, {c_[1]:.3f}, {c_[2]:.3f} Å)" for n_, c_ in ss["rv_cellph"]))
            st.button("Remove the unknown phases", key="rv_cell_clear", on_click=lambda: ss.pop("rv_cellph", None))
    _solve_ui(box, raws, nm)
    _lk1, _lk2 = box.columns([1, 2], vertical_alignment="center")
    _lk1.link_button("🔗 Open ICSD (licensed)", "https://icsd.fiz-karlsruhe.de/", width="stretch",
                     help="Opens ICSD in your browser. Search there with your own login, then download the CIF and upload it below.")
    _lk2.caption("Reticula never searches ICSD for you: each query is paid, so search there yourself and upload the CIF here "
                 "(also CSD, COD or any CIF file).")
    ups = box.file_uploader("…and / or your own CIF file(s)", type=["cif"], accept_multiple_files=True, key="rv_cifs")
    rad = c3.selectbox("Radiation", list(RV.RADIATION), key="rv_rad", help="Kα1 only = monochromated beam (e.g. STOE); Kα1 + Kα2 = Kβ filter (e.g. Rigaku, Bruker).")
    d1, d2, d3 = box.columns([3, 2, 2])
    rng = NE.slider_in(d1, "2θ range used (°)", float(lo - 2 if lo > 4 else 2.0), 120.0, (float(lo), float(max(hi, lo + 10))), 0.5, key="rv_rng")
    mode_l = d2.radio("Method", ["Rietveld (structure known)", "Rietveld + refine the atoms (structure)", "Le Bail (cell + peak shape only)"], key="rv_mode",
                      help="Rietveld uses the intensities of the CIF. 'Refine the atoms' also refines the atom positions (within the site symmetry) and one displacement "
                           "parameter per element, with soft bond restraints, and gives bond lengths, angles, a 3D view and a CIF. Le Bail lets every intensity free – "
                           "use it when the structure is not right.")
    hf = d3.checkbox("Replace Zr by Hf", subst.startswith("Hf"), key="rv_hf", help="For Hf analogues of the bundled Zr structures.")
    bgo, mono, offm, rscale = 6, 0.0, "sample displacement", 1.0
    from mfs import readers as R_
    glass_on, bg_pat = True, None                     # v2.18.2: glass-holder background (automatic) and an optional measured background
    inst_fix = None                                  # v2.18.2: {U..Y, fixed} from a measured standard (None = refined from the pattern)
    if LV.at_least("Advanced"):
        e1, e2 = box.columns(2)
        rscale = NE.slider_in(box, "Bond restraints (structure refinement)", 0.3, 3.0, 1.0, 0.1, key="rv_rscale",
                            help="< 1 tighter: keep the bond lengths of the starting structure; > 1 looser: let the data decide. 1 = default (0.06 → 0.025 Å).")
        bgo = NE.slider_in(e1, "Background polynomial order", 2, 14, 6, key="rv_bg")
        mono = e2.number_input("Monochromator 2θ (° – 0 = none, 26.6 = graphite)", 0.0, 60.0, 0.0, 0.1, key="rv_mono")
        offm = box.radio("Instrument offset", ["sample displacement (reflection geometry – default)", "zero shift (transmission / capillary / calibrated)", "both (only with a measured standard)"],
                         key="rv_off", horizontal=True, help="Zero shift and displacement are strongly correlated; refining both biases the cell (≈ 0.05 % on a real pattern).")
        glass_on = box.checkbox("Model the amorphous background of a glass sample holder (automatic)", True, key="rv_glass",
                                help="Adds broad humps (15–35° 2θ) to the background. Switch off for a zero-background holder.")
        bgf = box.file_uploader("Measured background (empty holder, .xy / .dat / .txt) – optional", key="rv_bgfile",
                                help="Pattern of the empty holder on the same instrument. It is scaled and subtracted in the refinement.")
        if bgf is not None:
            try:
                _d = R_.read_pxrd(io.BytesIO(bgf.getvalue()), bgf.name)
                bg_pat = (_d["x"].values, _d["y"].values)
            except Exception as e:
                box.error(f"Could not read the background file: {e}")
        prof_m = box.radio("Instrument profile", ["refined from this pattern (default)", "fixed from a standard measured on this instrument (Si SRM 640)"],
                           key="rv_prof", horizontal=True,
                           help="The peak shape (U, V, W, X, Y) is refined from the pattern by default. A standard (Si SRM 640 or LaB₆) measured on the same instrument "
                                "gives the instrument profile; it is then kept fixed, so the sample's peak widths are interpreted with the instrument's own contribution.")
        if prof_m.startswith("fixed"):
            std_up = box.file_uploader("Si standard pattern (.xy / .dat / .txt, same instrument and radiation)", key="rv_std")
            if std_up is not None:
                try:
                    inst_fix, inst_std = _std_profile(std_up)
                    box.caption(f"Instrument profile from the standard: U = {inst_std['U']:.4f}, V = {inst_std['V']:.4f}, W = {inst_std['W']:.4f}, "
                                f"X = {inst_std['X']:.4f}, Y = {inst_std['Y']:.4f} (°²/°). Standard fit: Rwp = {100 * inst_std['rwp']:.1f} %, cell a = {inst_std['cell_a']:.4f} Å (SRM 640: 5.4312 Å).")
                except Exception as e:
                    box.error(f"Could not fit the standard: {e}")
            else:
                box.caption("Upload the standard pattern to fix the profile; until then the profile is refined from the pattern.")
    phases_src = [(k.rsplit(" (", 1)[0], bundled[k]) for k in pick]
    for u in ups or []:
        phases_src.append((os.path.splitext(u.name)[0], u.getvalue()))
    po_txt = {}
    if LV.at_least("Advanced") and phases_src:
        with box.expander("📐 Preferred orientation (plate / needle crystals that lie in one direction)", expanded=False):
            st.caption("Give the axis (h k l) per phase, e.g. `0 0 1` for plates that lie flat; empty = random powder. The factor r is refined "
                       "(r < 1 plates, r > 1 needles) with the March–Dollase function.")
            for nm_, _src in phases_src:
                po_txt[nm_] = st.text_input(f"Axis of {nm_}", "", key=f"rv_po|{nm_}")

    def _axis(txt):
        try:
            v = [int(float(x)) for x in txt.replace(",", " ").split()]
            return tuple(v) if len(v) == 3 and any(v) else None
        except ValueError:
            return None
    store = ss.setdefault("rv_res", {})
    run = box.button("▶ Refine", key="rv_go", type="primary", width="stretch", disabled=not (phases_src or ss.get("rv_cellph")),
                     help="Runs the complete automatic strategy: start values, cell, peak shape, background, displacement." if phases_src else "Choose at least one phase.")
    auto = box.button("✨ Do everything automatically: identify the phases, refine the structure", key="rv_auto", width="stretch",
                      help="Searches the library for the phases (cell and zero shift included), then refines cell, peak shape, background, atom positions and displacement "
                           "parameters of the proposed phases. You get the figure, the statistics, the 3D structure and the CIF.")
    if run or auto:
        x, y = raws[nm]
        bar = box.progress(0.0, text="Starting …")
        try:
            if auto:
                bar.progress(0.02, text="Identifying the phases …")
                _m = [r_ for r_ in PID.match(x, y, RV.RADIATION[rad][0][0]) if r_["confidence"] != "possible"][:3] or PID.match(x, y, RV.RADIATION[rad][0][0])[:1]
                if not _m:
                    bar.empty()
                    box.warning("No phase of the library matches this pattern. Use the indexing (unknown phase) or upload a CIF.")
                    raise RuntimeError("no phase identified")
                phases_src = [(r_["name"], r_["file"]) for r_ in _m]
                mode_l = "Rietveld + refine the atoms (structure)"
                box.info("Phases found: " + ", ".join(f"{r_['name']} ({r_['confidence']}, cell {r_['cell_change_pct']:+.2f} %)" for r_ in _m))
            sub = {"Zr": "Hf"} if hf else None
            ph = [RV.Phase.from_cif(src, name, substitute=sub) for name, src in phases_src]
            ph += [RV.Phase.from_cell(c_, n_) for n_, c_ in ss.get("rv_cellph", [])]
            if ss.get("rv_cellph") and not mode_l.startswith("Le Bail"):
                box.info("A phase with only a cell (no structure) can only be fitted with Le Bail – the method was switched for this run.")
                mode_l = "Le Bail (cell + peak shape only)"
            for ph_, (name, _s) in zip(ph, phases_src):
                ph_.po_axis = _axis(po_txt.get(name, ""))
            res = RV.quick_refine(x, y, ph, rad, "lebail" if mode_l.startswith("Le Bail") else "rietveld", tt_range=rng, bg_order=bgo,
                                  mono_2theta=mono, offset={"s": "displacement", "z": "zero", "b": "both"}[offm[0]], progress=_rv_progress(bar),
                                  refine_structure=mode_l.startswith("Rietveld + refine"), restraint_scale=rscale, instrument=inst_fix, glass_bg=glass_on, bg_pattern=bg_pat)
            res["profile_fixed"] = bool(inst_fix)
            res["_phase_objs"] = ph                        # v2.18.4: for the 3D view, the CIF and the reflection list of any refinement
            store[nm] = res
            bar.empty()
        except Exception as e:
            bar.empty()
            box.error(f"Refinement failed: {e}")
    res = store.get(nm)
    if res is None:
        box.caption("Choose the phase(s) and press **Refine**. A pattern takes about a second.")
        return
    stt = res["stats"]
    rw = stt["Rwp"]
    verdict = "✅ good" if rw < 0.10 else ("🟡 acceptable" if rw < 0.20 else "🔴 poor")
    if res.get("profile_fixed"):
        box.caption("🔒 Peak profile (U, V, W, X, Y) fixed from the Si standard; the sample's cell and atoms are refined on top of it.")
    m = st.columns(2 + len(res["phases"]))
    m[0].metric("Rwp", f"{100 * rw:.1f} %", verdict)
    m[1].metric("GoF" if res["counts"] else "Rp", f"{math.sqrt(stt['chi2']):.2f}" if res["counts"] else f"{100 * stt['Rp']:.1f} %",
                help="GoF (= Rwp/Rexp) close to 1 means the fit is as good as the counting noise allows." if res["counts"] else "")
    for k, p in enumerate(res["phases"]):
        a = p["cell"][0]
        ea = (p["cell_esd"] or {}).get("a")
        m[2 + k].metric(p["name"][:22], f"{100 * p['weight_fraction']:.1f} wt %" if res["mode"] == "rietveld" else f"a = {a:.4f} Å",
                        f"a = {a:.4f}" + (f"({round(ea * 1e4):d})" if ea else "") + " Å" if res["mode"] == "rietveld" else None, delta_color="off")
    for n in res["notes"][:-1]:
        st.caption("⚠️ " + n)
    st.caption("ℹ️ " + res["notes"][-1])

    def draw(ax, res=res):
        RV.draw(ax, res)
    fig, ax = newfig(aspect=0.75)
    RV.draw(ax, res)
    show(fig, f"PXRD_Rietveld_{nm}" if res["mode"] == "rietveld" else f"PXRD_LeBail_{nm}",
         pd.DataFrame({"two_theta": res["x"], "observed": res["y"], "calculated": res["calc"], "background": res["bg"], "difference": res["resid"]}),
         draw, "Rietveld refinement" if res["mode"] == "rietveld" else "Le Bail fit", units={"x": "2theta"})
    tab = []
    for p in res["phases"]:
        row = {"phase": p["name"], "cell (Å)": ", ".join(f"{v:.4f}" for v in p["cell"][:3]) if p["system"] != "cubic" else f"{p['cell'][0]:.4f}",
               "esd a": (p["cell_esd"] or {}).get("a"), "V (Å³)": p["volume"], "density (g/cm³)": p["density"]}
        if res["mode"] == "rietveld":
            row.update({"weight %": 100 * p["weight_fraction"], "Bragg R (%)": 100 * p["Rb"], "B (Å²)": p["B"]})
            if p.get("po_r") is not None:
                row["texture r"] = p["po_r"]
        tab.append(row)
    st.dataframe(pd.DataFrame(tab).style.format(precision=4, na_rep="–"), hide_index=True, width="stretch")
    with LV.gate("Advanced", globals()):
        with st.expander("🔬 Refined profile parameters (± esd)"):
            st.dataframe(pd.DataFrame([{"parameter": v["label"], "value": v["value"], "± esd": v["esd"], "refined": v["refined"]} for v in res["global"].values()])
                         .style.format(precision=5, na_rep="–"), hide_index=True, width="stretch")
            if np.isfinite(res["size_nm_apparent"]) and res["size_nm_apparent"] < 1000:
                st.caption(f"Apparent size from the Lorentzian width: {res['size_nm_apparent']:.0f} nm, apparent strain {100 * res['strain_apparent']:.2f} % "
                           "(instrument broadening included – a lower limit for the size).")
    for _i, _p in enumerate(res["phases"]):                   # v2.16: refined structure (3D, tables, figure, CIF)
        if _p.get("structure"):
            from mfs import struct_ui as SU
            SU.render(res, _i, nm, {"show": show, "newfig": newfig, "dl": dl}, __version__)
    from mfs import refout as RO                                 # v2.18.4: 3D view, CIF, reflection list and .xy files for every refinement (not only the atom refinement)
    _pobjs = res.get("_phase_objs") or []
    for _i, _p in enumerate(res["phases"]):
        if not _p.get("structure"):
            RO.render(res, _i, _pobjs[_i] if _i < len(_pobjs) else None, nm, {"dl": dl}, __version__)
    _fn = re.sub(r"[^A-Za-z0-9_.-]+", "_", nm)
    dl(st, "⬇️ Everything to continue with (ZIP: CIF, .xy of the pattern, reflection lists, summary)", lambda: RO.bundle(res, _pobjs, nm, __version__), f"{_fn}_refinement_files.zip", "rv_zip",
       help="CIF per phase, observed / calculated / difference / background as .xy, the tick positions, the reflection lists (CSV and .hkl) and a JSON summary.")
    dl(st, "⬇️ Rietveld report (Excel: summary, phases, profile, pattern, notes)", lambda: RV.report_xlsx(res, nm), f"Rietveld_{_fn}.xlsx", "rv_xlsx")
    _dom = max(res["phases"], key=lambda p_: p_["weight_fraction"] if p_["weight_fraction"] == p_["weight_fraction"] else 0)
    for p in res["phases"]:
        vals = {f"{p['name']} cell a (Å)": p["cell"][0], f"{p['name']} volume (Å^3)": p["volume"], "Rwp (%)": 100 * rw}
        if p is _dom:
            vals["density_g_cm3"] = p["density"]
        if res["mode"] == "rietveld":
            vals[f"{p['name']} weight fraction (%)"] = 100 * p["weight_fraction"]
        store_put(nm, "pxrd", vals, "PXRD Rietveld" if res["mode"] == "rietveld" else "PXRD Le Bail")


def _band_tools(traces, kind, lo, hi, is_T):
    """v2.15: assign the bands of a spectrum from an editable list of typical ranges, and fit overlapping bands simultaneously."""
    from mfs import bands as BD, spectra as SP
    if not traces:
        return
    tech = "raman" if kind == "raman" else "ir"
    names = [t[0] for t in traces]
    c1, c2, c3 = st.columns([2, 1, 1])
    nm = c1.selectbox("Spectrum", names, key="bd_sample")
    prom = NE.slider_in(c2, "Peak prominence", 0.01, 0.5, 0.05, 0.01, key="bd_prom")
    nmax = c3.number_input("Max. bands", 1, 40, 12, key="bd_nmax")
    x, y = next((t[1], t[2]) for t in traces if t[0] == nm)
    pk = SP.pick_peaks(x, y, prom, int(nmax), (lo, hi), minima=is_T)
    ss = st.session_state
    ss.setdefault("band_user", [])
    with st.expander("✏️ My own assignments (take precedence; saved in the project)", expanded=False):
        own = st.data_editor(pd.DataFrame(ss["band_user"], columns=["from (cm⁻¹)", "to (cm⁻¹)", "assignment"]).astype({"assignment": str}),
                             num_rows="dynamic", key="bd_own", hide_index=True, width="stretch")
        ss["band_user"] = own.dropna(how="all").to_dict("records")
    rows = []
    for a in BD.assign(pk, ss["band_user"], tech):
        c = a["candidates"]
        rows.append({"band (cm⁻¹)": round(a["position"]), "assignment": c[0][0] if c else "–", "source": c[0][1] if c else "",
                     "other possibilities": "; ".join(t[0] for t in c[1:]) if len(c) > 1 else ""})
    if rows:
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        st.caption("Typical ranges from standard IR / Raman tables – indicative only: confirm with literature for your framework.")
    else:
        st.info("No bands found – lower the prominence.")
        return
    st.markdown("**Fit overlapping bands (pseudo-Voigt)**")
    f1, f2 = st.columns([3, 2])
    top = max(pk, key=lambda v: abs(float(np.interp(v, np.sort(x), y[np.argsort(x)])) - float(np.median(y))))      # strongest band → its neighbourhood
    near = sorted(v for v in pk if abs(v - top) <= 150)
    txt = f1.text_input("Band positions (cm⁻¹, comma separated)", ", ".join(f"{v:.0f}" for v in near), key="bd_pos")
    win = f2.text_input("Fit window (cm⁻¹, from–to)", f"{max(float(x.min()), min(near) - 60):.0f}-{min(float(x.max()), max(near) + 60):.0f}", key="bd_win")
    if st.button("▶ Fit bands", key="bd_go", type="primary"):
        try:
            cs = [float(v) for v in txt.replace(";", ",").split(",") if v.strip()]
            a_, b_ = sorted(float(v) for v in win.replace("–", "-").split("-")[:2])
            yy = -y if is_T else y
            r = BD.fit_bands(x, yy, cs, a_, b_)
            r["is_T"] = is_T
            ss["bd_res"] = (nm, r)
        except Exception as e:
            st.warning(f"Fit failed: {e}")
    got = ss.get("bd_res")
    if got and got[0] == nm:
        r = got[1]
        tab = pd.DataFrame([{"centre (cm⁻¹)": b["center"], "± ": b["center_err"], "FWHM (cm⁻¹)": b["fwhm"], "± FWHM": b["fwhm_err"], "area": b["area"],
                             "area %": 100 * b["area"] / sum(q["area"] for q in r["bands"]), "η (Lorentz share)": b["eta"]} for b in r["bands"]])
        st.dataframe(tab.style.format(precision=2), hide_index=True, width="stretch")
        st.caption(f"R² = {r['r2']:.4f}, rms = {r['rms']:.4g}. Areas are in the units of the plotted spectrum × cm⁻¹; compare them between spectra that are normalised the same way.")

        def draw(ax, r=r):
            BD.draw_fit(ax, r, invert_x=kind != "raman")
            ax.set_xlabel(P.LABELS[kind][0])
            ax.set_ylabel("−Transmittance (a.u.)" if r.get("is_T") else P.LABELS[kind][1])
        fig, ax = newfig(0.8); draw(ax)
        show(fig, ("FTIR" if kind != "raman" else "Raman") + "_bandfit", tab, draw, "Band fit", units={"x": "wavenumber"} if kind != "raman" else None)


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
    lo, hi = NE.slider_in(c1, "Range (cm⁻¹)", 100.0, 4000.0, (400.0, 2000.0) if "FTIR" in mode else (150.0, 1800.0), 10.0)
    conv = c2.checkbox("Convert transmittance → absorbance", False) if "transmittance" in mode else False
    hl = c2.text_input("Highlight bands (lo-hi:label; …)", "1700-1740:C=O; 1190-1230:S=O; 1030-1055:" if "FTIR" in mode else "")
    highlight = []
    for part in [p for p in hl.split(";") if p.strip()]:
        try:
            rng_, lab = part.split(":") if ":" in part else (part, "")
            a, b = [float(v) for v in rng_.split("-")]; highlight.append((a, b, lab.strip()))
        except Exception:
            pass
    multi = st.checkbox("The file(s) contain several spectra: one column per spectrum (first column = wavenumber)", False, key="ft_multi",
                        help="v2.18.2: each column after the first becomes its own spectrum; the column header is used as its name.")
    traces, rows = [], []
    for nm, f in samples:
        try:
            parts = R.read_ftir_multi(f, f.name) if multi else [(None, R.read_ftir(f, f.name))]
        except Exception as e:
            st.error(f"{f.name}: {e}"); continue
        for lab, d in parts:
            nm_ = nm if lab is None else f"{nm} · {lab}"
            y = A.transmittance_to_absorbance(d.y.values) if conv else d.y.values
            y = A.process(d.x.values, y, **opts)
            traces.append((nm_, d.x.values, y)); rows.append(pd.DataFrame({"sample": nm_, "wavenumber": d.x, "y_processed": y}))
    kind = "raman" if "Raman" in mode else ("ftir_A" if ("absorbance" in mode or conv) else "ftir_T")
    from mfs import spectra as SP
    is_T = kind == "ftir_T"
    with LV.gate("Advanced", globals()):
        with st.expander("🧪 Reference band, peak labels, zoomed inset, difference spectra", expanded=False):
            c1, c2, c3, c4 = st.columns(4)
            nb = c1.checkbox("Normalise to a reference band", False, key="ft_nb")
            dlo, dhi, dlab = (1540.0, 1620.0, "νas(COO⁻)") if kind != "raman" else (1630.0, 1690.0, "ν(C=C)")
            blo = c2.number_input("band from (cm⁻¹)", 50.0, 4000.0, dlo, key="ft_blo")
            bhi = c3.number_input("band to (cm⁻¹)", 50.0, 4000.0, dhi, key="ft_bhi")
            bmode = c4.radio("by", ["height", "area"], horizontal=True, key="ft_bm")
            c1, c2, c3, c4 = st.columns(4)
            pl = c1.checkbox("Label peaks", False, key="ft_pl")
            prom = NE.slider_in(c2, "min. prominence (fraction of range)", 0.01, 0.5, 0.08, 0.01, key="ft_prom")
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
    with LV.gate("Advanced", globals()):
        with st.expander("🏷️ Band assignment and multi-band fit", expanded=False):
            _band_tools(traces, kind, lo, hi, is_T)
    if len(traces) >= 2:
        with LV.gate("Advanced", globals()):
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
