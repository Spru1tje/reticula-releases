"""UV/Vis page (v2.9). Rendered from app.py with the app's helpers passed in (show, newfig, files_input …)."""
from __future__ import annotations

import os
import re

import numpy as np
import pandas as pd
import streamlit as st

from mfs import uvvis as U

YLAB = {"A": "Absorbance", "T": "Transmittance (%)", "R": "Reflectance (%)", "KM": "F(R) (Kubelka–Munk)",
        "dA": "dA/dλ (nm⁻¹)", "d2A": "d²A/dλ² (nm⁻²)"}


@st.cache_data(show_spinner=False, max_entries=64)
def _read(data: bytes, name: str):
    f = U.read_any(data, name)
    return f


def render(g):
    show, newfig, files_input, sample_table, store_put, EX = (g["show"], g["newfig"], g["files_input"], g["sample_table"],
                                                              g["store_put"], g["EX"])
    st.header("UV/Vis spectroscopy")
    st.caption("Implen NanoPhotometer / NPOS Excel export (.xlsx: results sheet + one tab per sample), Implen .ids (only if not "
               "encrypted – otherwise export to Excel in NPOS), Agilent Cary / Shimadzu / generic CSV-TXT, JCAMP-DX. "
               "Solutions (absorbance / %T) and solids (diffuse reflectance %R → Kubelka–Munk, Tauc).")
    files = files_input("UV/Vis files", ["xlsx", "xls", "csv", "txt", "dat", "ids", "jdx", "dx", "asc"],
                        "UVVis/*.*", "uv")
    if not files:
        return
    spectra, infos, results = [], [], []
    for f in files:
        try:
            uf = _read(f.getvalue(), f.name)
        except U.IDSEncrypted as e:
            st.warning(str(e)); continue
        except Exception as e:
            st.error(f"{f.name}: {e}"); continue
        for s in uf.spectra:
            s.source = f.name
        spectra.extend(uf.spectra)
        infos.append({"file": f.name, "instrument": uf.instrument or "–", "software": uf.software or "–",
                      "spectra": len(uf.spectra), "sheets": ", ".join(map(str, uf.sheets[:6])) + (" …" if len(uf.sheets) > 6 else ""),
                      "settings": "; ".join(f"{k}: {v}" for k, v in list(uf.meta.items())[:8])})
        if uf.results is not None:
            results.append((f.name, uf.results))
    if not spectra:
        st.info("No spectra found in the file(s)."); return
    with st.expander(f"📄 Files – {len(spectra)} spectra", expanded=False):
        st.dataframe(pd.DataFrame(infos).astype(str), hide_index=True, width="stretch")
        for fn, t in results:
            st.markdown(f"**Results table in {fn}** (as exported by the instrument)")
            st.dataframe(t.astype(str), hide_index=True, width="stretch")
    names = []
    for s in spectra:            # unique display names across files
        n, k = s.name, 2
        while n in names:
            n = f"{s.name} ({k})"; k += 1
        s.name = n; names.append(n)
    kinds = sorted({s.ytype for s in spectra})
    sel = g["sample_names_list"](names, "uv_samples", max_shown=12)
    tabs = st.tabs(["📈 Spectra", "📏 Calibration / assays", "⏱️ Kinetics", "🟫 Diffuse reflectance / Tauc", "∂ Derivatives"])

    # ------------------------------------------------------------------ processing shared by the tabs
    with tabs[0]:
        c1, c2, c3, c4 = st.columns(4)
        c1.caption("Spectra: use the 🧪 sample list above (show/hide, colour, rename, order).")
        yshow = c2.selectbox("Plot as", ["absorbance", "%T", "%R (as measured)", "Kubelka–Munk F(R)"],
                             index=2 if kinds == ["R"] else (1 if kinds == ["T"] else 0), key="uv_y")
        blank = c3.selectbox("Blank / baseline spectrum", ["none"] + names, key="uv_blank")
        xr = c4.slider("λ range (nm)", 150.0, 1100.0, (float(max(150, min(s.x.min() for s in spectra))),
                                                         float(min(1100, max(s.x.max() for s in spectra)))), 1.0, key="uv_xr")
        c1, c2, c3, c4 = st.columns(4)
        sm = c1.slider("Smoothing (SG window, points; 0 = off)", 0, 51, 0, 2, key="uv_sm")
        norm = c2.selectbox("Normalise", ["none", "max", "min-max", "at λ", "area"], key="uv_norm")
        nat = c2.number_input("λ for 'at λ' (nm)", 150.0, 1100.0, 280.0, 1.0, key="uv_nat") if norm == "at λ" else None
        off = g["offset_controls"](c3, "uv_off", 0.0, 0.0, 10.0, 0.05)
        zero = c3.checkbox("Zero at a non-absorbing λ", False, key="uv_zero")
        zl = c3.number_input("zero λ (nm)", 150.0, 1100.0, 750.0, 5.0, key="uv_zl") if zero else None
        pk = c4.checkbox("Mark λmax / peaks", True, key="uv_pk")
        prom = c4.slider("peak prominence (fraction)", 0.01, 0.5, 0.08, 0.01, key="uv_prom") if pk else 0.08
        path = c4.number_input("Path length (cm)", 0.001, 10.0, 1.0, 0.1, key="uv_path",
                               help="NanoPhotometer data are normalised to 10 mm unless the export says otherwise.")

    def y_of(s):
        if yshow == "absorbance":
            y = U.to_absorbance(s); yt = "A"
        elif yshow == "%T":
            y = U.to_percent_T(s); yt = "T"
        elif yshow.startswith("Kubelka"):
            y = U.kubelka_munk(s.y if s.ytype == "R" else U.to_percent_T(s)); yt = "KM"
        else:
            y = s.y.copy(); yt = s.ytype
        return y, yt

    lut = {s.name: s for s in spectra}
    proc, rows = [], []
    bl = lut.get(blank)
    for n0, n in sel:
        s = lut[n0]
        y, yt = y_of(s)
        if bl is not None and bl is not s:
            yb, _ = y_of(bl)
            y = U.subtract_blank(s.x, y, bl.x, yb)
        if sm >= 5:
            y = U.smooth(y, sm, 3)
        if zero and zl:
            y = y - U.value_at(s.x, y, zl, 2.0)
        y = U.normalise(s.x, y, norm, nat, xr)
        proc.append((n, s.x, y, yt))
        rows.append(pd.DataFrame({"sample": n, "wavelength_nm": s.x, "y": y}))
    g["set_even_offset"]("uv_off", proc, xlim=xr)
    ytype_plot = proc[0][3] if proc else "A"
    ylab = YLAB.get(ytype_plot, "Absorbance") + (" (norm.)" if norm != "none" else "")
    peak_rows = []

    with tabs[0]:
        if not proc:
            st.info("Select at least one spectrum.")
        else:
            def draw(ax, proc=proc):
                for i, (n, x, y, yt) in enumerate(proc):
                    m = (x >= xr[0]) & (x <= xr[1])
                    ax.plot(x[m], y[m] + i * off, lw=1, label=n)
                    if pk:
                        for lx, ly in U.peaks(x, y, xr, prom, 4, minima=(yt == "T")):
                            ax.annotate(f"{lx:.0f}", (lx, ly + i * off), xytext=(0, 3), textcoords="offset points",
                                        ha="center", va="bottom", fontsize="x-small")
                ax.set_xlim(*xr); ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel(ylab)
                if len(proc) <= 12:
                    ax.legend(frameon=False, fontsize="small")
            fig, ax = newfig(); draw(ax)
            show(fig, "UVVis_spectra", pd.concat(rows), draw, "UV/Vis spectra", units={"x": "wavelength"})
            for n, x, y, yt in proc:
                lm, am = U.lambda_max(x, y if yt != "T" else -y, xr)
                pks = U.peaks(x, y, xr, prom, 6, minima=(yt == "T"))
                peak_rows.append({"sample": n, "λmax (nm)": round(lm, 1), "value at λmax": round(float(np.interp(lm, x, y)), 4),
                                  "peaks (nm)": ", ".join(f"{p[0]:.0f}" for p in pks)})
            st.dataframe(pd.DataFrame(peak_rows), hide_index=True, width="stretch")
            if st.button("🗂️ Put λmax into the sample store", key="uv_store_lmax"):
                for r in peak_rows:
                    store_put(re.sub(r"\s*\(\d+\)$", "", r["sample"]), "uvvis", {"lambda_max_nm": r["λmax (nm)"],
                              "A_at_lambda_max": r["value at λmax"], "peaks_nm": r["peaks (nm)"]}, "UV/Vis")
                st.success(f"Stored λmax for {len(peak_rows)} sample(s).")

    # ------------------------------------------------------------------ calibration / assays
    with tabs[1]:
        c1, c2, c3 = st.columns(3)
        assay = c1.selectbox("Assay / method", ["custom"] + list(U.ASSAYS), key="uv_assay")
        a_def = U.ASSAYS.get(assay)
        lam = c2.number_input("Wavelength (nm)", 150.0, 1100.0, float(a_def[0]) if a_def else 280.0, 1.0, key="uv_cal_lam")
        hw = c2.number_input("± band for averaging (nm)", 0.0, 10.0, 0.0, 0.5, key="uv_cal_hw")
        model = c3.selectbox("Model", ["linear", "through origin", "quadratic"],
                             index=["linear", "through origin", "quadratic"].index(a_def[4]) if a_def else 0, key="uv_cal_model")
        unit = c3.text_input("Concentration unit", a_def[3] if a_def else "µM", key="uv_cal_unit")
        if a_def:
            st.caption(f"**{assay}** – λ = {a_def[0]} nm, standard: {a_def[1]}. {a_def[5]}")
        conc_guess = {}
        for n in names:
            m = re.search(r"(?<![\w.])(\d+(?:[.,]\d+)?)\s*(?:ug|µg|mg|mM|µM|uM|nM|M\b|g/L|ppm)", n, re.I)
            conc_guess[n] = float(m.group(1).replace(",", ".")) if m else None
        # NPOS results table: take 'Standard conc.' if present
        for fn, t in results:
            sc = next((c for c in t.columns if re.search(r"standard|std|conc", str(c), re.I)), None)
            nc = next((c for c in t.columns if re.search(r"sample", str(c), re.I)), None)
            if sc and nc:
                for _, r in t.iterrows():
                    v = U._num(r[sc])
                    if str(r[nc]) in conc_guess and np.isfinite(v):
                        conc_guess[str(r[nc])] = v
        tab = pd.DataFrame({"spectrum": names, f"A({lam:.0f} nm)": [round(U.value_at(lut[n].x, U.to_absorbance(lut[n]), lam, hw)
                                                                         - (U.value_at(bl.x, U.to_absorbance(bl), lam, hw) if bl is not None and bl is not lut[n] else 0), 5)
                                                                   for n in names],
                            "role": ["standard" if conc_guess[n] is not None else "unknown" for n in names],
                            f"c ({unit})": [conc_guess[n] for n in names], "dilution ×": 1.0})
        st.caption("Edit the roles / concentrations (blank-corrected absorbance with the blank chosen on the Spectra tab).")
        ed = st.data_editor(tab, hide_index=True, width="stretch", key="uv_cal_tab",
                            column_config={"role": st.column_config.SelectboxColumn(options=["standard", "unknown", "ignore"])},
                            disabled=["spectrum", f"A({lam:.0f} nm)"])
        stds = ed[(ed.role == "standard") & ed[f"c ({unit})"].notna()]
        acol = f"A({lam:.0f} nm)"
        if len(stds) >= (4 if model == "quadratic" else 3):
            try:
                cal = U.calibration(stds[f"c ({unit})"].astype(float), stds[acol].astype(float), model)
            except Exception as e:
                st.error(str(e)); cal = None
            if cal:
                unk = ed[ed.role == "unknown"].copy()
                cu, se = U.predict_conc(cal, unk[acol].astype(float).values) if len(unk) else (np.array([]), np.array([]))
                unk[f"c ({unit})"] = np.round(cu * unk["dilution ×"].astype(float).values, 5)
                unk[f"SE ({unit})"] = np.round(se * unk["dilution ×"].astype(float).values, 5)
                unk["in range"] = [(lo_ := min(cal["c"])) <= v / d <= max(cal["c"]) if np.isfinite(v) else False
                                   for v, d in zip(unk[f"c ({unit})"], unk["dilution ×"].astype(float))]
                mtxt = (f"A = {cal['slope']:.5g}·c + {cal['intercept']:.4g}" if model != "quadratic"
                        else f"A = {cal['quad']:.3g}·c² + {cal['slope']:.5g}·c + {cal['intercept']:.4g}")
                k = st.columns(5)
                k[0].metric("R²", f"{cal['R2']:.5f}"); k[1].metric("slope", f"{cal['slope']:.4g} ± {cal['slope_se']:.2g}")
                k[2].metric("s_y/x", f"{cal['s_yx']:.3g}"); k[3].metric(f"LOD ({unit})", f"{cal['LOD']:.3g}")
                k[4].metric(f"LOQ ({unit})", f"{cal['LOQ']:.3g}")
                st.caption(f"{mtxt}; n = {cal['n']}. LOD = 3.3·s_y/x / slope, LOQ = 10·s_y/x / slope (ICH Q2). "
                           + (f"ε = slope / l = {U.epsilon(cal, path):.4g} {unit}⁻¹ cm⁻¹." if model != "quadratic" else ""))

                def draw_cal(ax, cal=cal, unk=unk):
                    cc = np.asarray(cal["c"]); aa = np.asarray(cal["A"])
                    xx = np.linspace(0, max(cc.max(), np.nanmax(unk[f"c ({unit})"] / unk["dilution ×"].astype(float)) if len(unk) else 0) * 1.05, 200)
                    yy = (cal.get("quad", 0) * xx ** 2 + cal["slope"] * xx + cal["intercept"])
                    ax.plot(cc, aa, "o", ms=4, label="standards")
                    ax.plot(xx, yy, "-", lw=1, label=f"fit (R² = {cal['R2']:.4f})")
                    if len(unk):
                        ax.plot(unk[f"c ({unit})"] / unk["dilution ×"].astype(float), unk[acol], "s", ms=4, mfc="none", label="unknowns")
                    ax.axhline(cal["intercept"] + 3.3 * cal["s_yx"], ls=":", lw=0.7, c="0.5")
                    ax.set_xlabel(f"Concentration ({unit})"); ax.set_ylabel(f"A$_{{{lam:.0f}}}$")
                    ax.legend(frameon=False, fontsize="small")
                fig, ax = newfig(0.8); draw_cal(ax)
                show(fig, "UVVis_calibration", pd.concat([stds.assign(kind="standard"), unk.assign(kind="unknown")]), draw_cal,
                     "UV/Vis calibration")
                if len(unk):
                    st.dataframe(unk.astype(str), hide_index=True, width="stretch")
                    if st.button("🗂️ Put the concentrations into the sample store", key="uv_store_conc"):
                        for _, r in unk.iterrows():
                            store_put(re.sub(r"\s*\(\d+\)$", "", r["spectrum"]), "uvvis",
                                      {f"conc_{re.sub(r'[^A-Za-z0-9]+', '_', assay)}": r[f"c ({unit})"], "conc_unit": unit,
                                       "assay": assay, "lambda_nm": lam, "cal_R2": round(cal["R2"], 5)}, "UV/Vis")
                        st.success("Stored.")
                if assay.startswith(("OPA", "Ninhydrin")):
                    with st.expander("Degree of hydrolysis (Nielsen 2001)", expanded=False):
                        htot = st.number_input("h_tot (meq amino groups per g protein, e.g. 7.6 for whey, 8.0 BSA)", 0.1, 20.0, 8.0, 0.1, key="uv_htot")
                        pconc = st.number_input("protein in the assay (g/L)", 0.001, 100.0, 1.0, 0.1, key="uv_pconc")
                        if len(unk):
                            h = unk[f"c ({unit})"].astype(float) / (1.0 if unit == "mM" else 1000.0) / pconc     # meq/g
                            st.dataframe(pd.DataFrame({"spectrum": unk["spectrum"], "h (meq/g)": h.round(3),
                                                       "DH (%)": U.degree_of_hydrolysis(h, htot).round(1)}), hide_index=True)
        else:
            st.info(f"Mark at least {4 if model == 'quadratic' else 3} spectra as standards with a concentration.")

    # ------------------------------------------------------------------ kinetics
    with tabs[2]:
        c1, c2, c3, c4 = st.columns(4)
        klam = c1.number_input("λ (nm)", 150.0, 1100.0, 340.0, 1.0, key="uv_k_lam")
        kunit = c2.selectbox("time unit", ["min", "h", "s"], key="uv_k_unit")
        kfit = c3.checkbox("First-order fit A(t) = A∞ + (A₀ − A∞)e^(−kt)", True, key="uv_k_fit")
        kref = c4.number_input("reference λ to subtract (0 = none)", 0.0, 1100.0, 0.0, 5.0, key="uv_k_ref")
        kt = pd.DataFrame({"spectrum": names, f"t ({kunit})": [U.time_from_name(n, kunit) for n in names],
                           f"A({klam:.0f})": [U.value_at(lut[n].x, U.to_absorbance(lut[n]), klam)
                                              - (U.value_at(lut[n].x, U.to_absorbance(lut[n]), kref) if kref else 0) for n in names]})
        kt["use"] = kt[f"t ({kunit})"].notna()
        st.caption("Times are read from the names (e.g. 't40min', '2 h'); edit them here.")
        ke = st.data_editor(kt.round(5), hide_index=True, width="stretch", key="uv_k_tab", disabled=["spectrum", f"A({klam:.0f})"])
        kd = ke[ke.use & ke[f"t ({kunit})"].notna()].sort_values(f"t ({kunit})")
        if len(kd) >= 2:
            fit = None
            if kfit and len(kd) >= 4:
                try:
                    fit = U.first_order_fit(kd[f"t ({kunit})"].values, kd[f"A({klam:.0f})"].values)
                except Exception as e:
                    st.warning(f"Fit failed: {e}")

            def draw_k(ax, kd=kd, fit=fit):
                tt = kd[f"t ({kunit})"].values
                ax.plot(tt, kd[f"A({klam:.0f})"].values, "o", ms=4, label="data")
                if fit:
                    xx = np.linspace(0, tt.max(), 200)
                    ax.plot(xx, fit["fn"](xx), "-", lw=1, label=f"k = {fit['k']:.3g} {kunit}⁻¹")
                ax.set_xlabel(f"Time ({kunit})"); ax.set_ylabel(f"A$_{{{klam:.0f}}}$"); ax.legend(frameon=False, fontsize="small")
            fig, ax = newfig(0.75); draw_k(ax)
            show(fig, "UVVis_kinetics", kd, draw_k, "UV/Vis kinetics")
            if fit:
                k = st.columns(4)
                k[0].metric(f"k ({kunit}⁻¹)", f"{fit['k']:.4g}", f"± {fit['k_se']:.2g}", delta_color="off")
                k[1].metric(f"t½ ({kunit})", f"{fit['t_half']:.3g}"); k[2].metric("A₀ → A∞", f"{fit['A0']:.3f} → {fit['Ainf']:.3f}")
                k[3].metric("R²", f"{fit['R2']:.4f}")
            sname = st.text_input("Sample name for the Kinetics page", "UV-Vis kinetics", key="uv_k_name")
            if st.button("➡️ Send to the Kinetics page (as a sample/time/concentration table)", key="uv_k_send"):
                yv = kd[f"A({klam:.0f})"].values
                st.session_state["_uv_kin_csv"] = pd.DataFrame({"sample": sname, "time": kd[f"t ({kunit})"].values,
                                                                "concentration": (yv - yv.min()) if yv[0] < yv[-1] else yv}).to_csv(index=False).encode()
                st.session_state["_uv_kin_unit"] = kunit
                st.success("Sent – open Kinetics & adsorption and tick 'Use the table sent from UV/Vis'. "
                           "(Growing signals are sent as product formed, A − A_min.)")

    # ------------------------------------------------------------------ DRS / Tauc
    with tabs[3]:
        st.caption("Solids: measure diffuse reflectance (%R vs a white standard, BaSO₄/Spectralon); F(R) = (1−R)²/2R ∝ α/s. "
                   "Tauc: (F·hν)^n vs hν, n = 2 direct allowed, n = ½ indirect allowed; Eg = x-intercept of the linear edge. "
                   "ZrO₂ ≈ 5.0–5.2 eV; Zr/Hf-fumarate MOFs ≈ 3.9–4.1 eV (linker-based).")
        c1, c2, c3, c4 = st.columns(4)
        tsel = c1.multiselect("Spectra", names, default=[n for n in names if lut[n].ytype == "R"][:6] or names[:3], key="uv_t_sel")
        src = c2.radio("Signal", ["F(R) from %R", "absorbance (thin film / solution)"],
                       index=0 if all(lut[n].ytype == "R" for n in tsel or names[:1]) else 1, key="uv_t_src")
        tk = c3.radio("Transition", ["direct (n = 2)", "indirect (n = ½)"], key="uv_t_kind")
        auto = c4.checkbox("Automatic linear window", True, key="uv_t_auto")
        tnorm = c4.checkbox("Scale each curve to its maximum", True, key="uv_t_norm")
        win = None
        if not auto:
            w1, w2 = st.columns(2)
            win = (w1.number_input("fit from (eV)", 0.5, 8.0, 3.8, 0.01, key="uv_t_w1"), w2.number_input("fit to (eV)", 0.5, 8.0, 4.2, 0.01, key="uv_t_w2"))
        tr, eg_rows = [], []
        for n in tsel:
            s = lut[n]
            y = (s.y if s.ytype == "R" else U.to_percent_T(s)) if src.startswith("F(R)") else U.to_absorbance(s)
            if sm >= 5:
                y = U.smooth(y, sm, 3)
            E, T = U.tauc(s.x, y, tk, from_reflectance=src.startswith("F(R)"))
            try:
                bg = U.band_gap(E, T, win)
            except Exception as e:
                bg = None; st.caption(f"{n}: {e}")
            tr.append((n, E, T, bg))
            if bg:
                eg_rows.append({"sample": n, "Eg (eV)": round(bg["Eg"], 3), "λ edge (nm)": round(U.HC / bg["Eg"], 1),
                                "window (eV)": f"{bg['window'][0]:.2f}–{bg['window'][1]:.2f}", "R²": round(bg["R2"], 4)})
        if tr:
            nexp = "2" if tk.startswith("direct") else "1/2"

            def draw_t(ax, tr=tr):
                for n, E, T, bg in tr:
                    k_ = (1.0 / np.nanmax(T)) if tnorm and np.nanmax(T) > 0 else 1.0
                    l, = ax.plot(E, T * k_, lw=1, label=n)
                    if bg:
                        xx = np.linspace(bg["Eg"], bg["window"][1] + 0.15, 30)
                        yy = (bg["slope"] * xx + bg["intercept"]) * k_
                        keep = yy <= np.nanmax(T * k_) * 1.05
                        ax.plot(xx[keep], yy[keep], "--", lw=0.8, c=l.get_color())
                        ax.annotate(f"{bg['Eg']:.2f} eV", (bg["Eg"], 0), xytext=(0, 6), textcoords="offset points",
                                    ha="center", fontsize="x-small", color=l.get_color())
                ax.set_xlabel("hν (eV)")
                ax.set_ylabel((f"(F(R)·hν)$^{{{nexp}}}$" if src.startswith("F(R)") else f"(A·hν)$^{{{nexp}}}$") + (" (norm.)" if tnorm else ""))
                ax.set_ylim(bottom=0); ax.legend(frameon=False, fontsize="small")
            fig, ax = newfig(); draw_t(ax)
            show(fig, "UVVis_Tauc", pd.concat([pd.DataFrame({"sample": n, "E_eV": E, "tauc": T}) for n, E, T, _ in tr]), draw_t, "Tauc plot")
            if eg_rows:
                st.dataframe(pd.DataFrame(eg_rows), hide_index=True, width="stretch")
                if st.button("🗂️ Put the band gaps into the sample store", key="uv_store_eg"):
                    for r in eg_rows:
                        store_put(re.sub(r"\s*%R$|\s*\(\d+\)$", "", r["sample"]).strip(), "uvvis",
                                  {("Eg_direct_eV" if tk.startswith("direct") else "Eg_indirect_eV"): r["Eg (eV)"]}, "UV/Vis")
                    st.success("Stored.")

    # ------------------------------------------------------------------ derivatives
    with tabs[4]:
        c1, c2, c3 = st.columns(3)
        dord = c1.radio("Order", ["1st", "2nd"], horizontal=True, key="uv_d_ord")
        dwin = c2.slider("SG window (points)", 5, 61, 15, 2, key="uv_d_win")
        dpo = c3.slider("SG polynomial", 2, 5, 3, key="uv_d_po")
        o = 1 if dord == "1st" else 2
        dtr = [(n, x, U.derivative(x, y, o, dwin, dpo)) for n, x, y, _ in proc]
        if dtr:
            def draw_d(ax, dtr=dtr):
                for n, x, d in dtr:
                    m = (x >= xr[0]) & (x <= xr[1])
                    ax.plot(x[m], d[m], lw=1, label=n)
                ax.axhline(0, lw=0.5, c="0.6"); ax.set_xlim(*xr)
                ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel(("d" if o == 1 else "d²") + "y/dλ" + ("" if o == 1 else "²"))
                if len(dtr) <= 12:
                    ax.legend(frameon=False, fontsize="small")
            fig, ax = newfig(); draw_d(ax)
            show(fig, f"UVVis_derivative_{o}", pd.concat([pd.DataFrame({"sample": n, "wavelength_nm": x, f"d{o}": d}) for n, x, d in dtr]),
                 draw_d, f"UV/Vis derivative ({dord})")
            zc = []
            for n, x, d in dtr:
                m = (x >= xr[0]) & (x <= xr[1])
                xs, ds = x[m], d[m]
                if o == 1:
                    yo = next(p_[2] for p_ in proc if p_[0] == n)[m]
                    mx = U.maxima_from_derivative(xs, yo, ds)
                    zc.append({"sample": n, "maxima from dA/dλ = 0 (nm)": ", ".join(f"{v:.1f}" for v in mx[:12])})
                else:
                    from scipy.signal import find_peaks
                    idx, _ = find_peaks(-ds, prominence=0.05 * (np.ptp(ds) or 1))
                    zc.append({"sample": n, "bands from d²A/dλ² minima (nm)": ", ".join(f"{xs[i]:.1f}" for i in idx[:12])})
            st.dataframe(pd.DataFrame(zc), hide_index=True, width="stretch")
