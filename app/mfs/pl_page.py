"""Photoluminescence page (v2.14): steady-state PL, temperature-dependent PL (thermal quenching, peak shift / width) and ratiometric
luminescence thermometry (Boltzmann band ratio, relative sensitivity, temperature uncertainty). Rendered from app.py with the app's
helpers passed in (show, newfig, files_input ...)."""
from __future__ import annotations

import hashlib
import os

import numpy as np
import pandas as pd
import streamlit as st

from mfs import pl as PL

EXAMPLES = {"Thermometry – Er³⁺ upconversion, 293–323 K": "PL/PL_thermometry_Er_upconversion_EXAMPLE.csv",
            "Thermal quenching – phosphor, 80–500 K": "PL/PL_thermal_quenching_phosphor_EXAMPLE.csv"}
XLAB = {"nm": "Wavelength (nm)", "eV": "Energy (eV)", "cm-1": "Wavenumber (cm$^{-1}$)"}


def _spectra_from(files):
    out = []
    for f in files:
        try:
            for s in PL.read_pl(f.getvalue(), f.name):
                s["file"] = f.name
                out.append(s)
        except Exception as e:
            st.error(f"{f.name}: {e}")
    return out


def _band_default(spectra, lo, hi):
    """Two strongest peaks of the mean spectrum -> two integration windows split at the valley between them."""
    from scipy.signal import find_peaks
    wl = spectra[0]["wl"]
    m = (wl >= lo) & (wl <= hi)
    y = np.mean([np.interp(wl, s["wl"], s["I"]) for s in spectra], axis=0)[m]
    x = wl[m]
    if len(x) < 10:
        return None
    pk, pr = find_peaks(y, prominence=0.05 * (y.max() - y.min() + 1e-12), distance=max(3, len(x) // 60))
    if len(pk) < 2:
        return None
    top = sorted(pk[np.argsort(pr["prominences"])[-2:]])
    a, b = int(top[0]), int(top[1])
    v = a + int(np.argmin(y[a:b + 1]))
    w1 = PL.peak_stats(x, y, x[max(0, a - 15)], x[min(len(x) - 1, a + 15)])["fwhm"]
    w2 = PL.peak_stats(x, y, x[max(0, b - 15)], x[min(len(x) - 1, b + 15)])["fwhm"]
    w1 = 8.0 if not np.isfinite(w1) else w1
    w2 = 8.0 if not np.isfinite(w2) else w2
    return (float(max(x[0], x[a] - 2.2 * w1)), float(x[v])), (float(x[v]), float(min(x[-1], x[b] + 2.2 * w2)))


def render(g):
    show, newfig, files_input, store_put, EX, dl, safe_name, C = (g["show"], g["newfig"], g["files_input"], g["store_put"], g["EX"],
                                                                 g["dl"], g["safe_name"], g["C"])
    from mfs import style as S
    st.header("Photoluminescence – steady-state and temperature-dependent")
    st.caption("Emission (or excitation) spectra as CSV / TXT / DAT / Excel from Horiba, Edinburgh, Ocean Optics, Hitachi … : two columns "
               "(wavelength, intensity) per file, or ONE table with the wavelength in the first column and one column per temperature "
               "(header like `293.15 K`, `25 °C`). Several files = several temperatures: the temperature is read from the column header or "
               "the file name (`..._300K.txt`, `..._25C.csv`) and can be edited below. Typical outputs: spectra coloured by temperature, "
               "emission map, thermal quenching (activation energy, T₅₀), peak shift / width, and ratiometric thermometry "
               "(Boltzmann fit, relative sensitivity Sᵣ, temperature uncertainty σT).")
    ex = st.radio("Bundled EXAMPLE (synthetic)", list(EXAMPLES), horizontal=True, key="pl_example")
    files = files_input("PL spectra", ["csv", "txt", "dat", "xlsx", "xls"], EXAMPLES[ex], "pl", multiple=True)
    if not files:
        return
    spectra = _spectra_from(files)
    if not spectra:
        return
    # ------------------------------------------------------------ temperatures (editable)
    c1, c2 = st.columns([1, 3])
    unit = c1.radio("Plain numbers in headers are", ["auto", "K", "°C"], horizontal=True, key="pl_unit")
    tk = PL.to_kelvin([PL.parse_temperature(s["T_text"]) for s in spectra], {"auto": "auto", "K": "K", "°C": "C"}[unit])
    tkey = hashlib.md5("|".join(s["file"] + s["name"] for s in spectra).encode()).hexdigest()[:8] + unit
    tab0 = pd.DataFrame({"use": True, "spectrum": [s["name"] for s in spectra], "T (K)": [None if t is None else round(t, 2) for t in tk]})
    if len(spectra) > 12:
        c2.caption(f"{len(spectra)} spectra found.")
    ed = st.data_editor(tab0, key="pl_tab|" + tkey, hide_index=True, width="stretch", height=min(300, 38 + 35 * len(spectra)),
                        column_config={"T (K)": st.column_config.NumberColumn("T (K)", format="%.2f", help="edit if the temperature was not recognised")})
    sel = [(s, t) for s, (_, r), t in zip(spectra, ed.iterrows(), ed["T (K)"]) if r["use"]]
    if not sel:
        st.info("No spectrum selected."); return
    for s, t in sel:
        s["T"] = None if (t is None or (isinstance(t, float) and np.isnan(t))) else float(t)
    spectra = [s for s, _ in sel]
    n_t = sum(1 for s in spectra if s["T"] is not None)
    # ------------------------------------------------------------ processing
    wl_lo, wl_hi = float(min(s["wl"].min() for s in spectra)), float(max(s["wl"].max() for s in spectra))
    p1, p2, p3, p4 = st.columns(4)
    axis = p1.radio("x axis", ["nm", "eV", "cm-1"], horizontal=True, key="pl_axis", format_func=lambda a: {"nm": "nm", "eV": "eV", "cm-1": "cm⁻¹"}[a])
    base = p2.selectbox("Baseline", ["none", "minimum", "region (mean)"], key="pl_base")
    norm = p3.selectbox("Normalisation", ["none (counts)", "max", "area", "to the lowest temperature"], key="pl_norm")
    sm = p4.slider("Smoothing (points)", 0, 31, 0, 2, key="pl_sm", help="Savitzky–Golay, order 3; 0 = none")
    rng = st.slider("Wavelength range (nm)", wl_lo, wl_hi, (wl_lo, wl_hi), 0.5, key="pl_rng")
    reg = None
    if base.startswith("region"):
        reg = st.slider("Baseline region (nm)", wl_lo, wl_hi, (wl_lo, min(wl_hi, wl_lo + 0.05 * (wl_hi - wl_lo))), 0.5, key="pl_reg")
    proc = []
    for s in spectra:
        m = (s["wl"] >= rng[0]) & (s["wl"] <= rng[1])
        wl, y = s["wl"][m], PL.smooth(s["I"], sm)[m]
        y = PL.baseline_subtract(wl, y, {"none": "none", "minimum": "minimum", "region (mean)": "region"}[base], reg)
        proc.append({**s, "wl_c": wl, "y_c": y})
    t_order = sorted([p for p in proc if p["T"] is not None], key=lambda p: p["T"])
    ref_val = None
    if norm.startswith("to the lowest") and t_order:
        ref_val = float(np.max(t_order[0]["y_c"])) or 1.0
    for p in proc:
        y = PL.normalise(p["y_c"], {"none (counts)": "none", "max": "max", "area": "area", "to the lowest temperature": "ref"}[norm], ref_val)
        p["x"], p["y"] = PL.to_axis(p["wl_c"], y, axis)
    ylab = "Intensity (counts)" if norm.startswith("none") else ("Normalised intensity" if norm != "area" else "Intensity (area = 1)")
    xlab = XLAB[axis]
    names = ["📈 Spectra"] + (["🌡️ Temperature dependence", "🔥 Thermometry (band ratio)"] if n_t >= 3 else [])
    tabs = st.tabs(names)
    sample = os.path.splitext(files[0].name)[0] if hasattr(files[0], "name") else "PL"

    # ============================================================ spectra
    with tabs[0]:
        views = ["Overlay", "Stacked (waterfall)"] + (["Emission map"] if n_t >= 3 else [])
        v1, v2 = st.columns([3, 1])
        view = v1.radio("View", views, horizontal=True, key="pl_view")
        cmap = v2.selectbox("Colour scale", PL.CMAPS, key="pl_cmap", help="viridis / cividis are colour-blind safe")
        shown = sorted(proc, key=lambda p: (p["T"] is None, p["T"] or 0))

        def draw_spec(ax, shown=shown, view=view, cmap=cmap, xlab=xlab, ylab=ylab):
            if view == "Emission map":
                PL.draw_map(ax, [{"T": p["T"], "x": p["x"], "y": p["y"]} for p in shown if p["T"] is not None], xlab, "magma", ylab)
            elif view.startswith("Stacked"):
                PL.draw_ridge(ax, shown, xlab, ylab + " + offset", cmap)
            else:
                PL.draw_spectra(ax, shown, xlab, ylab, cmap)
        fig, ax = newfig(0.8); draw_spec(ax)
        stem = "PL_" + {"Overlay": "spectra", "Emission map": "map"}.get(view, "waterfall")
        df = pd.DataFrame({"x": shown[0]["x"], **{(f"{p['T']:.2f} K" if p["T"] is not None else p["name"]): np.interp(shown[0]["x"], p["x"], p["y"]) for p in shown}})
        show(fig, stem, df, draw_spec, "PL " + view.lower())
        st.markdown("**Peaks** – strongest band of each spectrum (position by parabola through the maximum, FWHM at half height above the window minimum)")
        rows = []
        for p in shown:
            ps = PL.peak_stats(p["wl_c"], p["y_c"])
            rows.append({"spectrum": p["name"], "T (K)": p["T"], "peak (nm)": round(ps["pos"], 2), "peak (eV)": round(PL.HC / ps["pos"], 4) if ps["pos"] else None,
                         "FWHM (nm)": round(ps["fwhm"], 2) if np.isfinite(ps["fwhm"]) else None,
                         "FWHM (meV)": round(1000 * PL.HC * ps["fwhm"] / ps["pos"] ** 2, 1) if np.isfinite(ps["fwhm"]) else None,
                         "height": round(ps["height"], 1), "area (nm·counts)": round(ps["area"], 0)})
        pk = pd.DataFrame(rows)
        st.dataframe(pk, hide_index=True, width="stretch")
        dl(st, "CSV peaks", pk.to_csv(index=False).encode("utf-8-sig"), "PL_peaks.csv", "pl_peaks_csv")

    if n_t < 3:
        if n_t == 0 and len(spectra) > 1:
            st.info("No temperatures recognised – fill the *T (K)* column above to unlock the temperature-dependent analysis.")
        return
    T = np.array([p["T"] for p in t_order])

    # ============================================================ temperature dependence
    with tabs[1]:
        st.caption("Integrated intensity in a window → thermal quenching: I(T) = I₀ / (1 + A·exp(−Eₐ/k_BT)), Eₐ the activation energy, "
                   "T₅₀ = Eₐ/(k_B ln A) the temperature at which the intensity has halved. The Arrhenius plot ln(I₀/I − 1) vs 1/T "
                   "is a straight line of slope −Eₐ/k_B.")
        w = st.slider("Integration window (nm)", rng[0], rng[1], rng, 0.5, key="pl_win")
        I_T = np.array([PL.band_integral(p["wl_c"], p["y_c"], w[0], w[1], "linear") for p in t_order])
        pos = np.array([PL.peak_stats(p["wl_c"], p["y_c"], w[0], w[1])["pos"] for p in t_order])
        fw = np.array([PL.peak_stats(p["wl_c"], p["y_c"], w[0], w[1])["fwhm"] for p in t_order])
        fit = PL.fit_quench(T, I_T)
        okq = fit.get("status") == "ok"
        if okq:
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Eₐ", f"{fit['Ea_meV']:.0f} meV", help=f"{fit['Ea_eV'] * PL.EV_TO_CM:.0f} cm⁻¹")
            m2.metric("T₅₀", f"{fit['T50']:.0f} K" if fit["T50"] else "–", help="temperature at which the intensity is half of I₀")
            m3.metric("A", f"{fit['A']:.3g}")
            m4.metric("R²", f"{fit['R2']:.4f}")
        else:
            st.info(f"Thermal-quenching fit: {fit.get('status')}.")
        Tf = np.linspace(T.min(), T.max(), 200)
        In = I_T / (np.max(I_T) or 1.0)

        def d_int(ax, T=T, In=In, okq=okq, fit=fit, Tf=Tf, mx=np.max(I_T) or 1.0):
            PL.draw_trend(ax, T, In, "Integrated intensity (norm.)", fit_T=Tf if okq else None, fit_y=fit["f"](Tf) / mx if okq else None,
                          fit_label="quenching fit")

        def d_arr(ax, T=T, I_T=I_T, fit=fit):
            PL.draw_arrhenius(ax, T, I_T, fit)

        def d_pos(ax, T=T, pos=pos):
            lin = np.polyfit(T, pos, 1) if np.isfinite(pos).all() else None
            PL.draw_trend(ax, T, pos, "Peak position (nm)", fit_T=T if lin is not None else None,
                          fit_y=np.polyval(lin, T) if lin is not None else None, fit_label=f"{lin[0]:+.4f} nm K$^{{-1}}$" if lin is not None else None)

        def d_fwhm(ax, T=T, fw=fw):
            PL.draw_trend(ax, T, fw, "FWHM (nm)", legend=False)
        panels = [("Intensity vs T (+ quenching fit)", d_int)] + ([("Arrhenius plot", d_arr)] if okq else []) + \
                 [("Peak position vs T", d_pos), ("FWHM vs T", d_fwhm)]
        import matplotlib.pyplot as plt
        fig, axs = C.compose([f for _, f in panels], ncols=2, width=max(g["WIDTH"], 6.6), panel_aspect=0.78, label_style="A", label_inside=False, wspace=0.75, hspace=0.5)
        res = pd.DataFrame({"T (K)": T, "integrated intensity": I_T, "I / I_max": In, "peak (nm)": pos, "FWHM (nm)": fw})
        show(fig, "PL_temperature_dependence", res, None, "PL temperature dependence")
        pick = st.selectbox("Single panel (e.g. for the combined figure)", ["(none)"] + [n for n, _ in panels], key="pl_single1")
        if pick != "(none)":
            fn = dict(panels)[pick]
            fig2, ax2 = newfig(0.8); fn(ax2)
            show(fig2, "PL_" + safe_name(pick), res, fn, pick)
        st.dataframe(res.round(4), hide_index=True, width="stretch")
        dl(st, "CSV temperature dependence", res.to_csv(index=False).encode("utf-8-sig"), "PL_temperature_dependence.csv", "pl_td_csv")
        if okq:
            store_put(sample, "pl", {"Ea_meV": fit["Ea_meV"], "T50_K": fit["T50"], "A": fit["A"], "R2_quench": fit["R2"]}, "Photoluminescence page")

    # ============================================================ thermometry
    with tabs[2]:
        st.caption("Ratiometric thermometry with two thermally coupled bands: Δ = I₁/I₂ = α·exp(−ΔE/k_BT); "
                   "relative sensitivity Sᵣ = 100 %·|(1/Δ)·∂Δ/∂T| = 100 %·ΔE/(k_BT²); temperature uncertainty "
                   "σT = (1/Sᵣ)·√(1/I₁ + 1/I₂) for photon counts (the usual figures of merit of luminescent thermometers).")
        dflt = _band_default([{"wl": p["wl_c"], "I": p["y_c"]} for p in t_order], rng[0], rng[1])
        if not dflt:
            st.info("No two separate emission bands were found automatically (a single-band sample can only be used for the quenching analysis). "
                    "If your spectra have two bands, set the two windows by hand.")
        b1d, b2d = dflt if dflt else ((rng[0], (rng[0] + rng[1]) / 2), ((rng[0] + rng[1]) / 2, rng[1]))
        b1 = st.slider("Band 1 window (nm) – I₁ (e.g. ²H₁₁/₂ ≈ 520 nm)", rng[0], rng[1], tuple(round(v, 1) for v in b1d), 0.5, key="pl_b1")
        b2 = st.slider("Band 2 window (nm) – I₂ (e.g. ⁴S₃/₂ ≈ 540 nm)", rng[0], rng[1], tuple(round(v, 1) for v in b2d), 0.5, key="pl_b2")
        raw_counts = all(np.allclose(p["I"], np.round(p["I"])) for p in t_order) and norm.startswith("none")
        cc1, cc2 = st.columns(2)
        counts = cc1.checkbox("Intensities are photon counts (enables σT)", raw_counts, key="pl_counts",
                              help="σT assumes Poisson statistics of the integrated counts; for arbitrary units it is meaningless.")
        offs = cc2.checkbox("Add a constant C to the Boltzmann law (Δ = α·exp(−ΔE/k_BT) + C)", False, key="pl_off")
        I1 = np.array([PL.band_integral(p["wl_c"], p["y_c"], b1[0], b1[1], "linear") for p in t_order])
        I2 = np.array([PL.band_integral(p["wl_c"], p["y_c"], b2[0], b2[1], "linear") for p in t_order])
        if (I2 <= 0).any() or (I1 <= 0).any():
            st.warning("A band has no (positive) integrated intensity in its window – adjust the windows."); return
        D = I1 / I2
        bf = PL.fit_boltzmann(T, D, offset=offs)
        if bf.get("status") != "ok":
            st.warning(f"Boltzmann fit: {bf.get('status')}"); return
        Sr = PL.relative_sensitivity(T, bf["dE_eV"])
        Sa = PL.absolute_sensitivity(T, bf)
        sT = PL.temperature_uncertainty(Sr, I1, I2) if counts else None
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("ΔE", f"{bf['dE_cm1']:.0f} cm⁻¹", help=f"{1000 * bf['dE_eV']:.1f} meV")
        m2.metric("Sᵣ at " + f"{T.min():.0f} K", f"{Sr[0]:.3f} % K⁻¹")
        m3.metric("σT" + (f" at {T.min():.0f} K" if counts else ""), f"{sT[0]:.2f} K" if counts else "–")
        m4.metric("R²", f"{bf['R2']:.4f}")
        Tf = np.linspace(T.min(), T.max(), 200)
        mp = [{"T": p["T"], "x": p["wl_c"][(p["wl_c"] >= min(b1[0], b2[0]) - 15) & (p["wl_c"] <= max(b1[1], b2[1]) + 15)],
               "y": p["y_c"][(p["wl_c"] >= min(b1[0], b2[0]) - 15) & (p["wl_c"] <= max(b1[1], b2[1]) + 15)]} for p in t_order]

        def d_map(ax, mp=mp):
            PL.draw_map(ax, mp, "Wavelength (nm)", "magma", "Intensity (counts)" if norm.startswith("none") else "Intensity")

        def d_D(ax, T=T, D=D, bf=bf, Tf=Tf):
            PL.draw_trend(ax, T, D, "Δ = $I_1$/$I_2$", fit_T=Tf, fit_y=bf["f"](Tf),
                          fit_label=f"fit: ΔE = {bf['dE_cm1']:.0f} cm$^{{-1}}$ ($R^2$ = {bf['R2']:.3f})")

        def d_Sr(ax, T=T, Sr=Sr, Tf=Tf, bf=bf):
            PL.draw_trend(ax, T, Sr, r"$S_\mathrm{r}$ (% K$^{-1}$)", fit_T=Tf, fit_y=PL.relative_sensitivity(Tf, bf["dE_eV"]), legend=False)

        def d_sT(ax, T=T, sT=sT):
            ax.plot(T, sT, "o-", ms=3.2, lw=0.9, color="#0072B2")
            ax.set_xlabel("Temperature (K)"); ax.set_ylabel(r"$\sigma_T$ (K)")
        panels = [("Emission map", d_map), ("Δ vs T (Boltzmann fit)", d_D), ("Relative sensitivity Sr", d_Sr)] + ([("Temperature uncertainty σT", d_sT)] if counts else [])
        fig, axs = C.compose([f for _, f in panels], ncols=2, width=max(g["WIDTH"], 6.6), panel_aspect=0.78, label_style="A", label_inside=False, wspace=0.75, hspace=0.5)
        res = pd.DataFrame({"T (K)": T, "I1": I1, "I2": I2, "Δ = I1/I2": D, "Δ fit": bf["f"](T), "Sr (% per K)": Sr, "Sa (per K)": Sa})
        if counts:
            res["σT (K)"] = sT
        show(fig, "PL_thermometry", res, None, "PL thermometry")
        pick = st.selectbox("Single panel (e.g. for the combined figure)", ["(none)"] + [n for n, _ in panels], key="pl_single2")
        if pick != "(none)":
            fn = dict(panels)[pick]
            fig2, ax2 = newfig(0.8); fn(ax2)
            show(fig2, "PL_" + safe_name(pick), res, fn, pick)
        st.dataframe(res.round(5), hide_index=True, width="stretch")
        dl(st, "CSV thermometry", res.to_csv(index=False).encode("utf-8-sig"), "PL_thermometry.csv", "pl_th_csv")
        store_put(sample, "pl", {"dE_cm1": bf["dE_cm1"], "alpha": bf["alpha"], "Sr_Tmin_pct_per_K": float(Sr[0]), "R2_boltzmann": bf["R2"],
                                 **({"sigmaT_Tmin_K": float(sT[0])} if counts else {})}, "Photoluminescence page")
