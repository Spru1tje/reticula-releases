"""Photoluminescence page (v2.14): steady-state PL, temperature-dependent PL (thermal quenching, peak shift / width) and ratiometric
luminescence thermometry (Boltzmann band ratio, relative sensitivity, temperature uncertainty). Rendered from app.py with the app's
helpers passed in (show, newfig, files_input ...)."""
from __future__ import annotations

import hashlib
import os

import numpy as np
import pandas as pd
import streamlit as st
from mfs import levels as LV

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


# ------------------------------------------------------------------------------------------------ v2.15: robust band ratio
def _pick_cb(target):
    """Callback of the buttons under the graph: the range brushed on the graph becomes band 1 / band 2 / the fit range."""
    ss = st.session_state
    try:
        sel = ss["pl_pick"]["selection"]["pl_brush"]
        lo, hi = sorted(float(v) for v in next(iter(sel.values())))
    except Exception:
        ss["pl_pick_msg"] = "Drag a range on the graph first."
        return
    ss[target] = (round(lo, 1), round(hi, 1))
    if target == "pl_b1" or target == "pl_b2":
        ss["pl_method"] = METHODS[1]
    ss["pl_pick_msg"] = f"Range {lo:.1f}–{hi:.1f} nm set."


METHODS = ["Automatic – global fit of all spectra (recommended)", "Integration windows (choose them on the graph)", "Compare both"]


def _brush_picker(spectra, lo, hi, marks):
    """Mean spectrum with a drag-to-select brush. `marks` = [(lo, hi, label, colour)] shown as shaded ranges. The brushed range is stored by Streamlit
    under st.session_state['pl_pick'] (read by _pick_cb)."""
    import altair as alt
    grid = np.asarray(spectra[0]["x"], float)
    m = (grid >= lo) & (grid <= hi)
    ym = np.mean([np.interp(grid, np.asarray(q["x"], float), np.asarray(q["y"], float)) for q in spectra], axis=0)
    step = max(1, int(m.sum() // 1200))
    df = pd.DataFrame({"wavelength (nm)": grid[m][::step], "intensity": ym[m][::step]})
    brush = alt.selection_interval(encodings=["x"], name="pl_brush")
    base = alt.Chart(df).mark_line(color="#0F172A", strokeWidth=1.4).encode(x=alt.X("wavelength (nm):Q", scale=alt.Scale(zero=False)),
                                                                           y=alt.Y("intensity:Q", title="mean intensity"))
    layers = []
    for a, b, lab, col in marks:
        layers.append(alt.Chart(pd.DataFrame({"a": [a], "b": [b], "label": [lab]})).mark_rect(opacity=0.18, color=col).encode(x="a:Q", x2="b:Q"))
    ch = (alt.layer(*layers, base) if layers else base).properties(height=170).add_params(brush)
    return st.altair_chart(ch, on_select="rerun", key="pl_pick", width="stretch")


def _ratio_block(t_order, T, rng, norm, counts_default, g):
    """Returns (I1, I2, D, D_err, sT_override, info) for the thermometry tab, or None. D_err is None for the integration method."""
    show, newfig, dl = g["show"], g["newfig"], g["dl"]
    from mfs import plratio as PR
    ss = st.session_state
    c0, c1 = st.columns([3, 2])
    method = c0.radio("Method", METHODS, key="pl_method", horizontal=False,
                      help="The global fit describes ALL spectra with the same bands (also hidden shoulders and dirty neighbours) and gives the ratio "
                           "with an error bar and a quality flag per spectrum. Integration windows cut the overlap at a chosen valley.")
    raw_counts = all(np.allclose(p["I"], np.round(p["I"])) for p in t_order) and norm.startswith("none")
    counts = c1.checkbox("Intensities are photon counts (statistics-based weights and σT)", counts_default, key="pl_counts",
                         help="Raw counts: points are weighted by Poisson statistics. Arbitrary units: the noise is estimated from the data.")
    # raw (unsmoothed, un-normalised) data inside the chosen wavelength range: the fit has its own background
    spec = []
    for p in t_order:
        w = (p["wl"] >= rng[0]) & (p["wl"] <= rng[1])
        spec.append({"T": p["T"], "x": p["wl"][w], "y": p["I"][w]})
    use_fit = not method.startswith("Integration")
    fit = None
    bands_txt, extra_txt, bgo, model, add_extra = "", "", 1, "auto", True
    if use_fit and LV.at_least("Advanced"):
        with st.expander("⚙️ Fit settings (automatic by default)", expanded=False):
            a1, a2, a3 = st.columns(3)
            bands_txt = a1.text_input("Ratio bands: centres (nm), empty = automatic", "", key="pl_bands", help="e.g. 524, 545 – band 1 is the first, band 2 the second.")
            extra_txt = a2.text_input("Extra (nuisance) bands (nm), empty = automatic", "", key="pl_extra",
                                      help="Bands that are not part of the ratio but overlap with it (a shoulder, a third emission).")
            add_extra = a2.checkbox("Let the program add extra bands when needed", True, key="pl_addextra")
            bgo = a3.slider("Background polynomial order", 0, 3, 1, key="pl_bgo")
            model = {"automatic (BIC)": "auto", "shared shape": "none", "shape drifts with T": "drift"}[
                a3.selectbox("Peak model", ["automatic (BIC)", "shared shape", "shape drifts with T"], key="pl_model")]
    flo, fhi = (float(rng[0]), float(rng[1]))
    if use_fit and LV.at_least("Advanced"):
        fw = st.slider("Fit range (nm)", float(rng[0]), float(rng[1]), (float(rng[0]), float(rng[1])), 0.5, key="pl_fitwin")
        flo, fhi = fw
    # ---------- graph with brush
    marks = []
    dflt = _band_default([{"wl": q["x"], "I": q["y"]} for q in spec], rng[0], rng[1])
    b1d, b2d = dflt if dflt else ((rng[0], (rng[0] + rng[1]) / 2), ((rng[0] + rng[1]) / 2, rng[1]))
    b1d, b2d = tuple(round(float(v), 1) for v in b1d), tuple(round(float(v), 1) for v in b2d)
    b1 = ss.get("pl_b1", b1d)
    b2 = ss.get("pl_b2", b2d)
    if not use_fit or method == METHODS[2]:
        if b1:
            marks.append((b1[0], b1[1], "band 1", "#0072B2"))
        if b2:
            marks.append((b2[0], b2[1], "band 2", "#D55E00"))
    if use_fit:
        marks.append((flo, fhi, "fit range", "#10B981"))
    _brush_picker(spec, float(rng[0]), float(rng[1]), marks)
    bc = st.columns(4)
    bc[0].button("→ band 1", key="pl_pb1", on_click=_pick_cb, args=("pl_b1",), width="stretch", help="Use the range you dragged on the graph as band 1")
    bc[1].button("→ band 2", key="pl_pb2", on_click=_pick_cb, args=("pl_b2",), width="stretch")
    if use_fit and LV.at_least("Advanced"):
        bc[2].button("→ fit range", key="pl_pfw", on_click=_pick_cb, args=("pl_fitwin",), width="stretch")
    if ss.get("pl_pick_msg"):
        bc[3].caption(ss.pop("pl_pick_msg"))
    st.caption("Drag over the graph to select a range, then press a button. With the automatic method nothing has to be chosen.")
    # ---------- integration windows
    I1w = I2w = None
    if (not use_fit) or method == METHODS[2]:
        w1 = st.slider("Band 1 window (nm) – I₁", float(rng[0]), float(rng[1]), b1d, 0.5, key="pl_b1")
        w2 = st.slider("Band 2 window (nm) – I₂", float(rng[0]), float(rng[1]), b2d, 0.5, key="pl_b2")
        I1w = np.array([PL.band_integral(p["wl_c"], p["y_c"], w1[0], w1[1], "linear") for p in t_order])
        I2w = np.array([PL.band_integral(p["wl_c"], p["y_c"], w2[0], w2[1], "linear") for p in t_order])
        if (I2w <= 0).any() or (I1w <= 0).any():
            st.warning("A band has no (positive) integrated intensity in its window – adjust the windows.")
            if not use_fit:
                return None
            I1w = I2w = None
    if not use_fit:
        return I1w, I2w, I1w / I2w, None, None, {"method": "windows"}
    # ---------- global fit
    try:
        man = [float(v) for v in bands_txt.replace(";", ",").split(",") if v.strip()]
        ext = [float(v) for v in extra_txt.replace(";", ",").split(",") if v.strip()]
    except ValueError:
        st.warning("Band positions not understood – the automatic choice is used.")
        man, ext = [], []
    with st.spinner("Fitting all spectra …"):
        if len(man) >= 2:
            bands = [(c, 10.0) for c in man[:2]]
            auto = PR.auto_bands(np.asarray(spec[0]["x"], float), np.mean([np.interp(spec[0]["x"], q["x"], q["y"]) for q in spec], axis=0), flo, fhi, 2)
            if auto:
                bands = [min(auto, key=lambda b: abs(b[0] - c)) if abs(min(auto, key=lambda b: abs(b[0] - c))[0] - c) < 8 else (c, 10.0) for c in man[:2]]
                bands = [(c0_, w_) for (c0_, w_), c in zip(bands, man[:2])]
            fit = PR.auto_fit(spec, counts, flo, fhi, bgo, 3 if add_extra and not ext else 0, bands=bands)
            if fit is not None and ext:
                fit = PR.fit_global(spec, fit["bands"] + [(c, 10.0) for c in ext], (flo, fhi), bgo, model, counts)
                fit["bands"], fit["extra"] = fit.get("bands", bands), [(c, 10.0) for c in ext]
        else:
            fit = PR.auto_fit(spec, counts, flo, fhi, bgo, 3 if add_extra and not ext else 0)
            if fit is not None and ext:
                fit = PR.fit_global(spec, fit["bands"] + [(c, 10.0) for c in ext], (flo, fhi), bgo, model, counts)
        if fit is not None and model != "auto" and not ext:
            fit = PR.fit_global(spec, fit["bands"] + fit["extra"], (flo, fhi), bgo, model, counts) | {"bands": fit["bands"], "extra": fit["extra"]}
        if fit is None:
            st.info("The automatic search did not find two separate bands. Choose the two windows by hand (method: integration windows) "
                    "or give the band centres in the fit settings.")
            return None
        fit = PR.robust_ratio(spec, fit["bands"], bg_order=bgo, counts=counts, main=fit)
    return _show_fit(fit, spec, T, I1w, I2w, method, counts, g)


def _show_fit(fit, spec, T, I1w, I2w, method, counts, g):
    show, newfig, dl = g["show"], g["newfig"], g["dl"]
    rows = fit["rows"]
    A = np.array([r["areas"][0] for r in rows])
    B = np.array([r["areas"][1] for r in rows])
    D = np.array([r["ratio"] for r in rows])
    stat = np.array([r["ratio_err"] for r in rows])
    tot = fit["ratio_total"]
    bad = [r for r in rows if r["flag"] != "good"]
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Bands", f"{fit['n_bands']}" + (f" (incl. {len(fit['extra'])} extra)" if fit.get("extra") else ""),
              help="The first two are the ratio bands: " + ", ".join(f"{c:.1f} nm" for c, _ in fit["bands"][:2]))
    m2.metric("Model", fit["model"])
    m3.metric("Typical ratio error", f"{100 * float(np.nanmedian(tot / D)):.1f} %", help="Statistical error + the spread of refits with other windows and the other model.")
    m4.metric("Spectra to check", f"{len(bad)} of {len(rows)}", "all good" if not bad else "see the table", delta_color="off" if not bad else "inverse")
    if fit.get("extra"):
        st.caption("➕ Extra (nuisance) bands added automatically: " + ", ".join(f"{c:.0f} nm" for c, _ in fit["extra"])
                   + " – they absorb a shoulder / neighbouring emission so that it does not distort the ratio.")
    if bad:
        st.warning("The fit is not trustworthy for: " + ", ".join(f"{r['T']:.1f} K ({r['flag']})" for r in bad)
                   + " – the residual shows structure (dirty / wavy spectrum, an unmodelled band). Look at the fit of those spectra below.")
    # which temperature to look at
    ti = st.select_slider("Show the fit of the spectrum at", options=list(range(len(rows))), value=len(rows) // 2, key="pl_fit_i",
                          format_func=lambda i: f"{rows[i]['T']:.1f} K")
    r = rows[ti]

    def draw_fit(ax, r=r):
        ax.plot(r["x"], r["y"], "o", ms=2.2, mfc="none", mec="#222222", mew=0.4, label="measured")
        ax.plot(r["x"], r["fit"], "-", color="#D62728", lw=1.0, label="fit")
        cols = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]
        for k, comp in enumerate(r["components"]):
            ax.fill_between(r["x"], r["background"], r["background"] + comp, color=cols[k % 5], alpha=0.35, lw=0,
                            label=(f"band {k + 1}" if k < 2 else "extra"))
        ax.plot(r["x"], r["background"], "--", color="#7F7F7F", lw=0.6)
        res = r["y"] - r["fit"]
        off = float(r["y"].min()) - 0.12 * float(np.ptp(r["y"])) - float(np.max(np.abs(res)))
        ax.plot(r["x"], off + res, "-", color="#7F7F7F", lw=0.6, label="residual")
        ax.axhline(off, color="#7F7F7F", lw=0.3)
        ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Intensity (counts)" if counts else "Intensity")
        ax.text(0.02, 0.97, f"{r['T']:.1f} K,  R$^2$ = {r['r2']:.4f}", transform=ax.transAxes, va="top", fontsize=6)
        ax.legend(frameon=False, fontsize=6, loc="upper right")
    fig, ax = newfig(0.75); draw_fit(ax)
    show(fig, "PL_ratio_fit", pd.DataFrame({"wavelength": r["x"], "measured": r["y"], "fit": r["fit"], "background": r["background"]}), draw_fit,
         "Band fit at one temperature")
    q = pd.DataFrame({"T (K)": [x["T"] for x in rows], "I1 (area)": A, "I2 (area)": B, "ratio I1/I2": D, "± statistical": stat,
                      "± systematic": fit["ratio_sys"], "± total": tot, "R²": [x["r2"] for x in rows],
                      "χ² (red.)": [x["chi2_red"] for x in rows], "residual correlation (DW)": [x["dw"] for x in rows], "quality": [x["flag"] for x in rows]})
    with st.expander("Quality of every fit (ratio ± error, R², residual statistics)", expanded=bool(bad)):
        st.dataframe(q.style.format(precision=4), hide_index=True, width="stretch")
        st.caption("Quality: **good** = random residual; **check** = visible structure or a poor χ²; **poor** = the model does not describe this spectrum. "
                   "DW (Durbin–Watson) ≈ 2 for random residuals, < 1 for wavy ones. Systematic error = spread of refits with a narrower / wider fit range and the other peak model.")
    if method == METHODS[2] and I1w is not None:
        with st.expander("Comparison: fit vs integration windows", expanded=True):
            Dw = I1w / I2w
            dev = 100 * (Dw / D - 1)

            def draw_cmp(ax, T=T, D=D, tot=tot, Dw=Dw):
                ax.errorbar(T, D, yerr=tot, fmt="o", ms=3.5, capsize=2, lw=0.8, color="#0072B2", label="global fit")
                ax.plot(T, Dw, "s", ms=3.5, mfc="none", color="#D55E00", label="integration windows")
                ax.set_xlabel("Temperature (K)"); ax.set_ylabel("Δ = $I_1$/$I_2$"); ax.legend(frameon=False, fontsize=6)
            fig2, ax2 = newfig(0.7); draw_cmp(ax2)
            show(fig2, "PL_ratio_compare", pd.DataFrame({"T (K)": T, "fit": D, "fit error": tot, "windows": Dw}), draw_cmp, "Fit vs windows")
            st.caption(f"The windows deviate from the global fit by {np.nanmedian(np.abs(dev)):.1f} % (median, up to {np.nanmax(np.abs(dev)):.1f} %): "
                       "that is the bias of cutting overlapping bands at a valley.")
    return A, B, D, tot, None, {"method": "fit", "fit": fit}


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
        raw_counts = all(np.allclose(p["I"], np.round(p["I"])) for p in t_order) and norm.startswith("none")
        offs = st.checkbox("Add a constant C to the Boltzmann law (Δ = α·exp(−ΔE/k_BT) + C)", False, key="pl_off")
        got = _ratio_block(t_order, T, rng, norm, raw_counts, g)
        if got is None:
            return
        I1, I2, D, D_err, _, info = got
        counts = bool(st.session_state.get("pl_counts", raw_counts))
        if (D <= 0).any() or not np.isfinite(D).all():
            st.warning("The ratio is not positive for every spectrum – check the bands / windows."); return
        bf = PL.fit_boltzmann(T, D, offset=offs, sigma=D_err)
        if bf.get("status") != "ok":
            st.warning(f"Boltzmann fit: {bf.get('status')}"); return
        Sr = PL.relative_sensitivity(T, bf["dE_eV"])
        Sa = PL.absolute_sensitivity(T, bf)
        sT = PL.temperature_uncertainty_from_ratio(Sr, D, D_err) if D_err is not None else (PL.temperature_uncertainty(Sr, I1, I2) if counts else None)
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("ΔE", f"{bf['dE_cm1']:.0f} cm⁻¹", f"± {bf['dE_err_cm1']:.0f}" if np.isfinite(bf.get("dE_err_cm1", np.nan)) else None, delta_color="off",
                  help=f"{1000 * bf['dE_eV']:.1f} meV")
        m2.metric("Sᵣ at " + f"{T.min():.0f} K", f"{Sr[0]:.3f} % K⁻¹")
        m3.metric("σT" + (f" at {T.min():.0f} K" if sT is not None else ""), f"{sT[0]:.2f} K" if sT is not None else "–")
        m4.metric("R²", f"{bf['R2']:.4f}")
        Tf = np.linspace(T.min(), T.max(), 200)
        _cb = [c for c, _ in (info.get("fit") or {}).get("bands", [])[:2]] if info.get("method") == "fit" else []
        _lo, _hi = ((min(_cb) - 40, max(_cb) + 40) if _cb else (float(rng[0]), float(rng[1])))
        mp = [{"T": p["T"], "x": p["wl_c"][(p["wl_c"] >= _lo) & (p["wl_c"] <= _hi)], "y": p["y_c"][(p["wl_c"] >= _lo) & (p["wl_c"] <= _hi)]} for p in t_order]

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
        panels = [("Emission map", d_map), ("Δ vs T (Boltzmann fit)", d_D), ("Relative sensitivity Sr", d_Sr)] + ([("Temperature uncertainty σT", d_sT)] if sT is not None else [])
        fig, axs = C.compose([f for _, f in panels], ncols=2, width=max(g["WIDTH"], 6.6), panel_aspect=0.78, label_style="A", label_inside=False, wspace=0.75, hspace=0.5)
        res = pd.DataFrame({"T (K)": T, "I1": I1, "I2": I2, "Δ = I1/I2": D, "Δ fit": bf["f"](T), "Sr (% per K)": Sr, "Sa (per K)": Sa})
        if D_err is not None:
            res["Δ error"] = D_err
        if sT is not None:
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
                                 **({"sigmaT_Tmin_K": float(sT[0])} if sT is not None else {})}, "Photoluminescence page")
