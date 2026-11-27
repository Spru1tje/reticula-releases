"""v2.18.4: the thermometry tab of the Photoluminescence page - ratio (two-band fit, manifold fit, windows), quality report with a confidence verdict,
uncertainty budget, repeats / runs / hysteresis, laser-heating check, known-answer self-check and the analysis record.

Design: the confidence verdict sits at the TOP of the results (filled in last, from every check), so the user sees at once whether the numbers can be
trusted and why. Every check says what it means and what to do. Numbers are never hidden, but a 'low' confidence is printed next to them."""
from __future__ import annotations

import hashlib
import math

import numpy as np
import pandas as pd
import streamlit as st
from mfs import numentry as NE  # every slider takes a typed number too (2.19.6)

from mfs import levels as LV
from mfs import pl as PL


# ------------------------------------------------------------------------------------------------------------------------------ helpers
def clamp_range(key, lo, hi):
    """v2.18.5: a range slider keeps its value in the session; after loading data with another range (e.g. 500–580 nm, then 900–1250 nm) the stored value
    lies outside the new limits and Streamlit refuses it. The stored value is clipped to the new limits, or dropped when nothing of it is left."""
    ss = st.session_state
    v = ss.get(key)
    if v is None:
        return
    try:
        a, b = (float(v[0]), float(v[1])) if isinstance(v, (tuple, list)) else (float(v), float(v))
    except Exception:
        ss.pop(key, None)
        return
    a2, b2 = max(min(a, b), lo), min(max(a, b), hi)
    if not isinstance(v, (tuple, list)):
        if not (lo <= float(v) <= hi):
            ss.pop(key, None)
        return
    if b2 <= a2:
        ss.pop(key, None)
    elif (a2, b2) != (a, b):
        ss[key] = (a2, b2)


def memo(parts, fn):
    """Session-level cache of slow fits (keyed by the data and every setting)."""
    store = st.session_state.setdefault("_pl_memo", {})
    k = hashlib.md5(repr(parts).encode()).hexdigest()
    if k not in store:
        if len(store) > 24:
            store.clear()
        store[k] = fn()
    return store[k]


def data_hash(spec):
    h = hashlib.md5()
    for q in spec:
        h.update(np.asarray(q["x"], float).tobytes())
        h.update(np.asarray(q["y"], float).tobytes())
        h.update(str(q.get("T")).encode())
    return h.hexdigest()


def _cut(p, rng, key):
    m = (p["wl"] >= rng[0]) & (p["wl"] <= rng[1])
    return np.asarray(p[key], float)[m], m


# ------------------------------------------------------------------------------------------------------------------------------ the ratio
def ratio_block(t_order, T, rng, norm, counts_default, g, METHODS, brush_picker, pick_cb, band_default):
    """Returns (I1, I2, D, D_stat, info) or None. D_stat = statistical error of every ratio. info: method, fit (fit methods), windows, Dw, sDw,
    alt_sets [(label, D array)], spec, D_total."""
    from mfs import plqa as QA
    from mfs import plratio as PR
    ss = st.session_state
    c0, c1 = st.columns([3, 2])
    method = c0.radio("Method", METHODS, key="pl_method", horizontal=False,
                      help="Global fit: all spectra described by the same two bands (+ extra bands). Integration windows: the classic area in two windows. "
                           "Manifold fit: for crowded, Stark-split bands (lanthanides) – every line gets its own shape and the lines in two ranges are summed.")
    counts = c1.checkbox("Intensities are photon counts (statistics-based weights and σT)", counts_default, key="pl_counts",
                         help="Raw counts: points are weighted by Poisson statistics (of the RAW counts, background included). Arbitrary units: the noise is estimated from the data.")
    spec = []
    for p in t_order:
        y, m = _cut(p, rng, "I_corr" if "I_corr" in p else "I")
        spec.append({"T": p["T"], "x": p["wl"][m], "y": y, "var": np.asarray(p.get("I_raw", p["I"]), float)[m] if counts else None})
    use_fit = not method.startswith("Integration")
    manifold = method == METHODS[3]
    bg_alts = []
    plan_ = t_order[0].get("_plan") if t_order else None
    if plan_ is not None and plan_.get("method") not in (None, "none"):
        from mfs import plbackground as BG_
        for lab_, pl2_ in BG_.alternatives(t_order[0]["wl"], plan_, max_order=2):
            try:
                rr_ = BG_.apply_plan([{"wl": p["wl"], "I": p["I"]} for p in t_order], pl2_, plan_.get("_dark"), bool(ss.get("pl_spike", False)))
                bg_alts.append((f"background: {lab_}", [{"T": p["T"], "x": p["wl"][(p["wl"] >= rng[0]) & (p["wl"] <= rng[1])], "y": r_["y"][(p["wl"] >= rng[0]) & (p["wl"] <= rng[1])],
                                                         "var": s_["var"]} for p, r_, s_ in zip(t_order, rr_, spec)]))
            except Exception:
                continue
    bands_txt, extra_txt, bgo, model, add_extra = "", "", 1, "auto", True
    if use_fit and LV.at_least("Advanced"):
        with st.expander("⚙️ Fit settings (automatic by default)", expanded=False):
            a1, a2, a3 = st.columns(3)
            if not manifold:
                bands_txt = a1.text_input("Ratio bands: centres (nm), empty = automatic", "", key="pl_bands", help="e.g. 524, 545 – band 1 is the first, band 2 the second.")
                extra_txt = a2.text_input("Extra (nuisance) bands (nm), empty = automatic", "", key="pl_extra",
                                          help="Bands that are not part of the ratio but overlap with it (a shoulder, a third emission).")
                add_extra = a2.checkbox("Let the program add extra bands when needed", True, key="pl_addextra")
            else:
                a1.caption("Manifold fit: every line between the fit-range limits is found automatically (maxima and shoulders) and gets its own shape; "
                           "lines left in the residual are added while the information criterion improves.")
            bgo = NE.slider_in(a3, "Background polynomial order", 0, 3, 1, key="pl_bgo")
            model = {"automatic (BIC)": "auto", "shared shape": "none", "shape drifts with T": "drift"}[
                a3.selectbox("Peak model", ["automatic (BIC)", "shared shape", "shape drifts with T"], key="pl_model")]
    flo, fhi = (float(rng[0]), float(rng[1]))
    user_range = False
    if use_fit and LV.at_least("Advanced"):
        clamp_range("pl_fitwin", float(rng[0]), float(rng[1]))
        fw = NE.slider("Fit range (nm)", float(rng[0]), float(rng[1]), (float(rng[0]), float(rng[1])), 0.5, key="pl_fitwin")
        user_range = (abs(fw[0] - rng[0]) > 1e-6 or abs(fw[1] - rng[1]) > 1e-6)
        flo, fhi = fw
    # ---------- graph with brush
    dflt = band_default([{"wl": q["x"], "I": q["y"]} for q in spec], rng[0], rng[1])
    b1d, b2d = dflt if dflt else ((rng[0], (rng[0] + rng[1]) / 2), ((rng[0] + rng[1]) / 2, rng[1]))
    b1d, b2d = tuple(round(float(v), 1) for v in b1d), tuple(round(float(v), 1) for v in b2d)
    clamp_range("pl_b1", float(rng[0]), float(rng[1]))
    clamp_range("pl_b2", float(rng[0]), float(rng[1]))
    b1, b2 = ss.get("pl_b1", b1d), ss.get("pl_b2", b2d)
    show_windows = (not use_fit) or method == METHODS[2] or manifold
    marks = []
    if show_windows:
        marks += [(b1[0], b1[1], "range 1" if manifold else "band 1", "#0072B2"), (b2[0], b2[1], "range 2" if manifold else "band 2", "#D55E00")]
    if use_fit and (user_range or not manifold):
        marks.append((flo, fhi, "fit range", "#10B981"))
    brush_picker([{"x": q["x"], "y": q["y"]} for q in spec], float(rng[0]), float(rng[1]), marks)
    bc = st.columns(4)
    bc[0].button("→ range 1" if manifold else "→ band 1", key="pl_pb1", on_click=pick_cb, args=("pl_b1",), width="stretch", help="Use the range you dragged on the graph")
    bc[1].button("→ range 2" if manifold else "→ band 2", key="pl_pb2", on_click=pick_cb, args=("pl_b2",), width="stretch")
    if use_fit and LV.at_least("Advanced"):
        bc[2].button("→ fit range", key="pl_pfw", on_click=pick_cb, args=("pl_fitwin",), width="stretch")
    if ss.get("pl_pick_msg"):
        bc[3].caption(ss.pop("pl_pick_msg"))
    st.caption("Drag over the graph to select a range, then press a button." + (" Range 1 holds the lines of the upper (thermally populated) level, range 2 those of the lower one; "
               "the order does not matter – the program puts the ratio the right way round." if manifold else " With the automatic method nothing has to be chosen."))
    if show_windows:
        w1 = NE.slider(("Range 1 (nm) – lines summed into I₁" if manifold else "Band 1 window (nm) – I₁"), float(rng[0]), float(rng[1]), b1d, 0.5, key="pl_b1")
        w2 = NE.slider(("Range 2 (nm) – lines summed into I₂" if manifold else "Band 2 window (nm) – I₂"), float(rng[0]), float(rng[1]), b2d, 0.5, key="pl_b2")
    else:
        w1, w2 = tuple(b1), tuple(b2)
    # ---------- integration windows: always computed (method, cross-check, alternatives)
    # after the standard background subtraction the windows integrate above ZERO (a straight line through the window ends would cut into neighbouring lines);
    # without it, the classic straight line through the window ends is the local background
    wb = "none" if (plan_ is not None and plan_.get("method") not in (None, "none")) else "linear"
    I1w, I2w, Dw, sDw = QA.window_ratio(spec, w1, w2, wb)
    win_ok = bool(np.all(I1w > 0) and np.all(I2w > 0))
    win_alts = []
    for f_, lab_ in ((-0.15, "windows 15 % narrower"), (0.15, "windows 15 % wider")):
        ww1 = (w1[0] - f_ * (w1[1] - w1[0]) / 2, w1[1] + f_ * (w1[1] - w1[0]) / 2)
        ww2 = (w2[0] - f_ * (w2[1] - w2[0]) / 2, w2[1] + f_ * (w2[1] - w2[0]) / 2)
        win_alts.append((lab_, QA.window_ratio(spec, ww1, ww2, wb)[2]))
    for lab_, sp_ in bg_alts:
        win_alts.append((lab_, QA.window_ratio(sp_, w1, w2, wb)[2]))
    if not use_fit:
        if not win_ok:
            st.warning("A band has no (positive) integrated intensity in its window – adjust the windows.")
            return None
        return I1w, I2w, Dw, sDw, {"method": "windows", "win_baseline": wb, "windows": (w1, w2), "Dw": Dw, "sDw": sDw, "alt_sets": win_alts, "spec": spec, "D_total": sDw}
    # ---------- fits (cached: a manifold fit takes seconds to minutes)
    dh = data_hash(spec)
    if manifold:
        if not user_range:
            flo = max(float(rng[0]), min(w1[0], w2[0]) - 15.0)
            fhi = min(float(rng[1]), max(w1[1], w2[1]) + 15.0)

        def run():
            f = PR.fit_manifolds(spec, w1, w2, counts, flo, fhi, bgo)
            if f is None:
                return None
            if model != "auto":
                f = PR.with_model(spec, f, model, f["bands"], f["window"], bgo, counts, 120, f["groups"]) | {"bands": f["bands"], "extra": [], "groups": f["groups"], "mode": "manifold",
                                                                                                        "ranges": f["ranges"]}
            return PR.robust_ratio(spec, f["bands"], bg_order=bgo, counts=counts, main=f, bg_alternatives=bg_alts)
        with st.spinner("Fitting every line of the two manifolds in all spectra … (cached afterwards)"):
            fit = memo(("manifold", dh, tuple(w1), tuple(w2), flo, fhi, bgo, model, counts, len(bg_alts)), run)
        if fit is None:
            st.warning("No line was found inside range 1 or inside range 2. Drag the two ranges over the two groups of lines.")
            return None
    else:
        try:
            man = [float(v) for v in bands_txt.replace(";", ",").split(",") if v.strip()]
            ext = [float(v) for v in extra_txt.replace(";", ",").split(",") if v.strip()]
        except ValueError:
            st.warning("Band positions not understood – the automatic choice is used.")
            man, ext = [], []

        def run():
            if len(man) >= 2:
                bands = [(c, 10.0) for c in man[:2]]
                auto = PR.auto_bands(np.asarray(spec[0]["x"], float), np.mean([np.interp(spec[0]["x"], q["x"], q["y"]) for q in spec], axis=0), flo, fhi, 2)
                if auto:
                    bands = [(c, min(auto, key=lambda b: abs(b[0] - c))[1]) if abs(min(auto, key=lambda b: abs(b[0] - c))[0] - c) < 8 else (c, 10.0) for c in man[:2]]
                f = PR.auto_fit(spec, counts, flo, fhi, bgo, 3 if add_extra and not ext else 0, bands=bands)
            else:
                f = PR.auto_fit(spec, counts, flo, fhi, bgo, 3 if add_extra and not ext else 0)
            if f is None:
                return None
            if ext:
                bb = f["bands"]
                f = PR.with_model(spec, f, model, bb + [(c, 10.0) for c in ext], (flo, fhi), bgo, counts)
                f["bands"], f["extra"] = bb, [(c, 10.0) for c in ext]
            elif model != "auto":
                f = PR.with_model(spec, f, model, f["bands"] + f["extra"], (flo, fhi), bgo, counts) | {"bands": f["bands"], "extra": f["extra"]}
            return PR.robust_ratio(spec, f["bands"], bg_order=bgo, counts=counts, main=f, bg_alternatives=bg_alts)
        with st.spinner("Fitting all spectra … (cached afterwards)"):
            fit = memo(("fit", dh, flo, fhi, bgo, model, add_extra, bands_txt, extra_txt, counts, len(bg_alts)), run)
        if fit is None:
            st.info("The automatic search did not find two separate bands. Choose the two windows by hand (method: integration windows), "
                    "give the band centres in the fit settings, or use the manifold fit for crowded bands.")
            return None
        share = PR.flag_share(fit)
        if share > 0.3:
            st.warning(f"⚠️ {round(share * len(fit['rows']))} of {len(fit['rows'])} spectra are not described well by two bands (wavy residual). This is typical for a crowded "
                       "band system, for example the Stark-split lines of a lanthanide: the extra bands then absorb part of the ratio bands, and the ratio is not reliable.")
            st.button("🧩 Switch to the manifold fit (every line its own shape; choose two ranges)", key="pl_to_manifold",
                      on_click=lambda: st.session_state.update({"pl_method": METHODS[3]}))
    res = show_fit(fit, spec, T, Dw, sDw, (w1, w2), win_alts, method, counts, g, METHODS)
    window_fit_panel(spec, w1, w2, counts, res[2], res[3], g)
    return res


def window_fit_table(spec, w1, w2, counts, D, D_stat, weighting="sqrt"):
    """Experimental window fit (mfs/plwindows.py): per temperature the ratio from the two chosen windows, using only the data inside them,
    compared with the global-fit ratio D (error D_stat). Pure function (no Streamlit), so it can be tested on synthetic data."""
    from mfs import plwindows as PW
    rows = PW.window_ratio_series(spec, tuple(w1), tuple(w2), counts=counts, weighting=weighting)
    T_arr = [float(q["T"]) for q in spec]
    cmp_ = PW.comparison_table(rows, T_arr, D, D_stat)          # S-Figuryte-3: the comparison is computed in plwindows only
    df = pd.DataFrame({
        "T (K)": cmp_["T"],
        "window ratio": cmp_["window ratio"],
        "window error": cmp_["window error"],
        "global ratio": cmp_["global ratio"],
        "global error": cmp_["global error"],
        "difference": cmp_["difference"],
        "disagree (> 3 σ)": cmp_["disagree"],
        "flags": cmp_["flags"],
    })
    return {"df": df, "n_compared": cmp_["n_compared"], "n_disagree": cmp_["n_disagree"], "rows": rows}


def window_fit_panel(spec, w1, w2, counts, D, D_stat, g):
    """Expander under the fit results: the window fit next to the global fit. Experimental; the flags say why a value is missing."""
    show, newfig = g["show"], g["newfig"]
    with st.expander("🪟 Window fit (experimental): ratio from the two chosen windows only", expanded=False):
        st.caption("Each window is fitted on its own data: one pseudo-Voigt band with a fixed Lorentzian share (eta = 0.5), a straight background, "
                   "and the neighbour band inside the joint fit. The error is statistical only: a shape mismatch of the real bands and low counts "
                   "(below 100 counts per point) are NOT in it. A disagreement with the global fit above 3 combined errors is a warning.")
        from mfs import plwindows as PW
        # the user chooses the weighting of the window fit; sqrt is the default, Poisson is opt-in (plwindows: weighting="poisson")
        wchoice = st.radio("Weighting of the window fit", ["sqrt (default)", "Poisson (opt-in)"], horizontal=True, key="pl_winweight",
                           help="sqrt: the default weights. Poisson: deviance weighting of the raw counts; it removes most of the low-count bias. "
                                "Changes the window-fit numbers in this panel only.")
        weighting = "poisson" if wchoice.startswith("Poisson") else "sqrt"
        for msg in PW.input_problems(spec, w1, w2)[:5]:           # S-Figuryte-3: say why a value may be missing or unreliable
            st.warning(msg)
        try:
            t = window_fit_table(spec, w1, w2, counts, D, D_stat, weighting)
        except ValueError as exc:
            st.warning(f"The window fit cannot be used: {exc}")
            return
        df = t["df"]
        st.dataframe(df, hide_index=True, width="stretch")
        st.caption(f"{t['n_compared']} temperatures compared with the global fit; {t['n_disagree']} disagree by more than 3 combined errors.")

        def draw(ax):
            ax.errorbar(df["T (K)"], df["global ratio"], yerr=df["global error"], fmt="o", ms=3, color="#0072B2", label="global fit")
            ax.errorbar(df["T (K)"], df["window ratio"], yerr=df["window error"], fmt="s", ms=3, mfc="none", color="#D55E00", label="window fit")
            ax.set_xlabel("Temperature (K)")
            ax.set_ylabel("Δ = $I_1$/$I_2$")
            ax.legend(frameon=False, fontsize=6)
        fig, ax = newfig(0.7)
        draw(ax)
        show(fig, "PL_window_fit", df, draw, "Window fit vs global fit")


def show_fit(fit, spec, T, Dw, sDw, windows, win_alts, method, counts, g, METHODS):
    show, newfig = g["show"], g["newfig"]
    rows = fit["rows"]
    A = np.array([r["areas"][0] for r in rows])
    B = np.array([r["areas"][1] for r in rows])
    D = np.array([r["ratio"] for r in rows])
    stat = np.array([r["ratio_err"] for r in rows])
    tot = fit["ratio_total"]
    bad = [r for r in rows if r["flag"] != "good"]
    groups = fit.get("groups")
    m1, m2, m3, m4 = st.columns(4)
    if groups:
        m1.metric("Lines", f"{len(fit['bands'])}", f"{len(groups[0])} in range 1, {len(groups[1])} in range 2", delta_color="off")
    else:
        m1.metric("Bands", f"{fit['n_bands']}" + (f" (incl. {len(fit['extra'])} extra)" if fit.get("extra") else ""),
                  help="The first two are the ratio bands: " + ", ".join(f"{c:.1f} nm" for c, _ in fit["bands"][:2]))
    m2.metric("Model", fit["model"])
    m3.metric("Typical ratio error", f"{100 * float(np.nanmedian(tot / D)):.1f} %", help="Statistical error + the spread of refits with other fit ranges, the other model and the other reasonable background subtractions.")
    m4.metric("Spectra to check", f"{len(bad)} of {len(rows)}", "all good" if not bad else "see the table", delta_color="off" if not bad else "inverse")
    if fit.get("extra"):
        st.caption("➕ Extra (nuisance) bands added automatically: " + ", ".join(f"{c:.0f} nm" for c, _ in fit["extra"])
                   + " – they absorb a shoulder / neighbouring emission so that it does not distort the ratio.")
    if bad:
        st.warning("The fit is not trustworthy for: " + ", ".join(f"{r['T']:.1f} K ({r['flag']})" for r in bad)
                   + " – the residual shows structure (dirty / wavy spectrum, an unmodelled band). Look at the fit of those spectra below.")
    ti = NE.select_slider("Show the fit of the spectrum at", options=list(range(len(rows))), value=len(rows) // 2, key="pl_fit_i",
                          format_func=lambda i: f"{rows[i]['T']:.1f} K")
    r = rows[ti]

    def draw_fit(ax, r=r, groups=groups):
        ax.plot(r["x"], r["y"], "o", ms=2.2, mfc="none", mec="#222222", mew=0.4, label="measured")
        ax.plot(r["x"], r["fit"], "-", color="#D62728", lw=1.0, label="fit")
        cols = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]
        seen = set()
        for k, comp in enumerate(r["components"]):
            if groups:
                kind = "I₁ lines" if k in groups[0] else ("I₂ lines" if k in groups[1] else "other lines")
                col = {"I₁ lines": "#0072B2", "I₂ lines": "#D55E00", "other lines": "#9E9E9E"}[kind]
            else:
                kind = f"band {k + 1}" if k < 2 else "extra"
                col = cols[k % 5]
            ax.fill_between(r["x"], r["background"], r["background"] + comp, color=col, alpha=0.35, lw=0, label=None if kind in seen else kind)
            seen.add(kind)
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
                      "± systematic": fit["ratio_sys"], **({"  of which background choice": fit["ratio_sys_bg"]} if fit.get("bg_alt_labels") else {}), "± total": tot, "R²": [x["r2"] for x in rows],
                      "χ² (red.)": [x["chi2_red"] for x in rows], "residual correlation (DW)": [x["dw"] for x in rows], "quality": [x["flag"] for x in rows]})
    with st.expander("Quality of every fit (ratio ± error, R², residual statistics)", expanded=bool(bad)):
        st.dataframe(q.style.format(precision=4), hide_index=True, width="stretch")
        st.caption("Quality: **good** = random residual; **check** = visible structure or a poor χ²; **poor** = the model does not describe this spectrum. "
                   "DW (Durbin–Watson) ≈ 2 for random residuals, < 1 for wavy ones. Systematic error = spread of refits with a narrower / wider fit range, the other peak model"
                   + (" and the other background subtractions." if fit.get("bg_alt_labels") else "."))
    if method == METHODS[2]:
        with st.expander("Comparison: fit vs integration windows", expanded=True):
            dev = 100 * (Dw / D - 1)

            def draw_cmp(ax, T=T, D=D, tot=tot, Dw=Dw):
                ax.errorbar(T, D, yerr=tot, fmt="o", ms=3.5, capsize=2, lw=0.8, color="#0072B2", label="global fit")
                ax.plot(T, Dw, "s", ms=3.5, mfc="none", color="#D55E00", label="integration windows")
                ax.set_xlabel("Temperature (K)"); ax.set_ylabel("Δ = $I_1$/$I_2$"); ax.legend(frameon=False, fontsize=6)
            fig2, ax2 = newfig(0.7); draw_cmp(ax2)
            show(fig2, "PL_ratio_compare", pd.DataFrame({"T (K)": T, "fit": D, "fit error": tot, "windows": Dw}), draw_cmp, "Fit vs windows")
            st.caption(f"The windows deviate from the global fit by {np.nanmedian(np.abs(dev)):.1f} % (median, up to {np.nanmax(np.abs(dev)):.1f} %): "
                       "the bias of cutting overlapping bands at a valley. What matters for a thermometer is whether both give the same energy gap – see the quality report.")
    alts = []
    if fit.get("fit_alt_ratios") is not None:
        alts += list(zip(fit.get("fit_alt_labels") or [], list(fit["fit_alt_ratios"])))
    if fit.get("bg_alt_ratios") is not None:
        alts += list(zip(fit.get("bg_alt_labels") or [], list(fit["bg_alt_ratios"])))
    return A, B, D, stat, {"method": "manifold" if groups else "fit", "fit": fit, "windows": windows, "Dw": Dw, "sDw": sDw, "alt_sets": alts, "spec": spec, "D_total": tot,
                           "win_alts": win_alts}


# ------------------------------------------------------------------------------------------------------------------------------ results
def _centroid_cm1(spec, lo, hi):
    """Intensity-weighted mean wavenumber of the mean spectrum inside lo..hi (nm)."""
    x = np.asarray(spec[0]["x"], float)
    y = np.mean([np.interp(x, q["x"], q["y"]) for q in spec], axis=0)
    m = (x >= min(lo, hi)) & (x <= max(lo, hi)) & (y > 0)
    if m.sum() < 3:
        return float("nan")
    return float(np.sum(y[m] * 1e7 / x[m]) / np.sum(y[m]))


def _mle_window_block(t_order, info, swapped, counts):
    """v2.18.6: Poisson maximum-likelihood ratio fit on the RAW counts of the two integration windows (background counts from the
    background plan). Returns the fit dict of mfs/plmle.fit, or None when it does not apply (not counts, not windows, too few points)."""
    from mfs import plmle as PM
    wins = info.get("windows")
    if not counts or not wins or len(wins) != 2 or any("I_raw" not in p for p in t_order):
        return None
    (a1, b1), (a2, b2) = wins
    N1, N2, B1, B2, TT = [], [], [], [], []
    for p in t_order:
        wl = np.asarray(p["wl"], float)
        raw = np.asarray(p["I_raw"], float)
        base = np.asarray(p.get("baseline", np.zeros_like(raw)), float)
        m1 = (wl >= min(a1, b1)) & (wl <= max(a1, b1))
        m2 = (wl >= min(a2, b2)) & (wl <= max(a2, b2))
        N1.append(raw[m1].sum()); B1.append(base[m1].sum())
        N2.append(raw[m2].sum()); B2.append(base[m2].sum())
        TT.append(float(p["T"]))
    N1, N2, B1, B2, TT = (np.asarray(v, float) for v in (N1, N2, B1, B2, TT))
    if swapped:                      # the rising ratio is band 2 / band 1 in the model
        N1, N2, B1, B2 = N2, N1, B2, B1
    ok = np.isfinite(N1) & np.isfinite(N2) & (N1 > 0) & (N2 > 0)
    if ok.sum() < 4:
        return None
    r = PM.fit(TT[ok], N1[ok], N2[ok], B1[ok], B2[ok])
    r["n"] = int(ok.sum())
    return r


def results(T, I1, I2, D, D_stat, info, t_order, offs, counts, sample, g, version="", settings=None):
    """Boltzmann calibration, quality report (verdict on top), uncertainty budget, experiment checks and downloads."""
    from mfs import plqa as QA
    show, newfig, dl, C, safe_name, store_put = g["show"], g["newfig"], g["dl"], g["C"], g["safe_name"], g["store_put"]
    ss = st.session_state
    top = st.container()
    checks = []
    alt_sets = list(info.get("alt_sets") or [])
    Dw, sDw = info.get("Dw"), info.get("sDw")
    # ---- orientation: Δ must RISE with temperature (band 1 = upper level); otherwise the bands are swapped
    swapped = False
    if len(T) >= 3 and np.polyfit(1.0 / T, np.log(D), 1)[0] > 0:
        swapped = True
        D_stat = D_stat / D ** 2
        D = 1.0 / D
        I1, I2 = I2, I1
        alt_sets = [(lab, 1.0 / np.asarray(a, float)) for lab, a in alt_sets]
        if Dw is not None:
            sDw, Dw = sDw / Dw ** 2, 1.0 / Dw
    bf = PL.fit_boltzmann(T, D, offset=offs, sigma=D_stat)
    if bf.get("status") != "ok":
        st.warning(f"Boltzmann fit: {bf.get('status')}")
        return
    bf["_swapped"] = swapped
    Sr = PL.relative_sensitivity(T, bf["dE_eV"])
    Sa = PL.absolute_sensitivity(T, bf)
    # ---- inputs of the budget
    with st.expander("🌡️ Temperature sensor and correction uncertainties (for the uncertainty budget)", expanded=False):
        s1, s2, s3 = st.columns(3)
        sT_rand = s1.number_input("Random error of each temperature reading (K)", 0.0, 20.0, 0.0, 0.05, key="pl_sT_rand",
                                  help="Scatter of the sensor reading and of the temperature stability at each step.")
        sT_sys = s2.number_input("Sample-vs-sensor offset / calibration (K)", 0.0, 50.0, 0.0, 0.1, key="pl_sT_sys",
                                 help="A common error of all temperatures: sensor calibration, sample not at the sensor temperature, laser heating.")
        resp_rel = s3.number_input("Uncertainty of the spectral-response correction (%)", 0.0, 50.0, 0.0, 0.5, key="pl_resp_u") / 100
        st.caption("General practice: a calibrated thermocouple or Pt100 next to the sample has about 0.1–0.5 K random error; the sample can differ from the sensor by "
                   "several kelvin (stage, cuvette, laser heating). Zero means 'not counted'.")
    runs = QA.between_runs(T, D, [p.get("file", "") for p in t_order], offs)
    B = QA.budget(T, D, D_stat, bf, sT_rand, sT_sys, alt_sets, runs["sd_cm1"], runs["n_runs"], offs, resp_rel)
    u_dE = B["u_dE_cm1"]
    T_pred = QA.predict_T(D, bf["alpha"], bf["dE_eV"], bf.get("C", 0.0))
    # ---- checks
    sat = QA.saturation(t_order)
    bad_sat = [s for s in sat if s["status"] != "ok"]
    if bad_sat:
        worst = "fail" if any(s["status"] == "fail" for s in bad_sat) else "warn"
        checks.append(QA.check("Detector saturation", worst, f"{len(bad_sat)} spectrum(s): " + "; ".join(f"{s['name']}: {s['reason']}" for s in bad_sat[:3]),
                               "A saturated peak is too low, so its band's share is too small and the ratio is wrong. Lower the slit, the power or the integration time."))
    else:
        checks.append(QA.check("Detector saturation", "ok", "No clipped or flat-topped peak found."))
    _, differ = QA.settings_table(t_order)
    if differ:
        checks.append(QA.check("Same settings", "warn", "The settings differ inside the series: " + ", ".join(differ) + ".", "Compare only spectra measured with the same settings."))
    n_bad = sum(1 for r in (info.get("fit") or {}).get("rows", []) if r["flag"] != "good") if info.get("fit") else 0
    # spectroscopic gap: the two band centres (two-band fit) or the centroids of the two ranges / windows
    spec = info.get("spec") or []
    if info.get("method") == "fit":
        cs = [c for c, _ in info["fit"]["bands"][:2]]
        gap_spec = PL.EV_TO_CM * abs(PL.HC / cs[0] - PL.HC / cs[1]) if len(cs) == 2 else None
    elif spec and info.get("windows"):
        (a1, b1_), (a2, b2_) = info["windows"]
        gap_spec = abs(_centroid_cm1(spec, a1, b1_) - _centroid_cm1(spec, a2, b2_))
    else:
        gap_spec = None
    checks += QA.boltzmann_checks(T, D, D_stat, bf, gap_spec if gap_spec and np.isfinite(gap_spec) else None, n_bad, len(T) if info.get("fit") else None)
    # v2.20.4: model choice (AIC), spectral response correction, and validation against an independent thermometer (optional)
    checks.append(QA.model_check(T, D, D_stat))
    checks.append(QA.correction_check(bool(ss.get("pl_cal_applied"))))
    with st.expander("🎯 Validation against an independent thermometer (optional, for a paper)", expanded=False):
        ref_txt = st.text_area("Reference temperatures (K), one per spectrum, in the order of the table", key="pl_ref_T", height=80,
                               help="Read them from a calibrated thermocouple, Pt100 or a certified reference measured at the same time as each spectrum. "
                                    "Leave empty if you have none; the report then says so.")
        refs = []
        for tok in ref_txt.replace(";", ",").replace(chr(10), ",").split(","):
            tok = tok.strip()
            if tok:
                try:
                    refs.append(float(tok))
                except ValueError:
                    refs = None
                    break
        if ref_txt.strip() and refs is None:
            st.warning("The reference list contains a value that is not a number.")
    checks.append(QA.validation_check(T_pred, refs if refs else None))      # the ratio thermometer's readout vs the reference
    # v2.18.6: Poisson maximum-likelihood fit on the raw counts (integration windows), compared with the subtraction result
    mle = _mle_window_block(t_order, info, swapped, counts)
    if mle is not None:
        if mle.get("converged") and np.isfinite(mle.get("sigma_cm1", np.nan)):
            diff_m = abs(abs(mle["dE_cm1"]) - abs(bf["dE_cm1"]))
            comb_m = math.sqrt(u_dE ** 2 + mle["sigma_cm1"] ** 2)
            stt_m = "ok" if diff_m <= 2 * comb_m else "warn"
            checks.append(QA.check("Poisson maximum likelihood (counts)", stt_m,
                                   f"The maximum-likelihood fit on the raw counts gives ΔE = {abs(mle['dE_cm1']):.0f} ± {mle['sigma_cm1']:.0f} cm⁻¹ "
                                   f"({mle['n']} temperatures); the subtraction method gives {abs(bf['dE_cm1']):.0f} cm⁻¹ "
                                   f"(difference {diff_m:.0f} cm⁻¹, {diff_m / comb_m if comb_m else float('inf'):.1f}× the combined uncertainty).",
                                   "" if stt_m == "ok" else "The background subtraction biases the ratio: check the windows and the background plan. "
                                                            "The maximum-likelihood value keeps every point and uses the Poisson error of the raw counts."))
        else:
            checks.append(QA.check("Poisson maximum likelihood (counts)", "warn",
                                   "The maximum-likelihood fit did not converge to a reliable optimum.",
                                   "Check the windows; a window without signal above the background makes the fit unidentifiable."))
        if mle.get("converged"):
            # v2.18.6: editable figure like the other PL figures (labels and legend movable, data as CSV)
            def draw_mle(ax_):
                Tf_ = np.linspace(float(np.min(T)), float(np.max(T)), 200)
                ax_.errorbar(1000.0 / np.asarray(T, float), np.log(np.asarray(D, float)), yerr=np.asarray(D_stat, float) / np.asarray(D, float),
                             fmt="o", ms=3.5, color="#1F4E79", ecolor="#7F9FBF", capsize=2, label="ratio (subtraction)")
                ax_.plot(1000.0 / Tf_, np.log(bf["f"](Tf_)), color="#7F8C8D", lw=1.2, label="subtraction fit")
                ax_.plot(1000.0 / Tf_, np.log(mle["alpha"] * np.exp(-mle["dE_cm1"] / PL.EV_TO_CM / (PL.KB_EV * Tf_))), color="#0F766E", lw=1.4,
                         label="Poisson maximum likelihood")
                ax_.set_xlabel("1000 / T (K⁻¹)"); ax_.set_ylabel("ln Δ"); ax_.legend(fontsize=7, frameon=False)
            fig_m, ax_m = newfig(0.75); draw_mle(ax_m)
            show(fig_m, "PL_mle_check", pd.DataFrame({"T (K)": np.asarray(T, float), "ratio": np.asarray(D, float), "ratio error": np.asarray(D_stat, float),
                                                      "subtraction fit": bf["f"](np.asarray(T, float)),
                                                      "Poisson ML fit": mle["alpha"] * np.exp(-mle["dE_cm1"] / PL.EV_TO_CM / (PL.KB_EV * np.asarray(T, float)))}),
                 draw_mle, "Poisson maximum likelihood vs subtraction")
    if info.get("method") in ("fit", "manifold") and Dw is not None and np.all(Dw > 0) and np.all(np.isfinite(Dw)):
        bw = PL.fit_boltzmann(T, Dw, offset=offs, sigma=sDw if np.all(sDw > 0) else None)
        if bw.get("status") == "ok":
            diff = abs(abs(bw["dE_cm1"]) - abs(bf["dE_cm1"]))
            comb = math.sqrt(u_dE ** 2 + (bw.get("dE_err_cm1") or 0) ** 2)
            stt = "ok" if diff <= 2 * comb else "warn"
            checks.append(QA.check("Independent method", stt, f"The integration windows give ΔE = {abs(bw['dE_cm1']):.0f} cm⁻¹, the fit {abs(bf['dE_cm1']):.0f} cm⁻¹ "
                                   f"(difference {diff:.0f} cm⁻¹, {diff / comb if comb else float('inf'):.1f}× the combined uncertainty).",
                                   "" if stt == "ok" else "The two ways of measuring the ratio disagree: the overlap of the bands, or what the extra lines absorb, decides the result. "
                                                         "Check the windows / ranges and the fit of the worst spectra."))
    if abs(bf["dE_cm1"]) > 0:
        rel = u_dE / abs(bf["dE_cm1"])
        checks.append(QA.check("Uncertainty of ΔE", "ok" if rel < 0.1 else ("warn" if rel < 0.3 else "fail"),
                               f"ΔE is known to ± {u_dE:.0f} cm⁻¹ ({100 * rel:.0f} %) with every component of the budget included."))
    if sT_rand == 0 and sT_sys == 0:
        checks.append(QA.check("Temperature sensor", "info", "No sensor uncertainty was given: the calibration uncertainty only contains the spectra. Fill it in above."))
    if not ss.get("pl_cal_applied"):
        checks.append(QA.check("Spectral response", "info", "No spectral-response correction. ΔE and Sᵣ do not depend on it (a response that does not change during the series multiplies every ratio by "
                               "the same factor), but α and the absolute ratio are specific to this instrument: a calibration curve cannot be moved to another spectrometer."))
    order = [p.get("_order", i) for i, p in enumerate(t_order)]
    hy = QA.hysteresis(T, D, order, bf["dE_eV"], D_stat)
    if hy.get("status") == "ok":
        lim = max(1.0, 2.0 * float(np.nanmedian(B["T_res"])))
        checks.append(QA.check("Heating vs cooling", "ok" if hy["max_dT"] <= lim else "warn",
                               f"Heating and cooling branch differ by up to {hy['max_dT']:.1f} K in apparent temperature.",
                               "" if hy["max_dT"] <= lim else "Hysteresis: the sample was not in equilibrium (wait longer at each step), or it changed during the run (drying, degradation)."))
    if runs["n_runs"] >= 2:
        checks.append(QA.check("Reproducibility", "ok" if runs["sd_cm1"] <= max(2 * B["rows"][0]["u_dE_cm1"], 0.05 * abs(bf["dE_cm1"])) else "warn",
                               f"{runs['n_runs']} independent runs: ΔE scatters by {runs['sd_cm1']:.0f} cm⁻¹ between them."))
    sc = ss.get("pl_selfcheck")
    if sc and sc[0] == _fit_key(info, bf):
        r_ = sc[1]
        if r_.get("status") == "ok":
            honest = 0.5 <= (r_["z_sd"] if r_["z_sd"] == r_["z_sd"] else 1) <= 2.0 and r_["coverage"] >= 0.6
            unbiased = abs(r_["bias_pct"]) < max(3.0, 100 * 2 * np.std([x["dE_cm1"] for x in r_["runs"]]) / r_["truth_cm1"] / math.sqrt(len(r_["runs"])))
            checks.append(QA.check("Simulation self-check (consistency only)", "ok" if honest and unbiased else "warn",
                                   f"On {len(r_['runs'])} simulated series with ΔE = {r_['truth_cm1']:.0f} cm⁻¹ the analysis returned {np.mean([x['dE_cm1'] for x in r_['runs']]):.0f} cm⁻¹ "
                                   f"(bias {r_['bias_pct']:+.1f} %); {100 * r_['coverage']:.0f} % of the runs within 2 error bars.",
                                   "" if honest and unbiased else "The analysis is biased or its error bars are too small for data like yours: trust the budget, not the fit error."))
    # ---- verdict on top
    v = QA.verdict(checks)
    with top:
        n_w = sum(c["status"] == "warn" for c in checks)
        n_f = sum(c["status"] == "fail" for c in checks)
        msg = {"high": "✅ **Confidence: HIGH** – every check passed.",
               "moderate": f"⚠️ **Confidence: MODERATE** – {n_w} warning(s). The numbers are usable with the caveats below.",
               "low": f"⛔ **Confidence: LOW** – {n_f} check(s) failed. Do not report these numbers before the problem is solved."}[v]
        {"high": st.success, "moderate": st.warning, "low": st.error}[v](msg)
        with st.expander("🔍 Quality report (every check, why it matters, what to do)", expanded=v != "high"):
            for c in sorted(checks, key=lambda c: {"fail": 0, "warn": 1, "info": 2, "ok": 3}[c["status"]]):
                st.markdown(f"{QA.ICON[c['status']]} **{c['name']}** – {c['text']}")
                if c.get("detail"):
                    st.caption("　" + c["detail"])
    # ---- metrics
    flag = None if v != "low" else "⛔ low confidence"
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("ΔE", f"{abs(bf['dE_cm1']):.0f} ± {u_dE:.0f} cm⁻¹", flag or f"{1000 * abs(bf['dE_eV']):.1f} meV", delta_color="inverse" if flag else "off",
              help="± = the total of the uncertainty budget (statistics, processing alternatives, temperature sensor, runs).")
    m2.metric("Sᵣ at " + f"{T.min():.0f} K", f"{Sr[0]:.3f} ± {Sr[0] * B['u_Sr_rel']:.3f} % K⁻¹")
    m3.metric(f"σT resolution at {T.min():.0f} K", f"{B['T_res'][0]:.2f} K", help="The precision of ONE reading: the temperature change the noise of one ratio corresponds to.")
    m4.metric(f"Calibration uncertainty at {T.min():.0f} K", f"{B['T_cal'][0]:.2f} K", help="The accuracy of the calibration curve T(Δ) at this temperature (Monte Carlo over ratios and sensor).")
    if swapped:
        st.caption("↔️ The ratio was inverted (Δ = I₂/I₁) because I₁/I₂ fell with temperature: the thermally populated level is band / range 2.")
    Tf = np.linspace(T.min(), T.max(), 200)
    mp = []
    lo_, hi_ = _map_range(info)
    for p in t_order:
        m = (p["wl_c"] >= lo_) & (p["wl_c"] <= hi_)
        mp.append({"T": p["T"], "x": p["wl_c"][m], "y": p["y_c"][m]})

    def d_map(ax, mp=mp):
        PL.draw_map(ax, mp, "Wavelength (nm)", "magma", "Intensity")

    def d_D(ax, T=T, D=D, bf=bf, Tf=Tf, e=D_stat):
        ax.errorbar(T, D, yerr=e, fmt="none", ecolor="#0072B2", lw=0.6, capsize=1.5)
        PL.draw_trend(ax, T, D, "Δ" + (" = $I_2$/$I_1$" if swapped else " = $I_1$/$I_2$"), fit_T=Tf, fit_y=bf["f"](Tf),
                      fit_label=f"ΔE = {abs(bf['dE_cm1']):.0f} ± {u_dE:.0f} cm$^{{-1}}$")

    def d_Sr(ax, T=T, Sr=Sr, Tf=Tf, bf=bf, u=B["u_Sr_rel"]):
        s = PL.relative_sensitivity(Tf, bf["dE_eV"])
        ax.fill_between(Tf, s * (1 - u), s * (1 + u), color="#0072B2", alpha=0.15, lw=0)
        PL.draw_trend(ax, T, Sr, r"$S_\mathrm{r}$ (% K$^{-1}$)", fit_T=Tf, fit_y=s, legend=False)

    def d_sT(ax, T=T, B=B):
        ax.plot(T, B["T_res"], "o-", ms=3.2, lw=0.9, color="#0072B2", label="resolution (one reading)")
        ax.plot(T, B["T_cal"], "s--", ms=3.0, lw=0.9, color="#D55E00", label="calibration (curve)")
        ax.set_xlabel("Temperature (K)"); ax.set_ylabel(r"$\sigma_T$ (K)"); ax.legend(frameon=False, fontsize=6)
    panels = [("Emission map", d_map), ("Δ vs T (Boltzmann fit)", d_D), ("Relative sensitivity Sr", d_Sr), ("Temperature uncertainty σT", d_sT)]
    fig, axs = C.compose([f for _, f in panels], ncols=2, width=max(g["WIDTH"], 6.6), panel_aspect=0.78, label_style="A", label_inside=False, wspace=0.75, hspace=0.5)
    res = pd.DataFrame({"T (K)": T, "I1": I1, "I2": I2, "Δ": D, "Δ statistical error": D_stat, "Δ fit": bf["f"](T), "T from Δ (K)": T_pred, "T from Δ − T set (K)": T_pred - T,
                        "Sr (% per K)": Sr, "Sa (per K)": Sa, "σT resolution (K)": B["T_res"], "σT calibration (K)": B["T_cal"]})
    show(fig, "PL_thermometry", res, None, "PL thermometry")
    pick = st.selectbox("Single panel (e.g. for the combined figure)", ["(none)"] + [n for n, _ in panels], key="pl_single2")
    if pick != "(none)":
        fn = dict(panels)[pick]
        fig2, ax2 = newfig(0.8); fn(ax2)
        show(fig2, "PL_" + safe_name(pick), res, fn, pick)
    st.dataframe(res.round(5), hide_index=True, width="stretch")
    with st.expander("📊 Uncertainty budget of ΔE (what the ± is made of)", expanded=False):
        bt = pd.DataFrame(B["rows"]).rename(columns={"component": "component", "u_dE_cm1": "u(ΔE) (cm⁻¹)", "u_alpha_rel_pct": "u(α) (%)", "note": "how"})
        st.dataframe(bt.style.format(precision=2, na_rep="–"), hide_index=True, width="stretch")
        st.caption(f"Total (square root of the sum of squares): u(ΔE) = {u_dE:.1f} cm⁻¹, u(α) = {B['u_alpha_rel_pct']:.1f} %. The processing rows refit the Boltzmann law on the ratios "
                   "obtained with every reasonable alternative (fit range, peak model, background, windows) – the spread is what the CHOICES cost. "
                   "Resolution = precision of one reading (noise of one ratio / Sᵣ); calibration = accuracy of the curve T(Δ).")
    with st.expander("🔁 Repeats, independent runs and heating / cooling", expanded=False):
        from mfs import plextras as _PX
        rep = _PX.repeat_uncertainty(T, D, lambda t_: PL.relative_sensitivity(t_, bf["dE_eV"]))
        if rep:
            st.markdown("**Repeats at one temperature**")
            st.dataframe(pd.DataFrame([{"T (K)": r_["T"], "repeats": r_["n"], "mean Δ": r_["mean"], "std Δ": r_["std"], "std (%)": r_["rel_pct"],
                                        "σT from the scatter (K)": r_["sigma_T"]} for r_ in rep]).style.format(precision=4), hide_index=True, width="stretch")
            st.caption("σT = (std Δ / Δ) / Sᵣ: the experimental resolution, including everything the counting statistics do not (drift, positioning, laser noise).")
        if runs["runs"]:
            st.markdown("**Independent runs** (one per file with at least 3 temperatures)")
            st.dataframe(pd.DataFrame(runs["runs"]).rename(columns={"run": "file", "n": "spectra", "dE_cm1": "ΔE (cm⁻¹)"}).style.format(precision=2), hide_index=True, width="stretch")
            if runs["n_runs"] < 2:
                st.caption("Measure the series two or three times (separate files) to get the reproducibility into the budget.")
        if hy.get("status") == "ok":
            st.markdown("**Heating vs cooling** (measurement order of the files / scans)")
            st.dataframe(pd.DataFrame(hy["rows"]).style.format(precision=2), hide_index=True, width="stretch")
        else:
            st.caption(f"Heating vs cooling: {hy.get('status')} – measure the series up and back down to test for hysteresis.")
        if not rep and runs["n_runs"] < 2:
            st.info("No repeats found. General practice: measure each temperature at least twice and the whole series in two or three runs.")
    _laser_heating(info, bf, swapped, t_order, counts, B, g)
    _self_check_ui(info, bf, T, counts, t_order)
    # ---- record and store
    sv = {"dE_cm1": abs(bf["dE_cm1"]), "u_dE_cm1": u_dE, "alpha": bf["alpha"], "Sr_Tmin_pct_per_K": float(Sr[0]), "R2_boltzmann": bf["R2"],
          "sigmaT_resolution_Tmin_K": float(B["T_res"][0]), "sigmaT_calibration_Tmin_K": float(B["T_cal"][0]), "confidence": v}
    store_put(sample, "pl", sv, "Photoluminescence page")
    files = sorted({(p.get("file", ""), p.get("file_md5", "")) for p in t_order})
    rec = __import__("mfs.plqa", fromlist=["x"]).analysis_record(version, [{"file": a, "md5": b} for a, b in files], {**(settings or {}), "ratio_method": info.get("method"),
                                                                    "windows_or_ranges_nm": info.get("windows"), "swapped": swapped, "boltzmann_offset": offs, "counts": counts,
                                                                    "sensor_random_K": sT_rand, "sensor_offset_K": sT_sys, "response_uncertainty": resp_rel,
                                                                    "fit_bands_nm": (info.get("fit") or {}).get("bands"), "fit_groups": (info.get("fit") or {}).get("groups"),
                                                                    "fit_model": (info.get("fit") or {}).get("model")},
                                                                    {**sv, "table": res.to_dict("records"), "budget": B["rows"]}, checks)
    c1, c2, c3 = st.columns(3)
    dl(c1, "⬇️ CSV thermometry", res.to_csv(index=False).encode("utf-8-sig"), "PL_thermometry.csv", "pl_th_csv")
    dl(c2, "⬇️ Analysis record (JSON: settings, file fingerprints, results, checks – for the SI)", rec.encode("utf-8"), "PL_thermometry_record.json", "pl_th_json")
    # v2.18.6: analysis report (PDF): verdict, settings, equations, calibration figure, table, budget, checks, limitations
    try:
        from mfs import plreport as PRPT
        # S-Figuryte-3: the fitted band shapes per spectrum (only for the plain two-band fit; the record is text, the numbers are not changed)
        shapes = None
        try:
            if info.get("method") == "fit" and info.get("spec") and info.get("windows"):
                from mfs import plwindows as PW
                (w1_, w2_) = info["windows"]
                # same defaults as the ratio fit (free_eta=False, weighting="sqrt", warm_start=False); keep them in step if one changes
                # the weighting the user chose in the window-fit panel (session key pl_winweight); sqrt when not chosen
                wgt_ = "poisson" if str(ss.get("pl_winweight", "")).startswith("Poisson") else "sqrt"
                shapes = PW.window_fit_record(info["spec"], w1_, w2_, counts, False, wgt_)
        except Exception:
            shapes = None
        rep_settings = {"sample": sample, "ratio method": info.get("method"), "windows or ranges (nm)": info.get("windows"),
                        "fitted band shapes (window fit, per spectrum)": shapes,
                        "band swapped (Δ inverted)": swapped, "Boltzmann offset C": offs, "counts (Poisson weights)": counts,
                        "temperature random error (K)": sT_rand, "temperature offset (K)": sT_sys, "response uncertainty": resp_rel,
                        "fit bands (nm)": (info.get("fit") or {}).get("bands"), "fit model": (info.get("fit") or {}).get("model")}
        pdf_bytes = PRPT.build_pdf(sample=sample, version=version, T=T, D=D, D_stat=D_stat, bf=bf, Sr=Sr, table=res, budget_rows=B["rows"],
                                   checks=checks, verdict=v, settings=rep_settings, u_dE_cm1=u_dE, T_res=B["T_res"], T_cal=B["T_cal"],
                                   sigma_scaled=any(c.get("title") == "Ratio noise" for c in checks),
                                   model=QA.model_comparison(T, D, D_stat), refs=(refs or None), T_pred=T_pred, corrected=bool(ss.get("pl_cal_applied")),
                                   band_info={"method": info.get("method"), "windows": info.get("windows"), "fit": (info.get("fit") or {}).get("bands")})
        dl(c3, "📄 Report (PDF)", pdf_bytes, "PL_analysis_report.pdf", "pl_report_pdf")
    except Exception as _e:
        st.caption(f"Report (PDF) could not be built: {_e}")


def _map_range(info):
    if info.get("windows"):
        (a1, b1), (a2, b2) = info["windows"]
        return min(a1, a2) - 25, max(b1, b2) + 25
    return -1e9, 1e9


def _fit_key(info, bf):
    f = info.get("fit")
    if not f:
        return None
    return hashlib.md5(repr((info.get("method"), [r["ratio"] for r in f["rows"]], round(bf["dE_cm1"], 3))).encode()).hexdigest()


def _self_check_ui(info, bf, T, counts, t_order):
    from mfs import plqa as QA
    if not info.get("fit"):
        return
    with st.expander("🧪 Simulation self-check – consistency only, not accuracy", expanded=False):
        st.caption("Series like yours are simulated from the same model with a known energy gap and analysed again. This tests the statistics and the code path, not the model itself: "
                   "if the band model is wrong for your sample, this check still passes. "
                   "with the same model. An unbiased analysis returns the known gap; honest error bars contain it in about 95 % of the runs (within 2 error bars). "
                   "Takes about as long as six fits.")
        if st.button("▶ Run the self-check (6 simulated series)", key="pl_selfcheck_go"):
            ped = [float(np.mean(p.get("baseline", np.zeros(1)))) for p in t_order]
            bar = st.progress(0.0, text="Simulating and refitting …")
            dE = abs(bf["dE_eV"])
            r = QA.self_check(info["fit"], T, dE, bf["alpha"], counts, ped if counts else None, 6, 0,
                              lambda f: bar.progress(f, text="Simulating and refitting …"), inverse=bool(bf.get("_swapped")))
            bar.empty()
            st.session_state["pl_selfcheck"] = (_fit_key(info, bf), r)
            st.rerun()
        sc = st.session_state.get("pl_selfcheck")
        if sc and sc[0] == _fit_key(info, bf) and sc[1].get("status") == "ok":
            r = sc[1]
            st.dataframe(pd.DataFrame([{"run": k + 1, "ΔE found (cm⁻¹)": x["dE_cm1"], "± (cm⁻¹)": x["err_cm1"], "deviation / error": x["z"]} for k, x in enumerate(r["runs"])])
                         .style.format(precision=2), hide_index=True, width="stretch")
            st.caption(f"Known ΔE {r['truth_cm1']:.0f} cm⁻¹; mean found {np.mean([x['dE_cm1'] for x in r['runs']]):.0f} cm⁻¹ (bias {r['bias_pct']:+.1f} %); spread of deviation/error "
                       f"{r['z_sd']:.2f} (1 = honest error bars, > 1 = too small); {100 * r['coverage']:.0f} % within 2 error bars.")


def _laser_heating(info, bf, swapped, t_order, counts, B, g):
    """Spectra at ONE set temperature and several excitation powers: the temperature read from the ratio should not depend on the power."""
    from mfs import plextras as PX
    from mfs import plqa as QA
    from mfs import plratio as PR
    with st.expander("🔥 Laser-heating check (spectra at one set temperature, several powers)", expanded=False):
        st.caption("Upload spectra measured at ONE set temperature with different excitation powers (the power in the name, e.g. `..._250mW`). The ratio is measured with "
                   "exactly the same method and turned into a temperature with your calibration; the slope is the laser heating per unit power.")
        ups = st.file_uploader("Power-series spectra", type=["csv", "txt", "dat", "xlsx", "xls", "fs"], accept_multiple_files=True, key="pl_heat_files")
        if not ups:
            return
        sp = []
        for u in ups:
            try:
                for s in PL.read_pl(u.getvalue(), u.name):
                    s["P"] = PX.parse_power(s["name"]) or PX.parse_power(u.name)
                    sp.append(s)
            except Exception as e:
                st.warning(f"{u.name}: {e}")
        sp = [s for s in sp if s.get("P")]
        if len(sp) < 3:
            st.info("Need at least three spectra with a power in the name.")
            return
        plan = t_order[0].get("_plan")
        if plan is not None:
            from mfs import plbackground as BG
            rr = BG.apply_plan([{"wl": s["wl"], "I": s["I"]} for s in sp], plan, plan.get("_dark"))
            for s, r in zip(sp, rr):
                s["I_corr"] = r["y"]
        x0 = info["spec"][0]["x"]
        lo, hi = float(x0.min()), float(x0.max())
        spec = [{"T": 300.0 + 1e-3 * k, "x": s["wl"][(s["wl"] >= lo) & (s["wl"] <= hi)], "y": np.asarray(s.get("I_corr", s["I"]), float)[(s["wl"] >= lo) & (s["wl"] <= hi)],
                 "var": np.asarray(s["I"], float)[(s["wl"] >= lo) & (s["wl"] <= hi)] if counts else None} for k, s in enumerate(sp)]
        if info.get("fit"):
            f = info["fit"]
            bl = list(f["bands"]) + list(f.get("extra", []))
            r = PR.fit_global(spec, bl, f["window"], f.get("bg_order", 1), "none", counts, 60, f.get("groups"),           # v2.19.5: warm start, ratios only
                              start=PR._start_for(f.get("params"), f.get("fit_bands", bl), bl), errors=False)
            Dp = np.array([q["ratio"] for q in r["rows"]])
        else:
            Dp = QA.window_ratio(spec, *info["windows"], info.get("win_baseline", "none"))[2]
        if swapped:
            Dp = 1.0 / Dp
        lh = QA.laser_heating([s["P"] for s in sp], Dp, bf["alpha"], bf["dE_eV"], bf.get("C", 0.0))
        if lh.get("status") != "ok":
            st.warning(lh.get("status"))
            return
        dT = lh["heating_at_max"]
        lim = max(1.0, 2 * float(np.nanmedian(B["T_res"])))
        (st.success if abs(dT) <= lim else st.warning)(
            f"Heating {lh['slope']:.3g} ± {lh['slope_err']:.2g} K per W; at the highest power the sample is {dT:+.1f} K warmer than at zero power "
            f"(extrapolated T₀ = {lh['T0']:.1f} K)." + ("" if abs(dT) <= lim else " Lower the power density for the calibration, or correct the temperatures."))
        st.dataframe(pd.DataFrame(lh["rows"]).style.format(precision=3), hide_index=True, width="stretch")
