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
            raw = f.getvalue()
            md5 = hashlib.md5(raw).hexdigest()
            for s in PL.read_pl(raw, f.name):
                s["file"] = f.name
                s["file_md5"] = md5
                s["_order"] = len(out)                      # v2.18.4: measurement order (files, then scans) - for the heating / cooling check
                s["I_raw"] = np.asarray(s["I"], float).copy()   # v2.18.4: before any correction - for the saturation check and the Poisson weights
                out.append(s)
        except Exception as e:
            st.error(f"{f.name}: {e}")
    if len(out) > 1 and len({s["name"] for s in out}) == len(out):      # v2.18.5: one series in several files -> the part that differs (e.g. 20C … 50C)
        short = PL._short_names([s["name"] for s in out])
        if len(set(short)) == len(short):
            for s, sh in zip(out, short):
                s.setdefault("meta", {}).setdefault("full_name", s["name"])
                if PL.parse_temperature(sh)[0] is not None or PL.parse_temperature(s["T_text"])[0] is None:
                    s["T_text"] = sh
                s["name"] = sh
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
    if (target == "pl_b1" or target == "pl_b2") and ss.get("pl_method") not in (METHODS[2], METHODS[3]):
        ss["pl_method"] = METHODS[1]
    ss["pl_pick_msg"] = f"Range {lo:.1f}–{hi:.1f} nm set."


METHODS = ["Automatic – global fit of all spectra (recommended)", "Integration windows (choose them on the graph)", "Compare both",
           "Manifold fit – crowded / Stark-split bands (sum of the lines in two ranges)"]


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


def _extras_panel(proc, rng, g):
    """v2.15.1: CIE 1931 colour of the spectra and the power dependence (log-log slope)."""
    show, newfig, dl = g["show"], g["newfig"], g["dl"]
    from mfs import plextras as PX
    with st.expander("🎨 CIE 1931 colour (x, y, colour temperature)", expanded=False):
        cc1, cc2 = st.columns(2)
        photon = cc1.checkbox("Intensities are photon counts (convert to energy before weighting)", True, key="pl_cie_ph")
        wp = cc2.selectbox("White point for the dominant wavelength and purity", list(PX.WHITE), key="pl_cie_wp")
        if not st.session_state.get("pl_cal_applied"):
            st.caption("⚠️ Colour coordinates need spectra corrected for the spectral response of the spectrometer; without the correction they are instrument-specific.")
        rows, pts, cols = [], [], []
        pal = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7F7F7F", "#000000"]
        for k, q in enumerate(proc):
            c = PX.cie_xy(q["wl_c"], q["y_c"], photon)
            dw = PX.dominant_wavelength(c["x"], c["y"], PX.WHITE[wp])
            cov = (min(780.0, float(q["wl_c"].max())) - max(380.0, float(q["wl_c"].min()))) / 400.0
            rows.append({"spectrum": q["name"], "T (K)": q["T"], "x": c["x"], "y": c["y"], "CCT (K)": c["CCT"],
                         ("complementary λ (nm)" if dw["complementary"] else "dominant λ (nm)"): dw["dominant_nm"], "purity (%)": dw["purity_pct"],
                         "visible range covered (%)": round(100 * max(0.0, min(1.0, cov))),
                         "in sRGB gamut": ("yes" if PX.xy_to_srgb(c["x"], c["y"])[1] else "no (not showable)") if np.isfinite(c["x"]) else "–"})
            if np.isfinite(c["x"]):
                pts.append((c["x"], c["y"])); cols.append(pal[k % len(pal)])
        tab = pd.DataFrame(rows)
        st.dataframe(tab.style.format(precision=4, na_rep="–"), hide_index=True, width="stretch")
        st.caption("Only the 380–780 nm part of each spectrum counts. Colour-matching functions: analytic fit of the CIE 1931 2° observer "
                   "(Wyman et al. 2013, errors ≈ 0.003); the colour temperature (McCamy) is only meaningful close to the Planckian locus.")
        if pts:
            def draw(ax, pts=pts, cols=cols):
                PX.draw_cie(ax, pts, colors=cols)
            fig, ax = newfig(1.0); draw(ax)
            show(fig, "PL_CIE", tab, draw, "CIE 1931 chromaticity")
    with st.expander("🔦 Power dependence (slope of log I vs log P)", expanded=False):
        st.caption("Needs spectra at several excitation powers. The power is read from the name (`250mW`, `1.5W`) and can be edited; the intensity is the "
                   "integral over the chosen wavelength range. Slope n ≈ 1: linear process; n ≈ 2 / 3: two- / three-photon upconversion; "
                   "n < 1 at high power: saturation.")
        base = pd.DataFrame({"use": True, "spectrum": [q["name"] for q in proc], "power (W)": [PX.parse_power(q["name"]) for q in proc],
                             "integral": [PL.band_integral(q["wl_c"], q["y_c"], rng[0], rng[1], "linear") for q in proc]})
        ed = st.data_editor(base, key="pl_pw_tab|" + "|".join(base["spectrum"].astype(str))[:200], hide_index=True, width="stretch",
                            disabled=["spectrum", "integral"])
        sel = ed[ed["use"] & ed["power (W)"].notna()]
        if len(sel) < 3:
            st.info("Give at least three spectra a power value (name them like `..._250mW` or type it in the table).")
        else:
            f = PX.fit_power(sel["power (W)"].to_numpy(float), sel["integral"].to_numpy(float))
            if f.get("status") == "ok":
                m1, m2, m3 = st.columns(3)
                m1.metric("slope n", f"{f['n']:.2f}", f"± {f['n_err']:.2f}", delta_color="off")
                m2.metric("R²", f"{f['R2']:.4f}")
                if f.get("saturating"):
                    m3.metric("saturation", "yes", f"n {f['slope_low']:.2f} → {f['slope_high']:.2f}", delta_color="inverse")
                P_ = sel["power (W)"].to_numpy(float); I_ = sel["integral"].to_numpy(float)

                def draw_p(ax, P_=P_, I_=I_, f=f):
                    ax.loglog(P_, I_, "o", ms=4, color="#0072B2")
                    pp = np.geomspace(P_.min(), P_.max(), 50)
                    ax.loglog(pp, f["f"](pp), "-", color="#D55E00", lw=1.0, label=f"n = {f['n']:.2f} ± {f['n_err']:.2f}")
                    ax.set_xlabel("Excitation power (W)"); ax.set_ylabel("Integrated intensity"); ax.legend(frameon=False, fontsize=6)
                fig, ax = newfig(0.8); draw_p(ax)
                show(fig, "PL_power_dependence", sel, draw_p, "Power dependence")
            else:
                st.warning(f.get("status"))


def _data_checks(spectra):
    """v2.18.4: what the files say about the measurement, and the checks that do not need any analysis (saturation, settings that differ)."""
    from mfs import plqa as QA
    sat = QA.saturation(spectra)
    rows, differ = QA.settings_table(spectra)
    bad = [s for s in sat if s["status"] != "ok"]
    worst = "fail" if any(s["status"] == "fail" for s in bad) else ("warn" if bad or differ else "ok")
    label = {"ok": "✅ Measurement settings and data checks – nothing found",
             "warn": "⚠️ Measurement settings and data checks – look at this",
             "fail": "⛔ Measurement settings and data checks – saturated spectra"}[worst]
    with st.expander(label, expanded=worst == "fail"):
        for s in bad:
            (st.error if s["status"] == "fail" else st.warning)(f"{s['name']}: {s['reason']}.")
        if differ:
            st.warning("These settings differ inside the series: " + ", ".join(differ) + ". Intensities (and ratios, when the response changes with the setting) are only comparable "
                       "with the same settings – or use 'per second' normalisation when only the dwell time differs.")
        tab = pd.DataFrame(rows)
        tab["max (raw counts)"] = [round(s["max"]) if np.isfinite(s["max"]) else None for s in sat]
        st.dataframe(tab, hide_index=True, width="stretch")
        st.caption("Read from the files where they contain it (Edinburgh .FS: detector and excitation source; Fluoracle text export: excitation, dwell time, repeats; "
                   "Ocean Optics: integration time and scans). Saturation: a maximum at the converter limit of an array detector, or three or more points with exactly the maximum.")


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


def _parse_windows(txt):
    """'900-925; 1090-1250' -> [(900.0, 925.0), (1090.0, 1250.0)] (also with commas and 'to')."""
    import re
    out = []
    for part in re.split(r"[;,]\s*(?=\d)|\n", str(txt)):
        m = re.findall(r"[-+]?\d+(?:[.]\d+)?", part.replace("–", "-").replace(" to ", "-"))
        if len(m) >= 2:
            a, b = float(m[0]), abs(float(m[1])) if float(m[1]) < 0 and float(m[0]) > 0 else float(m[1])
            out.append((min(a, b), max(a, b)))
    return out


def _background_block(spectra, mode, g):
    """v2.18.4: standard background subtraction of the whole series. Returns (plan, [{'y','baseline','spikes'}]) - the plan is the same for every spectrum."""
    from mfs import plbackground as BG
    show, newfig, dl = g["show"], g["newfig"], g["dl"]
    grid = np.asarray(spectra[0]["wl"], float)
    ys = [np.interp(grid, s["wl"], s["I"]) for s in spectra]
    kind = BG.range_kind(grid)
    det = next((s["meta"].get("detector") for s in spectra if s.get("meta") and s["meta"].get("detector")), "")
    with st.expander(f"🧹 Background subtraction – {mode}" + (f" · {kind} spectrum" + (f" · {det}" if det else "")), expanded=True):
        c1, c2 = st.columns(2)
        spike = c1.checkbox("Remove single-point spikes (cosmic rays, hot pixels)", False, key="pl_spike",
                            help="A point that differs from the median of its neighbours by more than 8 noise levels is replaced. Array detectors (CCD, InGaAs); not needed for a PMT scan.")
        dark, wins, forced = None, None, {}
        if mode == "Dark / background file":
            df_ = c2.file_uploader("Dark / background spectrum (same detector and settings)", type=["csv", "txt", "dat", "xlsx", "xls", "fs"], key="pl_dark")
            if df_ is not None:
                try:
                    d0 = PL.read_pl(df_.getvalue(), df_.name)[0]
                    dark = (d0["wl"], d0["I"])
                    c2.caption(f"✅ {df_.name}: {len(d0['wl'])} points. It is scaled in the signal-free windows and subtracted from every spectrum.")
                except Exception as e:
                    st.warning(f"Dark file not read: {e}")
            else:
                c2.info("Upload the dark / background spectrum; until then the automatic background is used.")
        if mode == "Manual windows":
            auto_w = BG.find_free_windows(grid, ys)
            txt = c2.text_input("Signal-free windows (nm), e.g. 900-925; 1090-1250", "; ".join(f"{a:.0f}-{b:.0f}" for a, b in auto_w), key="pl_win_txt",
                                help="Regions without emission. The start value is what the automatic search found: edit it.")
            wins = _parse_windows(txt) or None
            m1, m2 = st.columns(2)
            meth = m1.selectbox("Method", ["polynomial in the windows", "rubber band (lower convex hull)", "ModPoly (iterative polynomial)", "arPLS (penalised least squares)", "SNIP (peak clipping)"], key="pl_bg_meth")
            order = m2.slider("Polynomial order (0 = constant pedestal)", 0, 4, 1, key="pl_bg_order")
            forced = {"polynomial in the windows": ("poly", order), "rubber band (lower convex hull)": ("rubber", None), "ModPoly (iterative polynomial)": ("modpoly", order),
                      "arPLS (penalised least squares)": ("arpls", None), "SNIP (peak clipping)": ("snip", None)}[meth]
        plan = BG.autoselect(grid, ys, dark=dark, windows=wins, max_order=2)
        plan["_dark"] = dark
        if forced:
            plan["method"], plan["order"] = forced[0], forced[1]
            if plan["method"] == "poly" and plan["order"] == 0:
                plan["method"] = "offset"
            plan["notes"] = ["Manual: " + BG.LABELS.get(plan["method"], plan["method"]) + (f", order {plan['order']}" if plan["order"] is not None else "")] + plan["notes"][-1:]
        res = BG.apply_plan(spectra, plan, dark, spike)
        ped = [float(np.mean(r["baseline"])) for r in res]
        chk = [BG.check_after(s["wl"], r["y"], plan) for s, r in zip(spectra, res)]
        m = st.columns(4)
        m[0].metric("Method", BG.LABELS.get(plan["method"], plan["method"]).split(" (")[0], f"order {plan['order']}" if plan.get("order") is not None else None, delta_color="off")
        m[1].metric("Signal-free windows", f"{len(plan['windows'])}", f"{sum(b - a for a, b in plan['windows']):.0f} nm in total" if plan["windows"] else "none", delta_color="off")
        m[2].metric("Mean background", f"{np.mean(ped):,.0f}".replace(",", " "), f"{min(ped):,.0f} – {max(ped):,.0f}".replace(",", " "), delta_color="off")
        bias = np.nanmean([abs(c["bias"]) for c in chk]) if plan["windows"] else float("nan")
        m[3].metric("Left in the free windows", f"{bias:.2f} σ" if bias == bias else "–", "should be 0 ± 1 noise level (σ)", delta_color="off")
        if bias == bias and bias > 1.0:
            st.warning(f"After the subtraction the signal-free windows are not flat around zero ({bias:.1f} σ left). Change the windows or the method, or load a dark / background file.")
        for n in plan["notes"]:
            st.caption("ℹ️ " + n)
        if sum(r["spikes"] for r in res):
            st.caption(f"🧽 {sum(r['spikes'] for r in res)} spike(s) removed.")
        order_ = sorted(range(len(spectra)), key=lambda i: (spectra[i].get("T") is None, spectra[i].get("T") or 0))

        def draw_bg(ax, spectra=spectra, res=res, plan=plan, order_=order_):
            from mfs import style as S
            cols = S.colours_for([spectra[i]["name"] for i in order_])
            for c, i in zip(cols, order_):
                ax.plot(spectra[i]["wl"], spectra[i]["I"], color=c, lw=0.7)
                ax.plot(spectra[i]["wl"], res[i]["baseline"], color="k", lw=0.8, ls="--")
            for a, b in plan["windows"]:
                ax.axvspan(a, b, color="#5EEAD4", alpha=0.22, lw=0)
            ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Intensity (counts)")
            lo = min(float(np.percentile(s["I"], 1)) for s in spectra); hi = max(float(np.percentile(s["I"], 60)) for s in spectra)
            ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.35 * (hi - lo))
            ax.set_title("raw spectra, background (dashed), signal-free windows (green)", fontsize=7)
        fig, ax = newfig(0.55); draw_bg(ax)
        dfb = pd.DataFrame({"wavelength_nm": grid, **{f"{s['name']} raw": np.interp(grid, s["wl"], s["I"]) for s in spectra},
                            **{f"{s['name']} background": np.interp(grid, s["wl"], r["baseline"]) for s, r in zip(spectra, res)}})
        show(fig, "PL_background", dfb, draw_bg, "Background of the PL spectra")
        dfc = pd.DataFrame({"wavelength_nm": grid, **{s["name"]: np.interp(grid, s["wl"], r["y"]) for s, r in zip(spectra, res)}})
        dl(st, "⬇️ Background-corrected spectra (CSV)", dfc.to_csv(index=False).encode("utf-8-sig"), "PL_background_corrected.csv", "pl_bg_csv")
        if plan["candidates"]:
            with st.expander("Why this method? (the candidates, ranked by objective numbers)", expanded=False):
                tab = pd.DataFrame([{"method": BG.LABELS.get(c["method"], c["method"]), "order": c["order"], "BIC (lower = better)": None if c["bic"] != c["bic"] else round(c["bic"], 1),
                                     "bias in windows (σ)": None if c["bias"] != c["bias"] else round(c["bias"], 2), "rms in windows (σ)": None if c["rms"] != c["rms"] else round(c["rms"], 2),
                                     "points below −3σ": round(c["below"], 4)} for c in plan["candidates"]])
                st.dataframe(tab, hide_index=True, width="stretch")
                st.caption("A polynomial of higher order must clearly earn its extra parameter (ΔBIC > 10). The flexible methods (rubber band, arPLS, SNIP, ModPoly) can follow or swallow broad emission and are only used when no "
                           "signal-free window exists (the rubber band is the cautious default: it never rises above the data). 'Bias' is the mean left in the signal-free windows after the subtraction, 'rms' their scatter, both in noise levels σ.")
    return plan, res


def render(g):
    show, newfig, files_input, store_put, EX, dl, safe_name, C = (g["show"], g["newfig"], g["files_input"], g["store_put"], g["EX"],
                                                                 g["dl"], g["safe_name"], g["C"])
    st.header("Photoluminescence – steady-state and temperature-dependent")
    st.caption("Edinburgh Instruments `.FS` files (single scans and 'Multiple Scans') are read directly. "
               "Emission (or excitation) spectra as CSV / TXT / DAT / Excel from Horiba, Edinburgh, Ocean Optics, Hitachi … : two columns "
               "(wavelength, intensity) per file, or ONE table with the wavelength in the first column and one column per temperature "
               "(header like `293.15 K`, `25 °C`). Several files = several temperatures: the temperature is read from the column header or "
               "the file name (`..._300K.txt`, `..._25C.csv`) and can be edited below. Typical outputs: spectra coloured by temperature, "
               "emission map, thermal quenching (activation energy, T₅₀), peak shift / width, and ratiometric thermometry "
               "(Boltzmann fit, relative sensitivity Sᵣ, temperature uncertainty σT).")
    kind = st.radio("What did you measure?", ["Emission / excitation spectra", "Decay curves (lifetimes)"], horizontal=True, key="pl_kind")
    if kind.startswith("Decay"):
        from mfs import pl_decay_page as DP
        DP.render(g)
        return
    ex = st.radio("Bundled EXAMPLE (synthetic)", list(EXAMPLES), horizontal=True, key="pl_example")
    files = files_input("PL spectra", ["csv", "txt", "dat", "xlsx", "xls", "fs"], EXAMPLES[ex], "pl", multiple=True)
    if not files:
        return
    spectra = _spectra_from(files)
    if not spectra:
        return
    st.session_state["pl_cal_applied"] = False
    with st.expander("🎚️ Spectral-response correction (optional)", expanded=False):
        from mfs import plextras as _PX
        st.caption("The detector and grating do not see every wavelength equally well. The correction multiplies every spectrum by a factor per wavelength. "
                   "It changes the absolute ratio (α) and the band shapes, not the energy gap or Sᵣ of a thermometer (a constant response multiplies every ratio by "
                   "the same number). Corrected data are no longer photon counts.")
        cm = st.radio("Source", ["none", "correction file (wavelength, factor)", "calibration lamp (measured + certified spectrum)"], horizontal=True, key="pl_cal_src")
        _cw = _cfac = None
        try:
            if cm.startswith("correction"):
                cf = st.file_uploader("Correction file (two columns: wavelength in nm, factor)", type=["txt", "csv", "dat"], key="pl_cal")
                if cf is not None:
                    _cw, _cfac = _PX.read_correction(cf.getvalue())
            elif cm.startswith("calibration"):
                l1, l2 = st.columns(2)
                lm = l1.file_uploader("Lamp measured on YOUR spectrometer (same detector, grating, slits, filters)", type=["txt", "csv", "dat", "fs", "xlsx"], key="pl_lamp_meas")
                lc = l2.file_uploader("Certified lamp spectrum (from the lamp certificate)", type=["txt", "csv", "dat", "xlsx"], key="pl_lamp_cert")
                if lm is not None and lc is not None:
                    a_ = PL.read_pl(lm.getvalue(), lm.name)[0]
                    b_ = PL.read_pl(lc.getvalue(), lc.name)[0]
                    _cw, _cfac = _PX.correction_from_lamp(a_["wl"], a_["I"], b_["wl"], b_["I"])
            if _cw is not None:
                for s_ in spectra:
                    s_["I"] = _PX.apply_correction(s_["wl"], s_["I"], _cw, _cfac)
                cov = _PX.correction_coverage(min(s_["wl"].min() for s_ in spectra), max(s_["wl"].max() for s_ in spectra), _cw)
                st.session_state["pl_cal_applied"] = True
                (st.success if cov > 0.98 else st.warning)(f"Correction applied ({len(_cw)} points, {_cw.min():.0f}–{_cw.max():.0f} nm); it covers {100 * cov:.0f} % of the data range"
                                                          + ("." if cov > 0.98 else " – outside it the factor is held constant, which is NOT a correction."))
        except Exception as e:
            st.warning(f"Correction not applied: {e}")
    # ---- v2.18.4: excitation line and its second order (scattered laser light is not emission)
    _exc = next((s_.get("meta", {}).get("excitation_nm") for s_ in spectra if s_.get("meta", {}).get("excitation_nm")), None)
    with st.expander("✂️ Remove the excitation line (scattered light)", expanded=False):
        e1, e2, e3 = st.columns(3)
        # v2.18.5: never pre-filled from the file - the 'Fixed/Offset' value of an emission scan is the excitation monochromator position, which with an external laser
        # (or at the start of the scan range) is not the excitation wavelength; cutting there would remove real emission
        if _exc:
            e1.caption(f"The file gives {float(_exc):.0f} nm as the excitation / offset position – type it below only if it is the real excitation wavelength.")
        exc = e1.number_input("Excitation wavelength (nm, 0 = none)", 0.0, 3000.0, 0.0, 1.0, key="pl_exc",
                              help="Scattered excitation light and its second grating order (2 × λ) are cut out of every spectrum before anything else.")
        half = e2.number_input("Half-width of the cut (nm)", 0.5, 50.0, 5.0, 0.5, key="pl_exc_hw")
        sec = e3.checkbox("Also the second order (2 × λ)", True, key="pl_exc_2nd")
        if exc > 0:
            cuts = [(exc - half, exc + half)] + ([(2 * exc - half, 2 * exc + half)] if sec else [])
            n_cut = 0
            for s_ in spectra:
                keep = np.ones(len(s_["wl"]), bool)
                for a_, b_ in cuts:
                    keep &= ~((s_["wl"] >= a_) & (s_["wl"] <= b_))
                n_cut = max(n_cut, int((~keep).sum()))
                for k_ in ("wl", "I", "I_raw"):
                    s_[k_] = np.asarray(s_[k_])[keep]
            st.caption(f"{n_cut} points cut per spectrum (" + ", ".join(f"{a_:.0f}–{b_:.0f} nm" for a_, b_ in cuts) + ").")
    # ------------------------------------------------------------ temperatures (editable)
    c1, c2 = st.columns([1, 3])
    unit = c1.radio("Plain numbers in headers are", ["auto", "K", "°C"], horizontal=True, key="pl_unit")
    tk = PL.to_kelvin([PL.parse_temperature(s["T_text"]) for s in spectra], {"auto": "auto", "K": "K", "°C": "C"}[unit])
    tkey = hashlib.md5("|".join(s["file"] + s["name"] for s in spectra).encode()).hexdigest()[:8] + unit
    # v2.18.5: window and range widgets keep values from the PREVIOUS data set otherwise (e.g. a 505-580 nm window from the thermometry example integrated
    # nothing of the 450 nm band of the quenching example -> negative intensities). When the file set changes, they return to the full range.
    if st.session_state.get("_pl_src") != tkey:
        for _k in ("pl_win", "pl_b1", "pl_b2", "pl_fitwin", "pl_reg", "pl_rng", "pl_pick_msg"):
            st.session_state.pop(_k, None)
        st.session_state["_pl_src"] = tkey
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
    BG_OPTS = ["Automatic (recommended)", "Manual windows", "Dark / background file", "Constant: minimum", "Constant: mean of a region", "None"]
    base = p2.selectbox("Background", BG_OPTS, key="pl_bg_mode",
                        help="Automatic: the signal-free windows of the whole series are found and the pedestal / baseline is fitted there (order chosen by the information criterion). "
                             "Manual: you give the windows. Dark file: a measured background spectrum is subtracted.")
    has_dwell = all((s.get("meta") or {}).get("dwell_s") for s in spectra)
    norm = p3.selectbox("Normalisation", ["none (counts)", "max", "area", "to the lowest temperature"] + (["per second (dwell time × scans)"] if has_dwell else []), key="pl_norm",
                        help="Per second: divides by the dwell / integration time (and the number of scans) read from the files, so series measured with different settings can be compared.")
    sm = p4.slider("Smoothing (points)", 0, 31, 0, 2, key="pl_sm", help="Savitzky–Golay, order 3; 0 = none")
    from mfs.pl_thermo import clamp_range
    for _k in ("pl_rng", "pl_reg", "pl_win"):                       # v2.18.5: stored ranges from earlier data must fit the new data
        clamp_range(_k, wl_lo, wl_hi)
    rng = st.slider("Wavelength range (nm)", wl_lo, wl_hi, (wl_lo, wl_hi), 0.5, key="pl_rng")
    reg = None
    if base.startswith("Constant: mean"):
        reg = st.slider("Baseline region (nm)", wl_lo, wl_hi, (wl_lo, min(wl_hi, wl_lo + 0.05 * (wl_hi - wl_lo))), 0.5, key="pl_reg")
    # ---- v2.18.4: background subtraction on the WHOLE spectra of the series (same plan for every spectrum), before the range is cut
    plan, bg_res = None, None
    for s in spectra:
        s["I_corr"], s["baseline"] = s["I"], np.zeros_like(s["I"])
    if base in BG_OPTS[:3]:
        plan, bg_res = _background_block(spectra, base, g)
        if bg_res is not None:
            for s, r in zip(spectra, bg_res):
                s["I_corr"], s["baseline"], s["_plan"] = r["y"], r["baseline"], plan
    proc = []
    for s in spectra:
        m = (s["wl"] >= rng[0]) & (s["wl"] <= rng[1])
        wl, y = s["wl"][m], PL.smooth(s["I_corr"], sm)[m]
        y = PL.baseline_subtract(wl, y, {"Constant: minimum": "minimum", "Constant: mean of a region": "region"}.get(base, "none"), reg)
        proc.append({**s, "wl_c": wl, "y_c": y})
    t_order = sorted([p for p in proc if p["T"] is not None], key=lambda p: p["T"])
    ref_val = None
    if norm.startswith("to the lowest") and t_order:
        ref_val = float(np.max(t_order[0]["y_c"])) or 1.0
    for p in proc:
        if norm.startswith("per second"):
            m_ = p.get("meta") or {}
            y = p["y_c"] / (float(m_.get("dwell_s") or 1.0) * float(m_.get("scans") or 1))
        else:
            y = PL.normalise(p["y_c"], {"none (counts)": "none", "max": "max", "area": "area", "to the lowest temperature": "ref"}[norm], ref_val)
        p["x"], p["y"] = PL.to_axis(p["wl_c"], y, axis)
    ylab = "Intensity (counts)" if norm.startswith("none") else ("Intensity (counts s$^{-1}$)" if norm.startswith("per second") else ("Normalised intensity" if norm != "area" else "Intensity (area = 1)"))
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
        from mfs import pllines as LN
        ions = st.multiselect("Label the transitions of (typical ranges; the host shifts them a little)", list(LN.LINES), key="pl_ions",
                              help="Light bands with the usual emission ranges of these ions on the overlay (nm axis). A guide only: nothing is assigned automatically.")

        def draw_spec(ax, shown=shown, view=view, cmap=cmap, xlab=xlab, ylab=ylab, ions=ions, axis=axis):
            if view == "Emission map":
                PL.draw_map(ax, [{"T": p["T"], "x": p["x"], "y": p["y"]} for p in shown if p["T"] is not None], xlab, "magma", ylab)
            elif view.startswith("Stacked"):
                PL.draw_ridge(ax, shown, xlab, ylab + " + offset", cmap)
            else:
                PL.draw_spectra(ax, shown, xlab, ylab, cmap)
                if ions and axis == "nm":
                    LN.draw_labels(ax, ions, float(min(p["x"].min() for p in shown)), float(max(p["x"].max() for p in shown)), ax.get_ylim()[1] * 0.98)
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
        _data_checks(spectra)
        with LV.gate("Advanced", globals()):
            _extras_panel(proc, rng, g)

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
        clamp_range("pl_win", float(rng[0]), float(rng[1]))
        w = st.slider("Integration window (nm)", rng[0], rng[1], rng, 0.5, key="pl_win")
        I_T = np.array([PL.band_integral(p["wl_c"], p["y_c"], w[0], w[1], "linear") for p in t_order])
        pos = np.array([PL.peak_stats(p["wl_c"], p["y_c"], w[0], w[1])["pos"] for p in t_order])
        fw = np.array([PL.peak_stats(p["wl_c"], p["y_c"], w[0], w[1])["fwhm"] for p in t_order])
        fit = PL.fit_quench(T, I_T)
        okq = fit.get("status") == "ok"
        if okq:
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Eₐ", f"{fit['Ea_meV']:.0f} ± {fit['Ea_err_meV']:.0f} meV" if np.isfinite(fit.get("Ea_err_meV", np.nan)) else f"{fit['Ea_meV']:.0f} meV",
                      help=f"{fit['Ea_eV'] * PL.EV_TO_CM:.0f} cm⁻¹")
            m2.metric("T₅₀", (f"{fit['T50']:.0f} ± {fit['T50_err']:.0f} K" if np.isfinite(fit.get("T50_err", np.nan)) else f"{fit['T50']:.0f} K") if fit["T50"] else "–",
                      help="temperature at which the intensity is half of I₀")
            m3.metric("A", f"{fit['A']:.3g}")
            m4.metric("R²", f"{fit['R2']:.4f}")
            if fit.get("few_points_in_drop", 9) < 3:
                st.warning(f"Only {fit['few_points_in_drop']} temperature(s) lie on the drop of the intensity (between 90 % and 10 % of the maximum): Eₐ and T₅₀ are poorly "
                           "determined. Measure more temperatures where the intensity falls.")
            if T.max() < (fit["T50"] or 0):
                st.caption("ℹ️ T₅₀ lies above the measured range: it is an extrapolation.")
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
        from mfs import pl_thermo as TH
        st.caption("Ratiometric thermometry with two thermally coupled bands: Δ = I₁/I₂ = α·exp(−ΔE/k_BT); "
                   "relative sensitivity Sᵣ = 100 %·|(1/Δ)·∂Δ/∂T| = 100 %·ΔE/(k_BT²); temperature uncertainty σT = (σΔ/Δ)/Sᵣ. "
                   "The confidence verdict above the results collects every check (data, fit, physics, uncertainty, experiment).")
        raw_counts = all(np.allclose(p.get("I_raw", p["I"]), np.round(p.get("I_raw", p["I"]))) for p in t_order) and not st.session_state.get("pl_cal_applied")
        offs = st.checkbox("Add a constant C to the Boltzmann law (Δ = α·exp(−ΔE/k_BT) + C)", False, key="pl_off")
        got = TH.ratio_block(t_order, T, rng, norm, raw_counts, g, METHODS, _brush_picker, _pick_cb, _band_default)
        if got is None:
            return
        I1, I2, D, D_stat, info = got
        counts = bool(st.session_state.get("pl_counts", raw_counts))
        if (D <= 0).any() or not np.isfinite(D).all():
            st.warning("The ratio is not positive for every spectrum – check the bands / windows."); return
        settings = {"background": {k: (plan or {}).get(k) for k in ("method", "order", "windows")} if plan else base, "wavelength_range_nm": list(rng), "smoothing_points": sm,
                    "normalisation": norm, "excitation_cut_nm": st.session_state.get("pl_exc", 0.0), "spectral_correction": bool(st.session_state.get("pl_cal_applied"))}
        TH.results(T, I1, I2, D, D_stat, info, t_order, offs, counts, sample, g, g.get("__version__", ""), settings)
