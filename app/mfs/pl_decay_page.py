"""v2.18.4: the 'Decay curves (lifetimes)' part of the Photoluminescence page."""
from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st

from mfs import levels as LV
from mfs import pl as PL
from mfs import pldecay as DC


def _example():
    """Synthetic decays at five temperatures (bi-exponential, the long lifetime quenches) - for trying the page."""
    rng = np.random.default_rng(4)
    t = np.arange(0, 3000, 4.0)
    out = []
    for T in (293, 313, 333, 353, 373):
        tau = 400 / (1 + 2e4 * np.exp(-0.35 / (8.617e-5 * T)))
        tt = np.clip(t - 40, 0, None)
        y = (6000 * np.exp(-tt / tau) + 1500 * np.exp(-tt / 60)) * (t >= 40) + 15
        out.append({"name": f"{T} K", "T_text": f"{T} K", "t": t, "y": rng.poisson(y).astype(float), "unit": "µs"})
    return out


def render(g):
    show, dl = g["show"], g["dl"]
    st.caption("Decay curves as text / CSV / Excel: the first column is the time (unit from the header – `Time (ns)`, `µs`, `ms`; otherwise ns), then one column per decay "
               "(the header names it; a temperature in the name such as `300K` or `25C` is read). From time-correlated single-photon counting (TCSPC) or multichannel scaling (MCS). "
               "Export from the instrument software as text.")
    ex = st.checkbox("Use a synthetic EXAMPLE (5 temperatures, two lifetimes)", False, key="dc_ex")
    ups = st.file_uploader("Decay files", type=["csv", "txt", "dat", "xlsx", "xls"], accept_multiple_files=True, key="dc_files")
    curves = []
    if ex and not ups:
        curves = _example()
    for u in ups or []:
        try:
            for c in DC.read_decay(u.getvalue(), u.name):
                c["file"] = u.name
                curves.append(c)
        except Exception as e:
            st.error(f"{u.name}: {e}")
    if not curves:
        st.info("Upload decay curves, or tick the example.")
        return
    unit = curves[0]["unit"]
    c1, c2, c3 = st.columns(3)
    model_l = c1.radio("Model", ["Tail fit (exponentials after the peak)", "Rise and decay (fed level, upconversion)", "Reconvolution with a measured IRF"], key="dc_model",
                       help="Tail fit: the excitation pulse is short compared with the lifetime. Rise and decay: the emitting level is filled with a delay. "
                            "Reconvolution: lifetimes close to the instrument response need the measured instrument response function (IRF).")
    model = {"T": "tail", "R": "rise"}.get(model_l[0], "reconv")
    counts = c2.checkbox("Photon counts (Poisson weights)", all(np.allclose(c["y"], np.round(c["y"])) for c in curves), key="dc_counts")
    max_exp = c3.slider("Up to … exponentials", 1, 3, 3, key="dc_nexp", help="The simplest model that is not clearly beaten is kept (information criterion, separated lifetimes, every amplitude > 2 %).")
    irf = None
    if model == "reconv":
        iu = st.file_uploader("IRF (scatterer measured with the same settings, same time grid)", type=["csv", "txt", "dat"], key="dc_irf")
        if iu is None:
            st.info("Upload the instrument response (IRF) for the reconvolution fit.")
            return
        try:
            ic = DC.read_decay(iu.getvalue(), iu.name)[0]
            irf = ic
        except Exception as e:
            st.error(f"IRF: {e}")
            return
    tmin = float(min(c["t"].min() for c in curves))
    tmax = float(max(c["t"].max() for c in curves))
    rng = st.slider(f"Time range of the fit ({unit}) – empty start = automatic", tmin, tmax, (tmin, tmax), key="dc_rng")
    auto_start = st.checkbox("Start the tail fit automatically just after the peak", True, key="dc_auto") if model == "tail" else False
    rows, fits = [], []
    for c in curves:
        ir = None if irf is None else np.interp(c["t"], irf["t"], irf["y"], left=0.0, right=0.0)
        b, _ = DC.fit_auto(c["t"], c["y"], model, ir, None if (auto_start and model == "tail") else rng[0], rng[1], counts, max_exp)
        fits.append(b)
        T = PL.to_kelvin([PL.parse_temperature(c["T_text"])], "auto")[0]
        if b.get("status") != "ok":
            rows.append({"decay": c["name"], "T (K)": T, "status": b.get("status")})
            continue
        r = {"decay": c["name"], "T (K)": T, "model": f"{b['n_exp']} exp." if model != "rise" else "rise + decay"}
        for k, (tau, e) in enumerate(zip(b["taus"], b["tau_err"])):
            r[f"τ{k + 1} ({unit})"] = tau
            r[f"± τ{k + 1}"] = e
            if model != "rise":
                r[f"A{k + 1} (%)"] = 100 * b["frac_amp"][k]
        if model == "rise":
            r[f"rise ({unit})"], r["± rise"] = b["rise"], b["rise_err"]
        r[f"τ amplitude-avg ({unit})"] = b["tau_amp"]
        r[f"τ intensity-avg ({unit})"] = b["tau_int"]
        r["χ² (red.)"], r["DW"], r["quality"] = b["chi2_red"], b["dw"], b["flag"]
        rows.append(r)
    tab = pd.DataFrame(rows)
    bad = [r for r in rows if r.get("quality") not in ("good",)]
    if bad:
        st.warning(f"{len(bad)} of {len(rows)} decays are not described well (structured residual or χ² far from 1): " + ", ".join(str(r["decay"]) for r in bad[:6])
                   + ". Possible reasons: another lifetime component, a rise time (try 'Rise and decay'), the instrument response (try reconvolution), pile-up (count rate above "
                     "about 5 % of the repetition rate), or a wrong start of the fit.")
    else:
        st.success("All decays are described well (random residuals, χ² close to 1).")
    sel = st.selectbox("Show the fit of", list(range(len(curves))), format_func=lambda i: curves[i]["name"], key="dc_pick")
    b = fits[sel]
    if b.get("status") == "ok":
        import matplotlib.pyplot as plt
        from mfs import style as S
        S.apply_style(g.get("FONT", "Arial") if isinstance(g, dict) else "Arial")
        fig = plt.figure(figsize=(g["WIDTH"], g["WIDTH"] * 0.8))
        gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
        ax, rax = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
        DC.draw_decay(ax, rax, b, unit)
        ax.tick_params(labelbottom=False)
        txt = ", ".join(f"τ{k + 1} = {t_:.3g} ± {e_:.2g} {unit}" for k, (t_, e_) in enumerate(zip(b["taus"], b["tau_err"])))
        ax.set_title(txt + (f", rise {b['rise']:.3g} {unit}" if b.get("rise") else ""), fontsize=6)
        show(fig, "PL_decay_" + "".join(ch if ch.isalnum() else "_" for ch in curves[sel]["name"])[:30],
             pd.DataFrame({"time": b["t"], "measured": b["y"], "fit": b["fit"], "weighted residual": b["wres"]}), None, "Decay fit")
    st.dataframe(tab.style.format(precision=4, na_rep="–"), hide_index=True, width="stretch")
    st.caption("Amplitude-weighted average τ = Σaᵢτᵢ/Σaᵢ (rates, energy-transfer efficiency); intensity-weighted τ = Σaᵢτᵢ²/Σaᵢτᵢ (what most instruments report as 'average lifetime'). "
               "State which one you report. Lifetime errors are from the fit (scaled by √χ² when χ² > 1).")
    dl(st, "⬇️ CSV lifetimes", tab.to_csv(index=False).encode("utf-8-sig"), "PL_lifetimes.csv", "dc_csv")
    # lifetime vs temperature
    Ts = [(r["T (K)"], r[f"τ intensity-avg ({unit})"]) for r in rows if r.get("T (K)") is not None and f"τ intensity-avg ({unit})" in r]
    if len(Ts) >= 4:
        T = np.array([a for a, _ in Ts], float)
        tau = np.array([b_ for _, b_ in Ts], float)
        o = np.argsort(T)
        T, tau = T[o], tau[o]
        f = DC.fit_lifetime_T(T, tau)
        st.markdown("**Lifetime vs temperature** (intensity-weighted average)")
        if f.get("status") == "ok":
            m1, m2, m3 = st.columns(3)
            m1.metric("Eₐ (quenching)", f"{f['Ea_meV']:.0f} ± {f['Ea_err_meV']:.0f} meV")
            m2.metric("τ₀", f"{f['I0']:.3g} {unit}")
            m3.metric(f"Sᵣ at {T.max():.0f} K", f"{f['Sr'][-1]:.2f} % K⁻¹")
        else:
            st.caption(f"τ(T) fit: {f.get('status')}")
        newfig = g["newfig"]

        def d_tau(ax, T=T, tau=tau, f=f):
            Tf = np.linspace(T.min(), T.max(), 200)
            PL.draw_trend(ax, T, tau, f"τ ({unit})", fit_T=Tf if f.get("status") == "ok" else None, fit_y=f["f"](Tf) if f.get("status") == "ok" else None,
                          fit_label="τ₀ / (1 + C exp(−Eₐ/k_BT))")
        fig2, ax2 = newfig(0.75); d_tau(ax2)
        show(fig2, "PL_lifetime_vs_T", pd.DataFrame({"T (K)": T, "tau": tau}), d_tau, "Lifetime vs temperature")
    elif LV.at_least("Advanced"):
        st.caption("Give four or more decays a temperature in the name to fit τ(T) (lifetime thermometry).")
