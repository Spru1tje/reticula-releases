"""v2.18.6 / v2.20.4: analysis report (PDF) of the Photoluminescence thermometry tab.

One document per analysis: verdict, a methods section with the analysis chain and figures (v2.20.4), the choice between the
Boltzmann law and the law with an offset (AIC), the optional validation against an independent reference thermometer, the
settings used, the equations, the calibration figure, the thermometry table, the uncertainty budget, every quality check with
its advice, and the limitations. Built with the same page helper as the defect report (mfs/report.py), so the look matches."""
from __future__ import annotations

import datetime
import io

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch

from mfs.report import _Page

KB_EV = 8.617333262e-5
EV_CM = 8065.544
NAVY = "#1F4E79"
RED = "#C0392B"
PALE = "#EAF2F8"


def _checks_text(checks):
    rows = []
    for c in checks or []:
        rows.append((str(c.get("status", "")).upper(), str(c.get("title", "")), str(c.get("detail", "")), str(c.get("advice", "") or "")))
    return rows


def _dejavu(fn):
    """The figures use superscripts and Greek letters: draw them with DejaVu Sans wherever they are made (also outside the PDF)."""
    import functools

    @functools.wraps(fn)
    def wrapper(*a, **k):
        with plt.rc_context({"font.family": "DejaVu Sans"}):
            return fn(*a, **k)
    return wrapper


@_dejavu
def _calibration_figure(T, D, D_stat, bf, offs=False):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.0, 3.3), gridspec_kw={"width_ratios": [1.3, 1.0]})
    x = 1000.0 / np.asarray(T, float)
    y = np.log(np.asarray(D, float))
    err = np.asarray(D_stat, float) / np.asarray(D, float)
    a1.errorbar(x, y, yerr=err, fmt="o", ms=3.5, color=NAVY, ecolor="#7F9FBF", capsize=2, label="data (ln Δ)")
    xf = np.linspace(x.min(), x.max(), 200)
    Tf = 1000.0 / xf
    a1.plot(xf, np.log(bf["f"](Tf)), color=RED, lw=1.2, label="Boltzmann fit")
    a1.set_xlabel("1000 / T (K⁻¹)"); a1.set_ylabel("ln Δ"); a1.legend(fontsize=7, frameon=False)
    a1.set_title("Boltzmann calibration", fontsize=9)
    res = (np.asarray(D, float) - bf["f"](np.asarray(T, float))) / np.asarray(D, float)
    a2.axhline(0, color="black", lw=0.6)
    a2.plot(np.asarray(T, float), 100 * res, "o", ms=3.5, color=NAVY)
    a2.set_xlabel("T (K)"); a2.set_ylabel("residual (%)"); a2.set_title("Relative residuals", fontsize=9)
    fig.tight_layout()
    return fig


@_dejavu
def _flow_figure():
    """The analysis chain as six boxes with arrows (the methods figure)."""
    fig, ax = plt.subplots(figsize=(8.0, 1.8))
    ax.set_xlim(0, 10); ax.set_ylim(0, 2); ax.axis("off")
    NL = chr(10)
    steps = [("1  Spectra", "one spectrum per" + NL + "temperature" + NL + "(optional correction)"),
             ("2  Band ratio", "Δ = I₁ / I₂" + NL + "in windows or" + NL + "fitted bands"),
             ("3  Global fit", "ln Δ vs 1/T" + NL + "AIC: offset" + NL + "needed?"),
             ("4  Temperature", "T from the" + NL + "calibration;" + NL + "S_r and σT"),
             ("5  Uncertainty", "budget: counting," + NL + "fit, sensor," + NL + "correction"),
             ("6  Checks", "saturation, physics," + NL + "model, reference" + NL + "(if given)")]
    w, gap = 1.44, 0.26
    for i, (title, body) in enumerate(steps):
        x = 0.05 + i * (w + gap)
        ax.add_patch(FancyBboxPatch((x, 0.25), w, 1.6, boxstyle="round,pad=0.02,rounding_size=0.08", fc=PALE, ec=NAVY, lw=1.0))
        ax.text(x + w / 2, 1.58, title, ha="center", va="center", fontsize=7.5, fontweight="bold", color=NAVY)
        ax.text(x + w / 2, 0.85, body, ha="center", va="center", fontsize=5.8, color="#222222", linespacing=1.3)
        if i < len(steps) - 1:
            ax.annotate("", xy=(x + w + gap - 0.03, 1.05), xytext=(x + w + 0.03, 1.05),
                        arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.1))
    return fig


@_dejavu
def _model_figure(T, D, D_stat, mc):
    """Left: the same data with the two candidate laws. Right: their AIC, with the difference."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.0, 3.0), gridspec_kw={"width_ratios": [1.4, 1.0]})
    T = np.asarray(T, float)
    D = np.asarray(D, float)
    x = 1000.0 / T
    a1.errorbar(x, np.log(D), yerr=np.asarray(D_stat, float) / D, fmt="o", ms=3.5, color=NAVY, ecolor="#7F9FBF", capsize=2,
                label="data", zorder=3)
    xf = np.linspace(x.min(), x.max(), 200)
    Tf = 1000.0 / xf
    styles = {"Boltzmann": (NAVY, "-"), "Boltzmann + offset": (RED, "--")}
    models = (mc or {}).get("models", {})
    for name, (col, ls) in styles.items():
        m = models.get(name, {})
        pr = m.get("params")
        if not pr:
            continue
        c = pr[2] if len(pr) > 2 else 0.0
        yy = np.log(pr[0] * np.exp(-pr[1] / (KB_EV * Tf)) + c)
        a1.plot(xf, yy, color=col, ls=ls, lw=1.3, label=f"{name}: ΔE = {m['dE_cm1']:.0f} cm⁻¹")
    a1.set_xlabel("1000 / T (K⁻¹)"); a1.set_ylabel("ln Δ"); a1.legend(fontsize=6.5, frameon=False)
    a1.set_title("Two candidate laws, same data", fontsize=9)
    names = [n for n in styles if n in models and np.isfinite(models[n].get("aic", np.nan))]
    aics = [models[n]["aic"] for n in names]
    if names:
        a2.bar(range(len(names)), aics, color=[styles[n][0] for n in names], alpha=0.85)
        for i, v in enumerate(aics):
            a2.text(i, v, f"{v:.1f}", ha="center", va="bottom", fontsize=7)
        a2.set_xticks(range(len(names)))
        a2.set_xticklabels(["Boltzmann", "+ offset"][:len(names)], fontsize=7)
    a2.set_ylabel("AIC (lower = better)", fontsize=7)
    d = (mc or {}).get("delta_aic")
    a2.set_title("ΔAIC = " + (f"{d:.1f}" if d is not None and np.isfinite(d) else "–"), fontsize=9)
    fig.tight_layout()
    return fig


@_dejavu
def _validation_figure(T_meas, T_ref):
    """Residual plot: the ratio thermometer minus the independent reference, against the reference temperature.
    The 1 K and 3 K bands stay visible whatever the size of the deviation (a plot of T against T would hide them)."""
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    r = np.asarray(T_ref, float)
    d = np.asarray(T_meas, float) - r
    ymax = max(4.0, 1.25 * float(np.max(np.abs(d))))
    xlo, xhi = float(r.min()) - 5, float(r.max()) + 5
    ax.fill_between([xlo, xhi], -3, 3, color="#FDEBD0", label="±3 K (warn)", zorder=0)
    ax.fill_between([xlo, xhi], -1, 1, color="#D5F5E3", label="±1 K (ok)", zorder=1)
    ax.axhline(0, color="black", lw=0.8, zorder=2)
    ax.plot(r, d, "o", ms=4.5, color=NAVY, label="ratio − reference", zorder=3)
    ax.axhline(d.mean(), color=RED, lw=1.0, ls="--", label=f"bias {d.mean():+.2f} K", zorder=2)
    ax.set_xlim(xlo, xhi); ax.set_ylim(-ymax, ymax)
    ax.set_xlabel("reference temperature (K)"); ax.set_ylabel("deviation (K)")
    ax.set_title(f"RMS {np.sqrt(np.mean(d ** 2)):.2f} K, largest |Δ| {np.max(np.abs(d)):.2f} K over {len(d)} points", fontsize=9)
    ax.legend(fontsize=6.5, frameon=False, loc="upper right", ncol=2)
    fig.tight_layout()
    return fig


def build_pdf(**kw):
    """Return the PDF as bytes. The report uses DejaVu Sans locally (it has the superscripts and Greek letters the
    numbers need); the global matplotlib font of the program is not changed."""
    with plt.rc_context({"font.family": "DejaVu Sans"}):
        return _build_pdf(**kw)


def _methods(pg, *, band_info, corrected, mc, refs, T_pred, T_ref_ok, offs, counts):
    """Section 2: the analysis chain in plain words, the settings that were used, and the figure of the chain."""
    pg.text("2. Methods", size=10, bold=True)
    pg.text("The series of spectra was recorded at known temperatures. Every spectrum was processed in the same way:", size=8)
    steps = [
        "(1) the spectral response correction was " + ("applied to every spectrum." if corrected else "NOT applied (none was given)."),
        "(2) the two bands were measured " + (f"in the windows {band_info.get('windows')}" if band_info.get("windows")
                                            else (f"as fitted pseudo-Voigt bands {band_info.get('fit')}" if band_info.get("fit")
                                                  else "by the method chosen in the program")) + ";",
        "(3) the ratio Δ = I₁ / I₂ was formed for each spectrum; the fit uses the whole series at once (a global fit);",
        "(4) ln Δ was fitted against 1/T with the Boltzmann law" + (", with an added constant offset C" if offs else "")
        + " (a weighted straight line in ln Δ against 1/T);",
        "(5) every ratio was converted to a temperature with that calibration; the relative sensitivity S_r = 100·ΔE/(k_B·T²) "
        "and the temperature uncertainty σT = (σ_Δ/Δ)/(S_r/100) were computed from the uncertainty budget;",
        "(6) the quality checks of section 9 were run on the data and on the fit.",
    ]
    for s_ in steps:
        pg.text(s_, size=7.8, gap=0.004)
    if counts:
        pg.text("Counting statistics were used for the weights (Poisson errors of the counts).", size=7.5, gap=0.004)
    pg.figure(_flow_figure(), h=0.22)
    pg.text("Figure 1. The analysis chain. Each box is one step of the report; the checks at the end decide the verdict.",
            size=7, color="#444444", gap=0.006)


def _model_section(pg, *, mc, refs, T_pred, T, D, D_stat):
    """Section 3: which law describes the data (AIC) and, when a reference is given, how well the temperatures agree with it."""
    pg.text("3. Model choice and validation", size=10, bold=True)
    if not mc or mc.get("status") != "ok":
        pg.text("The model comparison could not be made (too few temperatures).", size=8)
    else:
        pg.text("Two laws were compared: the simple Boltzmann law and the same law with a constant offset C (a background that does not "
                "depend on temperature). The model with the lower AIC is preferred; an AIC difference above 2 is taken as a real "
                "preference, above 10 as a strong one. Result: " + mc.get("text", "") + ".", size=8)
        pg.figure(_model_figure(T, D, D_stat, mc), h=0.27)
        pg.text("Figure 2. Left: the same data with both candidate laws. Right: their AIC (lower is better).", size=7, color="#444444", gap=0.006)
    if refs is not None and len(refs) > 0:
        pg.text("Validation: the temperatures from the ratio calibration were compared with an independent reference thermometer measured "
                "at the same time as each spectrum. Bias = mean(T_ratio − T_ref); RMS = the root of the mean squared difference.", size=8)
        d = np.asarray(T_pred, float) - np.asarray(refs, float)
        pg.text(f"Result: bias {d.mean():+.2f} K, RMS {np.sqrt(np.mean(d ** 2)):.2f} K, largest |Δ| {np.max(np.abs(d)):.2f} K over {len(d)} points.",
                size=8)
        pg.figure(_validation_figure(T_pred, refs), h=0.25)
        pg.text("Figure 3. Deviation of the ratio-thermometer temperature from the reference, against the reference temperature; the bands are ±1 K and ±3 K.",
                size=7, color="#444444", gap=0.006)
    else:
        pg.text("Validation: no independent reference temperatures were given. The absolute accuracy of the temperatures is therefore "
                "NOT shown in this report; the precision (σT) is.", size=8, color="#7B241C")


def _build_pdf(*, sample, version, T, D, D_stat, bf, Sr, table, budget_rows, checks, verdict, settings, u_dE_cm1,
               T_res, T_cal, sigma_scaled=False, model=None, refs=None, T_pred=None, corrected=False, band_info=None):
    buf = io.BytesIO()
    refs = None if refs is None or len(refs) == 0 else np.asarray(refs, float)
    T_pred = np.asarray(T if T_pred is None else T_pred, float)
    with PdfPages(buf) as pdf:
        pg = _Page(pdf, f"Photoluminescence thermometry report – {sample}")
        pg.text(f"Generated {datetime.datetime.now():%Y-%m-%d %H:%M} by Reticula {version}. "
                "This report describes one analysis of one series of spectra; every number below is the output of that analysis.", size=8)
        dE = abs(bf["dE_cm1"])
        pg.text("1. Verdict", size=10, bold=True)
        pg.text(f"Confidence: {verdict}. The energy gap is ΔE = {dE:.0f} ± {u_dE_cm1:.0f} cm⁻¹ (combined standard uncertainty, all budget "
                f"components included; coverage factor k = 1). The Boltzmann law explains the ratio with R² = {bf['R2']:.4f}.")
        pg.text(f"Sensitivity at the lowest temperature: S_r = {Sr[0]:.3f} % K⁻¹. Temperature resolution at the lowest temperature: "
                f"σT,res = {T_res[0]:.2f} K; calibration uncertainty σT,cal = {T_cal[0]:.2f} K.")
        if sigma_scaled:
            pg.text("Note: the counting error was too small for the scatter of the ratio, so the uncertainty uses the measured scatter "
                    "(see the 'Ratio noise' check).", size=8)
        _methods(pg, band_info=band_info or {}, corrected=corrected, mc=model, refs=refs, T_pred=T_pred,
                 T_ref_ok=(refs is not None and len(refs) > 0), offs=bool((settings or {}).get("Boltzmann offset C")), counts=bool((settings or {}).get("counts (Poisson weights)")))
        _model_section(pg, mc=model, refs=refs, T_pred=T_pred, T=T, D=D, D_stat=D_stat)
        pg.text("4. Settings used", size=10, bold=True)
        for k, v in (settings or {}).items():
            pg.text(f"{k}: {v}", size=8, gap=0.004)
        pg.text("5. Equations", size=10, bold=True)
        pg.math(r"\Delta = \alpha\, e^{-\Delta E / k_B T}\quad\Leftrightarrow\quad \ln\Delta = \ln\alpha - \frac{\Delta E}{k_B}\,\frac{1}{T}")
        pg.math(r"S_r = \frac{100\,\Delta E}{k_B T^2}\ \ (\% \mathrm{K}^{-1}),\qquad \sigma_{T,\mathrm{res}} = \frac{\sigma_\Delta/\Delta}{S_r/100}")
        pg.text("The fit is a weighted straight line of ln Δ against 1/T (weights from the errors of Δ). The Boltzmann fit uses the "
                "parameters α and ΔE; the uncertainty of ΔE is the square root of its variance in the covariance matrix, scaled by √χ² "
                "only when the scatter is significant (Birge ratio).", size=8)
        pg.text("6. Calibration", size=10, bold=True)
        fig = _calibration_figure(T, D, D_stat, bf)
        pg.figure(fig, h=0.30)
        plt.close(fig)
        pg.text("7. Thermometry table", size=10, bold=True)
        if table is not None and len(table):
            cols = list(table.columns)
            pg.text("  ".join(str(c) for c in cols)[:200], size=6.5, gap=0.002, width=150)
            for _, row in table.head(40).iterrows():
                pg.text("  ".join(f"{row[c]:.4g}" if isinstance(row[c], (int, float, np.floating, np.integer)) else str(row[c]) for c in cols)[:200],
                        size=6.5, gap=0.001, width=150)
        pg.text("8. Uncertainty budget (ΔE)", size=10, bold=True)
        for r in budget_rows or []:
            pg.text(f"{r.get('component', r.get('name', ''))}: {r.get('value', '')}", size=7.5, gap=0.002)
        pg.text("9. Checks and advice", size=10, bold=True)
        for status, title, detail, advice in _checks_text(checks):
            pg.text(f"[{status}] {title}: {detail}", size=7.5, gap=0.002)
            if advice:
                pg.text(f"    advice: {advice}", size=7, gap=0.004, color="#444444")
        pg.text("10. Limitations", size=10, bold=True)
        pg.text("The report states what the program computed and checked. Points marked as general practice in the program are advice "
                "the program does not check. A high confidence verdict does not replace a reference temperature measurement on your own "
                "instrument. Band windows, background choices and the thermometer's temperature reading are the largest sources of error; "
                "see the budget.", size=8)
        pg.flush()
    return buf.getvalue()
