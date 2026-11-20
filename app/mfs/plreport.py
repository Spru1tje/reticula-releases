"""v2.18.6: analysis report (PDF) of the Photoluminescence thermometry tab.

One document per analysis: verdict, the settings used, the calibration equations and their numbers, the figure of the
Boltzmann calibration, the thermometry table, the uncertainty budget, every quality check with its advice, and the
limitations. Built with the same page helper as the defect report (mfs/report.py), so the look matches."""
from __future__ import annotations

import datetime
import io

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from mfs.report import _Page

KB_EV = 8.617333262e-5
EV_CM = 8065.544


def _checks_text(checks):
    rows = []
    for c in checks or []:
        rows.append((str(c.get("status", "")).upper(), str(c.get("title", "")), str(c.get("detail", "")), str(c.get("advice", "") or "")))
    return rows


def _calibration_figure(T, D, D_stat, bf, offs=False):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.0, 3.3), gridspec_kw={"width_ratios": [1.3, 1.0]})
    x = 1000.0 / np.asarray(T, float)
    y = np.log(np.asarray(D, float))
    err = np.asarray(D_stat, float) / np.asarray(D, float)
    a1.errorbar(x, y, yerr=err, fmt="o", ms=3.5, color="#1F4E79", ecolor="#7F9FBF", capsize=2, label="data (ln Δ)")
    xf = np.linspace(x.min(), x.max(), 200)
    Tf = 1000.0 / xf
    a1.plot(xf, np.log(bf["f"](Tf)), color="#C0392B", lw=1.2, label="Boltzmann fit")
    a1.set_xlabel("1000 / T (K⁻¹)"); a1.set_ylabel("ln Δ"); a1.legend(fontsize=7, frameon=False)
    a1.set_title("Boltzmann calibration", fontsize=9)
    res = (np.asarray(D, float) - bf["f"](np.asarray(T, float))) / np.asarray(D, float)
    a2.axhline(0, color="black", lw=0.6)
    a2.plot(np.asarray(T, float), 100 * res, "o", ms=3.5, color="#1F4E79")
    a2.set_xlabel("T (K)"); a2.set_ylabel("residual (%)"); a2.set_title("Relative residuals", fontsize=9)
    fig.tight_layout()
    return fig


def build_pdf(**kw):
    """Return the PDF as bytes. The report uses DejaVu Sans locally (it has the superscripts and Greek letters the
    numbers need); the global matplotlib font of the program is not changed."""
    with plt.rc_context({"font.family": "DejaVu Sans"}):
        return _build_pdf(**kw)


def _build_pdf(*, sample, version, T, D, D_stat, bf, Sr, table, budget_rows, checks, verdict, settings, u_dE_cm1,
               T_res, T_cal, sigma_scaled=False):
    buf = io.BytesIO()
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
        pg.text("2. Settings used", size=10, bold=True)
        for k, v in (settings or {}).items():
            pg.text(f"{k}: {v}", size=8, gap=0.004)
        pg.text("3. Equations", size=10, bold=True)
        pg.math(r"\Delta = \alpha\, e^{-\Delta E / k_B T}\quad\Leftrightarrow\quad \ln\Delta = \ln\alpha - \frac{\Delta E}{k_B}\,\frac{1}{T}")
        pg.math(r"S_r = \frac{100\,\Delta E}{k_B T^2}\ \ (\% \mathrm{K}^{-1}),\qquad \sigma_{T,\mathrm{res}} = \frac{\sigma_\Delta/\Delta}{S_r/100}")
        pg.text("The fit is a weighted straight line of ln Δ against 1/T (weights from the errors of Δ). The Boltzmann fit uses the "
                "parameters α and ΔE; the uncertainty of ΔE is the square root of its variance in the covariance matrix, scaled by √χ² "
                "only when the scatter is significant (Birge ratio).", size=8)
        pg.text("4. Calibration", size=10, bold=True)
        fig = _calibration_figure(T, D, D_stat, bf)
        pg.figure(fig, h=0.30)
        plt.close(fig)
        pg.text("5. Thermometry table", size=10, bold=True)
        if table is not None and len(table):
            cols = list(table.columns)
            pg.text("  ".join(str(c) for c in cols)[:200], size=6.5, gap=0.002, width=150)
            for _, row in table.head(40).iterrows():
                pg.text("  ".join(f"{row[c]:.4g}" if isinstance(row[c], (int, float, np.floating, np.integer)) else str(row[c]) for c in cols)[:200],
                        size=6.5, gap=0.001, width=150)
        pg.text("6. Uncertainty budget (ΔE)", size=10, bold=True)
        for r in budget_rows or []:
            pg.text(f"{r.get('component', r.get('name', ''))}: {r.get('value', '')}", size=7.5, gap=0.002)
        pg.text("7. Checks and advice", size=10, bold=True)
        for status, title, detail, advice in _checks_text(checks):
            pg.text(f"[{status}] {title}: {detail}", size=7.5, gap=0.002)
            if advice:
                pg.text(f"    advice: {advice}", size=7, gap=0.004, color="#444444")
        pg.text("8. Limitations", size=10, bold=True)
        pg.text("The report states what the program computed and checked. Points marked as general practice in the program are advice "
                "the program does not check. A high confidence verdict does not replace a reference temperature measurement on your own "
                "instrument. Band windows, background choices and the thermometer's temperature reading are the largest sources of error; "
                "see the budget.", size=8)
        pg.flush()
    return buf.getvalue()
