"""v2.19.5: makes the PNG figures of the 'faster manifold fit' section of the PL info document (mfs/assets/docs/pl_fig_speed_*.png).

Every number comes from tools/pl_speed_2195_data.json, measured on one computer with the real v2.19.2 / v2.19.4 code and the final 2.19.5
code (full page pipeline: line search, both shape models, robustness refits with other windows, the other model and two backgrounds),
plus a Monte Carlo of the ratio errors and the stuck-fit example (known answers: the scatter over noisy copies, the true ratio).
Run from the app folder:  python tools/make_pl_speed_figures.py"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
OUT = os.path.join(APP, "mfs", "assets", "docs")
DATA = os.path.join(HERE, "pl_speed_2195_data.json")
BLUE, RED, GREY, GREEN, LIGHT = "#1F4E79", "#C0392B", "#7F8C8D", "#0F766E", "#C8CED3"
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
                     "axes.spines.top": False, "axes.spines.right": False})

ORDER = ["test_series_6lines", "Er_example", "YbHo_composite", "heavy_13x800_9lines", "test_series_41spectra"]
NAMES = {"test_series_6lines": "Test series", "Er_example": "Er³⁺ example", "YbHo_composite": "Yb/Ho composite",
         "heavy_13x800_9lines": "Crowded series", "test_series_41spectra": "Long series"}
STEPS = ["v2.19.4", "s1_exactJ", "s2_warm", "s3_mixed", "s4_bounds", "s5_round2", "s6_round3", "s7_round4", "s8_round5"]


def label(d, k):
    m = d["datasets"][k]
    return f"{NAMES[k]}\n{m['n_spectra']} spectra × {m['n_points']} points, {m['n_lines']} lines"


def fmt(t):
    return f"{t:.2f} s" if t < 10 else f"{t:.0f} s"


def fig_versions(d):
    T = d["timings"]
    fig, ax = plt.subplots(figsize=(8.0, 4.0))
    y = np.arange(len(ORDER))[::-1]
    h = 0.26
    for j, (ver, col, name) in enumerate((("v2.19.2", LIGHT, "v2.19.2"), ("v2.19.4", GREY, "v2.19.4 (before)"), ("s8_round5", GREEN, "v2.19.5"))):
        vals = [T[ver][k] for k in ORDER]
        ax.barh(y + (1 - j) * h, vals, h, color=col, label=name)
        for yy, v in zip(y + (1 - j) * h, vals):
            ax.text(v * 1.08, yy, fmt(v), va="center", fontsize=6.8, color="black" if j < 2 else GREEN)
    ax.text(1.0, 1.0, "speed-up\nagainst v2.19.4", transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color=GREEN)
    for yy, k in zip(y, ORDER):
        f = T["v2.19.4"][k] / T["s8_round5"][k]
        ax.text(1.0, (yy - y.min() + 0.5) / (len(ORDER)), f"{f:.0f}×", transform=ax.transAxes, ha="right", va="center", fontsize=10,
                color=GREEN, fontweight="bold")
    ax.set_xscale("log")
    ax.set_yticks(y, [label(d, k) for k in ORDER], fontsize=7)
    ax.set_xlim(0.05, 3000)
    ax.set_ylim(y.min() - 0.5, y.max() + 0.5)
    ax.set_xlabel("time of the whole manifold analysis (s, log scale)")
    ax.set_title("Manifold fit: time per data set, same computer, same settings", loc="left")
    ax.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "pl_fig_speed_versions.png"), dpi=200)
    plt.close(fig)


def fig_steps(d):
    T = d["timings"]
    fig, ax = plt.subplots(figsize=(8.0, 3.3))
    x = np.arange(len(STEPS))
    for k, col, mk in (("heavy_13x800_9lines", BLUE, "o"), ("test_series_41spectra", RED, "s"), ("test_series_6lines", GREY, "^")):
        v = [T[s][k] for s in STEPS]
        ax.plot(x, v, marker=mk, color=col, lw=1.3, ms=4, label=NAMES[k] + f" ({d['datasets'][k]['n_spectra']} spectra)")
        ax.text(x[-1] + 0.15, v[-1], fmt(v[-1]), color=col, fontsize=7, va="center")
        ax.text(x[0] - 0.15, v[0], fmt(v[0]), color=col, fontsize=7, va="center", ha="right")
    ax.set_yscale("log")
    ax.set_xticks(x, ["start\n(v2.19.4)", "1", "2", "3", "4", "5", "6", "7", "8"])
    ax.set_xlim(-0.9, len(STEPS) - 0.2)
    ax.set_ylabel("time (s, log scale)")
    ax.set_xlabel("optimisation step (see the list in the text)")
    ax.set_title("Effect of every optimisation step (cumulative)", loc="left")
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "pl_fig_speed_steps.png"), dpi=200)
    plt.close(fig)


def fig_stuck(d):
    s = d["stuck"]
    x, y = np.array(s["x"]), np.array(s["y"])
    fo, fn = np.array(s["fit_old"]), np.array(s["fit_new"])
    fig = plt.figure(figsize=(8.0, 3.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.2, 1], width_ratios=[1.6, 1], hspace=0.08, wspace=0.28)
    a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[1, 0], sharex=a); c = fig.add_subplot(gs[:, 1])
    a.plot(x, y, ".", ms=1.8, color="black", label="data (spectrum at 293 K)")
    a.plot(x, fo, color=RED, lw=1.0, label=f"before: stuck, rss {s['rss_old']:.0f}")
    a.plot(x, fn, color=GREEN, lw=1.0, label=f"2.19.5: rss {s['rss_new']:.0f}")
    a.legend(fontsize=6.8, frameon=False, loc="upper left")
    a.set_ylabel("counts")
    a.set_title("Test series, window 500–575 nm, shared shape", loc="left")
    plt.setp(a.get_xticklabels(), visible=False)
    b.plot(x, y - fo, color=RED, lw=0.7)
    b.plot(x, y - fn, color=GREEN, lw=0.7)
    b.axhline(0, color="black", lw=0.5)
    b.set_xlabel("wavelength (nm)"); b.set_ylabel("residual")
    T = np.array(s["T"])
    c.plot(T, s["truth"], "k-", lw=1, label="true ratio")
    c.plot(T, s["ratio_old"], "o", color=RED, ms=3.5, mfc="none", label="before")
    c.plot(T, s["ratio_new"], "s", color=GREEN, ms=3, label="2.19.5")
    c.set_xlabel("T (K)"); c.set_ylabel("manifold ratio I₁/I₂")
    c.legend(fontsize=6.8, frameon=False)
    c.set_title("Ratio of every spectrum", loc="left")
    fig.savefig(os.path.join(OUT, "pl_fig_speed_stuck.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_errors(d):
    fig, axs = plt.subplots(1, 3, figsize=(8.0, 2.9), gridspec_kw={"width_ratios": [1, 1.15, 1.15]})
    m = d["mc_gauss"]
    T = 293 + 10 * np.arange(len(m["mc_sd"]))
    sc = d.get("stuck_copies", {})
    a0 = axs[0]
    a0.plot(T, 1e4 * np.array(m["mc_sd_old"]), "o-", color=RED, ms=3.5, lw=1, label="before 2.19.5")
    a0.plot(T, 1e4 * np.array(m["mc_sd"]), "s-", color=GREEN, ms=3.2, lw=1.2, label="2.19.5")
    a0.set_ylim(0, None)
    a0.set_xlabel("T (K)"); a0.set_ylabel("true scatter of the ratio (×10⁻⁴)")
    a0.set_title("A  Scatter of the ratio", loc="left")
    if sc:
        a0.text(0.35, 0.30, f"old fit stuck in\n{sc['old_stuck']} of {sc['n']} copies", transform=a0.transAxes, fontsize=6.8, color=RED)
    a0.legend(fontsize=6.6, frameon=False, loc="upper left")
    for ax, key, title in ((axs[1], "mc_gauss", "B  Gaussian (η on its bound 0)"), (axs[2], "mc_pv", "C  Pseudo-Voigt (η = 0.4)")):
        mm = d[key]
        sd = np.array(mm["mc_sd"])
        p = mm["rel_precision_of_sd"]
        ax.axhspan(1 - 2 * p, 1 + 2 * p, color=GREEN, alpha=0.10, lw=0)
        ax.axhline(1, color="black", lw=0.6)
        ax.plot(T, np.array(mm["err_old"]) / np.array(mm["mc_sd_old"]), "o-", color=RED, ms=3.5, lw=1, label="before 2.19.5")
        ax.plot(T, np.array(mm["err_all_params"]) / sd, "^-", color=GREY, ms=3.5, lw=1, label="exact Jacobian, all parameters")
        ax.plot(T, np.array(mm["err_new"]) / sd, "s-", color=GREEN, ms=3.2, lw=1.2, label="2.19.5")
        ax.set_title(f"{title}, {mm['n']} copies", loc="left", fontsize=7.8)
        ax.set_xlabel("T (K)")
        ax.set_yscale("log")
        ax.set_ylim(0.45, 6)
        ax.set_yticks([0.5, 0.7, 1, 1.5, 2, 3, 5], ["0.5", "0.7", "1", "1.5", "2", "3", "5"])
        ax.minorticks_off()
    axs[1].set_ylabel("reported error / true scatter")
    axs[2].legend(fontsize=6.4, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "pl_fig_speed_errors.png"), dpi=200)
    plt.close(fig)


def fig_table(d):
    T, R = d["timings"], d["results"]
    rows = []
    for k in ORDER:
        m, r = d["datasets"][k], R[k]
        rows.append([NAMES[k], f"{m['n_spectra']} × {m['n_points']}", str(m["n_lines"]), fmt(T["v2.19.2"][k]), fmt(T["v2.19.4"][k]), fmt(T["s8_round5"][k]),
                     f"{T['v2.19.4'][k] / T['s8_round5'][k]:.0f}×", f"{r['dE_old']:.1f} / {r['dE_new']:.1f}",
                     f"{r['model_old'].replace(' shape', '')} / {r['model_new'].replace(' shape', '')}"])
    cols = ["data set", "spectra ×\npoints", "lines", "v2.19.2", "v2.19.4", "v2.19.5", "speed-up\n(v2.19.4)", "ΔE (cm⁻¹)\nold / new", "shape model\nold / new"]
    fig, ax = plt.subplots(figsize=(8.0, 1.9))
    ax.axis("off")
    tb = ax.table(cellText=rows, colLabels=cols, loc="center", cellLoc="center", colLoc="center")
    tb.auto_set_font_size(False)
    tb.set_fontsize(7)
    tb.scale(1, 1.5)
    widths = [0.15, 0.09, 0.05, 0.08, 0.08, 0.08, 0.08, 0.14, 0.17]
    for (i, j), cell in tb.get_celld().items():
        cell.set_width(widths[j])
    for (i, j), cell in tb.get_celld().items():
        cell.set_edgecolor("#D5DBDF")
        if i == 0:
            cell.set_facecolor("#EAF1F4"); cell.set_text_props(fontweight="bold", color=BLUE)
        if j == 5 and i > 0:
            cell.set_text_props(color=GREEN, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "pl_fig_speed_table.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def make():
    d = json.load(open(DATA))
    os.makedirs(OUT, exist_ok=True)
    fig_versions(d)
    fig_steps(d)
    fig_stuck(d)
    fig_errors(d)
    fig_table(d)
    return d


if __name__ == "__main__":
    d = make()
    T = d["timings"]
    for k in ORDER:
        print(f"{k}: v2.19.2 {T['v2.19.2'][k]:.2f} s, v2.19.4 {T['v2.19.4'][k]:.2f} s, v2.19.5 {T['s8_round5'][k]:.2f} s")
