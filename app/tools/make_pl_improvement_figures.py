"""v2.18.6: makes the PNG figures of the PL improvement section of the info document (mfs/assets/docs/pl_*.png).

Every number in the figures comes from the synthetic tests of this release (known answer): the ratio fit (subtraction + ln-fit
against the Poisson maximum likelihood) and the background methods (rubber band, SNIP, arPLS, emission-masked polynomial).
Run from the app folder:  python tools/make_pl_improvement_figures.py"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
sys.path.insert(0, APP)
from mfs import plmle as PM                       # noqa: E402
from mfs import plbackground as BG                # noqa: E402

OUT = os.path.join(APP, "mfs", "assets", "docs")
KB = 8.617333262e-5
EV_CM = 8065.544
BLUE, RED, GREY, GREEN = "#1F4E79", "#C0392B", "#7F8C8D", "#0F766E"


def ratio_benchmark(n=300, dE=700.0, alpha=1.5, s2=2e4, B=50.0, seed=3):
    rng = np.random.default_rng(seed)
    T = np.arange(80.0, 401.0, 20.0)
    d_sub, e_sub, d_mle, e_mle = [], [], [], []
    for _ in range(n):
        D = alpha * np.exp(-dE / EV_CM / (KB * T))
        N1 = rng.poisson(B + D * s2).astype(float)
        N2 = rng.poisson(B + s2 * np.ones_like(T)).astype(float)
        S1, S2 = N1 - B, N2 - B
        ok = (S1 > 0) & (S2 > 0)
        lnD = np.log(S1[ok] / S2[ok]); sig = np.sqrt(1 / S1[ok] + 1 / S2[ok]); x = 1 / T[ok]; w = 1 / sig ** 2
        W = w.sum(); Sx = (w * x).sum(); Sy = (w * lnD).sum(); Sxx = (w * x * x).sum(); Sxy = (w * x * lnD).sum()
        b = (W * Sxy - Sx * Sy) / (W * Sxx - Sx ** 2)
        d_sub.append(-b * KB * EV_CM); e_sub.append(np.sqrt(W / (W * Sxx - Sx ** 2)) * KB * EV_CM)
        r = PM.fit(T, N1, N2, B, B)
        d_mle.append(r["dE_cm1"]); e_mle.append(r["sigma_cm1"])
    out = {}
    for name, d, e in (("subtraction", d_sub, e_sub), ("maximum likelihood", d_mle, e_mle)):
        d = np.array(d); e = np.array(e); z = (d - dE) / e
        out[name] = {"bias": 100 * (d.mean() - dE) / dE, "in1": 100 * np.mean(np.abs(z) <= 1), "in2": 100 * np.mean(np.abs(z) <= 2),
                     "z": z}
    return out


def background_benchmark(n=40, seed=5):
    rng = np.random.default_rng(seed)
    x = np.arange(500.0, 800.0, 0.5)

    def g(c, w, a):
        return a * np.exp(-0.5 * ((x - c) / w) ** 2)
    cases = {
        "A: sharp peaks\non a slope": (200 + 0.8 * (x - 500), g(560, 3, 800) + g(640, 3, 500)),
        "B: broad band\nover the range": (200 + 0.8 * (x - 500), g(650, 40, 1000)),
        "C: high noise\n(few counts)": (20.0 + 0.02 * (x - 500), g(650, 15, 60)),
        "D: concave\nbackground": (300 - 0.006 * (x - 650) ** 2, g(600, 4, 400)),
    }
    methods = {"rubber band": lambda y: BG.rubberband(x, y),
               "masked polynomial": lambda y: BG.masked_poly(x, y, max_order=1)[0]}
    res = {}
    for cname, (bg, pk) in cases.items():
        true = np.trapezoid(pk, x)
        res[cname] = {}
        for mname, fn in methods.items():
            v = []
            for _ in range(n):
                y = rng.poisson(bg + pk).astype(float)
                v.append(np.trapezoid(y - fn(y), x) / true - 1)
            v = np.array(v)
            res[cname][mname] = (100 * v.mean(), 100 * v.std(ddof=1))
    return res


def make():
    os.makedirs(OUT, exist_ok=True)
    rb = ratio_benchmark()
    # Figure 1: ratio fit, coverage and bias
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.0, 3.2), gridspec_kw={"width_ratios": [1.0, 1.2]})
    names = list(rb.keys())
    x = np.arange(2)
    in1 = [rb[k]["in1"] for k in names]; in2 = [rb[k]["in2"] for k in names]
    w = 0.36
    a1.bar(x - w / 2, in1, w, color=GREY, label="within 1σ")
    a1.bar(x + w / 2, in2, w, color=BLUE, label="within 2σ")
    a1.axhline(68.3, color=GREY, ls="--", lw=0.8); a1.axhline(95.4, color=BLUE, ls="--", lw=0.8)
    a1.set_xticks(x, ["subtraction\n+ ln-fit", "Poisson\nmaximum likelihood"], fontsize=7.5)
    a1.set_ylabel("% of 300 runs"); a1.set_ylim(0, 105); a1.legend(fontsize=7, frameon=False, loc="lower right")
    a1.set_title("Coverage of the 1σ and 2σ intervals", fontsize=8.5)
    bias = [rb[k]["bias"] for k in names]
    a2.bar(x, bias, 0.5, color=[GREY, GREEN])
    a2.axhline(0, color="black", lw=0.6)
    for i, b in enumerate(bias):
        a2.text(i, b / 2 if abs(b) > 0.5 else b + 0.12, f"{b:+.2f} %", ha="center", fontsize=8, color="white" if abs(b) > 0.5 else "black")
    a2.set_xticks(x, ["subtraction\n+ ln-fit", "Poisson\nmaximum likelihood"], fontsize=7.5)
    a2.set_ylabel("bias of ΔE (%)"); a2.set_title("Bias of the energy gap (true 700 cm⁻¹)", fontsize=8.5)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "pl_fig_ratio.png"), dpi=200); plt.close(fig)

    # Figure 2: background methods, area error per scenario
    bb = background_benchmark()
    fig, ax = plt.subplots(figsize=(8.0, 3.3))
    cases = list(bb.keys())
    xs = np.arange(len(cases))
    for j, (mname, col) in enumerate((("rubber band", GREY), ("masked polynomial", GREEN))):
        means = [bb[c][mname][0] for c in cases]; sds = [bb[c][mname][1] for c in cases]
        ax.bar(xs + (j - 0.5) * 0.36, means, 0.36, yerr=sds, color=col, capsize=2, label=mname)
        for i, m in enumerate(means):
            ax.text(xs[i] + (j - 0.5) * 0.36, m + (4 if m >= 0 else -10), f"{m:+.0f}%", ha="center", fontsize=7)
    ax.axhline(0, color="black", lw=0.6); ax.axhspan(-2, 2, color=GREEN, alpha=0.08)
    ax.set_xticks(xs, cases, fontsize=7.5); ax.set_ylabel("area error (%)")
    ax.set_title("Background methods: area error of the band (shaded: ±2 %)", fontsize=8.5)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "pl_fig_background.png"), dpi=200); plt.close(fig)
    return rb, bb


if __name__ == "__main__":
    rb, bb = make()
    for k in rb:
        print(f"{k}: bias {rb[k]['bias']:+.2f}%  within 1σ {rb[k]['in1']:.0f}%  within 2σ {rb[k]['in2']:.0f}%")
    for c in bb:
        print(c.replace("\n", " "), {m: f"{v[0]:+.1f}% (sd {v[1]:.1f})" for m, v in bb[c].items()})
