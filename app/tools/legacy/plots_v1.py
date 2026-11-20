"""Plotting functions. Every `draw_*` draws into a given Axes (so the figure composer can reuse it);
every `fig_*` returns a complete journal-styled Figure."""
from __future__ import annotations
import math
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, NullLocator
from matplotlib.patches import Rectangle
from . import style
from .style import colours_for, new_figure, SINGLE, DOUBLE, panel_label, no_y_ticks

LABELS = {
    "pxrd": (r"2$\theta$ (°)", "Intensity (a.u.)"),
    "ftir_T": (r"Wavenumber (cm$^{-1}$)", "Transmittance (a.u.)"),
    "ftir_A": (r"Wavenumber (cm$^{-1}$)", "Absorbance (a.u.)"),
    "raman": (r"Raman shift (cm$^{-1}$)", "Intensity (a.u.)"),
    "nmr": (r"$\delta$ (ppm)", ""),
    "iso": (r"Relative pressure, $p/p_0$", r"Quantity adsorbed (cm$^3$ g$^{-1}$ STP)"),
    "psd": ("Pore width (nm)", r"d$V$/d$w$ (cm$^3$ g$^{-1}$ nm$^{-1}$)"),
    "tga": ("Temperature (°C)", "Weight (%)"),
    "dtg": ("Temperature (°C)", r"DTG (% °C$^{-1}$)"),
    "xps": ("Binding energy (eV)", "Intensity (a.u.)"),
}


def _finish(ax, xlabel, ylabel, invert=False, xlim=None, ylim=None):
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    if invert and not ax.xaxis_inverted():
        ax.invert_xaxis()


# --------------------------------------------------------------------------- stacked spectra
def draw_stacked(ax, traces, kind="pxrd", offset=1.1, label_side="right", invert=None, xlim=None,
                 reference=None, sticks=None, highlight=None, lw=0.9, label_fs=None, colors=None,
                 ylabel=None, show_y=False, annotations=None):
    """Offset/stacked spectra. traces: list of (name, x, y); first item is drawn at the bottom.
    Each trace is min-max scaled to 1 inside the plotted window and shifted by `offset`.
    Labels sit above the trace at the right-hand end of the visible window."""
    if invert is None:
        invert = kind in ("ftir_T", "ftir_A", "nmr", "xps")
    cols = list(colors) if colors else colours_for([t[0] for t in traces])
    label_fs = label_fs or plt.rcParams["font.size"] - 1
    entries = list(traces)
    if reference is not None:
        entries = [reference] + entries
        cols = ["#7F7F7F"] + cols
    allx = np.concatenate([np.asarray(t[1], float) for t in entries])
    lo, hi = (min(xlim), max(xlim)) if xlim else (np.nanmin(allx), np.nanmax(allx))
    span = hi - lo
    # label anchor: right end of the visible window (high x if normal, low x if inverted axis)
    x_lab = (hi - 0.02 * span) if not invert else (lo + 0.02 * span)
    seg_lo, seg_hi = (lo + 0.6 * span, hi) if not invert else (lo, lo + 0.4 * span)
    base = 0.0
    for k, (nm, x, y) in enumerate(entries):
        x = np.asarray(x, float); y = np.asarray(y, float)
        m = (x >= lo) & (x <= hi) & np.isfinite(y)
        if not m.any():
            continue
        xv, yv = x[m], y[m]
        yy = (yv - np.nanmin(yv)) / ((np.nanmax(yv) - np.nanmin(yv)) or 1)
        yfull = (y - np.nanmin(yv)) / ((np.nanmax(yv) - np.nanmin(yv)) or 1)
        yfull = np.where(m, yfull, np.nan)
        ax.plot(x, yfull + base, color=cols[k], lw=lw, solid_joinstyle="round", solid_capstyle="round")
        seg = (xv >= seg_lo) & (xv <= seg_hi)
        ytop = np.nanmax(yy[seg]) if seg.any() else 0.3
        ax.text(x_lab, base + min(max(ytop, 0.12), 0.85) + 0.05, nm, color=cols[k], fontsize=label_fs,
                ha="right", va="bottom")
        base += offset
    top = base - offset + 1.12
    if sticks is not None and len(sticks):
        st = sticks[(sticks["two_theta"] >= lo) & (sticks["two_theta"] <= hi)]
        ax.vlines(st["two_theta"], -0.14, -0.14 + 0.12 * st["I"] / st["I"].max(), color="#7F7F7F", lw=0.6)
    if highlight:
        top += 0.18
        for (a, b, lab) in highlight:
            ax.axvspan(a, b, color="#EBEBEB", lw=0, zorder=0)
            if lab:
                ax.text((a + b) / 2, top - 0.02, lab, ha="center", va="top", fontsize=label_fs - 0.5, color="#555555")
    if annotations:
        top += 0.2
        for (xa, lab) in annotations:
            ax.text(xa, top - 0.04, lab, ha="center", va="top", fontsize=label_fs - 0.5, color="#444444")
            ax.plot([xa, xa], [top - 0.22, top - 0.27], color="#888888", lw=0.5)
    xl, yl = LABELS.get(kind, ("x", "y"))
    ax.set_xlabel(xl)
    lab_y = ylabel if ylabel is not None else yl
    if lab_y:
        ax.set_ylabel(lab_y)
    ax.set_xlim((hi, lo) if invert else (lo, hi))
    ax.set_ylim(-0.17 if sticks is not None else -0.06, top + 0.04)
    if not show_y:
        no_y_ticks(ax)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    return ax


def _fix_labels_inverted(ax):
    pass


def fig_stacked(traces, kind="pxrd", width=SINGLE, aspect=0.95, **kw):
    fig, ax = new_figure(width, aspect)
    draw_stacked(ax, traces, kind, **kw)
    return fig


# --------------------------------------------------------------------------- overlay (non-offset) curves
def draw_overlay(ax, traces, xlabel, ylabel, invert=False, xlim=None, ylim=None, legend=True, lw=1.0, markers=False, logx=False):
    cols = colours_for([t[0] for t in traces])
    for c, (nm, x, y) in zip(cols, traces):
        ax.plot(x, y, color=c, lw=lw, label=nm, marker="o" if markers else None, ms=2.5)
    if logx:
        ax.set_xscale("log")
    _finish(ax, xlabel, ylabel, invert, xlim, ylim)
    if legend:
        ax.legend(loc="best")
    return ax


# --------------------------------------------------------------------------- N2 sorption
def draw_isotherms(ax, isos, logx=False, legend=True, ms=3.0, ylim=None):
    """isos: list of (name, DataFrame[p_rel, q, branch])"""
    cols = colours_for([n for n, _ in isos])
    for c, (nm, df) in zip(cols, isos):
        a = df[df.branch == "ads"].sort_values("p_rel")
        d = df[df.branch == "des"].sort_values("p_rel")
        ax.plot(a.p_rel, a.q, "-", color=c, lw=0.8)
        ax.plot(a.p_rel, a.q, "o", color=c, ms=ms, mfc=c, mec=c, mew=0.6, label=nm)
        if len(d):
            ax.plot(d.p_rel, d.q, "-", color=c, lw=0.8)
            ax.plot(d.p_rel, d.q, "o", color=c, ms=ms, mfc="white", mec=c, mew=0.6)
    if logx:
        ax.set_xscale("log")
        ax.set_xlim(left=max(1e-6, min(df.p_rel[df.p_rel > 0].min() for _, df in isos) / 2), right=1.0)
    else:
        ax.set_xlim(0, 1.0)
        ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.set_ylim(bottom=0, top=ylim)
    _finish(ax, *LABELS["iso"])
    if legend:
        from matplotlib.lines import Line2D
        h, l = ax.get_legend_handles_labels()
        # matplotlib fills legend columns top-to-bottom: put the longest labels in the left column so the
        # right column never runs past the axes frame (e.g. "Sample C")
        vis = [len(re.sub(r"[$_{}\\]|mathrm", "", t)) for t in l]
        nleft = math.ceil(len(l) / 2)
        longest = set(sorted(range(len(l)), key=lambda i: -vis[i])[:nleft])
        order = [i for i in range(len(l)) if i in longest] + [i for i in range(len(l)) if i not in longest]
        h, l = [h[i] for i in order], [l[i] for i in order]
        ymax = max(df.q.max() for _, df in isos)
        ax.set_ylim(0, ylim or ymax * (1.18 + 0.075 * math.ceil(len(isos) / 2)))
        leg1 = ax.legend(h, l, loc="upper left", ncol=2, fontsize=plt.rcParams["font.size"] - 1.5, handletextpad=0.2,
                         columnspacing=0.6, borderaxespad=0.3)
        ax.add_artist(leg1)
        hb = [Line2D([], [], marker="o", ls="", color="k", ms=ms, label="adsorption"),
              Line2D([], [], marker="o", ls="", mfc="white", mec="k", ms=ms, label="desorption")]
        ax.legend(handles=hb, loc="lower right", fontsize=plt.rcParams["font.size"] - 1.5, handletextpad=0.3)
    return ax


def draw_psd(ax, psds, width_unit="nm", kind="dvdw", legend=False, xlim=(0.5, 5), offset=0.0, fill=True):
    cols = colours_for([n for n, _ in psds])
    base = 0
    for c, (nm, df) in zip(cols, psds):
        w = df["w"].values / (10 if df.attrs.get("width_unit", width_unit) == "Å" else 1)
        y = df[kind].values if kind in df else df.iloc[:, 1].values
        if kind == "dvdw" and df.attrs.get("width_unit", width_unit) == "Å":
            y = y * 10  # per Å -> per nm
        ax.plot(w, y + base, color=c, lw=1.0, label=nm)
        if fill:
            ax.fill_between(w, base, y + base, color=c, alpha=0.12, lw=0)
        base += offset
    ylab = LABELS["psd"][1] if kind == "dvdw" else (r"d$V$/dlog($w$) (cm$^3$ g$^{-1}$)" if kind == "dvdlogw" else r"Cumulative pore volume (cm$^3$ g$^{-1}$)")
    _finish(ax, LABELS["psd"][0], ylab, xlim=xlim)
    ax.set_ylim(bottom=0)
    if legend:
        ax.legend()
    return ax


def draw_bet(ax, bet, name="", color="k"):
    p, y = bet["plot_p"], bet["plot_y"]
    sel = (p >= bet["p_lo"]) & (p <= bet["p_hi"])
    ax.plot(p[~sel], y[~sel], "o", mfc="white", mec="#999999", ms=3)
    ax.plot(p[sel], y[sel], "o", color=color, ms=3)
    xx = np.linspace(bet["p_lo"], bet["p_hi"], 20)
    ax.plot(xx, bet["slope"] * xx + bet["intercept"], "-", color=color, lw=0.8)
    ax.set_xlabel(r"$p/p_0$")
    ax.set_ylabel(r"$(p/p_0)\,/\,[n(1-p/p_0)]$ (g cm$^{-3}$)")
    from matplotlib.ticker import ScalarFormatter
    for axis in (ax.xaxis, ax.yaxis):
        f = ScalarFormatter(useMathText=True); f.set_powerlimits((-2, 3)); axis.set_major_formatter(f)
    ax.text(0.04, 0.96, f"{name}\n$S_\\mathrm{{BET}}$ = {bet['area']:.0f} m$^2$ g$^{{-1}}$\n$C$ = {bet['C']:.0f}, $R^2$ = {bet['R2']:.4f}\n"
            f"$p/p_0$ {bet['p_lo']:.3g}–{bet['p_hi']:.3g}", transform=ax.transAxes, va="top", fontsize=plt.rcParams["font.size"] - 1)
    return ax


def draw_rouquerol(ax, bet, color="k"):
    ax.plot(bet["rq_p"], bet["rq"], "o-", color=color, ms=2.5, lw=0.7)
    ax.axvspan(bet["p_lo"], bet["p_hi"], color="#E0E0E0", lw=0, zorder=0)
    ax.set_xscale("log")
    ax.set_xlabel(r"$p/p_0$"); ax.set_ylabel(r"$n(1-p/p_0)$ (cm$^3$ g$^{-1}$)")
    return ax


# --------------------------------------------------------------------------- TGA
def draw_tga(ax, tgas, dtg_ax=None, normalise_end=False, xlim=None, legend=True, annotate=None):
    """tgas: list of (name, analysis-dict from analysis.tga_analysis)."""
    cols = colours_for([n for n, _ in tgas])
    for c, (nm, a) in zip(cols, tgas):
        w = a["w"] / a["W_end"] * 100 if normalise_end else a["w"]
        ax.plot(a["T"], w, color=c, lw=1.0, label=nm)
        if dtg_ax is not None:
            dtg_ax.plot(a["T"], a["dtg"], color=c, lw=0.8, ls="--")
    _finish(ax, LABELS["tga"][0], "Weight (% of residue)" if normalise_end else LABELS["tga"][1], xlim=xlim or (min(a["T"].min() for _, a in tgas), max(a["T"].max() for _, a in tgas)))
    if annotate:
        for (T0, T1, txt) in annotate:
            ax.axvspan(T0, T1, color="#EEEEEE", lw=0, zorder=0)
            ax.text((T0 + T1) / 2, ax.get_ylim()[0] + 0.03 * np.diff(ax.get_ylim())[0], txt, ha="center", va="bottom",
                    fontsize=plt.rcParams["font.size"] - 1.5, color="#555555")
    if dtg_ax is not None:
        dtg_ax.set_ylabel(LABELS["dtg"][1])
        dtg_ax.set_ylim(bottom=min(0, dtg_ax.get_ylim()[0]))
    if legend:
        ax.legend(loc="lower left")
    return ax


def fig_tga(tgas, width=SINGLE, dtg="panel", **kw):
    if dtg == "panel":
        style.apply_style()
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(width, width * 1.15), sharex=True, gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.06})
        draw_tga(a1, tgas, **kw)
        cols = colours_for([n for n, _ in tgas])
        for c, (nm, a) in zip(cols, tgas):
            a2.plot(a["T"], a["dtg"], color=c, lw=0.9)
        a1.set_xlabel("")
        a2.set_xlabel(LABELS["dtg"][0]); a2.set_ylabel(LABELS["dtg"][1])
        a2.set_ylim(bottom=0)
        return fig
    fig, ax = new_figure(width, 0.8)
    tw = ax.twinx() if dtg == "twin" else None
    if tw is not None:
        tw.tick_params(direction="in")
        ax.tick_params(right=False, which="both")
    draw_tga(ax, tgas, tw, **kw)
    return fig


# --------------------------------------------------------------------------- kinetics
def draw_kinetics(ax, sets, ylabel=r"$C/C_0$", time_label="Time (h)", show_fit=True, normalise=True, legend=True, annotate_k=True, k_unit="h$^{-1}$"):
    """sets: list of (name, t, c, c_err|None, fit|None)"""
    cols = colours_for([s[0] for s in sets])
    markers = ["o", "s", "^", "D", "v", "p", "h"]
    for k, (c, (nm, t, cc, err, fit)) in enumerate(zip(cols, sets)):
        c0 = fit["C0"] if (normalise and fit) else (cc[0] if normalise else 1)
        y = np.asarray(cc) / c0
        e = (np.asarray(err) / c0) if err is not None else None
        lab = nm
        if fit and annotate_k:
            lab = f"{nm} ($k$ = {fit['k']:.3g} ± {fit['k_err']:.1g} {k_unit})"
        ax.errorbar(t, y, yerr=e, fmt=markers[k % len(markers)], color=c, ms=3.5, mfc=c if k % 2 == 0 else "white", mec=c,
                    mew=0.8, elinewidth=0.6, capsize=1.8, label=lab, zorder=3)
        if fit and show_fit:
            tt = np.linspace(0, max(t) * 1.03, 200)
            ax.plot(tt, fit["model"](tt) / c0, "-", color=c, lw=0.9, zorder=2)
    ax.set_xlim(left=0)
    if normalise:
        ax.set_ylim(0, 1.42 if legend else 1.08)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    _finish(ax, time_label, ylabel)
    if legend:
        ax.legend(loc="upper center", ncol=2, fontsize=plt.rcParams["font.size"] - 1.8, handletextpad=0.3,
                  columnspacing=0.8, borderaxespad=0.25)
    return ax


def draw_k_bars(ax, names, k, kerr, ylabel=r"$k_\mathrm{obs}$ (10$^{-3}$ h$^{-1}$)", scale=1e3):
    cols = colours_for(names)
    xs = np.arange(len(names))
    ax.bar(xs, np.asarray(k) * scale, yerr=np.asarray(kerr) * scale, color=cols, width=0.65, edgecolor="black", linewidth=0.6,
           error_kw=dict(lw=0.7, capsize=2.5, capthick=0.7))
    ax.set_xticks(xs, names, rotation=30, ha="right")
    ax.tick_params(axis="x", which="minor", bottom=False, top=False)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)
    return ax


def draw_langmuir(ax, sets, xlabel=r"$C_\mathrm{e}$ (mM)", ylabel=r"$q_\mathrm{e}$ (µmol µmol$^{-1}_\mathrm{MOF}$)"):
    cols = colours_for([s[0] for s in sets])
    for c, (nm, ce, qe, fit) in zip(cols, sets):
        lab = nm + (f" ($q_\\mathrm{{max}}$ = {fit['q_max']:.2f}, $K$ = {fit['K']:.2g} mM$^{{-1}}$)" if fit else "")
        ax.plot(ce, qe, "o", color=c, ms=3.5, label=lab)
        if fit:
            xx = np.linspace(0, max(ce) * 1.05, 200)
            ax.plot(xx, fit["model"](xx), "-", color=c, lw=0.9)
    ax.set_xlim(left=0); ax.set_ylim(bottom=0)
    _finish(ax, xlabel, ylabel)
    ax.legend(loc="lower right", fontsize=plt.rcParams["font.size"] - 1.5)
    return ax


# --------------------------------------------------------------------------- bars (ICP / XRF / zeta)
def draw_grouped_bars(ax, df, value="mean", err="sd", group="sample", series="element", ylabel="wt%",
                      colors_by_group=True, legend=True, hatch_series=True):
    groups = list(dict.fromkeys(df[group]))
    sers = list(dict.fromkeys(df[series]))
    n = len(sers)
    w = 0.8 / n
    gcols = colours_for(groups)
    hatches = ["", "////", "....", "xxxx", "\\\\\\\\"]
    for j, s in enumerate(sers):
        sub = df[df[series] == s].set_index(group).reindex(groups)
        xs = np.arange(len(groups)) + (j - (n - 1) / 2) * w
        if colors_by_group and n <= 1:
            col = gcols
            ax.bar(xs, sub[value], w * 0.95, yerr=sub[err] if err in sub else None, color=col, edgecolor="black", lw=0.6,
                   error_kw=dict(lw=0.7, capsize=2.2, capthick=0.7), label=s)
        else:
            col = style.OKABE_ITO[(j + 1) % len(style.OKABE_ITO)]
            ax.bar(xs, sub[value], w * 0.95, yerr=sub[err] if err in sub else None, color=col, edgecolor="black", lw=0.6,
                   hatch=hatches[j % len(hatches)] if hatch_series and False else None,
                   error_kw=dict(lw=0.7, capsize=2.2, capthick=0.7), label=s)
    ax.set_xticks(np.arange(len(groups)), groups, rotation=25, ha="right")
    ax.tick_params(axis="x", which="minor", bottom=False, top=False)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)
    if legend and n > 1:
        ax.legend(loc="upper left", ncol=min(n, 3))
    return ax


# --------------------------------------------------------------------------- SEM
def nice_bar(length_nm_max):
    nice = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]
    cand = [n for n in nice if n <= length_nm_max]
    return cand[-1] if cand else nice[0]


import matplotlib.patheffects as pe
_OUTLINE = [pe.withStroke(linewidth=1.6, foreground="black")]


def draw_sem(ax, img, pixel_nm=None, crop_row=None, bar_nm=None, label=None, bar_color="white", loc="lower right", fs=None):
    im = img[:crop_row] if crop_row else img
    ax.imshow(im, cmap="gray", interpolation="lanczos", aspect="equal")
    ax.set_xticks([]); ax.set_yticks([])
    ax.tick_params(which="both", top=False, right=False, left=False, bottom=False)
    for s in ax.spines.values():
        s.set_linewidth(0.5)
    h, w = im.shape[:2]
    fs = fs or plt.rcParams["font.size"]
    if pixel_nm:
        L = bar_nm or nice_bar(0.22 * w * pixel_nm)
        px = L / pixel_nm
        txt = f"{L/1000:g} µm" if L >= 1000 else f"{L:g} nm"
        x1 = w * 0.95; x0 = x1 - px; y0 = h * 0.93
        ax.add_patch(Rectangle((x0, y0), px, h * 0.018, color=bar_color, lw=0))
        ax.text((x0 + x1) / 2, y0 - h * 0.015, txt, color=bar_color, ha="center", va="bottom", fontsize=fs, fontweight="bold",
                path_effects=_OUTLINE)
    if label:
        ax.text(0.03, 0.97, label, transform=ax.transAxes, color="white", va="top", ha="left", fontsize=fs + 1, fontweight="bold",
                path_effects=_OUTLINE)
    return ax


def fig_sem_panel(images, ncols=2, width=DOUBLE, labels=None, gap=0.02):
    """images: list of dict(image, pixel_nm, crop_row, title)"""
    n = len(images); nrows = math.ceil(n / ncols)
    ar = np.mean([(im["image"][:im.get("crop_row") or None].shape[0] / im["image"].shape[1]) for im in images])
    style.apply_style()
    fig, axs = plt.subplots(nrows, ncols, figsize=(width, width / ncols * ar * nrows * 1.0 + 0.1))
    axs = np.atleast_1d(axs).ravel()
    for k, ax in enumerate(axs):
        if k >= n:
            ax.axis("off"); continue
        d = images[k]
        lab = (labels[k] if labels else f"({chr(97+k)})")
        draw_sem(ax, d["image"], d.get("pixel_nm"), d.get("crop_row"), d.get("bar_nm"), lab)
        if d.get("title"):
            ax.text(0.97, 0.97, d["title"], transform=ax.transAxes, color="white", ha="right", va="top",
                    fontsize=plt.rcParams["font.size"], path_effects=_OUTLINE)
    fig.subplots_adjust(wspace=gap, hspace=gap, left=0, right=1, top=1, bottom=0)
    return fig


def draw_edx_map(ax, img, cmap_color="#D55E00", pixel_nm=None, label=None, element=None):
    from matplotlib.colors import LinearSegmentedColormap
    cm = LinearSegmentedColormap.from_list("el", ["#000000", cmap_color])
    ax.imshow(img, cmap=cm, interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    ax.tick_params(which="both", top=False, right=False, left=False, bottom=False)
    if element:
        ax.text(0.97, 0.97, element, transform=ax.transAxes, color="white", ha="right", va="top", fontweight="bold")
    if pixel_nm:
        draw_sem(ax, np.zeros((1, 1)), None)
    return ax


# --------------------------------------------------------------------------- SDS-PAGE
def draw_bands(ax, df, time_col="time", value_col="band_pct", group_col="band_label", time_label="Time (h)", ylabel="Band intensity (%)"):
    groups = list(dict.fromkeys(df[group_col]))
    cols = [style.OKABE_ITO[(i + 1) % 9] for i in range(len(groups))]
    mk = ["o", "s", "^", "D", "v"]
    for k, (g, c) in enumerate(zip(groups, cols)):
        s = df[df[group_col] == g].groupby(time_col)[value_col].agg(["mean", "std"]).reset_index()
        ax.errorbar(s[time_col], s["mean"], yerr=s["std"].fillna(0), fmt=mk[k % 5] + "-", color=c, ms=3.5, lw=0.9,
                    elinewidth=0.6, capsize=2, label=str(g))
    ax.set_xlim(left=0); ax.set_ylim(0, None)
    _finish(ax, time_label, ylabel)
    ax.legend(loc="best", fontsize=plt.rcParams["font.size"] - 1.5)
    return ax


# --------------------------------------------------------------------------- XPS
def shirley(x, y, tol=1e-6, maxit=50):
    x = np.asarray(x, float); y = np.asarray(y, float)
    o = np.argsort(x); x, y = x[o], y[o]
    yl, yr = y[0], y[-1]
    B = np.full_like(y, yl)
    for _ in range(maxit):
        Q = np.cumsum((y - B)[::-1])[::-1]  # integral from x to high BE
        Bn = yl + (yr - yl) * (Q[0] - Q) / (Q[0] if Q[0] else 1)
        if np.max(np.abs(Bn - B)) < tol * max(abs(yr), 1):
            B = Bn; break
        B = Bn
    inv = np.empty_like(o); inv[o] = np.arange(len(o))
    return B[inv]


def draw_xps(ax, traces, region=None, background=True, offset=0.9):
    cols = colours_for([t[0] for t in traces])
    base = 0
    for c, (nm, x, y) in zip(cols, traces):
        x = np.asarray(x); y = np.asarray(y)
        if region:
            m = (x >= region[0]) & (x <= region[1]); x, y = x[m], y[m]
        b = shirley(x, y) if background else np.zeros_like(y)
        s = (y - b).max() or 1
        ax.plot(x, (y - b) / s + base, color=c, lw=0.9)
        if len(traces) > 1:
            ax.text(x.max() - 0.02 * (x.max() - x.min()), base + 0.06, nm, color=c, ha="left", va="bottom", fontsize=plt.rcParams["font.size"] - 1.5)
        base += offset * 1.25
    _finish(ax, *LABELS["xps"], invert=True)
    ax.set_ylim(-0.08, base + 0.1)
    no_y_ticks(ax)
    return ax


# --------------------------------------------------------------------------- defect comparison
def fig_defect_comparison(res, width=DOUBLE):
    """Side-by-side dot chart: L, x, m, Cl and OS per method with the consensus band."""
    ms = [m for m in res["methods"] if m.status in ("ok", "warn")]
    cons = res["consensus"]
    style.apply_style()
    fig, axs = plt.subplots(1, 4, figsize=(width, max(1.9, 0.26 * len(ms) + 0.9)), sharey=True,
                            gridspec_kw={"wspace": 0.12})
    names = [SHORT.get(m.key, m.title[:30]) for m in ms]
    y = np.arange(len(ms))[::-1]
    quant = [("L", "Linkers per node, $L$"), ("x", "Missing linkers, $x$"), ("m", "Caps per node, $m$"), ("OS", "Open sites per node")]
    for ax, (key, lab) in zip(axs, quant):
        vals = [getattr(m, key) for m in ms]
        for yi, v, m in zip(y, vals, ms):
            if v is None:
                continue
            col = "#7F7F7F" if m.semi_quant else ("#D55E00" if m.status == "warn" else "#0072B2")
            ax.plot(v, yi, "o", color=col, ms=4.5, mfc=col if not m.semi_quant else "white", mew=0.9, zorder=3)
        if cons and key in cons:
            ax.axvline(cons[key], color="k", lw=0.8, ls="--", zorder=1)
            if key == "L":
                ax.axvspan(cons["L_min"], cons["L_max"], color="#DDDDDD", lw=0, zorder=0)
            if key == "x":
                ax.axvspan(res["input"].L0() - cons["L_max"], res["input"].L0() - cons["L_min"], color="#DDDDDD", lw=0, zorder=0)
        ax.set_xlabel(lab)
        ax.tick_params(axis="y", which="both", left=False, right=False)
        ax.xaxis.set_major_locator(MaxNLocator(3))
        ax.margins(x=0.14)
        ax.yaxis.set_minor_locator(NullLocator())
    axs[0].set_yticks(y, names)
    axs[0].set_ylim(-0.7, len(ms) - 0.3)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color="#0072B2", label="quantitative"),
         Line2D([], [], marker="o", ls="", color="#D55E00", label="warning"),
         Line2D([], [], marker="o", ls="", mfc="white", color="#7F7F7F", label="semi-quantitative"),
         Line2D([], [], ls="--", color="k", label="consensus (median)")]
    fig.legend(handles=h, loc="upper center", ncol=4, bbox_to_anchor=(0.55, 1.02 + 0.02 * (8 - len(ms)) / 8), frameon=False)
    return fig


SHORT = {"a": "(a) TGA plateau M-400", "b1": "(b) TGA M-200 + NMR", "b2": "(b') missing-cluster model",
         "c": "(c) TGA plateau + NMR", "d": "(d) qNMR + TGA + Cl (De Roo)", "e1": "(e) ICP : qNMR ratio",
         "e2": "(e') ICP mass balance", "f": "(f) pore volume (semi-q.)", "g": "(g) acid-base titration"}


def _short(title):
    t = title.split(")", 1)[-1].strip()
    t = t.replace("Experimental minimal formula: ", "")
    return (t[:38] + "…") if len(t) > 39 else t
