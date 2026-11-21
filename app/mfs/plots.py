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
                 ylabel=None, show_y=False, annotations=None, scale="each"):
    """Offset/stacked spectra. scale='each': every trace min-max scaled to 1 in the window; 'common': one factor for
    all traces (keeps relative intensities, e.g. after normalising to a reference band). traces: list of (name, x, y); first item is drawn at the bottom.
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
    rng_c = None
    if scale == "common":
        rs = []
        for (_, x, y) in entries:
            x = np.asarray(x, float); y = np.asarray(y, float); m = (x >= lo) & (x <= hi) & np.isfinite(y)
            if m.any():
                rs.append(np.nanmax(y[m]) - np.nanmin(y[m]))
        rng_c = float(np.median(rs)) if rs else 1.0      # median: one very strong trace does not flatten the others
    ax._mfs_stack = []
    for k, (nm, x, y) in enumerate(entries):
        x = np.asarray(x, float); y = np.asarray(y, float)
        m = (x >= lo) & (x <= hi) & np.isfinite(y)
        if not m.any():
            continue
        xv, yv = x[m], y[m]
        r_ = rng_c if rng_c else ((np.nanmax(yv) - np.nanmin(yv)) or 1)
        if invert and kind == "ftir_T" and rng_c:      # transmittance: align the tops (baseline 100 %)
            y0 = np.nanmax(yv) - r_
        else:
            y0 = np.nanmin(yv)
        yy = (yv - y0) / r_
        yfull = (y - y0) / r_
        ax._mfs_stack.append((nm, y0, r_, base, k))
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
_ENDS = []
ISO_DEFAULT_STYLE = {"marker": "o", "ms": 3.0, "ls": "-", "lw": 0.8, "des_marker": None, "des_fill": "white"}


def _thin(p, logx, n=34):
    """Indices of points to mark: at most ~n markers spread evenly along p/p0 (or log p/p0)."""
    p = np.asarray(p, float)
    if len(p) <= n:
        return np.arange(len(p))
    x = np.log10(np.clip(p, 1e-9, None)) if logx else p
    span = (x.max() - x.min()) or 1.0
    keep, last = [], -np.inf
    for i in np.argsort(x):
        if x[i] - last >= span / n:
            keep.append(i); last = x[i]
    keep += [int(np.argmax(x))]
    return np.array(sorted(set(keep)))


def _free_inset(ax, isos, logx, w=0.42, h=0.36, avoid=None):
    """Inset position (axes fraction) with the fewest isotherm points/lines under it."""
    y0, y1 = ax.get_ylim()
    pts = []
    for _, df in isos:
        p = np.asarray(df.p_rel, float); q = np.asarray(df.q, float)
        o = np.argsort(p); p, q = p[o], q[o]
        if logx:
            lo = np.log10(max(1e-7, p[p > 0].min() / 2)); xf = (np.log10(np.clip(p, 1e-9, None)) - lo) / (0 - lo)
        else:
            xf = p
        gx = np.linspace(0, 1, 120)
        pts.append(np.c_[gx, (np.interp(gx, xf, q) - y0) / (y1 - y0)])
    P_ = np.vstack(pts) if pts else np.zeros((0, 2))
    best, bc = None, None
    for x0 in (0.52, 0.36, 0.2, 0.08):
        for yb in (0.08, 0.3, 0.55):
            box = (x0, yb, w, h)
            if avoid is not None and not (x0 + w <= avoid[0] or x0 >= avoid[0] + avoid[2] or yb + h <= avoid[1] or yb >= avoid[1] + avoid[3]):
                continue
            inside = ((P_[:, 0] >= x0 - 0.02) & (P_[:, 0] <= x0 + w + 0.02) & (P_[:, 1] >= yb - 0.03) & (P_[:, 1] <= yb + h + 0.03)).sum()
            if bc is None or inside < bc:
                best, bc = box, inside
    return list(best)


def _end_labels(ax, ends, fs):
    """Curve names right of the plot at their last uptake, pushed apart so they never overlap (thin leader if moved)."""
    if not ends:
        return
    y0, y1 = ax.get_ylim()
    minsep = (y1 - y0) * 0.065
    ends = sorted(ends, key=lambda e: e[1])
    pos = [e[1] for e in ends]
    for _ in range(200):
        moved = False
        for k in range(1, len(pos)):
            if pos[k] - pos[k - 1] < minsep - 1e-12:
                mid = (pos[k] + pos[k - 1]) / 2; pos[k - 1] = mid - minsep / 2; pos[k] = mid + minsep / 2; moved = True
        if not moved:
            break
    for (nm, y, c), yp in zip(ends, pos):
        moved = abs(yp - y) > minsep * 0.15
        ax.annotate(nm, xy=(1.0, y), xycoords=("axes fraction", "data"), xytext=(1.06 if moved else 1.02, yp),
                    textcoords=("axes fraction", "data"), ha="left", va="center", fontsize=fs, color=c, annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color=c, lw=0.5, shrinkA=0, shrinkB=0) if moved else None)


def draw_isotherms(ax, isos, logx=False, legend=True, ms=3.0, ylim=None, styles=None, sort=False,
                   branch_key=True, legend_loc="lower right", connect=True, hysteresis=False, direct_labels=False,
                   inset_psd=None, inset_kind="dvdw", inset_xlim=(0.3, 2.0), thin=False):
    """N2 isotherms in the usual MOF style: adsorption = filled symbols, desorption = open symbols, lines joining
    the points in MEASUREMENT order (ads up to p/p0 max, then des back down), p/p0 0-1 (or log for micropore filling).
    styles: {sample: {ads_color, des_color, marker, des_marker, ms, ls, lw}} (any key optional).
    One single legend: one entry per sample + two key entries (filled = adsorption, open = desorption)."""
    from matplotlib.lines import Line2D
    styles = styles or {}
    cols = colours_for([n for n, _ in isos])
    global _ENDS
    _ENDS = []
    for c, (nm, df) in zip(cols, isos):
        st = {**ISO_DEFAULT_STYLE, "ms": ms, **(styles.get(nm) or {})}
        ca = st.get("ads_color") or c
        cd = st.get("des_color") or ca
        mk, mkd = st["marker"], st.get("des_marker") or st["marker"]
        d0 = df.sort_values("order") if ("order" in df and not sort) else df.sort_values("p_rel")
        a = d0[d0.branch == "ads"]
        d = d0[d0.branch == "des"]
        if sort:
            d = d.sort_values("p_rel", ascending=False)
        ls = st["ls"] if connect else "none"
        me = _thin(a.p_rel.values, logx) if thin else None
        ax.plot(a.p_rel, a.q, ls=ls, marker=mk, color=ca, ms=st["ms"], mfc=ca, mec=ca, mew=0.6, lw=st["lw"], label=nm,
                markevery=None if me is None else list(me))
        if len(d):
            # start the desorption line at the last adsorption point so the loop is closed (as in MicroActive)
            dp = np.r_[a.p_rel.values[-1:], d.p_rel.values] if len(a) else d.p_rel.values
            dq = np.r_[a.q.values[-1:], d.q.values] if len(a) else d.q.values
            ax.plot(dp, dq, ls=ls, color=cd, lw=st["lw"], marker="none", label="_nolegend_")
            md = _thin(d.p_rel.values, logx) if thin else None
            ax.plot(d.p_rel, d.q, ls="none", marker=mkd, color=cd, ms=st["ms"], mfc=st.get("des_fill") or "white",
                    mec=cd, mew=0.6, label="_nolegend_", markevery=None if md is None else list(md))
            if hysteresis and len(a) > 3 and len(d) > 3:
                pa, qa = a.p_rel.values, a.q.values
                o = np.argsort(d.p_rel.values); pd_, qd = d.p_rel.values[o], d.q.values[o]
                lo_, hi_ = max(pa.min(), pd_.min()), min(pa.max(), pd_.max())
                grid = np.linspace(lo_, hi_, 300)
                oa = np.argsort(pa)
                qa_i = np.interp(grid, pa[oa], qa[oa]); qd_i = np.interp(grid, pd_, qd)
                gap = qd_i - qa_i
                big = gap > 0.015 * max(qa.max(), 1e-9)
                if big.any():
                    ax.fill_between(grid, qa_i, qd_i, where=big, color=ca, alpha=0.13, lw=0, zorder=0, label="_nolegend_")
        ends_append = (nm, float(a.q.values[np.argmax(a.p_rel.values)]) if len(a) else float(d0.q.values[-1]), ca)
        _ENDS.append(ends_append)
    if logx:
        ax.set_xscale("log")
        ax.set_xlim(left=max(1e-7, min(df.p_rel[df.p_rel > 0].min() for _, df in isos) / 2), right=1.0)
    else:
        ax.set_xlim(0, 1.0)
        ax.xaxis.set_major_locator(MaxNLocator(5))
    ymax = max(df.q.max() for _, df in isos) if isos else 1
    ax.set_ylim(0, ylim or ymax * 1.08)
    _finish(ax, *LABELS["iso"])
    if logx:
        ax.set_xlabel(r"Relative pressure, $p/p_0$ (log scale)")
    if branch_key and any((df.branch == "des").any() for _, df in isos):
        # key entries as invisible proxies so they live in the SAME legend and can be renamed/hidden
        ax.plot([], [], ls="none", marker="o", color="#404040", ms=ms, label="adsorption")
        ax.plot([], [], ls="none", marker="o", mfc="white", mec="#404040", color="#404040", ms=ms, label="desorption")
    if inset_psd:
        keyed = direct_labels and branch_key and any((df.branch == "des").any() for _, df in isos)
        ia = ax.inset_axes(_free_inset(ax, isos, logx, avoid=(0.6, 0.0, 0.4, 0.22) if keyed else None))
        draw_psd(ia, inset_psd, kind=inset_kind, xlim=inset_xlim, legend=False, fill=True)
        ia.set_xlabel(ia.get_xlabel(), fontsize=plt.rcParams["font.size"] - 2.5, labelpad=1)
        ia.set_ylabel(r"d$V$/d$w$" if inset_kind == "dvdw" else ia.get_ylabel(), fontsize=plt.rcParams["font.size"] - 2.5, labelpad=1)
        ia.tick_params(labelsize=plt.rcParams["font.size"] - 3, length=2, pad=1)
        ia.set_yticks([])
    if direct_labels:
        fs = plt.rcParams["font.size"] - 1.5
        _end_labels(ax, list(_ENDS), fs)
        legend = False
        if branch_key and any((df.branch == "des").any() for _, df in isos):
            h, l = ax.get_legend_handles_labels()
            keep = [(hh, ll) for hh, ll in zip(h, l) if ll in ("adsorption", "desorption")]
            if keep:
                ax.legend([k[0] for k in keep], [k[1] for k in keep], loc="lower right", fontsize=fs, handletextpad=0.3,
                          borderaxespad=0.4, labelspacing=0.25, frameon=False)
    if legend:
        h, l = ax.get_legend_handles_labels()
        nsamp = len(isos)
        fs = plt.rcParams["font.size"] - 1.5
        if nsamp <= 3 and legend_loc == "lower right":
            # Type I isotherms leave the lower-right corner empty -> no extra headroom needed
            ax.legend(loc="lower right", fontsize=fs, handletextpad=0.3, borderaxespad=0.4, labelspacing=0.25)
        else:
            # many samples: two columns in head-room above the curves, longest labels in the left column
            vis = [len(re.sub(r"[$_{}\\]|mathrm", "", t)) for t in l]
            nleft = math.ceil(len(l) / 2)
            longest = set(sorted(range(len(l)), key=lambda i: -vis[i])[:nleft])
            order = [i for i in range(len(l)) if i in longest] + [i for i in range(len(l)) if i not in longest]
            h, l = [h[i] for i in order], [l[i] for i in order]
            ax.set_ylim(0, ylim or ymax * (1.12 + 0.085 * nleft))
            ax.legend(h, l, loc="upper left", ncol=2, fontsize=fs, handletextpad=0.3, columnspacing=0.6,
                      borderaxespad=0.3, labelspacing=0.25)
    return ax


def draw_psd(ax, psds, width_unit="nm", kind="dvdw", legend=False, xlim=(0.5, 5), offset=0.0, fill=True, cumulative=False):
    """PSD curves; cumulative=True adds the cumulative pore volume (dashed) on a twin right axis."""
    cols = colours_for([n for n, _ in psds])
    base = 0
    ax2 = ax.twinx() if cumulative and kind != "cum" else None
    for c, (nm, df) in zip(cols, psds):
        w = df["w"].values / (10 if df.attrs.get("width_unit", width_unit) == "Å" else 1)
        y = df[kind].values if kind in df else df.iloc[:, 1].values
        if kind == "dvdw" and df.attrs.get("width_unit", width_unit) == "Å":
            y = y * 10  # per Å -> per nm
        ax.plot(w, y + base, color=c, lw=1.0, label=nm)
        if fill:
            ax.fill_between(w, base, y + base, color=c, alpha=0.12, lw=0)
        if ax2 is not None:
            if "cum" in df:
                cv = df["cum"].values
            else:
                o = np.argsort(w); dv = df["dvdw"].values[o] * (10 if df.attrs.get("width_unit", width_unit) == "Å" else 1)
                cv = np.empty(len(w)); cv[o] = np.r_[0, np.cumsum(0.5 * (dv[1:] + dv[:-1]) * np.diff(w[o]))]
            ax2.plot(w, cv, color=c, lw=0.9, ls="--", label="_nolegend_")
        base += offset
    if ax2 is not None:
        ax2.set_ylabel(r"Cumulative pore volume (cm$^3$ g$^{-1}$, dashed)")
        ax2.set_ylim(bottom=0)
        ax2.spines["right"].set_visible(True)
    ylab = LABELS["psd"][1] if kind == "dvdw" else (r"d$V$/dlog($w$) (cm$^3$ g$^{-1}$)" if kind == "dvdlogw" else r"Cumulative pore volume (cm$^3$ g$^{-1}$)")
    _finish(ax, LABELS["psd"][0], ylab, xlim=xlim)
    ax.set_ylim(bottom=0)
    if legend:
        ax.legend(loc="center right" if ax2 is not None else "best", framealpha=0.92, edgecolor="none")
    return ax


def _pow10(v):
    v = abs(float(v)) or 1.0
    return int(np.floor(np.log10(v)))


def draw_bet(ax, bet, name="", color="k", overlay=None, box=True, show_excluded=True):
    """Linearised BET plot (v2.13 redesign, t94u): the fitted Rouquerol range shaded, fitted points filled, the other
    points of the plotted window open, fit line, and an annotation box (S_BET, C, R², range, n_m). y in 10^k units in the
    axis label (no offset text). overlay: [(name, bet, colour)] → several samples in one panel (annotation as a table)."""
    from matplotlib.ticker import MaxNLocator
    items = overlay or [(name, bet, color)]
    if overlay:
        return _draw_bet_overlay(ax, items, box)
    ymax_all = 0.0
    for nm, b, col in items:
        p, y = np.asarray(b["plot_p"]), np.asarray(b["plot_y"])
        sel = (p >= b["p_lo"]) & (p <= b["p_hi"])
        win = p <= min(0.35, b["p_hi"] * 1.35)
        ymax_all = max(ymax_all, float(np.max(y[win])) if win.any() else float(np.max(y)))
    k = _pow10(ymax_all)
    sc = 10.0 ** k
    for nm, b, col in items:
        p, y = np.asarray(b["plot_p"]), np.asarray(b["plot_y"])
        sel = (p >= b["p_lo"]) & (p <= b["p_hi"])
        win = p <= min(0.35, b["p_hi"] * 1.35)
        if show_excluded:
            ax.plot(p[win & ~sel], y[win & ~sel] / sc, "o", mfc="white", mec="#9AA5B1", ms=3, mew=0.6, zorder=2)
        ax.plot(p[sel], y[sel] / sc, "o", color=col, ms=3.4, mew=0.6, zorder=3, label=nm if overlay else None)
        xx = np.linspace(0, b["p_hi"] * 1.08, 40)
        ax.plot(xx, (b["slope"] * xx + b["intercept"]) / sc, "-", color=col, lw=0.9, zorder=2)
        if not overlay:
            ax.axvspan(b["p_lo"], b["p_hi"], color=col if col not in ("k", "#000000") else "#0072B2", alpha=0.08, lw=0, zorder=0)
    ax.set_xlabel(r"Relative pressure, $p/p_0$")
    ax.set_ylabel(r"$(p/p_0)\,/\,[n(1-p/p_0)]$ " + (f"($10^{{{k}}}$ g cm$^{{-3}}$)" if k else "(g cm$^{-3}$)"))
    pmax = max(min(0.35, b["p_hi"] * 1.35) for _, b, _ in items)
    ytop = max(max((b["slope"] * b["p_hi"] * 1.08 + b["intercept"]) / sc, ymax_all / sc) for _, b, _ in items)
    ax.set_xlim(0, pmax * 1.02); ax.set_ylim(0, ytop * 1.06)
    ax.xaxis.set_major_locator(MaxNLocator(5)); ax.yaxis.set_major_locator(MaxNLocator(5))
    from matplotlib.ticker import ScalarFormatter
    fx = ScalarFormatter(useOffset=False); fx.set_scientific(False); ax.xaxis.set_major_formatter(fx)
    fy = ScalarFormatter(useOffset=False); fy.set_scientific(False); ax.yaxis.set_major_formatter(fy)
    fs = plt.rcParams["font.size"] - 1
    if box and not overlay:
        nm_ = (f"{name}\n" if name else "")
        txt = (nm_ + f"$S_\\mathrm{{BET}}$ = {bet['area']:.0f} m$^2$ g$^{{-1}}$\n$C$ = {bet['C']:.0f}   $R^2$ = {bet['R2']:.4f}\n"
               f"$p/p_0$ = {bet['p_lo']:.2g}–{bet['p_hi']:.2g} ({bet['n']} pts)")
        if bet.get("nm") is not None:
            txt += f"\n$n_\\mathrm{{m}}$ = {bet['nm']:.1f} cm$^3$ g$^{{-1}}$"
        ax.text(0.96, 0.04, txt, transform=ax.transAxes, va="bottom", ha="right", fontsize=fs, linespacing=1.25,
                multialignment="left",
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#C8D1DC", lw=0.6, alpha=0.92), zorder=5)
    elif overlay and box:
        rows = "\n".join(f"{nm}: {b['area']:.0f} m$^2$ g$^{{-1}}$ (C {b['C']:.0f})" for nm, b, _ in items)
        ax.text(0.96, 0.04, rows, transform=ax.transAxes, va="bottom", ha="right", fontsize=fs - 0.5, linespacing=1.3,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#C8D1DC", lw=0.6, alpha=0.92), zorder=5)
    return ax


def _draw_bet_overlay(ax, items, box=True):
    """Several samples: log–log axes (BET ranges and n_m differ by orders of magnitude between MOFs), fitted points and
    fit line per sample, S_BET/C table in the empty upper-left corner."""
    from matplotlib.ticker import LogLocator
    lo, hi = np.inf, 0.0
    for nm, b, col in items:
        p, y = np.asarray(b["plot_p"]), np.asarray(b["plot_y"])
        sel = (p >= b["p_lo"]) & (p <= b["p_hi"]) & (y > 0)
        ax.plot(p[sel], y[sel], "o", color=col, ms=3.2, mew=0.6, zorder=3, label=nm)
        xx = np.geomspace(b["p_lo"] / 1.25, b["p_hi"] * 1.15, 50)
        yy = b["slope"] * xx + b["intercept"]
        ok = yy > 0
        ax.plot(xx[ok], yy[ok], "-", color=col, lw=0.9, zorder=2)
        if sel.any():
            lo, hi = min(lo, y[sel].min()), max(hi, y[sel].max())
    ax.set_xscale("log"); ax.set_yscale("log")
    if np.isfinite(lo) and hi > 0:
        ax.set_ylim(lo / 2.5, hi * 10 ** (0.27 * len(items) + 0.4))
    ax.xaxis.set_major_locator(LogLocator(numticks=6))
    ax.set_xlabel(r"Relative pressure, $p/p_0$")
    ax.set_ylabel(r"$(p/p_0)\,/\,[n(1-p/p_0)]$ (g cm$^{-3}$)")
    if box:   # S_BET and C in the legend (upper-left is empty on log–log BET axes)
        h = [l for l in ax.lines if l.get_label() in [nm for nm, _, _ in items]]
        lab = [f"{nm}: {b['area']:.0f} m$^2$ g$^{{-1}}$ (C {b['C']:.0f})" for nm, b, _ in items]
        ax.legend(h, lab, loc="upper left", fontsize=plt.rcParams["font.size"] - 2, frameon=False, handletextpad=0.3,
                  labelspacing=0.2, handlelength=1.0)
    return ax


def draw_rouquerol(ax, bet, color="k", name="", overlay=None):
    """Rouquerol consistency plot n(1−p/p0) vs p/p0 (log): the chosen BET range shaded and labelled, its points filled,
    the maximum (upper limit of the valid range) marked. overlay: [(name, bet, colour)]."""
    items = overlay or [(name, bet, color)]
    for nm, b, col in items:
        p, r = np.asarray(b["rq_p"]), np.asarray(b["rq"])
        sel = (p >= b["p_lo"]) & (p <= b["p_hi"])
        ax.plot(p, r, "-", color=col, lw=0.8, zorder=2, label=nm if overlay else None)
        ax.plot(p[~sel], r[~sel], "o", mfc="white", mec=col, ms=2.6, mew=0.6, zorder=3)
        ax.plot(p[sel], r[sel], "o", color=col, ms=3.0, zorder=3)
        im = int(np.argmax(r))
        ax.plot([p[im]], [r[im]], marker="v", color="#D55E00", ms=5, zorder=4, ls="none")
        if not overlay:
            ax.axvspan(b["p_lo"], b["p_hi"], color="#0072B2", alpha=0.10, lw=0, zorder=0)
            ax.text(np.sqrt(b["p_lo"] * b["p_hi"]), 0.03, "BET range", transform=ax.get_xaxis_transform(), ha="center",
                    va="bottom", fontsize=plt.rcParams["font.size"] - 2, color="#0B4F7C")
    ax.set_xscale("log")
    ax.set_xlabel(r"Relative pressure, $p/p_0$"); ax.set_ylabel(r"$n(1-p/p_0)$ (cm$^3$ g$^{-1}$ STP)")
    ax.set_ylim(bottom=0)
    if overlay:
        ax.set_ylim(0, ax.get_ylim()[1] * (1.12 + 0.07 * math.ceil(len(items) / 2)))
        ax.legend(fontsize=plt.rcParams["font.size"] - 2, frameon=False, loc="upper left", ncol=2, columnspacing=0.8,
                  handlelength=1.2, handletextpad=0.3, labelspacing=0.2)
    return ax


def draw_tplot(ax, tp, trng, color="k", name=""):
    """t-plot (Harkins–Jura): data open, fit window shaded, fit line extrapolated to t = 0 (intercept = micropores)."""
    ax.plot(tp["t"], tp["q"], "o", ms=3, color=color, mfc="white", mew=0.7)
    tt = np.linspace(0, max(np.max(tp["t"]), trng[1]), 50)
    ax.plot(tt, tp["slope"] * tt + tp["intercept"], "-", color="#D55E00", lw=0.9)
    ax.axvspan(*trng, color="#D55E00", alpha=0.08, lw=0)
    ax.set_xlabel(r"Statistical thickness, $t$ (Å)"); ax.set_ylabel(r"Quantity adsorbed (cm$^3$ g$^{-1}$ STP)")
    ax.set_xlim(0, None); ax.set_ylim(0, None)
    ax.text(0.96, 0.05, (f"{name}\n" if name else "") + f"$V_\\mathrm{{micro}}$ = {tp['V_micro']:.3f} cm$^3$ g$^{{-1}}$\n"
            f"$S_\\mathrm{{ext}}$ = {tp['S_ext']:.0f} m$^2$ g$^{{-1}}$", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=plt.rcParams["font.size"] - 1, bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#C8D1DC", lw=0.6))
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


SCALE_UNITS = {"pm": 1e-3, "Å": 0.1, "nm": 1.0, "µm": 1e3, "mm": 1e6}
SCALEBAR_DEFAULT = {"mode": "auto", "value": None, "unit": "auto", "length_px": None, "length_frac": None,
                    "position": "lower right", "color": "white", "thickness": 0.018, "font_size": 0.0, "box": False,
                    "box_alpha": 0.65, "margin": 0.04, "show_label": True, "label": "", "bold": True, "outline": True}


def format_length(nm, unit="auto", sig=3):
    """1500 nm -> '1.5 µm'; unit forces nm/µm/mm/Å/pm. Rounded to `sig` significant figures."""
    if unit in (None, "", "auto"):
        unit = "mm" if nm >= 1e6 else ("µm" if nm >= 1000 else ("nm" if nm >= 1 else "Å"))
    v = nm / SCALE_UNITS[unit]
    if v == 0:
        return f"0 {unit}"
    v = float(f"{v:.{sig}g}")
    return f"{v:g} {unit}"


def scalebar_geometry(pixel_nm, width_px, bar=None):
    """Resolve the bar settings into (length_px, length_nm, label).
    mode 'value' : the user gives the number + unit -> length in px follows from the pixel size (stays physical)
    mode 'length': the user gives the length (px or fraction of image width) -> label = length x pixel size
    mode 'auto'  : a round number close to 20 % of the width."""
    b = {**SCALEBAR_DEFAULT, **(bar or {})}
    if not pixel_nm or pixel_nm <= 0:
        return None
    if b["mode"] == "value" and b.get("value"):
        unit = b["unit"] if b["unit"] not in ("auto", None) else "nm"
        L_nm = float(b["value"]) * SCALE_UNITS[unit]
        px = L_nm / pixel_nm
        label = format_length(L_nm, unit)
    elif b["mode"] == "length" and (b.get("length_px") or b.get("length_frac")):
        px = float(b["length_px"] or b["length_frac"] * width_px)
        L_nm = px * pixel_nm
        label = format_length(L_nm, b["unit"])
    else:
        L_nm = nice_bar(0.22 * width_px * pixel_nm)
        px = L_nm / pixel_nm
        label = format_length(L_nm, b["unit"])
    if b.get("label"):
        label = b["label"]
    return px, L_nm, label


def draw_scalebar(ax, h, w, pixel_nm, bar=None, fs=None):
    b = {**SCALEBAR_DEFAULT, **(bar or {})}
    g = scalebar_geometry(pixel_nm, w, b)
    if g is None:
        return None
    px, L_nm, label = g
    fs = b["font_size"] or fs or plt.rcParams["font.size"]
    col = b["color"]
    other = "black" if col == "white" else "white"
    th = b["thickness"] * h if b["thickness"] < 1 else b["thickness"]      # fraction of height or px
    m = b["margin"]
    pos = b["position"]
    right, bottom = "right" in pos, "lower" in pos
    x0 = (w * (1 - m) - px) if right else w * m
    # text height in data px (approx from font size and axes size)
    fig = ax.figure
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    px_per_pt = h / max(bbox.height * 72, 1)
    th_txt = fs * 1.0 * px_per_pt if b["show_label"] else 0
    gap = 0.25 * fs * px_per_pt
    if bottom:
        ybar = h * (1 - m) - th
        ytxt = ybar - gap
        va = "bottom"
    else:
        ytxt = h * m
        ybar = ytxt + th_txt + gap
        va = "top"
    artists = []
    if b["box"]:
        pad = 0.25 * fs * px_per_pt
        top = min(ybar, ytxt - th_txt) - pad if bottom else ytxt - pad
        bot = ybar + th + pad
        tw = max(px, len(label) * fs * 0.6 * px_per_pt) if b["show_label"] else px
        cx = x0 + px / 2
        artists.append(ax.add_patch(Rectangle((cx - tw / 2 - pad, top), tw + 2 * pad, bot - top, color=other,
                                              alpha=b["box_alpha"], lw=0, zorder=4)))
    artists.append(ax.add_patch(Rectangle((x0, ybar), px, th, color=col, lw=0, zorder=5)))
    if b["show_label"]:
        pe_ = [pe.withStroke(linewidth=1.6, foreground=other)] if (b["outline"] and not b["box"]) else None
        # label aligned to the bar end on the image-edge side, so it can never run past the panel edge
        ha_ = b.get("label_align") or ("right" if right else "left")
        xt = {"right": x0 + px, "left": x0, "center": x0 + px / 2}[ha_]
        t = ax.text(xt, ytxt if bottom else ytxt, label, color=col, ha=ha_, va=va if bottom else "top",
                    fontsize=fs, fontweight="bold" if b["bold"] else "normal", path_effects=pe_, zorder=6)
        t.set_gid("mfs_note")
        artists.append(t)
    return {"length_px": px, "length_nm": L_nm, "label": label}


def draw_sem(ax, img, pixel_nm=None, crop_row=None, bar_nm=None, label=None, bar_color="white", loc="lower right", fs=None,
             bar=None):
    """SEM/TEM image with a physically correct scale bar. `bar` = scale-bar settings (see SCALEBAR_DEFAULT);
    the old keywords bar_nm / bar_color / loc still work."""
    im = img[:crop_row] if crop_row else img
    ax.imshow(im, cmap="gray", interpolation="lanczos", aspect="equal")
    ax.set_xticks([]); ax.set_yticks([])
    ax.tick_params(which="both", top=False, right=False, left=False, bottom=False)
    for s in ax.spines.values():
        s.set_linewidth(0.5)
    h, w = im.shape[:2]
    fs = fs or plt.rcParams["font.size"]
    b = dict(bar or {})
    if bar_nm and not b.get("mode"):
        b.update(mode="value", value=bar_nm, unit="nm")
    b.setdefault("color", bar_color); b.setdefault("position", loc)
    info = draw_scalebar(ax, h, w, pixel_nm, b, fs) if pixel_nm else None
    ax._mfs_scalebar = info
    if label:
        t = ax.text(0.03, 0.97, label, transform=ax.transAxes, color="white", va="top", ha="left", fontsize=fs + 1, fontweight="bold",
                    path_effects=_OUTLINE)
        t.set_gid("mfs_note")
    return ax


def draw_mag_label(ax, text, position="lower left", color="white", fs=None):
    """Magnification text (e.g. '×10 000') in a corner of an SEM panel, outlined like the scale-bar label."""
    right, bottom = "right" in position, "lower" in position
    other = "black" if color == "white" else "white"
    t = ax.text(0.97 if right else 0.03, 0.04 if bottom else 0.96, text, transform=ax.transAxes, color=color,
                ha="right" if right else "left", va="bottom" if bottom else "top", fontsize=fs or plt.rcParams["font.size"],
                fontweight="bold", path_effects=[pe.withStroke(linewidth=1.6, foreground=other)], zorder=6)
    t.set_gid("mfs_note")
    return t


def panel_font_size(panel_w_in, panel_h_in, frac=0.038, min_pt=4.5, max_pt=9.0):
    """Annotation font size for one image panel: a fraction of the panel's shorter side, clipped to [min_pt, max_pt]."""
    return float(np.clip(min(panel_w_in, panel_h_in) * 72 * frac, min_pt, max_pt))


def _inside(ax, artist, pad_px, renderer):
    bb = artist.get_window_extent(renderer); ab = ax.get_window_extent(renderer)
    return bb.x0 >= ab.x0 + pad_px - 0.5 and bb.x1 <= ab.x1 - pad_px + 0.5 and bb.y0 >= ab.y0 + pad_px - 0.5 and bb.y1 <= ab.y1 - pad_px + 0.5


def fit_annotations(fig, axs, pad_pt=2.0, min_pt=4.0, max_iter=16):
    """Shrink the scale-bar / magnification / panel-label text of every image panel (uniformly) until each text bbox
    lies inside its own axes with `pad_pt` padding. Returns the final font size."""
    fig.canvas.draw(); from .render import renderer_for; r = renderer_for(fig); pad = pad_pt * fig.dpi / 72
    texts = [(ax, t) for ax in axs for t in ax.texts if t.get_gid() == "mfs_note" and t.get_visible()]
    if not texts:
        return None
    fs = max(t.get_fontsize() for _, t in texts)
    for _ in range(max_iter):
        bad = [(ax, t) for ax, t in texts if not _inside(ax, t, pad, r)]
        if not bad or fs <= min_pt:
            break
        fs = max(min_pt, fs - 0.5)
        for _, t in texts:
            t.set_fontsize(fs)
        fig.canvas.draw(); from .render import renderer_for; r = renderer_for(fig)
    return fs


def fig_sem_panel(images, ncols=2, width=DOUBLE, labels=None, gap=0.02, font_size=None, min_pt=4.5, max_pt=9.0):
    """images: list of dict(image, pixel_nm, crop_row, title). Scale-bar label, bar thickness, magnification and panel
    labels are sized relative to each panel (one uniform font for all panels: the smallest panel decides, or
    `font_size` to override) and kept inside the panel bounds with padding."""
    n = len(images); nrows = math.ceil(n / ncols)
    ar = np.mean([(im["image"][:im.get("crop_row") or None].shape[0] / im["image"].shape[1]) for im in images])
    style.apply_style()
    fig, axs = plt.subplots(nrows, ncols, figsize=(width, width / ncols * ar * nrows * 1.0 + 0.1))
    axs = np.atleast_1d(axs).ravel()
    pw = width / ncols
    fs = float(font_size) if font_size else panel_font_size(pw, pw * ar, min_pt=min_pt, max_pt=max_pt)
    th_pt = float(np.clip(pw * ar * 72 * 0.018, 1.0, 3.5))               # bar thickness in points
    for k, ax in enumerate(axs):
        if k >= n:
            ax.axis("off"); continue
        d = images[k]
        lab = (labels[k] if labels else f"({chr(97+k)})")
        im_ = d["image"][:d.get("crop_row") or None]
        b = dict(d.get("bar") or {})
        b.setdefault("font_size", 0.0)
        if not b.get("font_size"):
            b["font_size"] = fs
        if b.get("thickness", SCALEBAR_DEFAULT["thickness"]) < 1:           # relative thickness -> points of this panel
            b["thickness"] = th_pt / 72 / (pw * ar) * im_.shape[0]
        draw_sem(ax, d["image"], d.get("pixel_nm"), d.get("crop_row"), d.get("bar_nm"), lab, bar=b, fs=fs)
        for t in ax.texts:
            if t.get_text() == lab:
                t.set_fontsize(fs)
        if d.get("mag_text"):
            draw_mag_label(ax, d["mag_text"], d.get("mag_pos") or "lower left", (d.get("bar") or {}).get("color", "white"), fs=fs)
        if d.get("title"):
            t = ax.text(0.97, 0.97, d["title"], transform=ax.transAxes, color="white", ha="right", va="top",
                        fontsize=fs, path_effects=_OUTLINE)
            t.set_gid("mfs_note")
    fig.subplots_adjust(wspace=gap, hspace=gap, left=0, right=1, top=1, bottom=0)
    fig._mfs_panel_fs = fit_annotations(fig, [a for a in axs[:n]], min_pt=min(4.0, fs))
    return fig


def draw_edx_map(ax, img, cmap_color="#D55E00", pixel_nm=None, label=None, element=None, bar=None):
    from matplotlib.colors import LinearSegmentedColormap
    cm = LinearSegmentedColormap.from_list("el", ["#000000", cmap_color])
    ax.imshow(img, cmap=cm, interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    ax.tick_params(which="both", top=False, right=False, left=False, bottom=False)
    if element:
        ax.text(0.97, 0.97, element, transform=ax.transAxes, color="white", ha="right", va="top", fontweight="bold")
    if pixel_nm:
        draw_scalebar(ax, img.shape[0], img.shape[1], pixel_nm, bar)
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
