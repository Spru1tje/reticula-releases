"""Multi-panel figure composer: arrange any draw_* functions on a grid with (a), (b), (c) labels.

v2.14: the panels of DIFFERENT techniques (UV/Vis, FTIR, PXRD ...) go into one figure like Fig. 2 of a paper: panel letters
inside or outside the axes, an optional caption under the figure, an automatic layout, and the settings made under each
single figure (axis titles, ranges, units, log, increments) travel with the panel."""
from __future__ import annotations
import math
import textwrap
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from . import style

LABEL_STYLES = ["(a)", "a", "A", "a)", "(A)"]


def auto_columns(n):
    """Sensible number of columns: 1, 2 and 3 panels in one row (like A | B | C), 4 as 2x2, then 3 per row."""
    return 1 if n <= 1 else (n if n <= 3 else (2 if n == 4 else 3))


def panel_label(k, label_style="(a)"):
    letter = chr(97 + k % 26)
    return {"(a)": f"({letter})", "a": letter, "A": letter.upper(), "a)": f"{letter})", "(A)": f"({letter.upper()})"}.get(label_style, f"({letter})")


def caption_lines(text, width_in, font_size):
    """Caption wrapped to the figure width (about 0.52 em per character)."""
    n = max(20, int(width_in * 72 / (0.52 * font_size)))
    out = []
    for para in str(text).split("\n"):
        out += textwrap.wrap(para, n) or [""]
    return out


def _group_boxes(fig):
    """Tight box (inches) of every panel INCLUDING its colour bar / inset axes / panel letter, keyed by its cell range in the grid."""
    rend = fig.canvas.get_renderer()
    inv = fig.dpi_scale_trans.inverted()
    out = {}
    for ax in fig.axes:
        sp = ax.get_subplotspec()
        if sp is None:
            continue
        top = sp.get_topmost_subplotspec()
        key = (top.rowspan.start, top.rowspan.stop, top.colspan.start, top.colspan.stop)
        try:
            bb = ax.get_tightbbox(rend).transformed(inv)
        except Exception:
            continue
        if key in out:
            o = out[key]
            out[key] = (min(o[0], bb.x0), min(o[1], bb.y0), max(o[2], bb.x1), max(o[3], bb.y1))
        else:
            out[key] = (bb.x0, bb.y0, bb.x1, bb.y1)
    return out


def resolve_overlaps(fig, gs, pad_in=0.06, max_iter=8):
    """Widens the gaps of the grid until neighbouring panels (with their colour bars, axis titles and letters) no longer overlap.
    Returns the final (wspace, hspace)."""
    fw, fh = fig.get_size_inches()
    ws, hs = gs.get_subplot_params(fig).wspace, gs.get_subplot_params(fig).hspace
    for _ in range(max_iter):
        fig.canvas.draw()
        g = _group_boxes(fig)
        need_w = need_h = 0.0
        keys = list(g)
        for a in keys:
            for b in keys:
                if a == b:
                    continue
                ra, rb = range(a[0], a[1]), range(b[0], b[1])
                ca, cb = range(a[2], a[3]), range(b[2], b[3])
                if b[2] == a[3] and set(ra) & set(rb):                    # b is right of a
                    need_w = max(need_w, g[a][2] - g[b][0] + pad_in)
                if b[0] == a[1] and set(ca) & set(cb):                    # b is below a
                    need_h = max(need_h, g[b][3] - g[a][1] + pad_in)
        if need_w <= 0 and need_h <= 0:
            break
        sp = gs.get_subplot_params(fig)
        ncols, nrows = gs.ncols, gs.nrows
        cell_w = fw * (sp.right - sp.left) / (ncols + sp.wspace * (ncols - 1)) if ncols > 1 else fw
        cell_h = fh * (sp.top - sp.bottom) / (nrows + sp.hspace * (nrows - 1)) if nrows > 1 else fh
        if need_w > 0 and ncols > 1:
            ws = min(1.6, ws + need_w / max(cell_w, 0.2))
        if need_h > 0 and nrows > 1:
            hs = min(1.6, hs + need_h / max(cell_h, 0.2))
        gs.update(wspace=ws, hspace=hs)
    return ws, hs


def compose(panels, ncols=2, width=style.DOUBLE, panel_aspect=0.78, label_style="(a)", wspace=0.42, hspace=0.42,
            label_x=-0.2, label_y=1.03, spans=None, font_size=None, label_inside=False, caption=None,
            auto_space=True, width_ratios=None, height_ratios=None, label_size=None, pad_in=0.06):
    """panels: list of callables f(ax) (or (callable, kwargs)). spans: optional list of (rowspan, colspan).
    label_inside: panel letters in the top-left corner inside the axes (bold, as in many papers) instead of above-left outside.
    caption: text printed under the figure (wrapped to the width).
    v2.15: auto_space = the gaps between the panels grow until nothing overlaps any more (colour bars, axis titles, panel letters of the
    neighbour); width_ratios / height_ratios = relative widths of the columns / heights of the rows; label_size = size of the panel letters (pt)."""
    style.apply_style(font_size) if font_size else style.apply_style()
    n = len(panels)
    nrows = math.ceil(n / ncols) if not spans else None
    if spans:
        # simple flow layout with spans
        nrows = 0; r = c = 0; slots = []
        for (rs, cs) in spans:
            if c + cs > ncols:
                r += 1; c = 0
            slots.append((r, c, rs, cs)); c += cs
            nrows = max(nrows, r + rs)
    else:
        slots = [(i // ncols, i % ncols, 1, 1) for i in range(n)]
    pw = width / ncols
    body_h = pw * panel_aspect * nrows + 0.15 * nrows
    cap_lines = caption_lines(caption, width, plt.rcParams["font.size"]) if caption else []
    cap_h = (len(cap_lines) * plt.rcParams["font.size"] * 1.35 / 72 + 0.12) if cap_lines else 0.0
    fig = plt.figure(figsize=(width, body_h + cap_h))
    wr = list(width_ratios)[:ncols] if width_ratios else None
    hr = list(height_ratios)[:nrows] if height_ratios else None
    gs = GridSpec(nrows, ncols, figure=fig, wspace=wspace, hspace=hspace, bottom=(cap_h / (body_h + cap_h) + 0.0), top=1.0,
                  width_ratios=wr if wr and len(wr) == ncols else None, height_ratios=hr if hr and len(hr) == nrows else None)
    axes = []
    for k, (p, (r, c, rs, cs)) in enumerate(zip(panels, slots)):
        ax = fig.add_subplot(gs[r:r + rs, c:c + cs])
        fn, kw = (p if isinstance(p, tuple) else (p, {}))
        fn(ax, **kw)
        lab = panel_label(k, label_style)
        if label_inside:
            ax.text(0.03, 0.97, lab, transform=ax.transAxes, fontsize=label_size or plt.rcParams["font.size"] + 2, fontweight="bold", va="top", ha="left",
                    bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.75), zorder=10)
        else:
            ax.text(label_x / cs, label_y, lab, transform=ax.transAxes, fontsize=label_size or plt.rcParams["font.size"] + 1,
                    fontweight="bold", va="bottom", ha="left")
        axes.append(ax)
    if auto_space and n > 1:
        resolve_overlaps(fig, gs, pad_in)
    if cap_lines:
        fig.text(0.01, (cap_h - 0.06) / (body_h + cap_h), "\n".join(cap_lines), ha="left", va="top", fontsize=plt.rcParams["font.size"],
                 linespacing=1.25)
    return fig, axes
