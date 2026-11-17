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


def compose(panels, ncols=2, width=style.DOUBLE, panel_aspect=0.78, label_style="(a)", wspace=0.42, hspace=0.42,
            label_x=-0.2, label_y=1.03, spans=None, font_size=None, label_inside=False, caption=None):
    """panels: list of callables f(ax) (or (callable, kwargs)). spans: optional list of (rowspan, colspan).
    label_inside: panel letters in the top-left corner inside the axes (bold, as in many papers) instead of above-left outside.
    caption: text printed under the figure (wrapped to the width)."""
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
    gs = GridSpec(nrows, ncols, figure=fig, wspace=wspace, hspace=hspace, bottom=(cap_h / (body_h + cap_h) + 0.0), top=1.0)
    axes = []
    for k, (p, (r, c, rs, cs)) in enumerate(zip(panels, slots)):
        ax = fig.add_subplot(gs[r:r + rs, c:c + cs])
        fn, kw = (p if isinstance(p, tuple) else (p, {}))
        fn(ax, **kw)
        lab = panel_label(k, label_style)
        if label_inside:
            ax.text(0.03, 0.97, lab, transform=ax.transAxes, fontsize=plt.rcParams["font.size"] + 2, fontweight="bold", va="top", ha="left",
                    bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.75), zorder=10)
        else:
            ax.text(label_x / cs, label_y, lab, transform=ax.transAxes, fontsize=plt.rcParams["font.size"] + 1,
                    fontweight="bold", va="bottom", ha="left")
        axes.append(ax)
    if cap_lines:
        fig.text(0.01, (cap_h - 0.06) / (body_h + cap_h), "\n".join(cap_lines), ha="left", va="top", fontsize=plt.rcParams["font.size"],
                 linespacing=1.25)
    return fig, axes
