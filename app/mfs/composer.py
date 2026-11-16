"""Multi-panel figure composer: arrange any draw_* functions on a grid with (a), (b), (c) labels."""
from __future__ import annotations
import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from . import style


def compose(panels, ncols=2, width=style.DOUBLE, panel_aspect=0.78, label_style="(a)", wspace=0.42, hspace=0.42,
            label_x=-0.2, label_y=1.03, spans=None, font_size=None):
    """panels: list of callables f(ax) (or (callable, kwargs)). spans: optional list of (rowspan, colspan)."""
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
    fig = plt.figure(figsize=(width, pw * panel_aspect * nrows + 0.15 * nrows))
    gs = GridSpec(nrows, ncols, figure=fig, wspace=wspace, hspace=hspace)
    axes = []
    for k, (p, (r, c, rs, cs)) in enumerate(zip(panels, slots)):
        ax = fig.add_subplot(gs[r:r + rs, c:c + cs])
        fn, kw = (p if isinstance(p, tuple) else (p, {}))
        fn(ax, **kw)
        letter = chr(97 + k)
        lab = {"(a)": f"({letter})", "a": letter, "A": letter.upper(), "a)": f"{letter})"}.get(label_style, f"({letter})")
        ax.text(label_x / cs, label_y, lab, transform=ax.transAxes, fontsize=plt.rcParams["font.size"] + 1,
                fontweight="bold", va="bottom", ha="left")
        axes.append(ax)
    return fig, axes
