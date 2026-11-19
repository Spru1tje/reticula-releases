"""v2.18.3: renders an info document (see helpdocs.py) as a PDF in the Reticula colours (dark blue page, side ribbon, light formula cards).
A document is {"title": str, "subtitle": str, "sections": [(kind, text), ...]} with kind
  h = heading, p = paragraph, b = bullet, n = numbered step (numbering restarts at every heading), f = formula (matplotlib mathtext), k = key point (highlighted box)."""
from __future__ import annotations

import io
import textwrap

import matplotlib
from matplotlib.figure import Figure
from matplotlib.patches import FancyBboxPatch, Rectangle

BG, RIB, TEAL, INK, HEAD, ACC, MUTED = "#0F2A52", "#17386B", "#0F766E", "#E2E8F0", "#FFFFFF", "#5EEAD4", "#94A3B8"
TAGLINE = "Reticula  ·  From raw data to figures, the easy way."
W_CHARS = 82


def render_pdf(doc: dict, png_dpi: int = 0):
    """One document -> PDF bytes (pages are added as needed). With png_dpi > 0 a list of PNG bytes (one per page) is returned instead (for checks).
    DejaVu Sans is used (it has the sub/superscript and Greek characters of the texts; the figure style's Arial has not)."""
    with matplotlib.rc_context({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "mathtext.fontset": "dejavusans"}):
        return _render(doc, png_dpi)


def _render(doc: dict, png_dpi: int = 0):
    from matplotlib.backends.backend_pdf import PdfPages
    pages: list[Figure] = []
    st = {"fig": None, "y": 0.0}

    def new_page(first=False):
        fig = Figure(figsize=(8.27, 11.69))
        fig.add_artist(Rectangle((0, 0), 1, 1, transform=fig.transFigure, color=BG, zorder=0))
        fig.add_artist(Rectangle((0, 0), 0.045, 1, transform=fig.transFigure, color=RIB, zorder=1))
        fig.add_artist(Rectangle((0.045, 0), 0.006, 1, transform=fig.transFigure, color=TEAL, zorder=1))
        fig.text(0.075, 0.025, TAGLINE, color=MUTED, fontsize=8)
        fig.text(0.93, 0.025, str(len(pages) + 1), color=MUTED, fontsize=8, ha="right")
        y = 0.955
        if first:
            for line in textwrap.wrap(doc["title"], 48):
                fig.text(0.075, y, line, color=HEAD, fontsize=17, fontweight="bold", va="top")
                y -= 0.032
            if doc.get("subtitle"):
                fig.text(0.075, y, doc["subtitle"], color=ACC, fontsize=10, va="top")
                y -= 0.026
            fig.add_artist(Rectangle((0.075, y), 0.2, 0.004, transform=fig.transFigure, color=ACC))
            y -= 0.028
        pages.append(fig)
        st["fig"], st["y"] = fig, y

    new_page(True)
    step = 0
    for kind, text in doc["sections"]:
        fig = st["fig"]
        need = {"h": 0.07, "f": 0.08, "k": 0.08}.get(kind, 0.045)
        if st["y"] - need < 0.06:
            new_page()
            fig = st["fig"]
        y = st["y"]
        if kind == "h":
            step = 0
            y -= 0.008
            fig.text(0.075, y, text, color=ACC, fontsize=12.5, fontweight="bold", va="top")
            y -= 0.03
        elif kind == "f":
            hh = 0.078 if (r"\sum" in text or r"\int" in text or r"\sqrt{\dfrac" in text) else 0.058
            fig.add_artist(FancyBboxPatch((0.12, y - hh - 0.002), 0.76, hh, boxstyle="round,pad=0.002,rounding_size=0.01",
                                          transform=fig.transFigure, fc="#F1F5F9", ec=ACC, lw=1))
            fig.text(0.5, y - hh / 2 - 0.002, text, fontsize=13, ha="center", va="center", color=BG)
            y -= hh + 0.018
        elif kind == "k":
            lines = textwrap.wrap(text, W_CHARS - 8)
            h = 0.018 * len(lines) + 0.016
            if y - h < 0.06:
                new_page(); fig = st["fig"]; y = st["y"]
            fig.add_artist(FancyBboxPatch((0.075, y - h), 0.82, h, boxstyle="round,pad=0.002,rounding_size=0.008",
                                          transform=fig.transFigure, fc=RIB, ec=ACC, lw=1))
            yy = y - 0.01
            for l in lines:
                fig.text(0.092, yy, l, fontsize=9.4, va="top", color=HEAD); yy -= 0.018
            y -= h + 0.012
        else:
            width = {"p": W_CHARS, "b": W_CHARS - 4, "n": W_CHARS - 6}[kind]
            if kind == "n":
                step += 1
            lines = textwrap.wrap(text, width)
            for k, l in enumerate(lines):
                if y < 0.06:
                    st["y"] = y; new_page(); fig = st["fig"]; y = st["y"]
                if k == 0 and kind == "b":
                    fig.text(0.08, y, "•", color=ACC, fontsize=10, va="top")
                if k == 0 and kind == "n":
                    fig.text(0.078, y, f"{step}", color=ACC, fontsize=9.6, va="top", fontweight="bold")
                fig.text(0.075 if kind == "p" else (0.095 if kind == "b" else 0.1), y, l, fontsize=9.6, va="top", color=INK)
                y -= 0.0185
            y -= 0.008
        st["y"] = y
    if png_dpi:
        outs = []
        for f in pages:
            b = io.BytesIO(); f.savefig(b, format="png", dpi=png_dpi); outs.append(b.getvalue())
        return outs
    buf = io.BytesIO()
    with PdfPages(buf) as pdf:
        for f in pages:
            pdf.savefig(f)
    return buf.getvalue()
