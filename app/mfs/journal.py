"""v2.15: journal check of a finished figure + one-click fix, and the one-Excel export of all figure data.

check(fig, width_in, ...) -> list of issues (each: key, level 'error'|'warn'|'ok', text) for the figure AS PRINTED at `width_in`
(the figure is scaled to that width, so a 7 in wide figure printed in a 3.3 in column has its fonts shrunk by 3.3/7).
fix(fig, width_in, ...) -> list of what was changed: sets the figure to the journal width (fonts and lines then print 1:1), raises
every font / line below the minimum, and swaps colours that colour-blind readers cannot tell apart for the Okabe-Ito palette.
"""
from __future__ import annotations

import io
import re

import numpy as np
from matplotlib import colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.text import Text

OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]
# journal presets: final width (in) of one column / two columns, minimum font (pt), minimum line (pt), minimum dpi
JOURNALS = {
    "ACS (single column 3.25 in)": (3.25, 5.0, 0.5, 600),
    "ACS (double column 7 in)": (7.0, 5.0, 0.5, 600),
    "RSC (single column 3.27 in)": (3.27, 6.0, 0.5, 600),
    "Elsevier (single column 3.54 in)": (3.54, 6.0, 0.5, 600),
    "Nature (single column 3.5 in)": (3.5, 5.0, 0.25, 600),
    "Wiley (single column 3.15 in)": (3.15, 6.0, 0.5, 600),
    "Custom": (3.5, 6.0, 0.5, 600),
}

# colour-vision deficiency simulation (Machado et al. 2009, severity 1.0), applied to linear RGB
_CVD = {
    "deuteranopia": np.array([[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]]),
    "protanopia": np.array([[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]]),
    "tritanopia": np.array([[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]]),
}


def _lin(c):
    c = np.asarray(c, float)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _lab(rgb_lin):
    M = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]])
    xyz = M @ np.clip(rgb_lin, 0, 1) / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > 216 / 24389, np.cbrt(xyz), (24389 / 27 * xyz + 16) / 116)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def delta_e(c1, c2, cvd: str | None = None) -> float:
    """CIE76 colour difference of two colours as seen by a normal or colour-blind (`cvd`) reader."""
    a, b = _lin(mcolors.to_rgb(c1)), _lin(mcolors.to_rgb(c2))
    if cvd:
        a, b = _CVD[cvd] @ a, _CVD[cvd] @ b
    return float(np.linalg.norm(_lab(a) - _lab(b)))


def _is_grey(c, tol=0.08):
    r, g, b = mcolors.to_rgb(c)
    return max(r, g, b) - min(r, g, b) < tol


def _texts(fig):
    for t in fig.findobj(Text):
        try:
            if t.get_visible() and t.get_text().strip():
                yield t
        except Exception:
            continue


def _lines(fig):
    for ln in fig.findobj(Line2D):
        try:
            if not ln.get_visible() or ln.get_linestyle() in ("None", "", " "):
                continue
            if len(ln.get_xdata()) < 2:
                continue
            yield ln
        except Exception:
            continue


def _series_colours(fig) -> list[str]:
    """The distinct non-grey colours of the plotted data (lines, scatter points, bars) in order of first appearance."""
    out: list[str] = []

    def add(c):
        try:
            h = mcolors.to_hex(c)
        except Exception:
            return
        if not _is_grey(h) and h not in out:
            out.append(h)
    for ln in _lines(fig):
        add(ln.get_color())
    for ax in fig.axes:
        for coll in ax.collections:
            try:
                fc = coll.get_facecolor()
                if len(fc):
                    add(fc[0])
            except Exception:
                pass
        for p in ax.patches:
            try:
                add(p.get_facecolor())
            except Exception:
                pass
    return out


def confusable_pairs(colours, min_de: float = 12.0):
    """[(colour a, colour b, kind, deltaE)] for every pair that is hard to tell apart for a colour-blind reader (kind = deuteranopia …)."""
    out = []
    for i in range(len(colours)):
        for j in range(i + 1, len(colours)):
            for k in _CVD:
                d = delta_e(colours[i], colours[j], k)
                if d < min_de and delta_e(colours[i], colours[j]) >= min_de:
                    out.append((colours[i], colours[j], k, round(d, 1)))
                    break
    return out


def check(fig, width_in: float, min_font: float = 6.0, min_line: float = 0.5, png_dpi: int | None = None, min_dpi: int = 600) -> list[dict]:
    res: list[dict] = []
    w0 = float(fig.get_figwidth())
    s = width_in / w0 if w0 else 1.0
    res.append({"key": "width", "level": "ok" if abs(s - 1) < 0.02 else "warn",
                "text": f"Figure is {w0:.2f} in wide; final width {width_in:.2f} in" + ("" if abs(s - 1) < 0.02 else f" → everything prints at {s * 100:.0f} % of the size on screen.")})
    sizes = [(t.get_fontsize() * s, t) for t in _texts(fig)]
    small = [(sz, t) for sz, t in sizes if sz < min_font - 1e-6]
    if small:
        lo = min(sz for sz, _ in small)
        res.append({"key": "font", "level": "error", "text": f"{len(small)} text(s) print smaller than {min_font:g} pt (smallest {lo:.1f} pt)."})
    else:
        res.append({"key": "font", "level": "ok", "text": f"Smallest text prints at {min((sz for sz, _ in sizes), default=0):.1f} pt (≥ {min_font:g} pt)."})
    lws = [ln.get_linewidth() * s for ln in _lines(fig)]
    for ax in fig.axes:
        lws += [sp.get_linewidth() * s for sp in ax.spines.values() if sp.get_visible()]
    thin = [x for x in lws if x < min_line - 1e-6]
    if thin:
        res.append({"key": "line", "level": "error", "text": f"{len(thin)} line(s) print thinner than {min_line:g} pt (thinnest {min(thin):.2f} pt)."})
    else:
        res.append({"key": "line", "level": "ok", "text": f"Thinnest line prints at {min(lws, default=0):.2f} pt (≥ {min_line:g} pt)."})
    cols = _series_colours(fig)
    pairs = confusable_pairs(cols)
    if pairs:
        res.append({"key": "colour", "level": "warn", "text": f"{len(pairs)} colour pair(s) are hard to tell apart with colour-blindness (e.g. {pairs[0][0]} / {pairs[0][1]}, {pairs[0][2]})."})
    else:
        res.append({"key": "colour", "level": "ok", "text": f"{len(cols)} data colour(s): distinguishable with the common types of colour-blindness."})
    if png_dpi is not None:
        res.append({"key": "dpi", "level": "ok" if png_dpi >= min_dpi else "warn",
                    "text": f"PNG/TIFF resolution {png_dpi} dpi" + ("" if png_dpi >= min_dpi else f" – journals ask ≥ {min_dpi} dpi for line art (sidebar → Export dpi).")})
    return res


def fix(fig, width_in: float, min_font: float = 6.0, min_line: float = 0.5, recolour: bool = True) -> list[str]:
    """Makes the figure fit the journal: width = final width (then pt = pt), fonts and lines raised to the minimum, colour-blind safe colours.
    Returns a list of what was changed."""
    done: list[str] = []
    w0, h0 = fig.get_size_inches()
    if abs(w0 - width_in) > 0.02 * width_in:
        fig.set_size_inches(width_in, h0 * width_in / w0)
        done.append(f"width {w0:.2f} → {width_in:.2f} in (height keeps the aspect ratio)")
    n = 0
    for t in _texts(fig):
        if t.get_fontsize() < min_font - 1e-6:
            t.set_fontsize(min_font)
            n += 1
    if n:
        done.append(f"{n} text(s) raised to {min_font:g} pt")
    n = 0
    for ln in _lines(fig):
        if ln.get_linewidth() < min_line - 1e-6:
            ln.set_linewidth(min_line)
            n += 1
    for ax in fig.axes:
        for sp in ax.spines.values():
            if sp.get_visible() and sp.get_linewidth() < min_line - 1e-6:
                sp.set_linewidth(min_line)
                n += 1
    if n:
        done.append(f"{n} line(s) raised to {min_line:g} pt")
    if recolour:
        cols = _series_colours(fig)
        if confusable_pairs(cols):
            mp = {c: OKABE_ITO[i % len(OKABE_ITO)] for i, c in enumerate(cols)}
            for ln in _lines(fig):
                h = mcolors.to_hex(ln.get_color())
                if h in mp:
                    ln.set_color(mp[h])
            for ax in fig.axes:
                for coll in ax.collections:
                    try:
                        fc = coll.get_facecolor()
                        if len(fc) and mcolors.to_hex(fc[0]) in mp:
                            coll.set_facecolor(mp[mcolors.to_hex(fc[0])])
                    except Exception:
                        pass
                for p in ax.patches:
                    h = mcolors.to_hex(p.get_facecolor())
                    if h in mp:
                        p.set_facecolor(mp[h])
            done.append("data colours replaced by the colour-blind safe Okabe-Ito palette")
    return done


def sheet_name(name: str, used: set) -> str:
    base = re.sub(r"[\[\]\*\?:/\\]", "_", str(name))[:28] or "figure"
    nm, k = base, 2
    while nm.lower() in used:
        nm = f"{base[:25]}_{k}"
        k += 1
    used.add(nm.lower())
    return nm


def figure_data_xlsx(figs: dict, version: str = "") -> bytes:
    """One Excel file with one sheet per figure (the plotted data) and an 'Index' sheet (figure, page, rows, columns, sheet).
    `figs` = {name: {"data": DataFrame | None, "page": str, "title": str}}."""
    import pandas as pd
    buf = io.BytesIO()
    used: set = {"index"}
    index = []
    with pd.ExcelWriter(buf, engine="openpyxl") as xw:
        for k, f in figs.items():
            d = f.get("data")
            if d is None or not hasattr(d, "to_excel") or not len(d):
                index.append({"figure": k, "page": f.get("page", ""), "title": f.get("title") or "", "sheet": "(no data table)", "rows": 0, "columns": ""})
                continue
            sn = sheet_name(k, used)
            d.to_excel(xw, sheet_name=sn, index=False)
            ws = xw.sheets[sn]
            for i, c in enumerate(d.columns, 1):
                ws.column_dimensions[ws.cell(1, i).column_letter].width = min(40, max(10, len(str(c)) + 2))
            ws.freeze_panes = "A2"
            index.append({"figure": k, "page": f.get("page", ""), "title": f.get("title") or "", "sheet": sn, "rows": int(len(d)),
                          "columns": ", ".join(map(str, d.columns))[:250]})
        pd.DataFrame(index).to_excel(xw, sheet_name="Index", index=False)
        ws = xw.sheets["Index"]
        for col, w in zip("ABCDEF", (28, 18, 34, 28, 8, 60)):
            ws.column_dimensions[col].width = w
        xw.book.move_sheet("Index", offset=-len(xw.book.sheetnames) + 1)
        xw.book.properties.creator = f"Reticula {version}".strip()
    return buf.getvalue()
