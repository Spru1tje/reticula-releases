"""Journal (ACS / Wiley-VCH) figure style, fixed colour-blind-safe palette and export helpers."""
from __future__ import annotations
import io, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# ---- widths (inches) -------------------------------------------------------
WIDTHS = {
    "ACS single column (3.25 in)": 3.25,
    "ACS double column (7.0 in)": 7.0,
    "Wiley single column (3.35 in / 85 mm)": 3.35,
    "Wiley double column (6.89 in / 175 mm)": 6.89,
    "1.5 column (5.0 in)": 5.0,
    # v2.11 journal presets
    "RSC single column (3.27 in / 8.3 cm)": 3.27,
    "RSC double column (6.73 in / 17.1 cm)": 6.73,
    "Elsevier single column (3.54 in / 90 mm)": 3.54,
    "Elsevier double column (7.48 in / 190 mm)": 7.48,
    "Nature single column (3.50 in / 89 mm)": 3.50,
    "Nature double column (7.20 in / 183 mm)": 7.20,
    "Science single column (2.24 in / 5.7 cm)": 2.24,
    "Science double column (4.76 in / 12.1 cm)": 4.76,
    "Slide (10 in, 16:9)": 10.0,
}
SINGLE, DOUBLE = 3.25, 7.0

# ---- Okabe-Ito palette (colour-blind safe) ---------------------------------
OKABE_ITO = ["#000000", "#0072B2", "#E69F00", "#009E73", "#D55E00",
             "#CC79A7", "#56B4E9", "#F0E442", "#7F7F7F"]
# Name -> colour rules (regex on the normalised sample name, most specific first). v2.13: there are NO rules by
# default - samples simply take the palette colours in list order, and users add their own rules in the sidebar.
# EXAMPLE_SAMPLE_COLOURS is the rule set for the bundled MOF-801 example data (one click in the sidebar loads it).
EXAMPLE_SAMPLE_COLOURS = [
    (r"zro2\s*@|@.*zro2|composite|zro2.?sulfo", "#D55E00"),   # composite  - vermillion
    (r"\bpm\b|physical\s*mix|phys\.?\s*mix|mixture", "#CC79A7"),  # physical mixture - purple
    (r"sulfo|so3h|sulfonat", "#0072B2"),                      # Sulfo-MOF-801 - blue
    (r"^zro2|zro2\s*(np|nano|colloid|only)?$|zirconia", "#E69F00"),  # ZrO2 - orange
    (r"mof.?801|parent|pristine|hf.?mof|^uio|^nu", "#000000"),  # parent - black
    (r"sim|calc|reference|ref\b|ccdc|cod", "#7F7F7F"),        # simulated - grey
]
_extra = ["#009E73", "#56B4E9", "#F0E442", "#882255", "#44AA99", "#999933"]
DEFAULT_SAMPLE_COLOURS: list = []
SAMPLE_COLOURS: list = list(DEFAULT_SAMPLE_COLOURS)
_POOL = [c for c in OKABE_ITO + _extra if c != "#7F7F7F"]      # palette order used when no rule matches (grey = references)


# Per-series colour overrides set by the user (sample name -> hex). Checked before the fixed mapping, so a
# colour chosen once is used for that sample in every figure.
SERIES_OVERRIDES: dict[str, str] = {}


def set_series_colours(d: dict):
    SERIES_OVERRIDES.clear(); SERIES_OVERRIDES.update({k: v for k, v in (d or {}).items() if v})


def _rule_match(pat, text) -> bool:
    """A user-typed rule may be an invalid regex while it is being edited - then it simply does not match."""
    try:
        return bool(pat) and re.search(pat, text) is not None
    except re.error:
        return False


def colour_for(name: str, index: int = 0) -> str:
    if name in SERIES_OVERRIDES:
        return SERIES_OVERRIDES[name]
    n = re.sub(r"[$_{}\\]|mathrm", "", (name or "")).strip().lower()
    if PALETTES.get(GLOBAL.get("palette")):                 # v2.11: a chosen palette replaces the sample-family colours
        p = palette(); return p[index % len(p)]
    for pat, col in SAMPLE_COLOURS:
        if _rule_match(pat, n):
            return col
    if not SAMPLE_COLOURS:                                  # no rules: plain palette order
        return _POOL[index % len(_POOL)]
    return _extra[index % len(_extra)]


def colours_for(names) -> list[str]:
    """User overrides first, then the name rules (if any), then palette order. A second member of the same family
    (or an unknown name) gets the next unused colour-blind-safe colour - without stealing the colour of a later
    sample that owns that family colour."""
    names = list(names)
    out = [None] * len(names)
    used = set()
    if PALETTES.get(GLOBAL.get("palette")):                 # v2.11 global palette: overrides first, then palette order
        p = palette(); free = [c for c in p]
        for i, nm in enumerate(names):
            if nm in SERIES_OVERRIDES:
                out[i] = SERIES_OVERRIDES[nm]
        j = 0
        for i in range(len(names)):
            if out[i] is None:
                out[i] = p[j % len(p)]; j += 1
        return out
    for i, nm in enumerate(names):
        if nm in SERIES_OVERRIDES:
            out[i] = SERIES_OVERRIDES[nm]; used.add(out[i])
    for i, nm in enumerate(names):
        if out[i] is not None:
            continue
        clean = re.sub(r"[$_{}\\]|mathrm", "", (nm or "")).lower()
        if any(_rule_match(p, clean) for p, _ in SAMPLE_COLOURS):
            c = colour_for(nm, 0)
            if c not in used:
                out[i] = c; used.add(c)
    pool = list(_POOL)
    # colours owned by a rule are kept free for the sample that matches it (nothing is reserved without rules)
    reserved = {colour_for(n, 0) for n in names if n not in SERIES_OVERRIDES} if SAMPLE_COLOURS else set()
    for i, nm in enumerate(names):
        if out[i] is None:
            c = next((p for p in pool if p not in used and p not in reserved), None) or next((p for p in pool if p not in used), "#444444")
            out[i] = c; used.add(c)
    return out


SERIF_FAMILIES = ("Times New Roman", "Liberation Serif", "Nimbus Roman", "DejaVu Serif", "Cambria", "Caladea", "Georgia", "Palatino Linotype",
                  "Book Antiqua", "Garamond", "Constantia")


def _font_family(pref=None):
    avail = {f.name for f in font_manager.fontManager.ttflist}
    if pref and pref != "auto":
        if pref in avail:
            return pref
        alias = {"Arial": ["Liberation Sans", "Nimbus Sans"], "Helvetica": ["Nimbus Sans", "Liberation Sans"],
                 "Times New Roman": ["Liberation Serif", "Nimbus Roman", "DejaVu Serif"], "Calibri": ["Carlito"],
                 "Cambria": ["Caladea"]}.get(pref, [])
        for f in alias:
            if f in avail:
                return f
    for f in ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"]:
        if f in avail:
            return f
    return "sans-serif"


# ---- v2.11 global figure style (every page). Defaults = the v2.10 look, so nothing changes unless the user changes it.
GLOBAL_DEFAULTS = {"font_size": 8.0, "family": "auto", "line_width": 1.0, "axes_width": 0.75, "marker_size": 3.5,
                   "tick_dir": "in", "minor": True, "palette": "Okabe–Ito (sample families)", "legend_frame": False}
GLOBAL: dict = dict(GLOBAL_DEFAULTS)
_FONT_BASE = ["auto", "Arial", "Helvetica", "Liberation Sans", "DejaVu Sans", "Times New Roman", "Calibri", "Cambria"]
_FONT_EXTRA = ["Arial Narrow", "Verdana", "Tahoma", "Segoe UI", "Trebuchet MS", "Century Gothic", "Candara", "Corbel", "Georgia", "Palatino Linotype",
               "Book Antiqua", "Garamond", "Constantia", "Liberation Serif", "DejaVu Serif", "Nimbus Sans", "Nimbus Roman", "Carlito", "Caladea",
               "Courier New", "Consolas", "DejaVu Sans Mono", "Liberation Mono", "Comic Sans MS"]


def available_fonts(extra=True):
    """Classic journal fonts always, plus any of the extra families that are installed on this computer."""
    avail = {f.name for f in font_manager.fontManager.ttflist}
    return _FONT_BASE + ([f for f in _FONT_EXTRA if f in avail and f not in _FONT_BASE] if extra else [])


FONT_FAMILIES = available_fonts()
PALETTES = {
    "Okabe–Ito (sample families)": None,       # v2.10 behaviour: fixed sample-family colours, then Okabe–Ito
    "Okabe–Ito": ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"],
    "Tol bright": ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"],
    "Tol muted": ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"],
    "Tol vibrant": ["#0077BB", "#33BBEE", "#009988", "#EE7733", "#CC3311", "#EE3377", "#BBBBBB"],
    "viridis (ordered series)": ["#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725"],
    "matplotlib tab10": ["#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"],
    "greyscale": ["#000000", "#555555", "#888888", "#AAAAAA", "#333333", "#777777"],
}
# journal presets: width key + style values (sizes from the author guidelines: final-size fonts 5-8 pt (Nature 5-7), lines >= 0.5 pt)
JOURNAL_PRESETS = {
    "ACS (J. Am. Chem. Soc., Chem. Mater. …)": dict(width="ACS single column (3.25 in)", font_size=8.0, family="Arial", line_width=1.0, axes_width=0.75),
    "Wiley-VCH (Angew., Adv. Mater. …)": dict(width="Wiley single column (3.35 in / 85 mm)", font_size=8.0, family="Arial", line_width=1.0, axes_width=0.75),
    "RSC (J. Mater. Chem., Chem. Sci. …)": dict(width="RSC single column (3.27 in / 8.3 cm)", font_size=7.0, family="Arial", line_width=0.9, axes_width=0.6),
    "Elsevier": dict(width="Elsevier single column (3.54 in / 90 mm)", font_size=8.0, family="Arial", line_width=1.0, axes_width=0.75),
    "Nature / Springer Nature": dict(width="Nature single column (3.50 in / 89 mm)", font_size=7.0, family="Helvetica", line_width=0.8, axes_width=0.6),
    "Science": dict(width="Science single column (2.24 in / 5.7 cm)", font_size=7.0, family="Helvetica", line_width=0.8, axes_width=0.6),
    "Thesis (A4, 1 column)": dict(width="1.5 column (5.0 in)", font_size=10.0, family="auto", line_width=1.2, axes_width=0.9),
    "Slides / poster": dict(width="Slide (10 in, 16:9)", font_size=16.0, family="auto", line_width=2.2, axes_width=1.5, marker_size=7.0),
}


def set_global(**kw):
    """Update the global style (unknown keys ignored) and re-apply it."""
    for k, v in kw.items():
        if k in GLOBAL_DEFAULTS and v is not None:
            GLOBAL[k] = v
    apply_style()


def palette():
    """The active palette as a list of colours (Okabe–Ito pool for the default sample-family mode)."""
    p = PALETTES.get(GLOBAL.get("palette"))
    return list(p) if p else [c for c in OKABE_ITO + _extra if c != "#7F7F7F"]


def apply_style(font_size: float | None = None):
    G = GLOBAL
    font_size = float(G["font_size"] if font_size is None else font_size)
    fam = _font_family(G.get("family"))
    lw, aw, tdir = float(G["line_width"]), float(G["axes_width"]), G.get("tick_dir", "in")
    k = aw / 0.75                                            # tick sizes/widths scale with the axes line width
    serif = fam in SERIF_FAMILIES
    plt.rcParams.update({
        "font.family": "serif" if serif else "sans-serif",
        "font.sans-serif": [fam, "Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": font_size,
        "axes.labelsize": font_size,
        "axes.titlesize": font_size,
        "xtick.labelsize": font_size - 1,
        "ytick.labelsize": font_size - 1,
        "legend.fontsize": font_size - 1,
        "legend.frameon": bool(G.get("legend_frame", False)),
        "legend.handlelength": 1.5,
        "axes.linewidth": aw,
        "lines.linewidth": lw,
        "lines.markersize": float(G["marker_size"]),
        "xtick.direction": tdir, "ytick.direction": tdir,
        "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 3.0 * k, "ytick.major.size": 3.0 * k,
        "xtick.minor.size": 1.6 * k, "ytick.minor.size": 1.6 * k,
        "xtick.major.width": aw, "ytick.major.width": aw,
        "xtick.minor.width": 0.6 * k, "ytick.minor.width": 0.6 * k,
        "xtick.minor.visible": bool(G.get("minor", True)), "ytick.minor.visible": bool(G.get("minor", True)),
        "axes.prop_cycle": (matplotlib.cycler(color=palette()) if PALETTES.get(G.get("palette"))
                            else matplotlib.rcParamsDefault["axes.prop_cycle"]),
        "axes.grid": False,
        "axes.unicode_minus": True,
        "mathtext.default": "regular",
        "mathtext.fontset": "custom",
        "mathtext.rm": fam, "mathtext.it": f"{fam}:italic", "mathtext.bf": f"{fam}:bold", "mathtext.sf": fam,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "svg.fonttype": "none", "svg.hashsalt": "mof-figure-studio",       # editable text in Illustrator/Inkscape
        "pdf.fonttype": 42, "ps.fonttype": 42,  # embed TrueType
        "figure.dpi": 150,
        "errorbar.capsize": 2.0,
    })
    if serif:
        plt.rcParams["font.serif"] = [fam, "Times New Roman", "Liberation Serif", "DejaVu Serif"]


apply_style()


def new_figure(width_in: float = SINGLE, aspect: float = 0.75, nrows=1, ncols=1, **kw):
    apply_style(kw.pop("font_size", None))
    fig, ax = plt.subplots(nrows, ncols, figsize=(width_in, width_in * aspect), **kw)
    return fig, ax


def panel_label(ax, letter: str, x=-0.18, y=1.02, bold=True, paren=True):
    txt = f"({letter})" if paren else letter
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=plt.rcParams["font.size"] + 1,
            fontweight="bold" if bold else "normal", va="bottom", ha="left")


def no_y_ticks(ax):
    """For offset/stacked spectra: arbitrary-unit axes carry no tick labels (journal convention)."""
    ax.set_yticks([])
    ax.tick_params(axis="y", which="both", left=False, right=False)


def export(fig, fmt: str, dpi: int = 600) -> bytes:
    buf = io.BytesIO()
    kw = {}
    if fmt in ("png", "tiff", "tif"):
        kw["dpi"] = dpi
    if fmt in ("tiff", "tif"):
        kw["pil_kwargs"] = {"compression": "tiff_lzw"}
        fmt = "tiff"
    if fmt == "pdf":
        kw["metadata"] = {"CreationDate": None, "Creator": "Reticula"}
    elif fmt == "svg":
        kw["metadata"] = {"Date": None, "Creator": "Reticula"}
    fig.savefig(buf, format=fmt, **kw)
    return buf.getvalue()


def save_all(fig, stem: str, formats=("png", "svg", "pdf", "tiff"), dpi=600):
    paths = []
    for f in formats:
        p = f"{stem}.{f}"
        with open(p, "wb") as fh:
            fh.write(export(fig, f, dpi))
        paths.append(p)
    return paths
