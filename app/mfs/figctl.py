"""Shared figure customisation (v2): one settings dict -> applied to ANY finished matplotlib figure.

Everything works post hoc on the drawn Figure, so every technique page (and the composer) gets the same options:
title, axis titles, axis ranges, display units (via tick re-labelling - the data are not touched), tick sides,
tick labels per side, direction, major/minor increments, a rebuilt legend (rename / reorder / hide entries,
position incl. outside, columns, font size, frame, marker scale) and JSON presets ("house style")."""
from __future__ import annotations
import copy, json, math, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend import Legend
from matplotlib.ticker import (AutoMinorLocator, FixedFormatter, FixedLocator, MaxNLocator, MultipleLocator,
                               NullLocator, NullFormatter, ScalarFormatter)

SIDES = ("bottom", "top", "left", "right")
TICK_PRESETS = {
    "All 4 sides": dict(bottom=True, top=True, left=True, right=True),
    "3 sides (no top)": dict(bottom=True, top=False, left=True, right=True),
    "3 sides (no right)": dict(bottom=True, top=True, left=True, right=False),
    "3 sides (no left)": dict(bottom=True, top=True, left=False, right=True),
    "3 sides (no bottom)": dict(bottom=False, top=True, left=True, right=True),
    "2 sides (bottom + left)": dict(bottom=True, top=False, left=True, right=False),
    "2 sides (bottom + right)": dict(bottom=True, top=False, left=False, right=True),
    "2 sides (top + bottom)": dict(bottom=True, top=True, left=False, right=False),
    "1 side (bottom)": dict(bottom=True, top=False, left=False, right=False),
    "1 side (left)": dict(bottom=False, top=False, left=True, right=False),
    "Custom": None,
}
LEGEND_LOCS = ["best", "upper right", "upper left", "lower left", "lower right", "right", "center left",
               "center right", "lower center", "upper center", "center", "outside right", "outside top",
               "outside bottom", "outside left"]


def default_settings():
    return {
        "version": 2,
        "title": "", "title_size": 0.0,
        "mathify": True,                       # turn ^{..} / _{..} / ^-1 into super/subscripts
        "font": {"label_size": 0.0, "tick_size": 0.0, "family": ""},   # 0 = keep; family "" = the figure style's font
        "ticks": {"preset": "All 4 sides", "sides": dict(bottom=True, top=True, left=True, right=True),
                  "labels": dict(bottom=True, top=False, left=True, right=False),
                  "direction": "in", "minor": True, "major_len": 3.5, "minor_len": 2.0, "width": 0.6,
                  "frame": "box"},     # box | match ticks (hide spines without ticks)
        "panels": {},   # "0": {xlabel, ylabel, xmin, xmax, ymin, ymax, x_unit, y_unit, xmajor, ymajor, xminor, yminor, xlog, ylog}
        "legend": {"custom": False, "show": True, "loc": "best", "ncol": 1, "fontsize": 0.0, "frame": False,
                   "markerscale": 1.0, "title": "", "panel": 0, "entries": [], "handlelength": 1.5,
                   "title_size": 0.0, "family": ""},
        "units_params": {"wavelength": 1.5406, "m0_mg": 10.0},
    }


def merge(base, over):
    """Deep-merge preset `over` into settings `base` (unknown keys kept, missing keys default)."""
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


HOUSE_KEYS = ("mathify", "font", "ticks", "title_size")
HOUSE_LEGEND_KEYS = ("loc", "ncol", "fontsize", "frame", "markerscale", "handlelength", "custom", "show", "title_size", "family")


def house_style(cfg):
    """Only the parts that make sense across figures (no axis labels, ranges or legend entry names)."""
    h = {k: copy.deepcopy(cfg[k]) for k in HOUSE_KEYS if k in cfg}
    h["legend"] = {k: cfg["legend"][k] for k in HOUSE_LEGEND_KEYS if k in cfg.get("legend", {})}
    h["kind"] = "Reticula house style"
    return h


def to_json(cfg) -> str:
    return json.dumps(cfg, indent=2, ensure_ascii=False)


def from_json(s) -> dict:
    d = json.loads(s)
    return merge(default_settings(), d)


# ----------------------------------------------------------------------------------------------- units
def _deg_to_d(x, lam):
    x = np.asarray(x, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return lam / (2 * np.sin(np.radians(x) / 2))


def _d_to_deg(d, lam):
    d = np.asarray(d, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 2 * np.degrees(np.arcsin(np.clip(lam / (2 * d), -1, 1)))


UNITS = {
    # family: [(name, axis label or None (=keep quantity, swap unit), unit text, fwd(x, P), inv(y, P))]
    "2theta": [("2θ (°)", None, "°", None, None),
               ("d-spacing (Å)", "$d$-spacing (Å)", "Å", lambda x, P: _deg_to_d(x, P["wavelength"]), lambda y, P: _d_to_deg(y, P["wavelength"])),
               ("q (Å⁻¹)", r"$q$ (Å$^{-1}$)", "Å$^{-1}$", lambda x, P: 4 * np.pi * np.sin(np.radians(np.asarray(x, float)) / 2) / P["wavelength"],
                lambda y, P: 2 * np.degrees(np.arcsin(np.clip(np.asarray(y, float) * P["wavelength"] / (4 * np.pi), -1, 1))))],
    "wavenumber": [("cm⁻¹", None, "cm$^{-1}$", None, None),
                   ("µm (wavelength)", "Wavelength (µm)", "µm", lambda x, P: 1e4 / np.asarray(x, float), lambda y, P: 1e4 / np.asarray(y, float)),
                   ("nm (wavelength)", "Wavelength (nm)", "nm", lambda x, P: 1e7 / np.asarray(x, float), lambda y, P: 1e7 / np.asarray(y, float))],
    "uptake": [("cm³ g⁻¹ STP", None, "cm$^3$ g$^{-1}$ STP", None, None),
               ("mmol g⁻¹", None, "mmol g$^{-1}$", lambda x, P: np.asarray(x, float) / 22.414, lambda y, P: np.asarray(y, float) * 22.414),
               ("mg g⁻¹ (N₂)", None, "mg g$^{-1}$", lambda x, P: np.asarray(x, float) / 22.414 * 28.0134, lambda y, P: np.asarray(y, float) / 28.0134 * 22.414)],
    "temperature": [("°C", None, "°C", None, None),
                    ("K", None, "K", lambda x, P: np.asarray(x, float) + 273.15, lambda y, P: np.asarray(y, float) - 273.15)],
    "mass_pct": [("wt%", None, "%", None, None),
                 ("mg (uses initial mass)", "Mass (mg)", "mg", lambda x, P: np.asarray(x, float) / 100 * P["m0_mg"], lambda y, P: np.asarray(y, float) * 100 / P["m0_mg"])],
    "time_h": [("h", None, "h", None, None),
               ("min", None, "min", lambda x, P: np.asarray(x, float) * 60, lambda y, P: np.asarray(y, float) / 60),
               ("s", None, "s", lambda x, P: np.asarray(x, float) * 3600, lambda y, P: np.asarray(y, float) / 3600),
               ("d", None, "d", lambda x, P: np.asarray(x, float) / 24, lambda y, P: np.asarray(y, float) * 24)],
    "pore_nm": [("nm", None, "nm", None, None),
                ("Å", None, "Å", lambda x, P: np.asarray(x, float) * 10, lambda y, P: np.asarray(y, float) / 10)],
    "conc_mM": [("mM", None, "mM", None, None),
                ("µM", None, "µM", lambda x, P: np.asarray(x, float) * 1000, lambda y, P: np.asarray(y, float) / 1000)],
}


def unit_names(family):
    return [u[0] for u in UNITS.get(family, [])]


def _unit(family, name):
    for u in UNITS.get(family, []):
        if u[0] == name:
            return u
    return None


def _swap_unit_in_label(label, new_unit):
    """'Quantity adsorbed (cm$^3$ g$^{-1}$ STP)' -> 'Quantity adsorbed (mmol g$^{-1}$)'."""
    m = re.match(r"^(.*?)(\s*\([^()]*\)\s*)$", label or "")
    base = m.group(1) if m else (label or "")
    return f"{base} ({new_unit})" if base else new_unit


# ----------------------------------------------------------------------------------------------- text
def mathify(s: str) -> str:
    """Free text -> mathtext for sub/superscripts: 'cm^-1' -> 'cm$^{-1}$', 'ZrO_2' -> 'ZrO$_{2}$', 'x^{2+}'.
    Text that already contains $...$ is left alone."""
    if not s or "$" in s:
        return s
    s = re.sub(r"\^\{([^}]*)\}", r"$^{\1}$", s)
    s = re.sub(r"_\{([^}]*)\}", r"$_{\1}$", s)
    s = re.sub(r"\^([-+−]?[0-9A-Za-z.]+)", lambda m: "$^{" + m.group(1).replace("−", "-") + "}$", s)
    s = re.sub(r"_([0-9A-Za-z])", r"$_{\1}$", s)
    return s.replace("$$", "")


def safe_text(s: str) -> str:
    """Never let a mistyped mathtext crash the figure: if it does not parse, show it literally."""
    if not s or "$" not in s:
        return s
    from matplotlib.mathtext import MathTextParser
    try:
        MathTextParser("path").parse(s)
        return s
    except Exception:
        return s.replace("$", r"\$")


def _txt(cfg, s):
    return safe_text(mathify(s) if cfg.get("mathify", True) else s)


# ----------------------------------------------------------------------------------------------- helpers
def data_axes(fig):
    """Axes that carry data (skip colour bars / empty helper axes), in drawing order."""
    out = []
    for ax in fig.axes:
        if ax.get_label() == "<colorbar>":
            continue
        if not ax.axison:
            # image panels (SEM) have axis off: still listed so titles/legends work, but ticks are skipped
            pass
        out.append(ax)
    return out


def is_twin(ax, fig):
    """A twinx/twiny axis drawn on top of another one (e.g. DTG on the TGA axis)."""
    for other in fig.axes:
        if other is ax:
            continue
        if np.allclose(other.get_position().bounds, ax.get_position().bounds) and fig.axes.index(other) < fig.axes.index(ax):
            return True
    return False


def describe(fig):
    """Current state of every panel (labels, limits, scale, tick visibility) - used as defaults/placeholders."""
    out = []
    for i, ax in enumerate(data_axes(fig)):
        out.append({"index": i, "xlabel": ax.get_xlabel(), "ylabel": ax.get_ylabel(), "title": ax.get_title(),
                    "xlim": tuple(float(v) for v in ax.get_xlim()), "ylim": tuple(float(v) for v in ax.get_ylim()),
                    "xscale": ax.get_xscale(), "yscale": ax.get_yscale(), "axison": bool(ax.axison),
                    "yticks_hidden": len(ax.get_yticks()) == 0, "xticks_hidden": len(ax.get_xticks()) == 0,
                    "twin": is_twin(ax, fig), "ylabel_right": ax.yaxis.get_label_position() == "right"})
    return out


def legend_entries(fig):
    """(label, handle) pairs from existing legends (incl. figure legends / proxies) + labelled artists."""
    seen, out = set(), []
    legs = [l for ax in fig.axes for l in ax.get_children() if isinstance(l, Legend)] + list(fig.legends)
    for lg in legs:
        hs = getattr(lg, "legend_handles", None) or getattr(lg, "legendHandles", [])
        for h, t in zip(hs, lg.get_texts()):
            lab = t.get_text()
            if lab and lab not in seen:
                seen.add(lab); out.append((lab, h))
    for ax in fig.axes:
        for h, lab in zip(*ax.get_legend_handles_labels()):
            if lab and not lab.startswith("_") and lab not in seen:
                seen.add(lab); out.append((lab, h))
    return out


def direct_labels(fig):
    """Text labels placed next to traces (stacked spectra) - they can be renamed like legend entries."""
    out = []
    for ax in fig.axes:
        for t in ax.texts:
            s = t.get_text()
            if s and not re.fullmatch(r"\(?[a-zA-Z]\)?", s) and len(s) < 60 and t.get_gid() != "mfs_note":
                out.append(s)
    return list(dict.fromkeys(out))


def nice_step(span, n=5):
    if not span or not np.isfinite(span) or span <= 0:
        return 1.0
    raw = span / n
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag:
            return m * mag
    return 10 * mag


def _fmt(v, step):
    dec = max(0, -int(math.floor(math.log10(step))) if step > 0 else 0)
    if step >= 1 and abs(step - round(step)) < 1e-9:
        dec = 0
    elif step > 0 and abs(step * 10 ** dec - round(step * 10 ** dec)) > 1e-6:
        dec += 1
    s = f"{v:.{dec}f}"
    return "0" if re.fullmatch(r"-?0\.?0*", s) else s.replace("-", "−")


# ----------------------------------------------------------------------------------------------- apply
def _axis_apply(ax, axis, p, cfg, minor_on):
    """Limits, units and tick increments for one axis ('x' or 'y') of one panel."""
    A = ax.xaxis if axis == "x" else ax.yaxis
    getlim, setlim = (ax.get_xlim, ax.set_xlim) if axis == "x" else (ax.get_ylim, ax.set_ylim)
    scale = ax.get_xscale() if axis == "x" else ax.get_yscale()
    if p.get(f"{axis}log") is not None and p.get(f"{axis}log") != (scale == "log"):
        old_lim = getlim()
        try:
            (ax.set_xscale if axis == "x" else ax.set_yscale)("log" if p.get(f"{axis}log") else "linear")
            A.get_majorticklocs()                # a log axis without positive data only fails here, not in set_*scale
            scale = "log" if p.get(f"{axis}log") else "linear"
        except Exception:                       # a log axis needs positive data: keep the axis as it was instead of failing
            (ax.set_xscale if axis == "x" else ax.set_yscale)(scale)
            setlim(*old_lim)
    hidden = len(A.get_majorticklocs()) == 0
    fam, uname = p.get(f"{axis}_family"), p.get(f"{axis}_unit")
    u = _unit(fam, uname) if fam and uname else None
    conv = u if (u and u[3] is not None) else None
    P = cfg.get("units_params", {})
    lo, hi = getlim()
    inverted = lo > hi
    vmin, vmax = p.get(f"{axis}min"), p.get(f"{axis}max")
    if conv:
        if vmin is not None or vmax is not None:
            cur = sorted([float(conv[3](lo, P)), float(conv[3](hi, P))])
            a = conv[4](vmin, P) if vmin is not None else None
            b = conv[4](vmax, P) if vmax is not None else None
            # map back; keep the original direction of the data axis
            dlo, dhi = sorted([v for v in (a if a is not None else lo, b if b is not None else hi)])
            setlim((dhi, dlo) if inverted else (dlo, dhi))
    else:
        if vmin is not None or vmax is not None:
            nlo = vmin if vmin is not None else min(lo, hi)
            nhi = vmax if vmax is not None else max(lo, hi)
            setlim((nhi, nlo) if inverted else (nlo, nhi))
    major, minor = p.get(f"{axis}major") or 0, p.get(f"{axis}minor") or 0
    if scale == "log":
        if not minor_on:
            A.set_minor_locator(NullLocator())
        return
    if hidden and not major:
        A.set_minor_locator(NullLocator())
        return
    lo, hi = getlim()
    if conv:
        f = lambda v: conv[3](v, P)
        g = lambda v: conv[4](v, P)
        with np.errstate(all="ignore"):
            d0, d1 = sorted([float(f(lo)), float(f(hi))])
        if not np.isfinite(d0):
            d0 = float(np.nanmin(f(np.linspace(min(lo, hi), max(lo, hi), 200)[1:])))
        if not np.isfinite(d1):
            d1 = float(np.nanmax(f(np.linspace(min(lo, hi), max(lo, hi), 200)[1:])))
        step = major or nice_step(d1 - d0)
        ticks = np.arange(math.ceil(d0 / step - 1e-9) * step, d1 + step * 1e-6, step)
        if len(ticks) > 60:
            step = nice_step(d1 - d0); ticks = np.arange(math.ceil(d0 / step) * step, d1 + step * 1e-6, step)
        pos = g(ticks)
        A.set_major_locator(FixedLocator(pos))
        A.set_major_formatter(FixedFormatter([_fmt(t, step) for t in ticks]))
        if minor_on:
            ms = minor or step / 2
            mt = np.arange(math.ceil(d0 / ms) * ms, d1 + ms * 1e-6, ms)
            A.set_minor_locator(FixedLocator(g(mt)) if len(mt) < 400 else NullLocator())
        else:
            A.set_minor_locator(NullLocator())
        lab = p.get(f"{axis}label") or (u[1] if u[1] else _swap_unit_in_label(ax.get_xlabel() if axis == "x" else ax.get_ylabel(), u[2]))
        (ax.set_xlabel if axis == "x" else ax.set_ylabel)(_txt(cfg, lab) if p.get(f"{axis}label") else lab)
        return
    if major:
        span = abs(hi - lo)
        if span / major > 200:
            major = nice_step(span)
        A.set_major_locator(MultipleLocator(major))
        if not isinstance(A.get_major_formatter(), (ScalarFormatter,)):
            A.set_major_formatter(ScalarFormatter())
    if minor_on:
        if minor:
            A.set_minor_locator(MultipleLocator(minor) if abs(hi - lo) / minor < 400 else AutoMinorLocator())
        else:
            A.set_minor_locator(AutoMinorLocator())
    else:
        A.set_minor_locator(NullLocator())


def apply_panel(ax, cfg, families=None, idx=0):
    """v2.14: the per-panel settings made under a single figure (axis titles, units, ranges, log, increments, title) on ONE axes - used
    when that figure becomes a panel of a combined figure. Ticks, fonts and legend follow the combined figure's own panel."""
    cfg = merge(default_settings(), cfg or {})
    p = dict(cfg["panels"].get(str(idx), {}))
    fam = (families or {}).get(idx) or (families or {}).get("all") or {}
    p.setdefault("x_family", fam.get("x")); p.setdefault("y_family", fam.get("y"))
    if p.get("xlabel") and not (p.get("x_unit") and _unit(p.get("x_family"), p.get("x_unit")) and _unit(p.get("x_family"), p.get("x_unit"))[3]):
        ax.set_xlabel(_txt(cfg, p["xlabel"]))
    if p.get("ylabel") and not (p.get("y_unit") and _unit(p.get("y_family"), p.get("y_unit")) and _unit(p.get("y_family"), p.get("y_unit"))[3]):
        ax.set_ylabel(_txt(cfg, p["ylabel"]))
    _axis_apply(ax, "x", p, cfg, cfg["ticks"].get("minor", True))
    _axis_apply(ax, "y", p, cfg, cfg["ticks"].get("minor", True))
    if cfg.get("title"):
        ax.set_title(_txt(cfg, cfg["title"]), fontsize=cfg.get("title_size") or None)
    return ax


def apply(fig, cfg, families=None):
    """Apply a settings dict to a finished figure. families: {panel_index or 'all': {'x': fam, 'y': fam}}."""
    cfg = merge(default_settings(), cfg or {})
    axes = data_axes(fig)
    T = cfg["ticks"]
    sides = dict(TICK_PRESETS.get(T.get("preset")) or T.get("sides") or {})
    for s in SIDES:
        sides.setdefault(s, True)
    labels = dict(T.get("labels") or {})
    fs = cfg["font"]
    multi = len([a for a in axes if a.axison and not is_twin(a, fig)]) > 1
    for i, ax in enumerate(axes):
        p = dict(cfg["panels"].get(str(i), {}))
        fam = (families or {}).get(i) or (families or {}).get("all") or {}
        p.setdefault("x_family", fam.get("x")); p.setdefault("y_family", fam.get("y"))
        if not ax.axison:
            continue
        twin = is_twin(ax, fig)
        if p.get("xlabel") and not (p.get("x_unit") and _unit(p.get("x_family"), p.get("x_unit")) and _unit(p.get("x_family"), p.get("x_unit"))[3]):
            ax.set_xlabel(_txt(cfg, p["xlabel"]))
        if p.get("ylabel") and not (p.get("y_unit") and _unit(p.get("y_family"), p.get("y_unit")) and _unit(p.get("y_family"), p.get("y_unit"))[3]):
            ax.set_ylabel(_txt(cfg, p["ylabel"]))
        _axis_apply(ax, "x", p, cfg, T.get("minor", True))
        _axis_apply(ax, "y", p, cfg, T.get("minor", True))
        kw = dict(direction=T.get("direction", "in"), width=T.get("width", 0.6))
        if twin:
            # twin y-axis (e.g. DTG): keeps its right-hand ticks and labels; only the style follows the preset
            ax.tick_params(axis="y", which="both", **kw)
            ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False, labeltop=False)
        else:
            right_twin = any(is_twin(o, fig) and np.allclose(o.get_position().bounds, ax.get_position().bounds) for o in axes if o is not ax)
            sd = dict(sides); lb = {s: labels.get(s, s in ("bottom", "left")) and sd[s] for s in SIDES}
            # panels that share an axis (e.g. TGA above DTG) keep their inner tick labels hidden
            if not ax.xaxis._major_tick_kw.get("label1On", True) and ax.get_shared_x_axes().get_siblings(ax) != [ax]:
                lb["bottom"] = False
            if not ax.yaxis._major_tick_kw.get("label1On", True) and ax.get_shared_y_axes().get_siblings(ax) != [ax]:
                lb["left"] = False
            if right_twin:
                sd["right"] = False; lb["right"] = False
            ax.tick_params(axis="x", which="both", bottom=sd["bottom"], top=sd["top"],
                           labelbottom=lb["bottom"], labeltop=lb["top"], **kw)
            ax.tick_params(axis="y", which="both", left=sd["left"], right=sd["right"],
                           labelleft=lb["left"], labelright=lb["right"], **kw)
            if lb["top"] and not lb["bottom"]:
                ax.xaxis.set_label_position("top")
            if lb["right"] and not lb["left"]:
                ax.yaxis.set_label_position("right")
            if T.get("frame") == "match ticks":
                for s in SIDES:
                    ax.spines[s].set_visible(bool(sd[s]) or (s == "right" and right_twin))
            else:
                for s in SIDES:
                    ax.spines[s].set_visible(True)
        ax.tick_params(which="major", length=T.get("major_len", 3.5))
        ax.tick_params(which="minor", length=T.get("minor_len", 2.0))
        if fs.get("tick_size"):
            ax.tick_params(which="both", labelsize=fs["tick_size"])
        if fs.get("label_size"):
            ax.xaxis.label.set_size(fs["label_size"]); ax.yaxis.label.set_size(fs["label_size"])
    # title
    if cfg.get("title"):
        t = _txt(cfg, cfg["title"])
        size = cfg.get("title_size") or plt.rcParams["font.size"] + 1
        if multi:
            fig.suptitle(t, fontsize=size)
        else:
            next(a for a in axes if a.axison or True).set_title(t, fontsize=size)
    _legend_apply(fig, cfg)
    _rename_direct(fig, cfg)
    _family_apply(fig, cfg)
    return fig


def _set_family(artists, fam):
    for t in artists:
        try:
            t.set_fontfamily(fam)
        except Exception:
            pass


def _family_apply(fig, cfg):
    """v2.14: one font family for everything in the figure (axis titles, ticks, texts, legend), and an own one for the legend."""
    from matplotlib.text import Text
    fam = (cfg.get("font") or {}).get("family") or ""
    if fam and fam != "auto":
        from . import style as _S
        fam = _S._font_family(fam)
        serif = fam in _S.SERIF_FAMILIES
        for ax in fig.axes:
            ax.get_xticklabels(); ax.get_yticklabels()                      # make sure the tick labels exist
            _set_family([t for t in ax.findobj(Text)], fam)
        _set_family([t for t in fig.texts], fam)
        for lg in list(fig.legends) + [l for ax in fig.axes for l in ax.get_children() if isinstance(l, Legend)]:
            _set_family(lg.get_texts() + [lg.get_title()], fam)
    L = cfg.get("legend") or {}
    lfam = L.get("family") or ""
    legs = list(fig.legends) + [l for ax in fig.axes for l in ax.get_children() if isinstance(l, Legend)]
    for lg in legs:
        if lfam and lfam != "auto":
            from . import style as _S
            _set_family(lg.get_texts() + [lg.get_title()], _S._font_family(lfam))
        if L.get("title_size"):
            lg.get_title().set_fontsize(float(L["title_size"]))


def _rename_direct(fig, cfg):
    L = cfg["legend"]
    if not L.get("entries"):
        return
    ren = {e["label"]: e for e in L["entries"]}
    for ax in fig.axes:
        for t in ax.texts:
            e = ren.get(t.get_text())
            if e:
                if not e.get("show", True):
                    t.set_visible(False)
                elif e.get("text") and e["text"] != e["label"]:
                    t.set_text(_txt(cfg, e["text"]))


def _needs_rebuild(L):
    """Settings that need a rebuilt legend (position, columns, marker size, panel, hiding, renaming). Font size, frame, title,
    title size and font family act on the legend that is already there - they work without ticking 'Customise the legend'."""
    return bool(L.get("custom") or L.get("show") is False or int(L.get("ncol", 1) or 1) != 1 or L.get("loc", "best") != "best"
                or abs(float(L.get("markerscale", 1) or 1) - 1.0) > 1e-9 or int(L.get("panel", 0) or 0) != 0
                or abs(float(L.get("handlelength", 1.5) or 1.5) - 1.5) > 1e-9
                or any((e.get("text") and e["text"] != e["label"]) or not e.get("show", True) for e in (L.get("entries") or [])))


def _existing_legends(fig):
    return list(fig.legends) + [l for ax in fig.axes for l in ax.get_children() if isinstance(l, Legend)]


def _direct_fontsize(fig, L):
    """The legend font size also sizes the sample labels written next to the traces (they are the legend of those figures)."""
    fs = L.get("fontsize")
    if not fs:
        return
    for ax in fig.axes:
        for t in ax.texts:
            s = t.get_text()
            if s and not re.fullmatch(r"\(?[a-zA-Z]\)?", s) and len(s) < 60 and t.get_gid() != "mfs_note":
                t.set_fontsize(float(fs))


def _legend_light(fig, cfg):
    L = cfg["legend"]
    for lg in _existing_legends(fig):
        if L.get("fontsize"):
            for t in lg.get_texts():
                t.set_fontsize(float(L["fontsize"]))
            if not L.get("title_size"):
                lg.get_title().set_fontsize(float(L["fontsize"]))
        if L.get("title"):
            lg.set_title(_txt(cfg, L["title"]))
        if L.get("frame"):
            lg.set_frame_on(True)
            fr = lg.get_frame(); fr.set_edgecolor("k"); fr.set_linewidth(0.6); fr.set_facecolor("white"); fr.set_alpha(1)
    _direct_fontsize(fig, L)


def _legend_apply(fig, cfg):
    L = cfg["legend"]
    if not _needs_rebuild(L):
        _legend_light(fig, cfg)
        return
    _direct_fontsize(fig, L)
    entries = legend_entries(fig)
    for ax in fig.axes:
        for lg in [c for c in ax.get_children() if isinstance(c, Legend)]:
            lg.remove()
    for lg in list(fig.legends):
        lg.remove()
    if not L.get("show", True) or not entries:
        return
    handles = dict(entries)
    order = L.get("entries") or [{"label": l, "text": l, "show": True} for l, _ in entries]
    known = {e["label"] for e in order}
    order = list(order) + [{"label": l, "text": l, "show": True} for l, _ in entries if l not in known]
    hs, ls = [], []
    for e in order:
        if e.get("show", True) and e["label"] in handles:
            hs.append(handles[e["label"]]); ls.append(_txt(cfg, e.get("text") or e["label"]))
    if not hs:
        return
    axes = [a for a in data_axes(fig) if a.axison] or fig.axes
    ax = axes[min(int(L.get("panel", 0) or 0), len(axes) - 1)]
    kw = dict(ncol=int(L.get("ncol", 1) or 1), frameon=bool(L.get("frame")), markerscale=float(L.get("markerscale", 1) or 1),
              fontsize=L.get("fontsize") or plt.rcParams["font.size"] - 1, handlelength=L.get("handlelength", 1.5),
              title=_txt(cfg, L.get("title")) or None, handletextpad=0.4, columnspacing=0.8, labelspacing=0.3)
    if L.get("frame"):
        kw.update(fancybox=False, edgecolor="k", framealpha=1)
    loc = L.get("loc", "best")
    if loc == "outside right":
        lg = ax.legend(hs, ls, loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0, **kw)
    elif loc == "outside left":
        lg = ax.legend(hs, ls, loc="center right", bbox_to_anchor=(-0.25, 0.5), borderaxespad=0, **kw)
    elif loc == "outside top":
        lg = ax.legend(hs, ls, loc="lower center", bbox_to_anchor=(0.5, 1.02), borderaxespad=0, **kw)
    elif loc == "outside bottom":
        lg = ax.legend(hs, ls, loc="upper center", bbox_to_anchor=(0.5, -0.22), borderaxespad=0, **kw)
    else:
        lg = ax.legend(hs, ls, loc=loc, **kw)
    if L.get("frame"):
        lg.get_frame().set_linewidth(0.6)
    return lg
