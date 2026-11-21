"""v2.12: movable / editable sample labels on every figure.

Movable items: the text labels written in data coordinates next to traces (stacked spectra, TGA, DSC …), other
annotation texts, and each legend (as a whole). Overrides are stored per figure stem in
st.session_state['lbl|<stem>'] = {item key: {"dx", "dy", "text"} | {"lx", "ly"} (legend, axes fraction)}.
A label's (dx, dy) is in DATA units relative to its automatic position, which itself sits on the trace's offset,
so a moved label follows its trace when the stacking offset changes. "Same place for every sample" stores, for every
label of the panel, the position that has the same distance to its own trace (baseline = minimum of the trace with the
label's colour) as the label that was moved. Overrides are applied to the figure before preview AND export."""
from __future__ import annotations
import base64, io, os, re
import numpy as np

PREFIX = "lbl|"
ORIG = "_mfs_orig_text"


def tag_originals(fig):
    """Remember every text's original string (before renaming by the customisation panel) = its stable key."""
    for ax in fig.axes:
        for t in ax.texts:
            if not hasattr(t, ORIG):
                setattr(t, ORIG, t.get_text())


def _movable_texts(fig):
    out = []
    for ai, ax in enumerate(fig.axes):
        seen = {}
        for t in ax.texts:
            s = getattr(t, ORIG, t.get_text())
            if not s or not t.get_visible() or re.fullmatch(r"\(?[a-zA-Z]\)?", s) or t.get_gid() in ("mfs_note", "mfs_fixed"):
                continue
            n = seen.get(s, 0); seen[s] = n + 1
            out.append((f"{ai}:{s}" + (f"#{n}" if n else ""), ai, ax, t))
    return out


def items(fig):
    """[(key, kind, ax_index, artist)] kind = 'data' (data-coordinate label), 'axes' (axes-fraction text), 'legend'."""
    res = []
    for k, ai, ax, t in _movable_texts(fig):
        res.append((k, "data" if t.get_transform() == ax.transData else ("axes" if t.get_transform() == ax.transAxes else "other"), ai, t))
    for ai, ax in enumerate(fig.axes):
        lg = ax.get_legend()
        if lg is not None and lg.get_visible():
            res.append((f"{ai}:__legend__", "legend", ai, lg))
    return [r for r in res if r[1] != "other"]


def _trace_base(ax, t):
    """Baseline (min y) of the line drawn in the label's colour, or None."""
    from matplotlib.colors import to_rgba
    try:
        c = to_rgba(t.get_color())
    except Exception:
        return None
    best = None
    for ln in ax.get_lines():
        try:
            if np.allclose(to_rgba(ln.get_color())[:3], c[:3], atol=0.02):
                y = np.asarray(ln.get_ydata(), float); y = y[np.isfinite(y)]
                if y.size > 2:
                    best = float(np.min(y)); break
        except Exception:
            continue
    return best


def apply(fig, ov):
    """Apply stored overrides to the figure (positions, texts, legend placement)."""
    if not ov:
        return fig
    for k, kind, ai, art in items(fig):
        o = ov.get(k)
        if not o:
            continue
        if kind == "legend":
            if "lx" in o:
                art.set_bbox_to_anchor((o["lx"], o["ly"]), transform=fig.axes[ai].transAxes); art._loc = 2   # upper left
            continue
        if o.get("text") is not None:
            art.set_text(o["text"])
        if o.get("dx") or o.get("dy"):
            if not hasattr(art, "_mfs_xy0"):
                art._mfs_xy0 = tuple(art.get_position())
            x0, y0 = art._mfs_xy0
            art.set_position((x0 + float(o.get("dx", 0)), y0 + float(o.get("dy", 0))))
    return fig


def _frac_to(ax, kind, fig, fx, fy):
    """Figure fraction -> data (kind 'data') or axes fraction ('axes', 'legend')."""
    disp = fig.transFigure.transform((fx, fy))
    return tuple((ax.transData if kind == "data" else ax.transAxes).inverted().transform(disp))


def handle(fig, ov, ev):
    """Apply one browser event {op: move|text|all|reset, id, fx, fy (new anchor, figure fraction), text} to `ov`."""
    ov = dict(ov or {})
    it = {k: (kind, ai, a) for k, kind, ai, a in items(fig)}
    op = ev.get("op")
    if op == "reset":
        return {}
    k = ev.get("id")
    if k not in it:
        return ov
    kind, ai, art = it[k]; ax = fig.axes[ai]
    if op == "text":
        o = dict(ov.get(k, {})); o["text"] = str(ev.get("text", ""))[:200]; ov[k] = o
        return ov
    if op == "move":
        x, y = _frac_to(ax, kind, fig, float(ev["fx"]), float(ev["fy"]))
        if kind == "legend":
            ov[k] = {"lx": x, "ly": y}; return ov
        x0, y0 = getattr(art, "_mfs_xy0", tuple(art.get_position()))
        o = dict(ov.get(k, {})); o["dx"], o["dy"] = x - x0, y - y0; ov[k] = o
        return ov
    if op == "all" and kind == "data":
        x0, y0 = getattr(art, "_mfs_xy0", tuple(art.get_position()))
        o = ov.get(k, {}); xs, ys = x0 + o.get("dx", 0), y0 + o.get("dy", 0)
        b0 = _trace_base(ax, art)
        for k2, (kd2, ai2, a2) in it.items():
            if ai2 != ai or kd2 != "data" or k2 == k:
                continue
            u0, v0 = getattr(a2, "_mfs_xy0", tuple(a2.get_position()))
            b2 = _trace_base(ax, a2)
            if b0 is not None and b2 is not None:      # same distance to its own trace offset
                tx, ty = xs, b2 + (ys - b0)
            else:                                       # no trace found: same displacement
                tx, ty = u0 + o.get("dx", 0), v0 + o.get("dy", 0)
            o2 = dict(ov.get(k2, {})); o2["dx"], o2["dy"] = tx - u0, ty - v0; ov[k2] = o2
    return ov


def payload(fig, dpi=110):
    """PNG of the whole figure (no tight crop, so figure fractions map 1:1) + the movable items' anchor boxes."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi)
    W, H = fig.get_size_inches() * dpi
    from .render import renderer_for
    r = renderer_for(fig)
    its = []
    for k, kind, ai, art in items(fig):
        try:
            bb = art.get_window_extent(r)
        except Exception:
            continue
        sc = dpi / fig.dpi
        if kind == "legend":
            ax_ = fig.axes[ai]
            anchor = (bb.x0 * sc, bb.y1 * sc)              # upper-left corner = legend anchor
            label = "Legend"
        else:
            anchor = tuple(np.asarray(fig.transFigure.inverted().transform(art.get_transform().transform(art.get_position()))) * (W, H))
            label = art.get_text()
        its.append({"id": k, "kind": kind, "label": label, "x0": bb.x0 * sc, "y0": H - bb.y1 * sc, "w": bb.width * sc, "h": bb.height * sc,
                    "ax": anchor[0], "ay": H - anchor[1]})
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode(), float(W), float(H), its


def component_dir():
    import sys
    cands = [os.path.join(os.path.dirname(os.path.abspath(__file__)), "components", "label_mover")]
    if getattr(sys, "_MEIPASS", None):
        cands.insert(0, os.path.join(sys._MEIPASS, "app", "mfs", "components", "label_mover"))
    for c in cands:
        if os.path.exists(os.path.join(c, "index.html")):
            return c
    return cands[-1]


_COMP = None


def editor(fig, key, component_key=None):
    """Draw the label editor; returns the latest browser event or None."""
    import streamlit.components.v1 as components
    global _COMP
    if _COMP is None:
        _COMP = components.declare_component("reticula_label_mover", path=component_dir())
    src, W, H, its = payload(fig)
    return _COMP(src=src, w=W, h=H, items=its, skey=key, default=None, key=component_key or ("lm_" + key))
