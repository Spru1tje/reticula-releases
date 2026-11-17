"""Streamlit glue for the v2.13 ¹H NMR page: drag-to-assign component payload/events (pure functions, testable)."""
from __future__ import annotations

import base64
import io
import os

import numpy as np

_COMP = None


def component_dir():
    import sys
    cands = [os.path.join(os.path.dirname(os.path.abspath(__file__)), "components", "nmr_assign")]
    if getattr(sys, "_MEIPASS", None):
        cands.insert(0, os.path.join(sys._MEIPASS, "app", "mfs", "components", "nmr_assign"))
    for c in cands:
        if os.path.exists(os.path.join(c, "index.html")):
            return c
    return cands[-1]


def payload(sp, mults, mol, groups, assignment, fixed=None, width=760, mol_size=(320, 260), xlim=None):
    """Images (data URLs) + chip and drop-zone pixel positions for the assignment component."""
    import matplotlib
    import matplotlib.pyplot as plt
    from mfs import nmrfig as F
    fixed = fixed or {}
    png, coords = F.mol_drawing(mol, groups, size=mol_size, notes=False)
    mol_src = "data:image/png;base64," + base64.b64encode(png).decode()
    dpi = 100
    fig, ax = plt.subplots(figsize=(width / dpi, 3.0), dpi=dpi)
    hi, lo = F.draw(ax, sp, mults, None, groups, xlim=xlim, show_structure=False, show_integrals=True, mark_unknown=False)
    try:
        fig.set_layout_engine("none")                      # the app style may enable tight/constrained layout at save
    except Exception:
        pass
    fig.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.17)
    fig.canvas.draw()
    W, H = fig.canvas.get_width_height()
    buf = io.BytesIO()
    with matplotlib.rc_context({"savefig.bbox": "standard"}):          # keep pixel ↔ data mapping exact
        fig.savefig(buf, format="png", dpi=dpi)
    plt.close(fig)
    spec_src = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()
    zones = []
    tr = ax.transData
    for i, m in enumerate(mults):
        if m.note and not m.overlap and not m.label:
            continue
        if not (lo <= m.centre <= hi):
            continue
        a_, b_ = (max(m.lines) + 0.01, min(m.lines) - 0.01) if m.lines else (m.hi, m.lo)
        (x0, y0), (x1, y1) = tr.transform((a_, 1.0)), tr.transform((b_, -0.05))
        xa, xb = sorted((x0, x1))
        if xb - xa < 10:
            c = (xa + xb) / 2; xa, xb = c - 5, c + 5
        zones.append(dict(mi=i, x0=float(xa), x1=float(xb), y0=float(H - max(y0, y1)), y1=float(H - min(y0, y1)),
                          short=f"{m.centre:.2f}", title=f"δ {m.centre:.2f} {m.mult} {m.nH if m.nH is not None else ''}H"
                          + (f" · {m.label}" if m.label else "")))
    chips = []
    allxy = np.array([coords[a] for a in coords]) if len(coords) else np.zeros((1, 2))
    cen = allxy.mean(axis=0)
    mw_, mh_ = mol_size

    def _out(a, dist):
        """Position `dist` px away from atom a, pointing away from its bonded heavy atoms (so the chip sits beside the atom, not on it)."""
        pos = np.asarray(coords[a], float)
        nb = [np.asarray(coords[n_.GetIdx()], float) for n_ in mol.GetAtomWithIdx(a).GetNeighbors() if n_.GetIdx() in coords]
        v = pos - (np.mean(nb, axis=0) if nb else cen)
        if np.linalg.norm(v) < 1.0:
            v = pos - cen
        if np.linalg.norm(v) < 1.0:
            v = np.array([1.0, -1.0])
        return pos + v / np.linalg.norm(v) * dist
    for g in groups:
        anchor = np.mean([coords[a] for a in g["atoms"]], axis=0)
        hetero = any(mol.GetAtomWithIdx(a).GetSymbol() != "C" for a in g["atoms"])
        xy = np.mean([_out(a, 30.0 if hetero else 24.0) for a in g["atoms"]], axis=0)
        sh = "–" if g["shift"] is None else f"{g['shift']:.2f}"
        on = ", ".join(f"{mol.GetAtomWithIdx(a).GetSymbol()}{a + 1}" for a in g["atoms"])
        chips.append(dict(gid=g["gid"], x=float(xy[0]), y=float(xy[1]), ax=float(anchor[0]), ay=float(anchor[1]),
                          text="H" + g.get("tag", ",".join(map(str, g["labels"]))), fixed=g["gid"] in fixed,
                          title=f"H{g.get('tag', '')}: {g['nH']} H on {on} · predicted {sh} ppm ({g['conf']})"))
    for _ in range(40):                                  # keep the chips apart from each other and inside the picture
        moved = False
        for i in range(len(chips)):
            for j in range(i + 1, len(chips)):
                dx, dy = chips[j]["x"] - chips[i]["x"], chips[j]["y"] - chips[i]["y"]
                d = (dx * dx + dy * dy) ** 0.5
                if d < 30.0:
                    ux, uy = (dx / d, dy / d) if d > 1e-6 else (1.0, 0.0)
                    push = (30.0 - d) / 2
                    chips[i]["x"] -= ux * push; chips[i]["y"] -= uy * push
                    chips[j]["x"] += ux * push; chips[j]["y"] += uy * push
                    moved = True
        for c in chips:
            c["x"] = float(min(max(c["x"], 16.0), mw_ - 16.0)); c["y"] = float(min(max(c["y"], 12.0), mh_ - 12.0))
        if not moved:
            break
    return dict(mol_src=mol_src, spec_src=spec_src, mw=mol_size[0], mh=mol_size[1], sw=W, zones=zones, groups=chips,
                assign={str(k): v for k, v in assignment.items()})


def handle(fixed, ev):
    """Apply a component event to the manual-corrections dict {gid: mult index | None}; returns the new dict."""
    fixed = dict(fixed or {})
    if not ev or not isinstance(ev, dict):
        return fixed
    op = ev.get("op")
    if op == "reset":
        return {}
    if op == "assign" and ev.get("gid") is not None:
        fixed[int(ev["gid"])] = int(ev["mi"])
    elif op == "unassign" and ev.get("gid") is not None:
        fixed[int(ev["gid"])] = None
    return fixed


def widget(sp, mults, mol, groups, assignment, fixed, key, xlim=None):
    import streamlit.components.v1 as components
    global _COMP
    if _COMP is None:
        _COMP = components.declare_component("reticula_nmr_assign", path=component_dir())
    p = payload(sp, mults, mol, groups, assignment, fixed, xlim=xlim)
    return _COMP(skey=key, default=None, key="nmra_" + key, **p)
