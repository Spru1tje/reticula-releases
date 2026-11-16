"""v2.11 unified sample list: one component for every multi-trace page (show/hide, colour, name, drag-to-reorder,
shift-click ranges, all/none/invert/reverse). The order drives the stacking and the legend. State lives in
st.session_state['_sl|<key>'] (tracked by the page-state module, so it is saved with the page and the save dialog).
The pure functions mirror the component's JS and are unit-tested."""
from __future__ import annotations
import os, sys

PREFIX = "_sl|"


def component_dir():
    cands = [os.path.join(os.path.dirname(os.path.abspath(__file__)), "components", "sample_list")]
    if getattr(sys, "_MEIPASS", None):
        cands.insert(0, os.path.join(sys._MEIPASS, "app", "mfs", "components", "sample_list"))
    for c in cands:
        if os.path.exists(os.path.join(c, "index.html")):
            return c
    return cands[-1]


# ---- pure logic (same rules as the JS)
def merge(state, ids, default_names, max_shown=None):
    """Keep order / names / visibility of files already in the list, append new files, drop removed ones.
    New entries are shown, except beyond `max_shown` (on a fresh list)."""
    state = [dict(o) for o in (state or []) if isinstance(o, dict) and "id" in o]
    known = {o["id"] for o in state}
    out = [o for o in state if o["id"] in set(ids)]
    fresh = not state
    for k, (i, nm) in enumerate(zip(ids, default_names)):
        if i not in known:
            out.append({"id": i, "name": nm, "show": not (fresh and max_shown is not None and k >= max_shown)})
    for o in out:
        o.setdefault("name", o["id"]); o["show"] = bool(o.get("show", True))
    return out


def set_range(items, a, b, show):
    items = [dict(o) for o in items]
    for k in range(min(a, b), max(a, b) + 1):
        items[k]["show"] = bool(show)
    return items


def move(items, a, b):
    items = [dict(o) for o in items]
    if 0 <= a < len(items) and 0 <= b < len(items) and a != b:
        o = items.pop(a); items.insert(b, o)
    return items


def all_(items):
    return [dict(o, show=True) for o in items]


def none_(items):
    return [dict(o, show=False) for o in items]


def invert(items):
    return [dict(o, show=not o["show"]) for o in items]


def visible(items):
    """[(id, unique name)] of the shown samples in list order (duplicate names get ' (2)', ...)."""
    seen, out = {}, []
    for o in items:
        if not o.get("show", True):
            continue
        n = str(o.get("name") or o["id"]); seen[n] = seen.get(n, 0) + 1
        out.append((o["id"], n if seen[n] == 1 else f"{n} ({seen[n]})"))
    return out


def colours(items, colour_for, colours_for):
    """Swatch colour per item: the colour the figure will use (visible ones in plotting order)."""
    vis = [o for o in items if o.get("show", True)]
    vc = dict(zip([id(o) for o in vis], colours_for([o["name"] for o in vis]))) if vis else {}
    return [vc.get(id(o)) or colour_for(o["name"], i) for i, o in enumerate(items)]


def plain(name):
    """Name without mathtext markup (ZrO$_2$ -> ZrO2), for captions."""
    import re
    return re.sub(r"[$_^{}\\]", "", str(name)).replace("mathrm", "")


def _hex(c):
    try:
        from matplotlib.colors import to_hex
        return to_hex(c)
    except Exception:
        return "#444444"


_COMP = None


def render(items, skey, component_key, applied_t=0):
    import streamlit.components.v1 as components
    global _COMP
    if _COMP is None:
        _COMP = components.declare_component("reticula_sample_list", path=component_dir())
    return _COMP(items=items, skey=skey, debounce=250, applied_t=applied_t or 0, default=None, key=component_key)


def sample_list(files_or_ids, key, default_names, label="🧪 Samples – show/hide, colour, rename, drag to reorder", expanded=True,
                container=None, max_shown=None):
    """Draw the list; returns the full item list (in order). `files_or_ids`: file objects (``.name``) or ids."""
    import streamlit as st
    from mfs import style as S
    ids = [getattr(f, "name", f) for f in files_or_ids]
    sk = PREFIX + key
    ss = st.session_state
    items = merge(ss.get(sk), ids, default_names, max_shown)
    ck = f"slc_{key}"
    v = ss.get(ck)                     # the component's latest value: apply it BEFORE drawing (args then match the browser)
    if isinstance(v, dict) and v.get("key") == key and v.get("t") and v.get("t") != ss.get(sk + "|t"):
        ss[sk + "|t"] = v["t"]
        old_cols = {o["id"]: _hex(c) for o, c in zip(items, colours(items, S.colour_for, S.colours_for))}
        new = merge([{"id": o["id"], "name": o["name"], "show": o["show"]} for o in v.get("items", []) if "id" in o],
                    ids, default_names)
        sc = dict(ss.get("series_colours") or {})
        changed_col = False
        for o in v.get("items", []):
            oc = old_cols.get(o.get("id"))
            if oc and o.get("colour") and str(o["colour"]).lower() != oc.lower():
                sc[o["name"]] = o["colour"]; changed_col = True
        if changed_col:
            ss["series_colours"] = sc; S.set_series_colours(sc)
        items = new
    box = (container or st).expander(label, expanded=expanded)
    with box:
        cols = colours(items, S.colour_for, S.colours_for)
        payload = [{"id": o["id"], "name": o["name"], "show": o["show"], "colour": _hex(c)} for o, c in zip(items, cols)]
        render(payload, key, component_key=ck, applied_t=ss.get(sk + "|t") or 0)
        ss[sk] = items
        vis = visible(items)
        st.caption(f"Plot order ({len(vis)} of {len(items)} shown; first = top of the stack and of the legend): "
                   + (" → ".join(plain(n) for _, n in vis) if vis else "nothing selected"))
    try:
        from mfs import pagestate as PS
        PS.track(sk)
    except Exception:
        pass
    return items
