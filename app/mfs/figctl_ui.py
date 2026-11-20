"""Streamlit control panel for mfs.figctl - the same panel under every figure of every page."""
from __future__ import annotations
import hashlib, json, re
import numpy as np
import streamlit as st
from . import figctl as F
from . import style as S

NOT_SERIES = re.compile(r"^(adsorption|desorption|quantitative|warning|semi-quantitative|consensus.*|fit|model|"
                        r"simulated.*|blank.*|ideal.*|_.*)$", re.I)


# ----------------------------------------------------------------------------- house style (sidebar)
def house():
    return st.session_state.get("house_style") or {}


def _clear_fc_keys():
    for k in list(st.session_state.keys()):
        if isinstance(k, str) and k.startswith("fc_"):
            del st.session_state[k]


def _load_house_cb():
    up = st.session_state.get("house_upload")
    if up is None:
        return
    try:
        d = json.loads(up.getvalue().decode("utf-8"))
        st.session_state["house_style"] = F.house_style(F.merge(F.default_settings(), d))
        st.session_state["house_msg"] = f"House style '{up.name}' loaded and applied to every figure."
        _clear_fc_keys()
    except Exception as e:
        st.session_state["house_msg"] = f"Could not read {up.name}: {e}"


def sidebar_house_style():
    with st.sidebar.expander("🎨 House style preset (JSON)"):
        st.caption("Tick, legend, font and frame settings shared by ALL figures. Save the panel of any figure as a house "
                   "style, or load a JSON saved earlier.")
        st.file_uploader("Load house style", type=["json"], key="house_upload", on_change=_load_house_cb)
        if st.session_state.get("house_msg"):
            st.info(st.session_state.pop("house_msg"))
        if house():
            st.download_button("Download current house style", F.to_json(house()).encode(), "house_style.json",
                               mime="application/json", on_click="ignore", key="house_dl")
            if st.button("Reset house style", key="house_reset"):
                st.session_state["house_style"] = {}
                _clear_fc_keys(); st.rerun()


def _base():
    return F.merge(F.default_settings(), house())


# ----------------------------------------------------------------------------- state <-> cfg
def _k(key, *parts):
    return "fc_" + key + "_" + "_".join(str(p) for p in parts)


def _get(key, name, default):
    return st.session_state.get(_k(key, name), default)


def cfg_to_state(key, cfg):
    """Write a settings dict into the widget keys of figure `key` (used when a preset is loaded)."""
    cfg = F.merge(F.default_settings(), cfg)
    ss = st.session_state
    ss[_k(key, "title")] = cfg.get("title", "")
    ss[_k(key, "title_size")] = float(cfg.get("title_size") or 0.0)
    ss[_k(key, "mathify")] = bool(cfg.get("mathify", True))
    T = cfg["ticks"]
    ss[_k(key, "preset")] = T.get("preset", "All 4 sides")
    for s in F.SIDES:
        ss[_k(key, "side", s)] = bool(T["sides"].get(s, True))
        ss[_k(key, "lab", s)] = bool(T["labels"].get(s, s in ("bottom", "left")))
    for n in ("direction", "minor", "major_len", "minor_len", "width", "frame"):
        ss[_k(key, n)] = T.get(n)
    ss[_k(key, "label_size")] = float(cfg["font"].get("label_size") or 0.0)
    ss[_k(key, "tick_size")] = float(cfg["font"].get("tick_size") or 0.0)
    ss[_k(key, "family")] = cfg["font"].get("family") or ""
    L = cfg["legend"]
    ss[_k(key, "leg", "title_size")] = float(L.get("title_size") or 0.0)
    ss[_k(key, "leg", "family")] = L.get("family") or ""
    for n in ("custom", "show", "loc", "ncol", "fontsize", "frame", "markerscale", "title", "panel"):
        ss[_k(key, "leg", n)] = L.get(n)
    if L.get("entries"):
        ss[_k(key, "leg_entries_json")] = L["entries"]
    for i, p in (cfg.get("panels") or {}).items():
        for n, v in p.items():
            if n.endswith("major") or n.endswith("minor"):
                ss[_k(key, "p", i, n, "n")] = float(v or 0.0); ss[_k(key, "p", i, n, "s")] = float(v or 0.0)
            else:
                ss[_k(key, "p", i, n)] = v
    ss[_k(key, "wavelength")] = float(cfg["units_params"].get("wavelength", 1.5406))
    ss[_k(key, "m0")] = float(cfg["units_params"].get("m0_mg", 10.0))


def _load_fig_cb(key):
    up = st.session_state.get(_k(key, "upload"))
    if up is None:
        return
    try:
        cfg = json.loads(up.getvalue().decode("utf-8"))
        cfg_to_state(key, cfg)
        st.session_state[_k(key, "msg")] = f"Settings from {up.name} applied."
    except Exception as e:
        st.session_state[_k(key, "msg")] = f"Could not read {up.name}: {e}"


def _sync(src, dst, lo=None, hi=None):
    v = st.session_state.get(src, 0.0)
    if lo is not None:
        v = float(min(max(v, lo), hi))
    st.session_state[dst] = v


def _inc(col, key, i, name, label, span, help_=None):
    """Slider + number input kept in sync (0 = automatic)."""
    span = float(abs(span)) if span and np.isfinite(span) and span > 0 else 1.0
    step = F.nice_step(span, 5) / 20
    hi = float(F.nice_step(span, 1))
    ks, kn = _k(key, "p", i, name, "s"), _k(key, "p", i, name, "n")
    st.session_state.setdefault(kn, 0.0)
    st.session_state.setdefault(ks, 0.0)
    st.session_state[ks] = float(min(max(st.session_state[ks], 0.0), hi))
    col.slider(label, 0.0, hi, step=step, key=ks, on_change=_sync, args=(ks, kn), help=help_ or "0 = automatic")
    col.number_input(label + " (exact)", 0.0, None, step=step, key=kn, on_change=_sync, args=(kn, ks, 0.0, hi),
                     format="%g", label_visibility="collapsed")
    return float(st.session_state.get(kn) or 0.0)


# ----------------------------------------------------------------------------- the panel
def panel(fig, key, units=None, series=None, expanded=False):
    """Draw the customisation panel for `fig` and return (settings, families).
    units: {'x': family, 'y': family} for all panels or {index: {...}}."""
    key = re.sub(r"[^A-Za-z0-9_]+", "_", key)
    base = _base()
    desc = F.describe(fig)
    data_panels = [d for d in desc if d["axison"]]
    fams = {}
    for d in desc:
        u = (units or {}).get(d["index"]) if units and d["index"] in units else (units if units and ("x" in units or "y" in units) else {})
        fams[d["index"]] = u or {}
    entries = F.legend_entries(fig)
    direct = F.direct_labels(fig)
    cfg = F.merge(F.default_settings(), base)
    with st.expander("🎛️ Customise this figure – title, axes, units, ticks, legend, colours, presets", expanded=expanded):
        if st.session_state.get(_k(key, "msg")):
            st.info(st.session_state.pop(_k(key, "msg")))
        from mfs import levels as LV
        adv = LV.at_least("Advanced")                  # v2.15: Standard shows the essentials; ticks and presets are Advanced
        _tabs = st.tabs(["🖋️ Title & axes", "🏷️ Legend", "🎨 Style & colours"] + (["📏 Axes & ticks", "💾 Preset"] if adv else []))
        t1, t3, t4 = _tabs[:3]
        t2, t5 = (_tabs[3], _tabs[4]) if adv else (None, None)
        # ---------------- title & axes
        with t1:
            c1, c2, c3 = st.columns([3, 1, 1])
            cfg["title"] = c1.text_input("Figure title (optional)", base.get("title", ""), key=_k(key, "title"))
            if adv:                      # Standard: title size, sub/superscript rule and font sizes stay at their defaults
                cfg["title_size"] = c2.number_input("Title size (pt, 0 = auto)", 0.0, 24.0, float(base.get("title_size") or 0.0), 0.5, key=_k(key, "title_size"))
                cfg["mathify"] = c3.checkbox("^ and _ → super/subscript", base.get("mathify", True), key=_k(key, "mathify"),
                                             help="cm^-1 → cm⁻¹, ZrO_2 → ZrO₂, x^{2+}. Text with $…$ is used as mathtext as is. Unicode (², ₂, ⁻¹, θ, µ) always works.")
            if adv:
                c1, c2, c3 = st.columns(3)
                cfg["font"]["label_size"] = c1.number_input("Axis-title size (pt, 0 = style)", 0.0, 24.0, float(base["font"].get("label_size") or 0.0), 0.5, key=_k(key, "label_size"))
                cfg["font"]["tick_size"] = c2.number_input("Tick-label size (pt, 0 = style)", 0.0, 24.0, float(base["font"].get("tick_size") or 0.0), 0.5, key=_k(key, "tick_size"))
                _fam_opts = [""] + [f for f in S.available_fonts() if f != "auto"]
                cfg["font"]["family"] = c3.selectbox("Font (this figure)", _fam_opts, index=_fam_opts.index(base["font"].get("family") or "") if (base["font"].get("family") or "") in _fam_opts else 0,
                                                     format_func=lambda f: f or "(figure style)", key=_k(key, "family"),
                                                     help="Overrides the font of the figure style for this figure only (axis titles, ticks, texts and legend).")
            if not data_panels:
                st.caption("Image figure – no axes to customise (use the scale-bar controls above).")
            pidx = 0
            if len(data_panels) > 1:
                pidx = st.selectbox("Panel", [d["index"] for d in data_panels],
                                    format_func=lambda i: f"panel {i + 1}: {desc[i]['ylabel'] or desc[i]['xlabel'] or ''}"[:60], key=_k(key, "panel_sel"))
            for d in data_panels:
                i = d["index"]
                if i != pidx:
                    # keep values of panels not shown
                    p = {n[len(_k(key, "p", i, "")):]: v for n, v in st.session_state.items()
                         if isinstance(n, str) and n.startswith(_k(key, "p", i, "")) and not n.endswith("_s")}
                    cfg["panels"][str(i)] = {n[:-2] if n.endswith("_n") else n: v for n, v in p.items()}
                    continue
                p = {}
                fam = fams.get(i, {})
                for ax_ in ("x", "y"):
                    c1, c2, c3, c4 = st.columns([3, 2, 1, 1])
                    cur_lab = d[f"{ax_}label"]
                    p[f"{ax_}label"] = c1.text_input(f"{ax_.upper()}-axis title", "", key=_k(key, "p", i, f"{ax_}label"),
                                                     placeholder=cur_lab or "(none)", help="Leave empty to keep the default title")
                    if fam.get(ax_):
                        names = F.unit_names(fam[ax_])
                        p[f"{ax_}_unit"] = c2.selectbox(f"{ax_.upper()} unit", names, key=_k(key, "p", i, f"{ax_}_unit"))
                        p[f"{ax_}_family"] = fam[ax_]
                    else:
                        c2.caption(f"{ax_.upper()} unit: as in the data")
                    lim = d[f"{ax_}lim"]
                    p[f"{ax_}min"] = c3.number_input(f"{ax_} min", value=None, key=_k(key, "p", i, f"{ax_}min"), placeholder=f"{min(lim):.4g}", format="%g")
                    p[f"{ax_}max"] = c4.number_input(f"{ax_} max", value=None, key=_k(key, "p", i, f"{ax_}max"), placeholder=f"{max(lim):.4g}", format="%g")
                if fam.get("x") == "2theta" and p.get("x_unit") not in (None, "2θ (°)"):
                    cfg["units_params"]["wavelength"] = st.number_input("Wavelength λ (Å) for d / q", 0.1, 3.0, 1.5406, format="%.4f", key=_k(key, "wavelength"))
                if fam.get("y") == "mass_pct" and p.get("y_unit", "").startswith("mg"):
                    cfg["units_params"]["m0_mg"] = st.number_input("Initial sample mass (mg)", 0.01, 1000.0, 10.0, key=_k(key, "m0"))
                c1, c2 = st.columns(2)
                p["xlog"] = c1.checkbox("log x", d["xscale"] == "log", key=_k(key, "p", i, "xlog"))
                p["ylog"] = c2.checkbox("log y", d["yscale"] == "log", key=_k(key, "p", i, "ylog"))
                st.caption("Ranges are entered in the units shown on the axis. Empty = automatic.")
                cfg["panels"][str(i)] = {**cfg["panels"].get(str(i), {}), **p}
        # ---------------- ticks
        if adv:
            with t2:
                T = cfg["ticks"]
                c1, c2, c3 = st.columns(3)
                presets = list(F.TICK_PRESETS)
                T["preset"] = c1.selectbox("Ticks on", presets, index=presets.index(base["ticks"].get("preset", "All 4 sides")), key=_k(key, "preset"))
                T["direction"] = c2.radio("Direction", ["in", "out", "inout"], index=["in", "out", "inout"].index(base["ticks"].get("direction", "in")),
                                          horizontal=True, key=_k(key, "direction"))
                T["frame"] = c3.radio("Frame", ["box", "match ticks"], index=["box", "match ticks"].index(base["ticks"].get("frame", "box")),
                                      horizontal=True, key=_k(key, "frame"), help="'match ticks' hides the axis lines on sides without ticks")
                sides = F.TICK_PRESETS.get(T["preset"]) or base["ticks"]["sides"]
                st.markdown("**Ticks / tick labels per side**")
                cols = st.columns(4)
                for c, s in zip(cols, F.SIDES):
                    T["sides"][s] = c.checkbox(f"ticks {s}", sides[s], key=_k(key, "side", s), disabled=T["preset"] != "Custom")
                    if T["preset"] != "Custom":
                        T["sides"][s] = sides[s]
                    T["labels"][s] = c.checkbox(f"labels {s}", base["ticks"]["labels"].get(s, s in ("bottom", "left")), key=_k(key, "lab", s))
                c1, c2, c3, c4 = st.columns(4)
                T["minor"] = c1.checkbox("Minor ticks", base["ticks"].get("minor", True), key=_k(key, "minor"))
                T["major_len"] = c2.number_input("Major length (pt)", 0.0, 12.0, float(base["ticks"].get("major_len", 3.5)), 0.5, key=_k(key, "major_len"))
                T["minor_len"] = c3.number_input("Minor length (pt)", 0.0, 8.0, float(base["ticks"].get("minor_len", 2.0)), 0.5, key=_k(key, "minor_len"))
                T["width"] = c4.number_input("Tick width (pt)", 0.1, 3.0, float(base["ticks"].get("width", 0.6)), 0.1, key=_k(key, "width"))
                for d in data_panels:
                    i = d["index"]
                    if len(data_panels) > 1:
                        st.markdown(f"**Increments – panel {i + 1}**")
                    else:
                        st.markdown("**Increments** (in the displayed unit; 0 = automatic)")
                    p = cfg["panels"].setdefault(str(i), {})
                    c1, c2, c3, c4 = st.columns(4)
                    spans = {}
                    for ax_ in ("x", "y"):
                        lo, hi = d[f"{ax_}lim"]
                        fam = fams.get(i, {}).get(ax_); un = p.get(f"{ax_}_unit")
                        u = F._unit(fam, un) if fam and un else None
                        if u and u[3]:
                            with np.errstate(all="ignore"):
                                a, b = u[3](lo, cfg["units_params"]), u[3](hi, cfg["units_params"])
                            spans[ax_] = abs(float(b) - float(a)) if np.isfinite(a) and np.isfinite(b) else abs(hi - lo)
                        else:
                            spans[ax_] = abs(hi - lo)
                    p["xmajor"] = _inc(c1, key, i, "xmajor", "x major", spans["x"])
                    p["xminor"] = _inc(c2, key, i, "xminor", "x minor", spans["x"] / 4)
                    p["ymajor"] = _inc(c3, key, i, "ymajor", "y major", spans["y"])
                    p["yminor"] = _inc(c4, key, i, "yminor", "y minor", spans["y"] / 4)
                    if d["xscale"] == "log" or d["yscale"] == "log":
                        st.caption("Log axes use decades; increments apply to linear axes only.")
                    if d["yticks_hidden"]:
                        st.caption("Stacked/offset plot: the y axis is in arbitrary units and has no ticks unless you set a y increment.")
        # ---------------- legend
        with t3:
            L = cfg["legend"]
            if not entries and not direct:
                st.caption("This figure has no legend entries.")
            L["custom"] = st.checkbox("Customise the legend (otherwise the default layout is kept)", bool(base["legend"].get("custom", False)),
                                      key=_k(key, "leg", "custom"))
            c1, c2, c3, c4 = st.columns(4)
            L["show"] = c1.checkbox("Show legend", base["legend"].get("show", True), key=_k(key, "leg", "show"))
            locs = F.LEGEND_LOCS
            L["loc"] = c2.selectbox("Position", locs, index=locs.index(base["legend"].get("loc", "best")) if base["legend"].get("loc", "best") in locs else 0,
                                    key=_k(key, "leg", "loc"))
            L["ncol"] = c3.number_input("Columns", 1, 8, int(base["legend"].get("ncol", 1) or 1), key=_k(key, "leg", "ncol"))
            L["fontsize"] = c4.number_input("Font size (pt, 0 = auto)", 0.0, 20.0, float(base["legend"].get("fontsize") or 0.0), 0.5, key=_k(key, "leg", "fontsize"))
            c1, c2, c3, c4 = st.columns(4)
            L["frame"] = c1.checkbox("Frame", bool(base["legend"].get("frame", False)), key=_k(key, "leg", "frame"))
            L["markerscale"] = c2.number_input("Marker size ×", 0.2, 5.0, float(base["legend"].get("markerscale", 1.0) or 1.0), 0.1, key=_k(key, "leg", "markerscale"))
            L["title"] = c3.text_input("Legend title", "", key=_k(key, "leg", "title"))
            c1, c2, c3, c4 = st.columns(4)
            L["title_size"] = c1.number_input("Legend-title size (pt, 0 = same as entries)", 0.0, 24.0, float(base["legend"].get("title_size") or 0.0), 0.5, key=_k(key, "leg", "title_size"))
            _lfam = [""] + [f for f in S.available_fonts() if f != "auto"]
            L["family"] = c2.selectbox("Legend font", _lfam, index=_lfam.index(base["legend"].get("family") or "") if (base["legend"].get("family") or "") in _lfam else 0,
                                       format_func=lambda f: f or "(same as figure)", key=_k(key, "leg", "family"))
            if len(data_panels) > 1:
                L["panel"] = c4.selectbox("Put legend on panel", list(range(len(data_panels))), format_func=lambda i: f"panel {i + 1}", key=_k(key, "leg", "panel"))
            labels = [l for l, _ in entries] + [t for t in direct if t not in dict(entries)]
            if labels:
                st.caption("Rename, hide and reorder entries (also renames labels written next to traces). Sort by 'order'.")
                import pandas as pd
                hk = hashlib.md5("|".join(labels).encode()).hexdigest()[:8]
                pre = {e["label"]: e for e in st.session_state.get(_k(key, "leg_entries_json"), [])}
                df = pd.DataFrame({"entry": labels, "new text": [pre.get(l, {}).get("text", l) for l in labels],
                                   "show": [pre.get(l, {}).get("show", True) for l in labels],
                                   "order": [pre.get(l, {}).get("order", i) for i, l in enumerate(labels)]})
                ed = st.data_editor(df, key=_k(key, "leg_ed", hk), hide_index=True, disabled=["entry"],
                                    column_config={"order": st.column_config.NumberColumn(step=1)})
                ed = ed.sort_values("order")
                L["entries"] = [{"label": r["entry"], "text": r["new text"] or r["entry"], "show": bool(r["show"]), "order": int(r["order"])}
                                for _, r in ed.iterrows()]
                if any(e["text"] != e["label"] or not e["show"] for e in L["entries"]) and not L["custom"]:
                    # renaming legend entries needs a rebuilt legend; direct labels are renamed regardless
                    if entries:
                        L["custom"] = True
        # ---------------- colours
        with t4:
            cfg["simple_black"] = st.checkbox("⚫ Simple mode – all black, samples told apart by symbol or line style (sleek journal look)", base.get("simple_black", False),
                                              key=_k(key, "simple_black"), help="Everything becomes black. Data points get one symbol per sample (● ■ ▲ ◆ ▼ …; the second group of a scatter is open). "
                                                                                 "Continuous curves (spectra, patterns, fits) get one line style per sample (solid, dashed, dotted, dash-dot) with no markers. "
                                                                                 "Bars get a hatch per series, filled areas a light grey. Switch it off to return to the colours.")
            names = [n for n in series if n] if series is not None else None
            if names is None:
                names = [l for l, _ in entries if not NOT_SERIES.match(l)] + [t for t in direct if not NOT_SERIES.match(t)]
            names = list(dict.fromkeys(names))
            if not names:
                st.caption("No sample series in this figure.")
            else:
                st.caption("Colour per sample – used for this sample in every figure (click ↺ to go back to the fixed palette).")
                cols = st.columns(min(4, len(names)))
                cur = st.session_state.setdefault("series_colours", {})
                defaults = S.colours_for(names)
                for j, (nm, dc) in enumerate(zip(names, defaults)):
                    c = cols[j % len(cols)]
                    ck = "sc_" + key + "_" + hashlib.md5(nm.encode()).hexdigest()[:10]
                    st.session_state[ck] = cur.get(nm, dc)   # one colour per sample, shared by every figure
                    c.color_picker(re.sub(r"[$_{}\\]|mathrm", "", nm)[:28], key=ck, on_change=_colour_cb, args=(nm, ck))
                if st.button("↺ reset colours of these samples", key=_k(key, "colreset")):
                    for nm in names:
                        cur.pop(nm, None)
                    S.set_series_colours(cur); st.rerun()
        # ---------------- preset
        if adv:
            with t5:
                c1, c2, c3 = st.columns(3)
                c1.download_button("Save these settings (JSON)", F.to_json(cfg).encode(), f"{key}_figure_settings.json",
                                   mime="application/json", on_click="ignore", key=_k(key, "dl"))
                c2.download_button("Save as house style (JSON)", F.to_json(F.house_style(cfg)).encode(), "house_style.json",
                                   mime="application/json", on_click="ignore", key=_k(key, "dlh"))
                if c3.button("Use as house style for all figures", key=_k(key, "mkhouse")):
                    st.session_state["house_style"] = F.house_style(cfg)
                    st.session_state["house_msg"] = "House style set from this figure."
                    st.rerun()
                st.file_uploader("Load settings for this figure (JSON)", type=["json"], key=_k(key, "upload"),
                                 on_change=_load_fig_cb, args=(key,))
    return cfg, fams


def _colour_cb(name, ck):
    cur = st.session_state.setdefault("series_colours", {})
    if ck not in st.session_state:       # the picker was not drawn in this run (the sample list changed): nothing to store
        return
    cur[name] = st.session_state[ck]
    S.set_series_colours(cur)
