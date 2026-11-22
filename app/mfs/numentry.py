"""Every slider of Reticula can also take a typed number (2.19.6).

slider(...) and select_slider(...) take the arguments of st.slider / st.select_slider and return the same value. A number box sits
next to the slider and shares its value: typing a number moves the slider, and dragging the slider updates the box. A range (a tuple)
gets two boxes. Typed values are clamped to the slider range and a range is kept in order, so the result always lies inside the range."""
from __future__ import annotations

import streamlit as st


def _k(key, part=""):
    return f"{key}__n{part}"


def _clamp(v, lo, hi):
    return min(max(float(v), float(lo)), float(hi))


def _pull(key, parts):
    """slider -> number boxes (on_change of the slider)."""
    v = st.session_state[key]
    if parts == 2:
        st.session_state[_k(key, "0")], st.session_state[_k(key, "1")] = v[0], v[1]
    else:
        st.session_state[_k(key, "0")] = v


def _push(key, parts, lo, hi):
    """number boxes -> slider (on_change of a number box): clamped to the range, a range kept in order."""
    if parts == 2:
        a = _clamp(st.session_state[_k(key, "0")], lo, hi)
        b = _clamp(st.session_state[_k(key, "1")], lo, hi)
        a, b = min(a, b), max(a, b)
        st.session_state[key] = (a, b)
        st.session_state[_k(key, "0")], st.session_state[_k(key, "1")] = a, b
    else:
        a = _clamp(st.session_state[_k(key, "0")], lo, hi)
        st.session_state[key] = a
        st.session_state[_k(key, "0")] = a


def _init(key, parts, value):
    """Sets the state of the widgets before they are created (no default is passed to them afterwards: no conflict in Streamlit)."""
    if key not in st.session_state:
        st.session_state[key] = (float(value[0]), float(value[1])) if parts == 2 else value
    if _k(key, "0") not in st.session_state:
        cur = st.session_state[key]
        if parts == 2:
            st.session_state[_k(key, "0")], st.session_state[_k(key, "1")] = cur[0], cur[1]
        else:
            st.session_state[_k(key, "0")] = cur


def _boxes(key, parts, lo, hi, step, fmt, label):
    """The number box (or two) next to a slider, in the column that is already open."""
    kw = dict(min_value=float(lo), max_value=float(hi), step=float(step) if step else None, format=fmt, label_visibility="collapsed")
    if parts == 2:
        c0, c1 = st.columns(2)
        with c0:
            st.number_input(f"{label} from", key=_k(key, "0"), on_change=_push, args=(key, 2, lo, hi), **kw)
        with c1:
            st.number_input(f"{label} to", key=_k(key, "1"), on_change=_push, args=(key, 2, lo, hi), **kw)
    else:
        st.number_input(label, key=_k(key, "0"), on_change=_push, args=(key, 1, lo, hi), **kw)


def _then(key, parts, cb, args, kwargs):
    """on_change of the slider: keeps the number boxes in step, then runs the caller's own on_change (if any)."""
    _pull(key, parts)
    if cb is not None:
        cb(*(args or ()), **(kwargs or {}))


def slider(label, min_value, max_value, value=None, step=None, format=None, key=None, help=None, **kw):
    """st.slider with a number box next to it (two boxes for a range). Same return value as st.slider."""
    key = key or f"nsl_{label}"
    cb, cb_args, cb_kw = kw.pop("on_change", None), kw.pop("args", None), kw.pop("kwargs", None)
    parts = 2 if isinstance(value, (tuple, list)) else 1
    lo, hi = min_value, max_value
    if value is None:
        value = (lo, hi) if parts == 2 else lo
    _init(key, parts, value)
    last = key + "__last"
    if st.session_state.get(last) != st.session_state[key]:       # the value was set by other code (e.g. a click in a chart): follow it
        _pull(key, parts)
    col_s, col_n = st.columns([3, 2], vertical_alignment="bottom")
    with col_s:
        out = st.slider(label, lo, hi, key=key, step=step, format=format, help=help, on_change=_then,
                        args=(key, parts, cb if callable(cb) else None, cb_args, cb_kw), **kw)
    with col_n:
        _boxes(key, parts, lo, hi, step, format, label)
    st.session_state[last] = out
    return out


def select_slider(label, options, value=None, format_func=None, key=None, help=None, **kw):
    """st.select_slider with a number box for the INDEX of the option (the slider keeps its labels, e.g. '293.0 K').
    Returns the option, as st.select_slider does."""
    opts = list(options)
    key = key or f"nss_{label}"
    cb, cb_args, cb_kw = kw.pop("on_change", None), kw.pop("args", None), kw.pop("kwargs", None)
    idx_default = opts.index(value) if value is not None else len(opts) // 2
    ik = key + "__idx"
    if ik not in st.session_state:
        st.session_state[ik] = idx_default
    if _k(ik, "0") not in st.session_state:
        st.session_state[_k(ik, "0")] = idx_default

    def _pull_idx():
        st.session_state[_k(ik, "0")] = opts.index(st.session_state[key])
        if callable(cb):
            cb(*(cb_args or ()), **(cb_kw or {}))

    def _push_idx():
        st.session_state[key] = opts[int(_clamp(st.session_state[_k(ik, "0")], 0, len(opts) - 1))]
        st.session_state[_k(ik, "0")] = opts.index(st.session_state[key])

    if key not in st.session_state:
        st.session_state[key] = opts[idx_default]
    col_s, col_n = st.columns([3, 2], vertical_alignment="bottom")
    with col_s:
        extra = {} if format_func is None else {"format_func": format_func}      # None would be called as the formatter
        out = st.select_slider(label, options=opts, key=key, help=help, on_change=_pull_idx, **extra, **kw)
    with col_n:
        st.number_input(f"{label} (index)", min_value=0, max_value=len(opts) - 1, step=1, key=_k(ik, "0"),
                        on_change=_push_idx, label_visibility="collapsed")
    return out
