"""Page-state persistence + unsaved-work guard (v2.9).

Streamlit forgets a widget's value as soon as the page that drew it stops rendering. Reticula draws one page per run
(sidebar radio), so every page switch used to wipe uploads, settings and results.  This module:

* records, at the end of every run, the values of the keyed widgets the current page drew (only value-type widgets:
  checkbox/toggle, radio, selectbox, multiselect, slider, number/text/date/colour inputs …; never buttons, uploaders,
  data editors or custom components) in a plain session dict ``_ps_wp`` and writes them back into ``st.session_state``
  at the start of the next run, so the widget finds its old value when the page is shown again;
* wraps ``file_uploader``: uploaded bytes are kept per widget key and handed back (as file-like objects with ``.name``)
  when the page is re-mounted with an empty uploader; removing the file with ✖ on the same page forgets it;
* tracks per page which values differ from the page's baseline (first-seen defaults or the last save) → "unsaved";
* intercepts navigation (radio ``on_change`` / :func:`request`) and asks Save / Don't save / Cancel in an ``st.dialog``.
  Save = new baseline + snapshot to ``<data dir>/workspace/<page>.json`` (+ uploaded files), restored on next launch;
  Don't save = revert the page to its last saved state (or to defaults if never saved); Cancel = stay.
"""
from __future__ import annotations

import hashlib
import io
import json
import os
import re
import datetime as _dt

import streamlit as st

_ALLOWED = {"bool_value", "int_value", "double_value", "string_value", "int_array_value", "double_array_value",
            "string_array_value"}
_SKIP_PREFIX = ("_", "FormSubmitter", "pk_", "psx_")
GLOBAL_KEYS = {"nav"}
WS_DIR = None            # set by configure()


def configure(ws_dir):
    global WS_DIR
    WS_DIR = ws_dir


def enabled():
    v = st.session_state.get("_ps_guard")
    return (os.environ.get("MFS_NAV_GUARD", "1") != "0") if v is None else bool(v)


# ----------------------------------------------------------------------------- helpers
def _simple(v):
    if isinstance(v, (bool, int, float, str)) or v is None:
        return True
    if isinstance(v, (_dt.date, _dt.time)):
        return True
    if isinstance(v, (list, tuple)):
        return all(isinstance(x, (bool, int, float, str, _dt.date, _dt.time)) for x in v)
    return False


def _plain_data(v, depth=0):
    """Scalars and (nested) lists / tuples / str-keyed dicts of them, e.g. the sample-list state (v2.11)."""
    if depth > 6:
        return False
    if _simple(v):
        return True
    if isinstance(v, (list, tuple)):
        return all(_plain_data(x, depth + 1) for x in v)
    if isinstance(v, dict):
        return all(isinstance(k, str) and _plain_data(x, depth + 1) for k, x in v.items())
    return False


def _jsonable(v):
    if isinstance(v, _dt.datetime):
        return {"__dt__": v.isoformat()}
    if isinstance(v, _dt.date):
        return {"__d__": v.isoformat()}
    if isinstance(v, _dt.time):
        return {"__t__": v.isoformat()}
    if isinstance(v, tuple):
        return {"__tuple__": [_jsonable(x) for x in v]}
    if isinstance(v, list):
        return [_jsonable(x) for x in v]
    return v


def _unjson(v):
    if isinstance(v, dict):
        if "__dt__" in v:
            return _dt.datetime.fromisoformat(v["__dt__"])
        if "__d__" in v:
            return _dt.date.fromisoformat(v["__d__"])
        if "__t__" in v:
            return _dt.time.fromisoformat(v["__t__"])
        if "__tuple__" in v:
            return tuple(_unjson(x) for x in v["__tuple__"])
    if isinstance(v, list):
        return [_unjson(x) for x in v]
    return v


def _key_live(k):
    """v2.11: is a widget with user key `k` already drawn in this run? (O(1); the full scan cost ~1.5 ms per widget)"""
    try:
        from streamlit.runtime.scriptrunner_utils.script_run_context import get_script_run_ctx
        ctx = get_script_run_ctx()
        wid = ctx.session_state._state._key_id_mapper._key_id_mapping.get(k)
        return wid is not None and wid in ctx.shared.widget_ids_this_run
    except Exception:
        return k in _rendered_widget_keys()


def _rendered_widget_keys():
    """{user_key: value_type} of keyed widgets drawn so far in this run (internal API, guarded)."""
    try:
        from streamlit.runtime.scriptrunner_utils.script_run_context import get_script_run_ctx
        ctx = get_script_run_ctx()
        state = ctx.session_state._state
        live = ctx.shared.widget_ids_this_run          # ThreadSafeSet: membership only
        idk = state._key_id_mapper._id_key_mapping
        meta = state._new_widget_state.widget_metadata
        out = {}
        for wid in list(idk):
            if wid not in live:
                continue
            k = idk.get(wid)
            m = meta.get(wid)
            if k and m is not None:
                out[k] = m.value_type
        return out
    except Exception:
        return {}


def _ss():
    ss = st.session_state
    for k, d in (("_ps_wp", {}), ("_ps_pages", {}), ("_ps_base", {}), ("_ps_up", {}), ("_ps_saved", {}),
                 ("_ps_extra", {}), ("_ps_live_cur", [])):
        if k not in ss:
            ss[k] = d
    return ss


def slug(page):
    return re.sub(r"[^A-Za-z0-9]+", "_", page).strip("_").lower() or "page"


# ----------------------------------------------------------------------------- uploads
class RestoredFile(io.BytesIO):
    """Stand-in for a Streamlit UploadedFile rebuilt from bytes kept in the session (or a saved snapshot)."""

    def __init__(self, name, type_, data):
        super().__init__(data)
        self.name, self.type, self.size = name, type_, len(data)
        self.file_id = "restored-" + hashlib.md5(name.encode() + data[:4096] + str(len(data)).encode()).hexdigest()[:12]
        self.restored = True

    def __eq__(self, o):
        return isinstance(o, RestoredFile) and o.file_id == self.file_id

    def __hash__(self):
        return hash(self.file_id)


def _fp(recs):
    return [(n, len(b), hashlib.md5(b).hexdigest()[:10]) for n, _, b in (recs or [])]


_ORIG_FU = None


def install_uploader_hook():
    """Patch DeltaGenerator.file_uploader once per process."""
    global _ORIG_FU
    from streamlit.delta_generator import DeltaGenerator
    if getattr(DeltaGenerator.file_uploader, "_ps_wrapped", False):
        return
    _ORIG_FU = DeltaGenerator.file_uploader

    def file_uploader(self, label, type=None, accept_multiple_files=False, key=None, *a, **kw):
        res = _ORIG_FU(self, label, type, accept_multiple_files, key, *a, **kw)
        if key is None:
            return res
        try:
            return _after_upload(self, key, res, accept_multiple_files)
        except Exception:
            return res

    file_uploader._ps_wrapped = True
    file_uploader.__doc__ = _ORIG_FU.__doc__
    DeltaGenerator.file_uploader = file_uploader


def _after_upload(dg, key, res, multiple):
    if str(key).startswith("_"):              # private uploaders (project open, …) are not page work
        return res
    ss = _ss()
    files = (list(res) if multiple else [res]) if res else []
    files = [f for f in files if f is not None]
    page = ss.get("_ps_page")
    if page:
        ss["_ps_pages"].setdefault(page, {}).setdefault("up", [])
        if key not in ss["_ps_pages"][page]["up"]:
            ss["_ps_pages"][page]["up"].append(key)
    live_prev = set(ss.get("_ps_live_prev", []))
    if files:
        ss["_ps_up"][key] = [(f.name, getattr(f, "type", "") or "", f.getvalue()) for f in files]
        ss["_ps_live_cur"].append(key)
        return res
    recs = ss["_ps_up"].get(key)
    if not recs:
        return res
    if key in live_prev:                     # the user removed the file(s) on this very page
        ss["_ps_up"].pop(key, None)
        return res
    objs = [RestoredFile(n, t, b) for n, t, b in recs]
    c1, c2 = dg.columns([5, 1])
    c1.caption("↩ kept from earlier: " + ", ".join(f"`{o.name}`" for o in objs) + " – drop new files above to replace")
    if c2.button("✖", key=f"_psforget_{key}", help="Forget the kept file(s)"):
        ss["_ps_up"].pop(key, None)
        st.rerun()
    return objs if multiple else objs[0]


# ----------------------------------------------------------------------------- per-run hooks
def begin():
    """Call at the very top of the script, before any widget."""
    ss = _ss()
    if "_ps_loaded" not in ss:
        ss["_ps_loaded"] = True
        _load_snapshots()
    to = ss.pop("_ps_nav_to", None)
    if to:
        ss["nav"] = to
    req = ss.pop("_ps_req", None)
    if req:
        _route(req)
    ss["_ps_rendering"] = False
    ss["_ps_live_prev"] = list(ss.get("_ps_live_cur", []))
    ss["_ps_live_cur"] = []
    for k, v in list(ss["_ps_wp"].items()):
        if k not in ss:
            try:
                ss[k] = v
            except Exception:
                ss["_ps_wp"].pop(k, None)


def mark_global():
    """Call after the sidebar: every keyed widget drawn so far is global (not page state)."""
    GLOBAL_KEYS.update(_rendered_widget_keys())


def track(key):
    """Count a non-widget session key (e.g. particle exclusions) as part of the current page's work."""
    ss = _ss()
    page = ss.get("_ps_page")
    if page:
        ex = ss["_ps_pages"].setdefault(page, {}).setdefault("extra", [])
        if key not in ex:
            ex.append(key)
            ss["_ps_base"].setdefault(page, {}).setdefault("extra:" + key, _sig(ss.get(key)))


def start_page(page):
    ss = _ss()
    ss["_ps_page"] = page
    ss["_ps_rendering"] = True


def end(page):
    """Call after the page has been drawn."""
    ss = _ss()
    ss["_ps_rendering"] = False
    rec = ss["_ps_pages"].setdefault(page, {})
    keys = rec.setdefault("keys", [])
    base = ss["_ps_base"].setdefault(page, {})
    for k, vt in _rendered_widget_keys().items():
        if k in GLOBAL_KEYS or k.startswith(_SKIP_PREFIX) or vt not in _ALLOWED or k not in ss:
            continue
        v = ss[k]
        if not _simple(v):
            continue
        if k not in keys:
            keys.append(k)
        base.setdefault("w:" + k, _sig(v))
        ss["_ps_wp"][k] = v
    for k in rec.get("up", []):
        base.setdefault("u:" + k, _sig(None))


def _sig(v):
    try:
        return hashlib.md5(repr(_jsonable(v)).encode()).hexdigest()[:12]
    except Exception:
        return repr(v)[:64]


def _current(page):
    ss = _ss()
    rec = ss["_ps_pages"].get(page, {})
    cur = {}
    for k in rec.get("keys", []):
        if k in ss:
            cur["w:" + k] = _sig(ss[k])
        elif k in ss["_ps_wp"]:
            cur["w:" + k] = _sig(ss["_ps_wp"][k])
    for k in rec.get("up", []):
        recs = ss["_ps_up"].get(k)
        cur["u:" + k] = _sig(_fp(recs) if recs else None)
    for k in rec.get("extra", []):
        cur["extra:" + k] = _sig(ss.get(k))
    return cur


def changes(page):
    """Names of the values that differ from the page's baseline."""
    base = _ss()["_ps_base"].get(page, {})
    out = []
    for k, v in _current(page).items():
        if k in base and base[k] != v:
            kind, name = k.split(":", 1)
            out.append(f"uploaded file(s) [{name}]" if kind == "u" else name)
    return out


def dirty(page):
    return bool(page) and page != "Home" and bool(changes(page))


# ----------------------------------------------------------------------------- save / discard
def save(page, to_disk=True):
    ss = _ss()
    ss["_ps_base"][page] = dict(_current(page))
    rec = ss["_ps_pages"].get(page, {})
    snap = {"page": page, "saved": _dt.datetime.now().isoformat(timespec="seconds"),
            "widgets": {k: (ss[k] if k in ss else ss["_ps_wp"].get(k)) for k in rec.get("keys", [])
                        if (k in ss or k in ss["_ps_wp"])},
            "uploads": {k: ss["_ps_up"][k] for k in rec.get("up", []) if ss["_ps_up"].get(k)},
            "extra": {k: ss.get(k) for k in rec.get("extra", []) if k in ss}}
    ss["_ps_saved"][page] = snap
    if to_disk and WS_DIR:
        try:
            _write_snapshot(snap)
        except Exception as e:      # never lose the session because the disk is full / read-only
            st.toast(f"Could not write the workspace file: {e}")
    return snap


def discard(page):
    """Revert the page to its last saved state, or to the defaults if it was never saved."""
    ss = _ss()
    rec = ss["_ps_pages"].get(page, {})
    snap = ss["_ps_saved"].get(page)
    for k in rec.get("keys", []):
        ss["_ps_wp"].pop(k, None)
        if k in ss:
            del ss[k]
    for k in rec.get("up", []):
        ss["_ps_up"].pop(k, None)
    for k in rec.get("extra", []):
        if k in ss:
            del ss[k]
    if snap:
        for k, v in snap["widgets"].items():
            ss["_ps_wp"][k] = v
        for k, v in snap["uploads"].items():
            ss["_ps_up"][k] = [tuple(x) for x in v]
        for k, v in snap.get("extra", {}).items():
            ss[k] = v
        ss["_ps_base"][page] = {}
        _rebase_from_snapshot(page, snap)
    else:
        ss["_ps_base"][page] = {}


def _rebase_from_snapshot(page, snap):
    b = _ss()["_ps_base"].setdefault(page, {})
    for k, v in snap["widgets"].items():
        b["w:" + k] = _sig(v)
    for k in _ss()["_ps_pages"].get(page, {}).get("up", []):
        recs = snap["uploads"].get(k)
        b["u:" + k] = _sig(_fp(recs) if recs else None)
    for k, v in snap.get("extra", {}).items():
        b["extra:" + k] = _sig(v)


def _write_snapshot(snap, ws_dir=None):
    WS = ws_dir or WS_DIR
    d = os.path.join(WS, slug(snap["page"]))
    os.makedirs(d, exist_ok=True)
    ups = {}
    for k, recs in snap["uploads"].items():
        ups[k] = []
        for i, (n, t, b) in enumerate(recs):
            fn = f"{slug(k)}__{i}__{re.sub(r'[^A-Za-z0-9._-]+', '_', n)}"
            with open(os.path.join(d, fn), "wb") as fh:
                fh.write(b)
            ups[k].append([n, t, fn])
    extra = {k: _jsonable(v) for k, v in snap.get("extra", {}).items() if _plain_data(v) or isinstance(v, dict)}
    out = {"page": snap["page"], "saved": snap["saved"], "widgets": {k: _jsonable(v) for k, v in snap["widgets"].items()},
           "uploads": ups, "extra": extra}
    tmp = os.path.join(WS, slug(snap["page"]) + ".json.tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1, default=str)
    os.replace(tmp, os.path.join(WS, slug(snap["page"]) + ".json"))


def current_snapshot(page, ss=None):
    """The page's present work (saved or not) in the snapshot format, without touching its baseline (v2.11 project files).
    `ss`: a plain dict copy of the session state (download callbacks run outside the script thread)."""
    ss = _ss() if ss is None else ss
    rec = ss["_ps_pages"].get(page, {})
    return {"page": page, "saved": _dt.datetime.now().isoformat(timespec="seconds"),
            "widgets": {k: (ss[k] if k in ss else ss["_ps_wp"].get(k)) for k in rec.get("keys", [])
                        if (k in ss or k in ss["_ps_wp"]) and _simple(ss[k] if k in ss else ss["_ps_wp"].get(k))},
            "uploads": {k: ss["_ps_up"][k] for k in rec.get("up", []) if ss["_ps_up"].get(k)},
            "extra": {k: ss.get(k) for k in rec.get("extra", []) if k in ss}}


def reset_all_pages():
    """Forget every page's widget values, uploads and extras in this session (before opening a project)."""
    ss = _ss()
    for page, rec in list(ss["_ps_pages"].items()):
        for k in rec.get("keys", []) + rec.get("extra", []):
            if k in ss:
                try:
                    del ss[k]
                except Exception:
                    pass
    for k in [k for k in list(ss.keys()) if str(k).startswith(AUTO_PREFIX)]:
        try:
            del ss[k]
        except Exception:
            pass
    for k in ("_ps_wp", "_ps_pages", "_ps_base", "_ps_up", "_ps_saved"):
        ss[k] = {}


def _load_snapshots(ws_dir=None, override=False):
    WS = ws_dir or WS_DIR
    if not WS or not os.path.isdir(WS):
        return
    ss = _ss()
    for fn in sorted(os.listdir(WS)):
        if not fn.endswith(".json"):
            continue
        try:
            js = json.load(open(os.path.join(WS, fn), encoding="utf-8"))
            page = js["page"]
            d = os.path.join(WS, slug(page))
            widgets = {k: _unjson(v) for k, v in js.get("widgets", {}).items()}
            uploads = {k: [(n, t, open(os.path.join(d, f), "rb").read()) for n, t, f in v
                           if os.path.exists(os.path.join(d, f))] for k, v in js.get("uploads", {}).items()}
            extra = {k: _unjson(v) for k, v in js.get("extra", {}).items()}
            ss["_ps_wp"].update(widgets)
            ss["_ps_up"].update({k: v for k, v in uploads.items() if v})
            for k, v in extra.items():
                if override:
                    ss[k] = v
                else:
                    ss.setdefault(k, v)
            rec = ss["_ps_pages"].setdefault(page, {})
            rec["keys"] = list(widgets)
            rec["up"] = list(uploads)
            rec["extra"] = list(extra)
            snap = {"page": page, "saved": js.get("saved"), "widgets": widgets, "uploads": uploads, "extra": extra}
            ss["_ps_saved"][page] = snap
            _rebase_from_snapshot(page, snap)
        except Exception:
            continue


def saved_pages():
    return {p: s.get("saved") for p, s in _ss()["_ps_saved"].items()}


def forget_all():
    ss = _ss()
    if WS_DIR and os.path.isdir(WS_DIR):
        import shutil
        shutil.rmtree(WS_DIR, ignore_errors=True)
    ss["_ps_saved"] = {}


# ----------------------------------------------------------------------------- navigation
def _route(target):
    ss = _ss()
    cur = ss.get("_ps_page")
    if enabled() and cur and cur != target and dirty(cur):
        ss["_ps_pending"] = target
    else:
        ss["nav"] = target


def on_nav():
    """on_change of the sidebar radio: undo the switch if the current page has unsaved work, and ask."""
    ss = _ss()
    new, cur = ss.get("nav"), ss.get("_ps_page")
    if enabled() and cur and new != cur and dirty(cur):
        ss["nav"] = cur
        ss["_ps_pending"] = new


def request(target):
    """Programmatic navigation (Home tiles, 'open page' buttons). Safe inside callbacks."""
    _ss()["_ps_req"] = target


def maybe_dialog():
    """Show the Save / Don't save / Cancel dialog while a held-back navigation is pending. Call after the sidebar."""
    ss = _ss()
    target = ss.pop("_ps_pending", None)
    if target:
        ss["_ps_dlg"] = (ss.get("_ps_page"), target)
    if ss.get("_ps_dlg"):
        _dialog(*ss["_ps_dlg"])


def _dlg_close():
    _ss()["_ps_dlg"] = None


def _dlg_save():
    ss = _ss()
    cur, target = ss.get("_ps_dlg") or (None, None)
    if cur:
        try:
            save(cur)
        except Exception as e:
            import traceback; ss["_ps_err"] = traceback.format_exc()
        ss["nav"] = target
    ss["_ps_dlg"] = None


def _dlg_discard():
    ss = _ss()
    cur, target = ss.get("_ps_dlg") or (None, None)
    if cur:
        discard(cur)
        ss["nav"] = target
    ss["_ps_dlg"] = None


@st.dialog("Unsaved work", on_dismiss=_dlg_close)
def _dialog(cur, target):
    if not _ss().get("_ps_dlg"):          # a button callback already decided (dialog-only rerun) → redraw the app
        st.rerun()
    ch = changes(cur)
    try:
        from mfs.ui import DISPLAY as _DSP
        cur_l, target_l = _DSP.get(cur, cur), _DSP.get(target, target)
    except Exception:
        cur_l, target_l = cur, target
    st.markdown(f"**Do you want to save your current work?**  \nYou changed **{len(ch)}** setting(s)/file(s) on "
                f"*{cur_l}* and are switching to *{target_l}*.")
    if ch:
        st.caption("Changed: " + ", ".join(f"`{c.split('|')[3] if c.startswith(AUTO_PREFIX) and c.count('|') >= 4 else c}`"
                                           for c in ch[:12]) + (" …" if len(ch) > 12 else ""))
    st.caption("Save keeps it in this session and in the workspace folder (restored on next start). Don't save reverts "
               "the page to its last saved state. Results already written to the sample store stay there.")
    c1, c2, c3 = st.columns(3)
    c1.button("💾 Save", type="primary", key="_psd_save", width="stretch", on_click=_dlg_save)
    c2.button("Don't save", key="_psd_discard", width="stretch", on_click=_dlg_discard)
    c3.button("Cancel", key="_psd_cancel", width="stretch", on_click=_dlg_close)


def status_badge(page):
    """Small sidebar line: unsaved / saved state of the current page, plus a manual Save button."""
    if not page or page == "Home":
        return
    ch = changes(page)
    s = _ss()["_ps_saved"].get(page)
    box = st.sidebar.container(key="ps_status")
    if ch:
        c1, c2 = box.columns([3, 2])
        c1.caption(f"● {len(ch)} unsaved change(s) on this page")
        if c2.button("💾 Save page", key="_ps_save_btn", width="stretch"):
            save(page); st.rerun()
    elif s:
        box.caption(f"✓ saved {str(s.get('saved', ''))[11:16] or ''}")


# ----------------------------------------------------------------------------- auto keys for un-keyed widgets
# Most Reticula widgets have no key, and Streamlit cannot hand an un-keyed widget its old value after a page switch.
# While a page is being drawn, value widgets without a key get a deterministic one built from the page, the widget
# type and ALL its parameters (label, options, min/max, default …) – exactly what identifies an un-keyed widget in
# Streamlit – so behaviour on the page is unchanged, but the value can be stored and restored by key.
_AUTO = ("checkbox", "toggle", "radio", "selectbox", "multiselect", "slider", "select_slider", "number_input",
         "text_input", "text_area", "date_input", "time_input", "color_picker", "pills", "segmented_control")
_NOHASH = {"help", "on_change", "args", "kwargs", "format_func", "label_visibility", "disabled", "placeholder"}
AUTO_PREFIX = "ak|"


def _auto_key(page, kind, args, kw):
    parts = [repr(a) for a in args] + [f"{k}={v!r}" for k, v in sorted(kw.items()) if k not in _NOHASH]
    h = hashlib.md5(("\x1f".join(parts)).encode("utf-8", "replace")).hexdigest()[:10]
    lab = re.sub(r"\s+", " ", str(args[0]) if args else "")[:40]
    base = f"{AUTO_PREFIX}{slug(page)}|{kind}|{lab}|{h}"
    k, n = base, 1
    while _key_live(k):
        n += 1
        k = f"{base}#{n}"
    return k


def _rendering():
    if st.session_state.get("_ps_rendering"):
        return True
    try:                                        # fragment rerun (e.g. the particle picker): the page is "rendering"
        from streamlit.runtime.scriptrunner_utils.script_run_context import get_script_run_ctx
        return bool(get_script_run_ctx().fragment_ids_this_run)
    except Exception:
        return False


def install_autokey_hook():
    import inspect
    from streamlit.delta_generator import DeltaGenerator
    for kind in _AUTO:
        orig = getattr(DeltaGenerator, kind, None)
        if orig is None or getattr(orig, "_ps_wrapped", False):
            continue
        try:
            names = list(inspect.signature(orig).parameters)
            kidx = names.index("key") - 1          # minus self
        except Exception:
            continue

        def make(orig=orig, kind=kind, kidx=kidx):
            def w(self, *args, **kw):
                try:
                    if kw.get("key") is None and len(args) <= kidx and _rendering():
                        kw.pop("key", None)
                        kw["key"] = _auto_key(st.session_state.get("_ps_page") or "", kind, args, kw)
                except Exception:
                    pass
                return orig(self, *args, **kw)
            w._ps_wrapped = True
            w.__doc__ = orig.__doc__
            w.__name__ = kind
            return w
        setattr(DeltaGenerator, kind, make())


# ----------------------------------------------------------------------------- v2.11 undo / redo (every page)
# After every run the page's undoable state (its value widgets, tracked extras such as sample-list order/visibility and
# particle exclusions, plus keys registered with track_undo, e.g. manual outline edits) is compared with the head of the
# page's history; a change pushes a new entry. Undo/Redo restore an entry from a button callback, i.e. before the widgets
# of the next run are drawn. Custom components that keep their own state get `_ps_rev` to adopt the restored value.
HIST_MAX = 60


def track_undo(key):
    """Make a non-widget session key part of the current page's undo history (not of its saved page state)."""
    ss = _ss()
    page = ss.get("_ps_page")
    if page:
        u = ss["_ps_pages"].setdefault(page, {}).setdefault("undo", [])
        if key not in u:
            u.append(key)


def _undo_state(page):
    import copy
    ss = _ss()
    rec = ss["_ps_pages"].get(page, {})
    w = {}
    for k in rec.get("keys", []):
        v = ss[k] if k in ss else ss["_ps_wp"].get(k, None)
        if (k in ss or k in ss["_ps_wp"]) and _simple(v):
            w[k] = v
    x = {}
    for k in list(rec.get("extra", [])) + list(rec.get("undo", [])):
        if k in ss:
            try:
                x[k] = copy.deepcopy(ss[k])
            except Exception:
                pass
    return {"w": w, "x": x}


def _hist(page):
    return _ss().setdefault("_ps_hist", {}).setdefault(page, {"stack": [], "pos": -1, "skip": False})


def _label(k):
    if k.startswith(AUTO_PREFIX) and k.count("|") >= 4:
        return k.split("|")[3] or k.split("|")[2]
    if k.startswith("_sl|"):
        return "sample list"
    if k.startswith("psx_") or k == "_ps_excl":
        return "particle exclusions"
    if k.startswith("pse_"):
        return "outline edits"
    if k.startswith("psc_"):
        return "click measurements"
    return k


def _diff(a, b):
    ks = []
    for part in ("w", "x"):
        for k in set(a[part]) | set(b[part]):
            if _sig(a[part].get(k)) != _sig(b[part].get(k)):
                ks.append(_label(k))
    out = []
    for k in ks:
        if k not in out:
            out.append(k)
    return out


def history_record(page):
    """Call after the page has been drawn (after end())."""
    if not page or page == "Home":
        return
    h = _hist(page)
    cur = _undo_state(page)
    sig = _sig(cur)
    if h["skip"] and h["stack"]:                    # the run right after an undo/redo: adopt what was actually drawn
        h["skip"] = False
        h["stack"][h["pos"]].update(sig=sig, state=cur)
        return
    h["skip"] = False
    if h["stack"] and h["stack"][h["pos"]]["sig"] == sig:
        return
    lab = _diff(h["stack"][h["pos"]]["state"], cur) if h["stack"] else []
    if h["stack"] and not lab:                      # only key sets changed (new widget appeared): no user edit
        h["stack"][h["pos"]].update(sig=sig, state=cur)
        return
    del h["stack"][h["pos"] + 1:]
    h["stack"].append({"sig": sig, "state": cur, "label": ", ".join(lab[:3]) + (" …" if len(lab) > 3 else "")})
    if len(h["stack"]) > HIST_MAX:
        del h["stack"][0]
    h["pos"] = len(h["stack"]) - 1


def _restore(page, state):
    import copy
    ss = _ss()
    now = _undo_state(page)
    for k, v in state["w"].items():
        if _sig(now["w"].get(k)) != _sig(v):
            try:
                ss[k] = v
                ss["_ps_wp"][k] = v
            except Exception:
                pass
    for k, v in state["x"].items():
        ss[k] = copy.deepcopy(v)
    for k in set(now["x"]) - set(state["x"]):       # registered later than this entry: back to "not set"
        if k in ss and not str(k).startswith("_ps_"):
            try:
                del ss[k]
            except Exception:
                pass
    ss["_ps_rev"] = int(ss.get("_ps_rev", 0)) + 1


def can_undo(page):
    h = _hist(page); return h["pos"] > 0


def can_redo(page):
    h = _hist(page); return 0 <= h["pos"] < len(h["stack"]) - 1


def undo(page):
    h = _hist(page)
    if h["pos"] > 0:
        h["pos"] -= 1; h["skip"] = True
        _restore(page, h["stack"][h["pos"]]["state"])


def redo(page):
    h = _hist(page)
    if h["pos"] < len(h["stack"]) - 1:
        h["pos"] += 1; h["skip"] = True
        _restore(page, h["stack"][h["pos"]]["state"])


def history_box(page):
    """Sidebar ↶ Undo / ↷ Redo for the current page."""
    if not page or page == "Home":
        return
    h = _hist(page)
    box = st.sidebar.container(key="ps_undo")
    c1, c2, c3 = box.columns([1, 1, 1.3], vertical_alignment="center")
    u_lab = h["stack"][h["pos"]]["label"] if can_undo(page) else ""
    r_lab = h["stack"][h["pos"] + 1]["label"] if can_redo(page) else ""
    c1.button("↶ Undo", key="_ps_undo", on_click=undo, args=(page,), disabled=not can_undo(page), width="stretch",
              help=("Undo: " + u_lab) if u_lab else "Nothing to undo on this page")
    c2.button("↷ Redo", key="_ps_redo", on_click=redo, args=(page,), disabled=not can_redo(page), width="stretch",
              help=("Redo: " + r_lab) if r_lab else "Nothing to redo")
    c3.caption(f"{max(h['pos'], 0)} step(s)" + (f" · {len(h['stack']) - 1 - h['pos']} redo" if can_redo(page) else ""))
