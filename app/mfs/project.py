"""v2.11 project files: the whole session in one `.reticula` file (a zip).

    manifest.json            app, version, created, current page, pages
    workspace/<page>.json    every page's widget values / extras (pagestate snapshot format, saved or not)
    workspace/<page>/...     the uploaded files of that page
    globals.json             figure style + colours + other sidebar settings, sample store, defect-calculator input
Opening a project replaces the session's page work (and the on-disk workspace, so it is also there after a restart)."""
from __future__ import annotations

import datetime as _dt
import io
import json
import os
import shutil
import tempfile
import zipfile
from pathlib import Path

import streamlit as st

from . import pagestate as PS

EXT = ".reticula"
FORMAT = 1
GLOBAL_KEYS = ("_lv", "ps_templates", "band_user", "nmr_lib_user", "series_colours", "si_collect", "sample_store", "def_in", "store_persist", "undo_limit")
GLOBAL_PREFIX = ("gs_",)            # global figure style widgets (M7)
GLOBAL_SKIP = {"gs_reset"}          # buttons: their state is read-only and must never be written back


def _global_ok(k: str) -> bool:
    return k not in GLOBAL_SKIP and (k in GLOBAL_KEYS or k.startswith(GLOBAL_PREFIX))
_OPENED_FROM_ARGV = set()           # per process: a project passed on the command line is opened once


def _jsonable_globals(ss):
    out = {}
    for k in list(ss.keys()):
        k = str(k)
        if _global_ok(k):
            v = ss[k]
            try:
                json.dumps(PS._jsonable(v), default=str)
            except Exception:
                continue
            out[k] = PS._jsonable(v)
    return out


def export_bytes(current_page: str | None = None, version: str = "", ss=None) -> bytes:
    ss = PS._ss() if ss is None else ss
    tmp = Path(tempfile.mkdtemp(prefix="reticula_prj_"))
    try:
        ws = tmp / "workspace"; ws.mkdir()
        pages = []
        for page in list(ss["_ps_pages"]):
            snap = PS.current_snapshot(page, ss)
            if snap["widgets"] or snap["uploads"] or snap["extra"]:
                PS._write_snapshot(snap, str(ws)); pages.append(page)
        (tmp / "globals.json").write_text(json.dumps(_jsonable_globals(ss), ensure_ascii=False, indent=1, default=str), encoding="utf-8")
        man = {"app": "Reticula", "format": FORMAT, "version": version, "created": _dt.datetime.now().isoformat(timespec="seconds"),
               "page": current_page, "pages": pages}
        (tmp / "manifest.json").write_text(json.dumps(man, indent=1), encoding="utf-8")
        bio = io.BytesIO()
        with zipfile.ZipFile(bio, "w", zipfile.ZIP_DEFLATED) as z:
            for p in sorted(tmp.rglob("*")):
                if p.is_file():
                    z.write(p, p.relative_to(tmp).as_posix())
        return bio.getvalue()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def read_manifest(data: bytes) -> dict:
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        man = json.loads(z.read("manifest.json"))
    if man.get("app") != "Reticula":
        raise ValueError("not a Reticula project file")
    return man


def request_open(data: bytes, name: str = "project"):
    """Queue a project to be opened at the top of the next run (widgets can only be reset before they are drawn)."""
    read_manifest(data)
    st.session_state["_proj_pending"] = (name, data)


def apply_pending(store_path: str | None = None):
    """Call right after pagestate.begin(), before any widget. Returns the manifest if a project was opened."""
    ss = st.session_state
    item = ss.pop("_proj_pending", None)
    if not item:
        return None
    name, data = item
    tmp = Path(tempfile.mkdtemp(prefix="reticula_open_"))
    try:
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            for m in z.namelist():                      # no path traversal
                if m.startswith("/") or ".." in Path(m).parts:
                    raise ValueError(f"unsafe member {m}")
            z.extractall(tmp)
        man = json.loads((tmp / "manifest.json").read_text(encoding="utf-8"))
        PS.reset_all_pages()
        PS._load_snapshots(str(tmp / "workspace"), override=True)
        for k, v in ss["_ps_wp"].items():               # what begin() does for a fresh run
            try:
                ss[k] = v
            except Exception:
                pass
        g = json.loads((tmp / "globals.json").read_text(encoding="utf-8")) if (tmp / "globals.json").exists() else {}
        for k, v in g.items():
            if not _global_ok(str(k)):
                continue
            try:
                ss[k] = PS._unjson(v)
            except Exception:
                pass
        if store_path and "sample_store" in g:
            try:
                Path(store_path).parent.mkdir(parents=True, exist_ok=True)
                Path(store_path).write_text(json.dumps(g["sample_store"], ensure_ascii=False, indent=1, default=str), encoding="utf-8")
            except Exception:
                pass
        if PS.WS_DIR:                                   # the opened project becomes the saved workspace (restored on restart)
            shutil.rmtree(PS.WS_DIR, ignore_errors=True)
            try:
                shutil.copytree(tmp / "workspace", PS.WS_DIR)
            except Exception:
                pass
        if man.get("page"):
            ss["nav"] = man["page"]
        ss["_proj_opened"] = {"name": name, **man}
        return man
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def open_from_argv():
    """A project path handed over by the launcher (double-click on a .reticula file) -> opened once per process."""
    p = os.environ.get("MFS_OPEN_PROJECT")
    if not p or p in _OPENED_FROM_ARGV or not os.path.isfile(p):
        return
    _OPENED_FROM_ARGV.add(p)
    try:
        request_open(Path(p).read_bytes(), Path(p).name)
    except Exception as e:
        st.session_state["_proj_err"] = f"{Path(p).name}: {e}"


def sidebar(current_page, version):
    """📁 Project box: save the whole session to a .reticula file / open one."""
    ss = st.session_state
    with st.sidebar.expander("📁 Project (.reticula)", expanded=False):
        op = ss.get("_proj_opened")
        if op:
            st.caption(f"Opened **{op['name']}** (saved {str(op.get('created', ''))[:16].replace('T', ' ')}, Reticula {op.get('version', '?')}, "
                       f"{len(op.get('pages', []))} page(s)).")
        if ss.get("_proj_err"):
            st.error(ss.pop("_proj_err"))
        nm = st.text_input("Project name", ss.get("_proj_name", "my_project"), key="_proj_name")
        PS._ss()
        frozen = {k: ss[k] for k in list(ss.keys())}      # the callback runs outside the script thread
        st.download_button("💾 Save project", data=lambda: export_bytes(current_page, version, frozen),
                           file_name=(nm or "project").strip().replace(" ", "_") + EXT, mime="application/zip",
                           key="_proj_save", width="stretch",
                           help="Every page's settings, uploaded files, sample store and figure style in one file.")
        up = st.file_uploader("Open a project", type=["reticula", "zip"], key="_proj_up")
        if up is not None:
            def _go(up=up):
                try:
                    request_open(up.getvalue(), up.name)
                except Exception as e:
                    ss["_proj_err"] = f"{up.name}: {e}"
            st.button(f"📂 Open {up.name} (replaces the current work)", key="_proj_open", on_click=_go, width="stretch")
