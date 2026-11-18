"""v2.16: one-click problem report. A small ZIP the user can send: program version and environment, the page, the settings in the page, the last
error with its traceback, the names of the loaded files and the figures made - NEVER the data itself (no spectra, images or tables) and no
personal paths (the user's home folder is replaced by ~)."""
from __future__ import annotations

import datetime
import io
import json
import os
import platform
import re
import sys
import traceback
import zipfile

_SKIP = re.compile(r"^(fc_|sc_|data_|exportbar_|add_|si_|db_|bp_|FormSubmitter|nav$|cmap$|house_|store_|_)")


def _anon(text: str) -> str:
    """Replaces the home folder and the user name in paths."""
    home = os.path.expanduser("~")
    out = str(text).replace(home, "~").replace(home.replace("\\", "/"), "~")
    user = os.path.basename(home)
    if user and len(user) > 2:
        out = re.sub(re.escape(user), "<user>", out, flags=re.I)
    return out


def record_error(session: dict, page: str, exc: BaseException):
    """Remembers the last error of a page (called from the page dispatcher)."""
    session["_last_error"] = {"page": page, "type": type(exc).__name__, "message": _anon(str(exc))[:800],
                              "traceback": _anon("".join(traceback.format_exception(type(exc), exc, exc.__traceback__)))[-6000:],
                              "time": datetime.datetime.now().isoformat(timespec="seconds")}


def _versions() -> dict:
    out = {}
    for m in ("streamlit", "numpy", "scipy", "pandas", "matplotlib", "pymatgen", "PIL", "openpyxl", "rdkit", "altair"):
        try:
            mod = __import__(m)
            out[m] = getattr(mod, "__version__", "?")
        except Exception:
            out[m] = "not installed"
    return out


def build(session, page: str, version: str, sources=None, note: str = "") -> bytes:
    """The report ZIP (bytes). `session` = st.session_state (any mapping)."""
    params = {}
    for k, v in list(session.items()):
        if not isinstance(k, str) or _SKIP.match(k):
            continue
        if isinstance(v, (int, float, bool)) or (isinstance(v, str) and len(v) <= 120):
            if k.startswith("ak|"):
                parts = k.split("|")
                k = parts[3] if len(parts) >= 4 else k
            params[k] = v
    figs = [{"figure": k, "page": r.get("page"), "title": r.get("title"), "rows": r.get("n_rows")} for k, r in (session.get("si_figs") or {}).items()]
    store = {sid: sorted((rec.get("sections") or {})) for sid, rec in (session.get("sample_store") or {}).items()}
    lv = session.get("_lv") or {}
    rep = {"program": "Figuryte", "version": version, "created": datetime.datetime.now().isoformat(timespec="seconds"),
           "page": page, "detail_level": lv.get(page) if lv.get(page) not in (None, "global") else lv.get("global"),
           "platform": platform.platform(), "python": sys.version.split()[0], "frozen": bool(getattr(sys, "frozen", False)),
           "libraries": _versions(), "loaded_files": [_anon(s) for s in (sources or [])][:60], "settings_on_this_page": params,
           "figures_made": figs, "sample_store_sections": store, "last_error": session.get("_last_error"),
           "user_note": note.strip()[:4000]}
    md = (f"# Figuryte problem report\n\nVersion {version} · page **{page}** · {rep['created']}\n\n"
          f"**What the user was doing:** {note.strip() or '(not given)'}\n\n")
    err = session.get("_last_error")
    if err:
        md += f"**Last error on '{err['page']}':** `{err['type']}: {err['message']}`\n\n```\n{err['traceback']}\n```\n\n"
    else:
        md += "No error was recorded in this session.\n\n"
    md += ("The ZIP contains the version, environment, the settings of the page, the file NAMES and the figures made - no data. "
           "Send it together with the file that triggered the problem if possible.\n")
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("report.json", json.dumps(rep, indent=1, ensure_ascii=False, default=str))
        z.writestr("report.md", md)
    return buf.getvalue()
