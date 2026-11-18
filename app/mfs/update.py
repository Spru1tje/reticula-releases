"""v2.11 update check (opt-in, asks once). Never downloads or installs anything by itself.

On start (at most once a day per process) the release feed is read in a background thread:
    default  https://api.github.com/repos/Spru1tje/reticula-releases/releases/latest   (public repo, no token)
    override env MFS_UPDATE_URL or the URL saved in the settings
If the release tag is newer than this version and not skipped, a non-blocking banner offers
Download (opens the installer link in the browser) / Later (this session) / Skip this version.
Offline, timeouts, rate limits or a malformed answer -> silently nothing."""
from __future__ import annotations

import json
import os
import platform
import re
import threading
import time
from pathlib import Path

DEFAULT_URL = "https://api.github.com/repos/Spru1tje/reticula-releases/releases/latest"
TIMEOUT = 6.0
_STATE = {"result": None, "error": None, "running": False, "t": 0.0, "url": None}
_LOCK = threading.Lock()


# ----------------------------------------------------------------------------- settings
def settings_path(data_dir) -> Path:
    return Path(data_dir) / "update.json"


def load_settings(data_dir) -> dict:
    try:
        return json.loads(settings_path(data_dir).read_text(encoding="utf-8"))
    except Exception:
        return {}


def save_settings(data_dir, d: dict):
    try:
        p = settings_path(data_dir); p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(d, indent=1), encoding="utf-8")
    except Exception:
        pass


def feed_url(cfg: dict) -> str:
    return os.environ.get("MFS_UPDATE_URL") or cfg.get("url") or DEFAULT_URL


# ----------------------------------------------------------------------------- versions / assets
def parse_version(s) -> tuple:
    m = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", str(s or ""))
    return tuple(int(x or 0) for x in m.groups()) if m else (0,)


def is_newer(latest, current) -> bool:
    return parse_version(latest) > parse_version(current)


_MAC_ARCH = {"arm64": ("arm64", "aarch64", "apple-silicon", "applesilicon"), "x86_64": ("x86_64", "x64", "amd64", "intel")}


def _mac_arch(machine: str | None) -> str:
    m = (machine or platform.machine() or "").lower()
    return "arm64" if m in ("arm64", "aarch64") else ("x86_64" if m in ("x86_64", "amd64", "i386", "x64") else m)


def pick_asset(rel: dict, system: str | None = None, machine: str | None = None) -> dict | None:
    """Windows: Reticula_Setup_<ver>.exe (single file). macOS (v2.12): the .dmg/.pkg named macOS_<machine>
    (arm64 / x86_64 from platform.machine()), then a universal one, then an untagged one; an Intel Mac is never offered
    an arm64-only image (an Apple-silicon Mac may take an Intel one, it runs under Rosetta). None -> release page."""
    system = system or platform.system()
    assets = rel.get("assets") or []
    names = [(a.get("name") or "", a) for a in assets]
    if system == "Windows":
        for n, a in names:
            if re.fullmatch(r"Reticula_Setup_[\w.\-]+\.exe", n, re.I):
                return a
        for n, a in names:
            if n.lower().endswith(".exe"):
                return a
    elif system == "Darwin":
        arch = _mac_arch(machine)
        imgs = [(n, a) for ext in (".dmg", ".pkg") for n, a in names if n.lower().endswith(ext)]
        imgs += [(n, a) for n, a in names if "mac" in n.lower() and n.lower().endswith(".zip")]

        def tag(n):
            low = n.lower()
            for k, al in _MAC_ARCH.items():
                if any(re.search(r"(^|[_\-. ])" + re.escape(x) + r"([_\-. ]|$)", low) for x in al):
                    return k
            return "universal" if "universal" in low else None
        for want in (arch, "universal", None):
            for n, a in imgs:
                if tag(n) == want:
                    return a
        if arch == "arm64":
            for n, a in imgs:
                if tag(n) == "x86_64":
                    return a
    return None


def summarize(rel: dict, current: str, system: str | None = None) -> dict | None:
    tag = rel.get("tag_name") or rel.get("name") or ""
    if not tag or rel.get("draft") or rel.get("prerelease"):
        return None
    ver = ".".join(str(x) for x in parse_version(tag))
    a = pick_asset(rel, system)
    return {"version": ver, "tag": tag, "newer": is_newer(tag, current), "url": (a or {}).get("browser_download_url") or rel.get("html_url"),
            "asset": (a or {}).get("name"), "size": (a or {}).get("size"), "page": rel.get("html_url"),
            "notes": (rel.get("body") or "")[:1200], "published": rel.get("published_at")}


# ----------------------------------------------------------------------------- fetch
def fetch(url: str, timeout: float = TIMEOUT) -> dict:
    import urllib.request
    req = urllib.request.Request(url, headers={"Accept": "application/vnd.github+json", "User-Agent": "Reticula-update-check",
                                               "X-GitHub-Api-Version": "2022-11-28"})
    with urllib.request.urlopen(req, timeout=timeout) as r:          # nosec - fixed https URL or the user's own
        return json.loads(r.read(2_000_000).decode("utf-8"))


def _worker(url, current):
    try:
        res = summarize(fetch(url), current)
        with _LOCK:
            _STATE.update(result=res, error=None)
    except Exception as e:                                           # offline, 403 rate limit, DNS, JSON … -> silent
        with _LOCK:
            _STATE.update(result=None, error=f"{type(e).__name__}: {e}"[:200])
    finally:
        with _LOCK:
            _STATE["running"] = False


def start_check(url, current, force=False, max_age=86400.0):
    """Start a background check unless one ran recently for the same URL. Returns True if one is running now."""
    with _LOCK:
        if _STATE["running"]:
            return True
        if not force and _STATE["url"] == url and time.time() - _STATE["t"] < max_age and (_STATE["result"] or _STATE["error"]):
            return False
        _STATE.update(running=True, t=time.time(), url=url, result=None, error=None)
    threading.Thread(target=_worker, args=(url, current), daemon=True, name="reticula-update").start()
    return True


def status():
    with _LOCK:
        return dict(_STATE)


# ----------------------------------------------------------------------------- UI
def banner(data_dir, current):
    """Top-of-page banner (non-blocking). Call once per run before the page."""
    import streamlit as st
    ss = st.session_state
    if os.environ.get("MFS_NO_UPDATE_CHECK") == "1":
        return
    cfg = load_settings(data_dir)
    if cfg.get("opt_in") is None:                                   # ask once
        if ss.get("_upd_asked"):
            return
        box = st.container(border=True, key="upd_ask")
        c1, c2, c3 = box.columns([5, 1.3, 1], vertical_alignment="center")
        c1.markdown("🔔 **Check for Figuryte updates when it starts?** It only reads the public release list on GitHub "
                    "(no data is sent, nothing is installed). You can change this later under ℹ️ About.")
        if c2.button("Yes, check", key="upd_yes", type="primary", width="stretch"):
            cfg["opt_in"] = True; save_settings(data_dir, cfg); ss["_upd_asked"] = True; st.rerun()
        if c3.button("No", key="upd_no", width="stretch"):
            cfg["opt_in"] = False; save_settings(data_dir, cfg); ss["_upd_asked"] = True; st.rerun()
        return
    if not cfg.get("opt_in"):
        return
    url = feed_url(cfg)
    running = start_check(url, current)
    s = status()
    if running or s["running"]:
        _poll()                                                     # reruns the app once when the answer is in
        return
    r = s["result"]
    if not r or not r["newer"] or r["version"] == cfg.get("skip") or ss.get("_upd_later") == r["version"]:
        return
    box = st.container(border=True, key="upd_banner")
    c1, c2, c3, c4 = box.columns([4.6, 1.4, 0.9, 1.5], vertical_alignment="center")
    c1.markdown(f"⬆️ **Figuryte {r['version']} is available** (you have {current})."
                + (f"  \n<span style='opacity:.7'>{r['asset']}" + (f" · {r['size'] / 1e6:.0f} MB" if r.get("size") else "") + "</span>"
                   if r.get("asset") else ""), unsafe_allow_html=True)
    c2.link_button("⬇️ Download", r["url"], type="primary", width="stretch",
                   help="Opens the download in your browser. Run the installer yourself when convenient; your data and settings are kept.")
    if c3.button("Later", key="upd_later", width="stretch"):
        ss["_upd_later"] = r["version"]; st.rerun()
    if c4.button("Skip this version", key="upd_skip", width="stretch"):
        cfg["skip"] = r["version"]; save_settings(data_dir, cfg); st.rerun()
    if r.get("notes"):
        with box.expander("What's new"):
            st.markdown(r["notes"])


def _poll():
    import streamlit as st

    @st.fragment(run_every=1.5)
    def _wait():
        if not status()["running"]:
            st.rerun()
    _wait()


def settings_box(data_dir, current):
    """Inside the About expander: opt in/out, URL, check now."""
    import streamlit as st
    cfg = load_settings(data_dir)
    def _cb():
        c = load_settings(data_dir); c["opt_in"] = bool(st.session_state.get("_upd_opt")); save_settings(data_dir, c)
    st.session_state["_upd_opt"] = bool(cfg.get("opt_in"))          # the file is the truth (the banner may have changed it)
    on = st.checkbox("Check for updates on start", key="_upd_opt", on_change=_cb)
    s = status()
    if s.get("result"):
        r = s["result"]
        st.caption(f"Latest release: {r['version']}" + (" – newer" if r["newer"] else " – you are up to date")
                   + (f" (skipped)" if r["version"] == cfg.get("skip") else ""))
    elif s.get("error") and on:
        st.caption("Last check failed (offline?) – will try again next start.")
    if st.button("Check now", key="_upd_now", disabled=not on):
        cfg.pop("skip", None); save_settings(data_dir, cfg); st.session_state.pop("_upd_later", None)
        start_check(feed_url(cfg), current, force=True); st.rerun()
    st.caption(f"Feed: {feed_url(cfg)}")
