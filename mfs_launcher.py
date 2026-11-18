"""Reticula - desktop launcher (Windows; macOS/Linux-safe: Windows-only calls are guarded).

Double-click -> a small control window opens, the Streamlit server starts in a hidden child process
(same .exe, '--serve' mode) and the default browser opens at http://127.0.0.1:<port>/ once /_stcore/health is ok.
Closing the control window (or 'Quit') stops the server. No console window."""
import os, sys, socket, subprocess, threading, time, webbrowser

APP_NAME = "Reticula"
TAGLINE = "Publication figures for framework materials"
FROZEN = getattr(sys, "frozen", False)
BASE = sys._MEIPASS if FROZEN else os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.environ.get("MFS_APP_DIR") or os.path.join(BASE, "app")
APP = os.path.join(APP_DIR, "app.py")
ICON = os.path.join(BASE, "reticula.ico")


def documents_dir():
    if os.name != "nt":
        return os.path.join(os.path.expanduser("~"), "Documents")
    try:  # the real Documents folder (also when redirected to OneDrive)
        import ctypes
        from ctypes import wintypes
        buf = ctypes.create_unicode_buffer(wintypes.MAX_PATH)
        if ctypes.windll.shell32.SHGetFolderPathW(None, 5, None, 0, buf) == 0 and buf.value:
            return buf.value
    except Exception:
        pass
    return os.path.join(os.path.expanduser("~"), "Documents")


# macOS: settings, workspace, logs and exports live in Application Support (not Documents).
if sys.platform == "darwin":
    DATA_DIR = os.path.join(os.path.expanduser("~"), "Library", "Application Support", APP_NAME)
else:
    DATA_DIR = os.path.join(documents_dir(), APP_NAME)
EXPORT_DIR = os.path.join(DATA_DIR, "exports")
if os.name == "nt":
    LOG_DIR = os.path.join(os.environ.get("LOCALAPPDATA") or DATA_DIR, APP_NAME)
elif sys.platform == "darwin":
    LOG_DIR = DATA_DIR
else:
    LOG_DIR = os.path.join(os.environ.get("XDG_STATE_HOME") or os.path.join(os.path.expanduser("~"), ".local", "state"), APP_NAME)
LOG = os.path.join(LOG_DIR, "server_log.txt")



def open_path(path):
    """Open a folder/file with the OS file manager. os.startfile is Windows-only."""
    if not path:
        return
    try:
        if os.name == "nt":
            os.startfile(path)  # noqa: E1101
        elif sys.platform == "darwin":
            subprocess.Popen(["open", path])
        else:
            subprocess.Popen(["xdg-open", path])
    except Exception:
        pass


HOST = "127.0.0.1"          # one address everywhere (server, health check, browser URL): no localhost/::1 mix-ups
PIDFILE = os.path.join(LOG_DIR, "server.pid")


def app_version():
    try:
        sys.path.insert(0, APP_DIR)
        from mfs import __version__
        return __version__
    except Exception:
        return "0"


def port_open(port):
    for fam, addr in ((socket.AF_INET, "127.0.0.1"), (socket.AF_INET6, "::1")):
        try:
            with socket.socket(fam) as s:
                s.settimeout(0.3)
                if s.connect_ex((addr, port)) == 0:
                    return True
        except OSError:
            pass
    return False


def http_get(url, timeout=2.0):
    import urllib.request
    req = urllib.request.Request(url, headers={"Cache-Control": "no-cache"})
    with urllib.request.build_opener(urllib.request.ProxyHandler({})).open(req, timeout=timeout) as r:   # never via a proxy
        return r.status, r.headers.get("Content-Type", ""), r.read()


def health(port):
    try:
        st, _, body = http_get(f"http://{HOST}:{port}/_stcore/health", 1.5)
        return st == 200 and body.strip() == b"ok"
    except Exception:
        return False


def frontend_check(port):
    """Index page + one JS bundle: returns (ok, message). A JS bundle served as text/plain = the Windows registry
    .js MIME mapping problem (white page)."""
    import re as _re
    try:
        st, ct, body = http_get(f"http://{HOST}:{port}/", 5)
        m = _re.search(rb'src="(\.?/?static/js/[^"]+\.js)"', body)
        if st != 200 or not m:
            return False, f"index.html incomplete (HTTP {st}, {len(body)} bytes, no JS bundle)"
        src = m.group(1).decode().lstrip(".").lstrip("/")
        st2, ct2, js = http_get(f"http://{HOST}:{port}/{src}", 10)
        if st2 != 200 or len(js) < 1000:
            return False, f"JS bundle missing ({src}: HTTP {st2}, {len(js)} bytes)"
        if "javascript" not in ct2:
            return False, f"JS bundle served as '{ct2}' (expected application/javascript)"
        return True, f"frontend ok ({src.rsplit('/', 1)[-1]}: {ct2})"
    except Exception as e:
        return False, f"frontend check failed: {e}"


def _kill_pid(pid):
    try:
        if os.name == "nt":
            subprocess.run(["taskkill", "/PID", str(pid), "/T", "/F"], creationflags=0x08000000, capture_output=True, timeout=15)
        else:
            os.kill(pid, 9)
        return True
    except Exception:
        return False


def stop_old_instances(log=print):
    """A Reticula server left over from an earlier start (PID file), or an older Reticula (2.2-2.6, no PID file) still
    listening on 8501-8510: stop it, so the browser never talks to an old frontend. Other programs are left alone."""
    try:
        if os.path.exists(PIDFILE):
            pid = int(open(PIDFILE, encoding="utf-8").read().split()[0])
            if pid != os.getpid() and _kill_pid(pid):
                log(f"stopped a previous Reticula server (PID {pid})")
            os.remove(PIDFILE)
    except Exception:
        pass
    if os.name != "nt":
        return
    try:
        out = subprocess.run(["netstat", "-ano", "-p", "TCP"], capture_output=True, text=True, creationflags=0x08000000, timeout=15).stdout
        pids = set()
        for line in out.splitlines():
            t = line.split()
            if len(t) >= 5 and t[3].upper() == "LISTENING" and t[1].rsplit(":", 1)[-1].isdigit() and 8501 <= int(t[1].rsplit(":", 1)[-1]) <= 8510:
                pids.add(int(t[4]))
        for pid in pids:
            if pid == os.getpid():
                continue
            name = subprocess.run(["tasklist", "/FI", f"PID eq {pid}", "/FO", "CSV", "/NH"], capture_output=True, text=True,
                                  creationflags=0x08000000, timeout=15).stdout
            if "reticula" in name.lower() or "reticula" in name.lower() or "mof_figure_studio" in name.lower() or "moffigure" in name.lower():
                if _kill_pid(pid):
                    log(f"stopped an older Reticula server (PID {pid})")
    except Exception as e:
        log(f"old-instance check skipped: {e}")


def fix_mimetypes():
    """Windows may map .js to text/plain (or .css to something odd) in HKCR; browsers then refuse the Streamlit bundle
    (strict MIME checking for module scripts) and show a WHITE PAGE. Read the registry first, then override."""
    import mimetypes
    mimetypes.init()
    for ext, t in ((".js", "application/javascript"), (".mjs", "application/javascript"), (".css", "text/css"),
                   (".json", "application/json"), (".map", "application/json"), (".svg", "image/svg+xml"),
                   (".woff", "font/woff"), (".woff2", "font/woff2"), (".ttf", "font/ttf"), (".wasm", "application/wasm"),
                   (".html", "text/html"), (".png", "image/png"), (".ico", "image/x-icon")):
        mimetypes.add_type(t, ext, strict=True)
        mimetypes.add_type(t, ext, strict=False)


def free_port(start=8501):
    for p in range(start, start + 100):
        if not port_open(p):
            return p
    return start



# ----------------------------------------------------------------------------- v2.11: loading screen (splash page)
STEPS = 8   # number of '[startup] i/N' steps the server child logs (see prewarm)

SPLASH_HTML = """<!doctype html><html><head><meta charset="utf-8"><title>Starting Reticula…</title>
<link rel="icon" href="data:image/png;base64,%(logo)s">
<style>
html,body{margin:0;height:100%%;font-family:"Segoe UI",system-ui,Arial,sans-serif;background:#0B1F3A;color:#E6EEF8}
.wrap{display:flex;flex-direction:column;align-items:center;justify-content:center;height:100%%;text-align:center}
img{width:112px;height:112px;margin-bottom:14px;filter:drop-shadow(0 6px 18px rgba(0,0,0,.35))}
h1{margin:0;font-size:30px;font-weight:700;letter-spacing:.3px}
.tag{color:#7DD3C8;font-size:14px;margin:4px 0 26px}
.msg{font-size:17px;margin-bottom:12px}
.bar{width:420px;max-width:80vw;height:8px;background:#16325A;border-radius:6px;overflow:hidden}
.fill{height:100%%;width:3%%;background:linear-gradient(90deg,#14B8A6,#5EEAD4);border-radius:6px;transition:width .4s}
.spin{display:inline-block;width:14px;height:14px;border:2px solid #5EEAD4;border-right-color:transparent;border-radius:50%%;
 animation:r 0.9s linear infinite;vertical-align:-2px;margin-right:8px}
@keyframes r{to{transform:rotate(360deg)}}
.stage{font-size:13px;color:#9FB4CF;margin-top:10px;min-height:18px}
pre{font-size:11px;color:#6F86A6;text-align:left;width:560px;max-width:90vw;white-space:pre-wrap;margin-top:22px;min-height:70px}
.foot{position:fixed;bottom:12px;font-size:11px;color:#6F86A6;width:100%%;text-align:center}
.err{color:#FCA5A5}
</style></head><body><div class="wrap">
<img src="data:image/png;base64,%(logo)s" alt="Reticula">
<h1>Reticula</h1><div class="tag">Publication figures for framework materials · version %(ver)s</div>
<div class="msg" id="msg"><span class="spin"></span>Starting Reticula…</div>
<div class="bar"><div class="fill" id="fill"></div></div>
<div class="stage" id="stage">starting the Python runtime…</div>
<pre id="log"></pre></div>
<div class="foot">This page opens the app by itself as soon as it is ready (the first start after installing can take up to a minute).
Keep the small Reticula window open — closing it stops the app.</div>
<script>
var misses=0, t0=Date.now();
function tick(){
 fetch("/status?"+Date.now(),{cache:"no-store"}).then(function(r){return r.json()}).then(function(s){
  misses=0;
  document.getElementById("fill").style.width=Math.max(3,Math.min(100,s.pct))+"%%";
  document.getElementById("stage").textContent=s.stage+"  ("+Math.round(s.elapsed)+" s)";
  document.getElementById("log").textContent=(s.log||[]).join("\\n");
  if(s.error){document.getElementById("msg").innerHTML='<span class="err">'+s.error+'</span>';return;}
  if(s.ready&&s.url){document.getElementById("msg").textContent="Ready – opening Reticula…";location.replace(s.url);return;}
  setTimeout(tick,400);
 }).catch(function(){ if(++misses>15){document.getElementById("msg").innerHTML=
   '<span class="err">The Reticula launcher is not running any more. Start Reticula again from the Start menu.</span>';return;}
   setTimeout(tick,700);});
}
tick();
</script></body></html>"""


class Splash:
    """Tiny local HTTP server (own port, works offline) serving the loading page + /status from the launcher state."""
    def __init__(self, ver, log_path, port_start=8590):
        import base64, http.server, socketserver
        self.state = {"ready": False, "url": None, "error": None, "t0": time.time()}
        try:
            logo = base64.b64encode(open(os.path.join(APP_DIR, "assets", "reticula_icon.png"), "rb").read()).decode()
        except Exception:
            logo = ""
        page = (SPLASH_HTML % {"logo": logo, "ver": ver}).encode("utf-8")
        me = self

        class H(http.server.BaseHTTPRequestHandler):
            def do_GET(self):
                if self.path.startswith("/status"):
                    import json
                    body = json.dumps(me.status()).encode(); ct = "application/json"
                elif self.path in ("/", "/index.html") or self.path.startswith("/?"):
                    body, ct = page, "text/html; charset=utf-8"
                else:
                    self.send_response(404); self.end_headers(); return
                self.send_response(200); self.send_header("Content-Type", ct); self.send_header("Cache-Control", "no-store")
                self.send_header("Content-Length", str(len(body))); self.end_headers(); self.wfile.write(body)

            def log_message(self, *a):
                pass

        class S(socketserver.ThreadingMixIn, http.server.HTTPServer):
            daemon_threads = True
            allow_reuse_address = False
        self.log_path = log_path
        self.httpd = None
        for p in range(port_start, port_start + 40):
            if port_open(p):
                continue
            try:
                self.httpd = S((HOST, p), H); self.port = p; break
            except OSError:
                continue
        if self.httpd is None:
            self.url = None; return
        self.url = f"http://{HOST}:{self.port}/"
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()

    def status(self):
        pct, stage, lines = 3.0, "starting the Python runtime…", []
        try:
            with open(self.log_path, encoding="utf-8", errors="replace") as fh:
                lines = [l.rstrip() for l in fh.readlines()[-200:] if l.strip()]
        except Exception:
            pass
        import re as _re
        for l in lines:
            m = _re.search(r"\[startup\] (\d+)/(\d+) (.*)", l)
            if m:
                pct = 5 + 85 * int(m.group(1)) / max(1, int(m.group(2))); stage = m.group(3)
        if self.state["ready"]:
            pct, stage = 100, "ready"
        return {"pct": pct, "stage": stage, "log": [l[:160] for l in lines[-6:]], "ready": self.state["ready"],
                "url": self.state["url"], "error": self.state["error"], "elapsed": time.time() - self.state["t0"]}


def prewarm(log=print):
    """v2.11: import what the first page needs while the loading screen is shown, so the first browser connection does not
    sit on a white page (modules stay in sys.modules for the Streamlit script runs). Logs '[startup] i/N …' for the splash."""
    import importlib
    t = time.time()
    groups = [("numerical libraries (numpy, pandas)", ["numpy", "pandas"]),
              ("plotting (matplotlib, fonts)", ["matplotlib", "matplotlib.pyplot"]),
              ("scientific routines (scipy)", ["scipy.signal", "scipy.optimize", "scipy.ndimage", "scipy.interpolate"]),
              ("image libraries (Pillow)", ["PIL.Image"]),
              ("web framework (Streamlit)", ["streamlit", "streamlit.web.cli", "streamlit.components.v1"]),
              ("Reticula modules", ["mfs", "mfs.style", "mfs.ui", "mfs.readers", "mfs.analysis", "mfs.plots", "mfs.figctl",
                                    "mfs.figctl_ui", "mfs.pagestate", "mfs.chem", "mfs.defects", "mfs.composer", "mfs.report"])]
    for i, (label, mods) in enumerate(groups, start=1):
        log(f"[startup] {i}/{STEPS} loading {label}…")
        for m in mods:
            try:
                importlib.import_module(m)
                if m == "matplotlib":
                    importlib.import_module("matplotlib").use("Agg")
            except Exception as e:
                log(f"   ({m}: {type(e).__name__}: {e})")
    try:      # first-start font cache: build it now instead of during the first figure
        from matplotlib import font_manager as _fm
        _fm.findfont("Arial")
    except Exception:
        pass
    if os.environ.get("MFS_SLOW_START"):          # tests/screenshots only: keep the loading screen up a little longer
        time.sleep(float(os.environ["MFS_SLOW_START"]))
    log(f"[startup] {len(groups) + 1}/{STEPS} starting the web server… (libraries loaded in {time.time() - t:.1f} s)")

# ----------------------------------------------------------------------------- child: the Streamlit server
def watch_parent(ppid):
    """Exit when the launcher window process goes away (also if it is killed)."""
    if os.name != "nt":                      # POSIX: re-parented to init/launchd when the launcher dies
        while True:
            time.sleep(2)
            if os.getppid() != ppid:
                os._exit(0)
    try:
        import ctypes
        k = ctypes.windll.kernel32
        h = k.OpenProcess(0x00100000, False, ppid)  # SYNCHRONIZE
        if h:
            k.WaitForSingleObject(h, 0xFFFFFFFF)
            os._exit(0)
    except Exception:
        pass


def serve(port, ppid):
    os.makedirs(LOG_DIR, exist_ok=True)
    log = open(LOG, "w", encoding="utf-8", buffering=1)
    sys.stdout = sys.stderr = log
    if ppid:
        threading.Thread(target=watch_parent, args=(ppid,), daemon=True).start()
    fix_mimetypes()
    try:
        open(PIDFILE, "w", encoding="utf-8").write(f"{os.getpid()} {port} {app_version()}")
    except Exception:
        pass
    for k in ("HTTP_PROXY", "HTTPS_PROXY", "http_proxy", "https_proxy"):
        os.environ.pop(k, None)
    os.environ["NO_PROXY"] = os.environ["no_proxy"] = "127.0.0.1,localhost,::1"
    os.environ["MFS_EXPORT_DIR"] = EXPORT_DIR
    os.environ.setdefault("MFS_STORE", os.path.join(DATA_DIR, "sample_store.json"))
    os.makedirs(EXPORT_DIR, exist_ok=True)
    os.makedirs(DATA_DIR, exist_ok=True)
    os.chdir(APP_DIR)
    if APP_DIR not in sys.path:
        sys.path.insert(0, APP_DIR)
    if os.environ.get("MFS_NO_PREWARM") != "1":
        prewarm(lambda m: print(m, flush=True))
    from streamlit.web import cli as stcli
    sys.argv = ["streamlit", "run", APP,
                "--global.developmentMode=false",
                "--server.port", str(port), "--server.address", HOST, "--server.headless", "true",
                "--browser.serverAddress", HOST, "--browser.serverPort", str(port),
                "--server.enableCORS", "true", "--server.enableXsrfProtection", "true",
                "--server.fileWatcherType", "none", "--server.runOnSave", "false", "--server.maxUploadSize", "500",
                "--browser.gatherUsageStats", "false", "--theme.base", "light", "--client.toolbarMode", "minimal"]
    import mimetypes
    print(f"{APP_NAME} {app_version()}: starting on {HOST}:{port}; .js -> {mimetypes.guess_type('a.js')[0]}", flush=True)
    sys.exit(stcli.main())


# ----------------------------------------------------------------------------- parent: control window
def project_arg(argv):
    """v2.11: a .reticula project passed on the command line (double-click / 'Open with') -> absolute path or None."""
    for a in argv[1:]:
        if not a.startswith("--") and a.lower().endswith((".reticula", ".reticula", ".zip")) and os.path.isfile(a):
            return os.path.abspath(a)
    return None


def main():
    os.makedirs(LOG_DIR, exist_ok=True)
    prj = project_arg(sys.argv)
    if prj:                                    # the server process inherits it; the app opens it once
        os.environ["MFS_OPEN_PROJECT"] = prj
    stop_old_instances()
    time.sleep(0.5)
    port = free_port()
    ver = app_version()
    base_url = f"http://{HOST}:{port}/"

    def fresh_url():          # version + start time: an old cached index.html/bundle is never reused
        return f"{base_url}?v={ver}-{int(time.time())}"
    url = fresh_url()
    cmd = [sys.executable, "--serve", str(port), str(os.getpid())] if FROZEN else \
          [sys.executable, os.path.abspath(__file__), "--serve", str(port), str(os.getpid())]
    flags = 0x08000000 if os.name == "nt" else 0  # CREATE_NO_WINDOW
    try:                      # the loading screen reads this log: never show the previous start's lines
        open(LOG, "w", encoding="utf-8").close()
    except Exception:
        pass
    child = subprocess.Popen(cmd, creationflags=flags, stdin=subprocess.DEVNULL,
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    state = {"ready": False, "opened": False, "fmsg": "", "t_ready": None}
    t0 = time.time()
    splash = Splash(ver, LOG)
    if os.environ.get("MFS_SPLASH_FILE"):          # tests: where the loading page is
        try:
            open(os.environ["MFS_SPLASH_FILE"], "w", encoding="utf-8").write(f"{splash.url}\n{base_url}\n")
        except Exception:
            pass
    no_gui = os.environ.get("MFS_NO_GUI") == "1"
    no_browser = os.environ.get("MFS_NO_BROWSER") == "1"
    if not no_browser and splash.url:               # v2.11: the loading screen appears at once (no white page)
        state["opened"] = True
        webbrowser.open(splash.url)

    def monitor():                                   # health polling in the background; the splash + window read `state`
        while not state["ready"]:
            if child.poll() is not None:
                splash.state["error"] = f"The server stopped unexpectedly (code {child.returncode}). Details: {LOG}"
                return
            if health(port):
                ok, state["fmsg"] = frontend_check(port)
                state["t_ready"] = time.time() - t0
                splash.state["url"] = fresh_url()
                splash.state["ready"] = state["ready"] = True
                print(f"{APP_NAME}: health ok after {state['t_ready']:.1f} s, {state['fmsg']}", flush=True)
                if not no_browser and not state["opened"]:
                    state["opened"] = True
                    webbrowser.open(fresh_url())
                return
            time.sleep(0.4)
    threading.Thread(target=monitor, daemon=True).start()

    def stop():
        try:
            child.terminate()
            child.wait(5)
        except Exception:
            try:
                child.kill()
            except Exception:
                pass

    try:
        if no_gui:
            raise ImportError
        import tkinter as tk
    except Exception:  # headless fallback: wait for the server, wait for the child
        try:
            child.wait()
        except KeyboardInterrupt:
            stop()
        return

    root = tk.Tk()
    root.title(APP_NAME)
    try:
        if sys.platform == "darwin":
            png = os.path.join(APP_DIR, "assets", "reticula_icon.png")
            if os.path.isfile(png):
                root.iconphoto(True, tk.PhotoImage(file=png))
        else:
            root.iconbitmap(default=ICON)
    except Exception:
        pass
    root.resizable(False, False)
    root.configure(padx=18, pady=14)
    root.configure(bg="#0B1F3A")
    try:
        logo = tk.PhotoImage(file=os.path.join(APP_DIR, "assets", "reticula_icon.png")).subsample(4, 4)
    except Exception:
        logo = None
    head = tk.Frame(root, bg="#0B1F3A"); head.pack(anchor="w", fill="x")
    if logo is not None:
        tk.Label(head, image=logo, bg="#0B1F3A").pack(side="left", padx=(0, 10))
        root._logo = logo
    tt = tk.Frame(head, bg="#0B1F3A"); tt.pack(side="left")
    tk.Label(tt, text=APP_NAME, font=("Segoe UI", 16, "bold"), fg="white", bg="#0B1F3A").pack(anchor="w")
    tk.Label(tt, text=TAGLINE, font=("Segoe UI", 9), fg="#7DD3C8", bg="#0B1F3A").pack(anchor="w")
    tk.Label(tt, text=f"Version {app_version()}", font=("Segoe UI", 8), fg="#E6EEF8", bg="#0B1F3A").pack(anchor="w")
    tk.Label(tt, text="Built by Stijn Germonpr\u00e9", font=("Segoe UI", 8), fg="#9FB4CF", bg="#0B1F3A").pack(anchor="w")
    msg = tk.StringVar(value="Starting the figure server…\nThe first start can take up to a minute.")
    tk.Label(root, textvariable=msg, justify="left", font=("Segoe UI", 10), fg="#E6EEF8", bg="#0B1F3A").pack(anchor="w", pady=(10, 10))
    row = tk.Frame(root, bg="#0B1F3A"); row.pack(fill="x")
    b_open = tk.Button(row, text="Open in browser", width=16, state="disabled", bg="#0F766E", fg="white",
                       activebackground="#115E59", activeforeground="white", relief="flat", command=lambda: webbrowser.open(fresh_url()))
    b_open.pack(side="left")
    b_reload = tk.Button(row, text="Reload app", width=11, state="disabled", command=lambda: webbrowser.open(fresh_url()))
    b_reload.pack(side="left", padx=(6, 0))
    tk.Button(row, text="Open exports folder", width=18,
              command=lambda: (os.makedirs(EXPORT_DIR, exist_ok=True), open_path(EXPORT_DIR))).pack(side="left", padx=6)
    tk.Button(row, text="Quit", width=10, command=lambda: (stop(), root.destroy())).pack(side="right")
    tk.Label(root, text="Keep this window open while you work. Closing it stops the app.",
             font=("Segoe UI", 8), fg="#9FB4CF", bg="#0B1F3A").pack(anchor="w", pady=(10, 0))
    status = tk.StringVar(value="status: starting server…")
    tk.Label(root, textvariable=status, font=("Segoe UI", 8), fg="#7DD3C8", bg="#0B1F3A", justify="left").pack(anchor="w")
    lrow = tk.Frame(root, bg="#0B1F3A"); lrow.pack(fill="x", pady=(2, 0))
    tk.Label(lrow, text=f"Log: {LOG}", font=("Segoe UI", 8), fg="#9FB4CF", bg="#0B1F3A").pack(side="left")
    tk.Button(lrow, text="Open log", width=9, command=lambda: open_path(LOG) if os.path.exists(LOG) else None).pack(side="right")
    try:
        from tkinter import ttk
        sty = ttk.Style(root)
        try:
            sty.theme_use("clam")
        except Exception:
            pass
        sty.configure("R.Horizontal.TProgressbar", troughcolor="#16325A", background="#14B8A6", bordercolor="#0B1F3A",
                      lightcolor="#14B8A6", darkcolor="#14B8A6")
        pbar = ttk.Progressbar(root, orient="horizontal", length=420, mode="determinate", maximum=100, style="R.Horizontal.TProgressbar")
        pbar.pack(anchor="w", pady=(8, 0))
    except Exception:
        pbar = None
    stage = tk.StringVar(value="Loading… starting the Python runtime")
    tk.Label(root, textvariable=stage, font=("Segoe UI", 8), fg="#E6EEF8", bg="#0B1F3A").pack(anchor="w")

    def poll():
        if child.poll() is not None and not state["ready"]:
            msg.set(f"The server stopped unexpectedly (code {child.returncode}).\nDetails: {LOG}")
            b_open.config(state="disabled")
            return
        if state["ready"]:
            if b_open["state"] == "disabled":
                status.set(f"status: server healthy ({state['t_ready']:.0f} s); {state['fmsg']}")
                msg.set(f"Running at {base_url}\nFigures saved with 'Save to folder' go to:\n{EXPORT_DIR}")
                b_open.config(state="normal"); b_reload.config(state="normal")
                stage.set("Ready ✓  (the app is open in your browser)")
                if pbar is not None:
                    pbar["value"] = 100
            return
        st_ = splash.status()
        if pbar is not None:
            pbar["value"] = st_["pct"]
        stage.set(f"Loading… {st_['stage']}")
        status.set(f"status: waiting for the server on {HOST}:{port} ({time.time() - t0:.0f} s)"
                   + ("" if port_open(port) else " – port not open yet"))
        root.after(400, poll)

    root.protocol("WM_DELETE_WINDOW", lambda: (stop(), root.destroy()))
    root.after(300, poll)
    try:
        root.mainloop()
    finally:
        stop()


def selftest(out):
    """Import every reader backend the app loads lazily and run the CIF simulation (support/diagnostics)."""
    import importlib, traceback, json
    res = {}
    for m in ["numpy", "scipy.optimize", "scipy.signal", "scipy.sparse.linalg", "scipy.spatial", "pandas", "pyarrow",
              "matplotlib.backends.backend_pdf", "matplotlib.backends.backend_svg", "PIL.Image", "PIL.ImageDraw", "PIL.ImageFont", "openpyxl", "xlrd",
              "tifffile", "nmrglue", "brukeropusreader", "ncempy.io", "h5py", "hdf5plugin", "pymatgen.core",
              "pymatgen.analysis.diffraction.xrd", "streamlit", "tkinter", "pytesseract", "streamlit_image_coordinates"]:
        try:
            importlib.import_module(m); res[m] = "ok"
        except Exception as e:
            res[m] = f"FAIL {type(e).__name__}: {e}"
    try:
        sys.path.insert(0, APP_DIR)
        from mfs import analysis as A, porosity
        df, sticks, formula = A.simulate_pxrd(os.path.join(APP_DIR, "references", "UiO-66_Zr_COD4512074.cif"), substitute={"Zr": "Hf"})
        res["simulate_pxrd"] = f"ok {formula} {len(sticks)} reflections"
    except Exception:
        res["simulate_pxrd"] = "FAIL " + traceback.format_exc()[-400:]
    try:      # SEM calibration: bundled Tesseract OCR on the built-in example (no metadata)
        from mfs import semcal, readers as R
        ex = os.path.join(APP_DIR, "example_data", "SEM_calibration", "JEOL-style_noMetadata_x25000.tif")
        d = R.read_sem(open(ex, "rb"), ex)
        c = semcal.calibrate(d["image"])
        res["sem_ocr"] = (f"ok tesseract={semcal.tesseract_cmd()} bar={c['bar'] and c['bar']['length_px']} px "
                          f"label={c['label'] and c['label']['shown']} mag={c['mag_ocr'] and c['mag_ocr']['value']} "
                          f"nm_per_px={c['nm_per_px']}" if c["label"] and c["mag_ocr"] else f"FAIL {c}")
    except Exception:
        res["sem_ocr"] = "FAIL " + traceback.format_exc()[-400:]
    try:      # DSC: bundled Mettler STARe example -> melting peak
        from mfs import dsc
        ex = os.path.join(APP_DIR, "example_data", "DSC", "MOF-PEG_Mettler-STARe_EXAMPLE.txt")
        d = dsc.read_dsc(ex)
        c = dict(d["curves"][0], id=0)
        r = [x for x in dsc.evaluate_all(c, "linear", d["meta"]["mass_mg"]) if x["kind"] == "endo"][0]
        res["dsc"] = (f"ok peak={r['T_peak']:.1f} onset={r['T_onset']:.1f} dH={r['dH']:.1f} J/g"
                      if abs(r["T_peak"] - 120) < 0.5 and abs(r["dH"] + 80) < 2 else f"FAIL {r['T_peak']} {r['dH']}")
    except Exception:
        res["dsc"] = "FAIL " + traceback.format_exc()[-400:]
    try:      # Batch A: DFT-kernel PSD on the bundled N2 example, kinetics CI, SI table (DOCX/XLSX)
        from mfs import psd, kinetics as K, si_export as SI, samplestore as SS, readers as R
        ex = os.path.join(APP_DIR, "example_data", "Sorption", "MOF801_3Flex_MicroActive.xlsx")
        iso = R.read_sorption(ex, ex)["iso"]; a_ = iso[iso.branch.astype(str).str.startswith("ads")]
        f = psd.dft_fit(a_.p_rel.values.astype(float), a_.q.values.astype(float))
        st_ = {}; SS.put(st_, "x", "defects", {"L": 5.0, "x": 1.0, "OS": 1.6})
        dx = SI.to_docx(SI.defect_rows(st_.get("sample_store", st_))); xx = SI.to_xlsx(SI.defect_rows(st_.get("sample_store", st_)))
        k = K.arrhenius([40, 50, 60], [1e-5, 2.3e-5, 5e-5])
        res["batchA"] = (f"ok dft_peak={f['w'][f['dVdw'].argmax()]:.2f} nm resid={f['residual_rel']:.3f} docx={len(dx)} xlsx={len(xx)} "
                         f"Ea={k['Ea_kJ']:.1f}") if dx[:2] == b"PK" and xx[:2] == b"PK" else "FAIL docx/xlsx"
    except Exception:
        res["batchA"] = "FAIL " + traceback.format_exc()[-400:]
    try:      # Batch B: XPS doublet fit, TGA steps / MS, DSC cp (sapphire), FTIR inset
        from mfs import xps as XP, tga as TG, dsc as DSC2, spectra as SP, readers as R, analysis as A
        E = os.path.join(APP_DIR, "example_data")
        d = R.read_xps(os.path.join(E, "XPS", "ZrO2@Sulfo-MOF801_S2p_CasaXPS.csv"), "S2p.csv")
        fx = XP.fit(d.x.values, d.y.values, [{"name": "S", "region": "S 2p", "center": 168.5}], "shirley")
        t = R.read_tga(os.path.join(E, "TGA", "MOF801_Netzsch.txt"), "MOF801_Netzsch.txt"); st_ = TG.steps(A.tga_analysis(t["T"], t["w"], 21), 1.0)
        ms = TG.read_ms(os.path.join(E, "TGA", "MOF801_QMS_ion_currents_EXAMPLE.csv"))
        P_ = os.path.join(E, "DSC", "Cp_sapphire"); rd = lambda f: DSC2.read_dsc(os.path.join(P_, f))["curves"][0]
        cp = DSC2.cp_from_sapphire(rd("Cp_sample_MOF-801_EXAMPLE.csv"), rd("Cp_sapphire_EXAMPLE.csv"), rd("Cp_blank_EXAMPLE.csv"), 9.80, 25.30)
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(); ax.plot([1, 2], [1, 2]); SP.add_inset(ax, lambda a: a.plot([1, 2], [1, 2]), (1.2, 1.5)); plt.close(fig)
        cp100 = float(__import__("numpy").interp(100, cp["T"], cp["Cp"]))
        ok = abs(fx["rows"][0]["area (total)"] - 2250) < 100 and len(st_) == 3 and "m/z 44" in ms and abs(cp100 - 0.74) < 0.02
        res["batchB"] = (f"ok S2p area={fx['rows'][0]['area (total)']:.0f} steps={len(st_)} ms={len(ms)} cp100={cp100:.3f}" if ok else
                         f"FAIL S2p={fx['rows'][0]['area (total)']:.0f} steps={len(st_)} cp100={cp100:.3f}")
    except Exception:
        res["batchB"] = "FAIL " + traceback.format_exc()[-400:]
    try:      # Batch C: particle sizing on the bundled SEM, EDX overlay + line scan, ratios with propagated σ
        import numpy as np
        from mfs import particles as PA, edx as EDX, ratios as RT, readers as R
        E = os.path.join(APP_DIR, "example_data")
        f = os.path.join(E, "SEM", "MOF-801_SEM.tif"); m = R.read_sem(open(f, "rb"), f)
        lab, fa, _t = PA.segment(m["image"][:960], dark_particles=False); rows = PA.measure(lab, 2.6, fa)
        fit = PA.lognormal_fit([r_["ECD (nm)"] for r_ in rows])
        from PIL import Image
        mp = {el: np.array(Image.open(os.path.join(E, "EDX", f"ZrO2_at_Sulfo-MOF801_EDXmap_{el}.png")).convert("L")) for el in ("Zr", "Hf", "S")}
        rgb = EDX.overlay(mp, {"Zr": "#E69F00", "Hf": "#56B4E9", "S": "#F0E442"}); d_, pr = EDX.line_scan(mp["Hf"], (10, 128), (300, 128))
        r = RT.zr_hf_s({"Zr": 9.337, "Hf": 45.962, "S": 0.773}, {"Zr": 0.47, "Hf": 2.3, "S": 0.04})
        ok = 10 <= fit["n"] <= 60 and 150 < fit["mean"] < 600 and rgb.shape == (256, 320, 3) and len(pr) > 200 and 0.3 < r["S per M6 node"][0] < 0.5
        res["batchC"] = (f"ok particles n={fit['n']} ECD={fit['mean']:.0f}±{fit['sd']:.0f} nm S/M6={r['S per M6 node'][0]:.3f}±{r['S per M6 node'][1]:.3f}"
                         if ok else f"FAIL n={fit['n']} mean={fit['mean']:.0f} S/M6={r['S per M6 node'][0]:.3f}")
    except Exception:
        res["batchC"] = "FAIL " + traceback.format_exc()[-400:]
    try:      # Batch D: batch processing, SI bundle, SQLite sample DB from the bundled pipeline overview, journal template
        import tempfile, zipfile, io as _io
        from mfs import batchproc as BP, si_bundle as SB, sampledb as DB, templates as TP, composer as C
        E = os.path.join(APP_DIR, "example_data")
        fs = BP.list_folder(os.path.join(E, "PXRD"), "PXRD"); z, log = BP.run("PXRD", [(f, f) for f in fs], formats=("pdf",), dpi=100)
        con = DB.connect(os.path.join(tempfile.mkdtemp(), "t.db"))
        n = DB.import_master(con, os.path.join(APP_DIR, "example_data", "SampleDB", "Sample_overview_EXAMPLE.xlsx"))
        si = SB.build({"x": {"page": "PXRD", "pdf": b"%PDF", "csv": b"a"}}, {}, "selftest")
        fig, _ = C.compose([lambda ax: ax.plot([0, 1])], ncols=1, width=TP.width("ACS", "single"), font_size=7)
        chk = TP.check(fig, "ACS", "single")
        ok = all(l_["status"] == "ok" for l_ in log) and n >= 8 and "SI_parameters.xlsx" in zipfile.ZipFile(_io.BytesIO(si)).namelist() and chk[0][0] == "ok"
        res["batchD"] = f"ok batch={len(log)} files db={n} samples si+template ok" if ok else f"FAIL batch={log} db={n} chk={chk}"
    except Exception:
        res["batchD"] = "FAIL " + traceback.format_exc()[-400:]
    try:      # Batch E: Structure (CIF) on the bundled UiO-66 / MOF-801 CIFs: parse, grow, ORTEP draw, SI tables, PXRD, pores
        import io as _io2
        from mfs import cif as CX, cifplot as CP, cifsi as CS, cifpxrd as CPX, poretools as PT
        import matplotlib.pyplot as _plt
        R_ = os.path.join(APP_DIR, "references")
        S8 = CX.read(os.path.join(R_, "MOF-801_Zr_COD4121459.cif"))[0]
        sites, bonds = S8.grow("cell", part="major")
        f_, ax_ = _plt.subplots(); CP.draw(ax_, S8, sites, bonds, style="ortep", polyhedra=True, view=CP.view_matrix(S8, "best", sites=sites[:60]))
        bb = _io2.BytesIO(); f_.savefig(bb, format="png"); _plt.close(f_)
        g = CX.geom_tables(S8); dx = CS.docx(CS.si_sections(S8, g, CX.hbonds(S8)))
        sim, stx = CPX.simulate(os.path.join(R_, "MOF-801_Zr_COD4121459.cif"), (3, 30))
        pr = PT.analyse(CX.read(os.path.join(R_, "UiO-66_Zr_COD4512074.cif"))[0], 0.6, n_mc=40)
        ok = len(sites) > 200 and bb.tell() > 5000 and len(dx) > 3000 and len(stx) > 10 and 6 < pr["Di (Å)"] < 10
        res["batchE"] = (f"ok MOF-801 {len(sites)} atoms, {len(stx)} hkl; UiO-66 Di={pr['Di (Å)']:.1f} Df={pr['Df (Å)']:.1f} A"
                         if ok else f"FAIL sites={len(sites)} png={bb.tell()} docx={len(dx)} hkl={len(stx)} Di={pr['Di (Å)']}")
    except Exception:
        res["batchE"] = "FAIL " + traceback.format_exc()[-400:]
    try:      # v2.8 particle watershed + panel annotations inside the panels
        import numpy as np
        from mfs import particles as PA2, plots as PL
        yy, xx = np.mgrid[:200, :300]; im_ = np.full((200, 300), 200, np.uint8)
        for cx in (90, 150, 210):
            im_[(xx - cx) ** 2 + (yy - 100) ** 2 < 32 ** 2] = 40
        lab_, f_, t_, inf_ = PA2.segment_ws(im_, True, exclude_holes=False)
        n_ = len([i for i in np.unique(lab_) if i])
        fig_ = PL.fig_sem_panel([{"image": im_, "pixel_nm": 5.0, "mag_text": "x1000"}] * 4, 4, 3.25)
        fs_ = fig_._mfs_panel_fs
        res["particles_panels"] = f"ok watershed {n_}/3 discs, panel font {fs_:.1f} pt" if n_ == 3 and fs_ else f"FAIL n={n_} fs={fs_}"
    except Exception:
        res["particles_panels"] = "FAIL " + traceback.format_exc()[-400:]
    try:                                       # v2.9: UV/Vis, Rigaku .rasx, picker logic, page-state module, credits
        import io as _io, zipfile as _zf, glob as _gl
        import numpy as np
        from mfs import uvvis as UV, rigaku as RG, particle_picker as PKr, pagestate as PSr, AUTHOR
        exd = os.path.join(APP_DIR, "example_data", "UVVis")
        uf = UV.read_any(open(os.path.join(exd, "NanoPhotometer_BCA_EXAMPLE.xlsx"), "rb").read(), "bca.xlsx")
        drs = UV.read_any(open(os.path.join(exd, "DRS_reflectance_EXAMPLE.csv"), "rb").read(), "drs.csv")
        E_, T_ = UV.tauc(drs.spectra[0].x, drs.spectra[0].y, "direct", True); eg = UV.band_gap(E_, T_)["Eg"]
        b_ = _io.BytesIO()
        with _zf.ZipFile(b_, "w") as z:
            z.writestr("Data0/Profile0.txt", "\n".join(f"{5 + i * 0.02:.2f}\t{100 + i}\t2.0" for i in range(200)))
            z.writestr("Data0/MesurementConditions0.xml", "<M><ScanInformation><AxisName>TwoThetaTheta</AxisName></ScanInformation>"
                       "<XrayGenerator><TargetName>Cu</TargetName><KAlpha1>1.540593</KAlpha1></XrayGenerator></M>")
        sc = RG.read_any(b_.getvalue(), "t.rasx")[0]
        ok = (len(uf.spectra) == 11 and uf.instrument.startswith("Implen") and 4.8 < eg < 5.4 and sc["is_2theta_theta"]
              and sc["y"][0] == 200.0 and PKr.toggle([1], 1) == [] and AUTHOR.endswith("pr\u00e9"))
        res["v29"] = (f"ok uvvis {len(uf.spectra)} spectra ({uf.instrument}), ZrO2 Eg={eg:.2f} eV, rasx {len(sc['x'])} pts "
                      f"I=counts*att, author={AUTHOR}" if ok else f"FAIL uv={len(uf.spectra)} eg={eg} rasx={sc['y'][:2]}")
    except Exception:
        res["v29"] = "FAIL " + traceback.format_exc()[-400:]
    try:                                       # v2.10: adaptive particle segmentation + learned shape prior
        import numpy as np
        from mfs import particles as PA3, shapeprior as SPr
        yy, xx = np.mgrid[:300, :480]
        bg = 225 - 40 * xx / 480; ab = np.zeros((300, 480))
        for k, (cx, cy) in enumerate([(60 + 70 * i, 70 + 80 * j) for i in range(6) for j in range(3)]):
            u, v = (xx - cx) * 0.94 + (yy - cy) * 0.34, -(xx - cx) * 0.34 + (yy - cy) * 0.94
            ab += (0.55 if k % 2 else 0.22) * ((np.abs(u) <= 13) & (np.abs(v) <= 13))
        im3 = np.clip(bg * np.exp(-ab) + np.random.default_rng(0).normal(0, 5, ab.shape), 0, 255).astype(np.uint8)
        l3, f3, t3, i3 = PA3.segment_adaptive(im3, shape_prior=True)
        n3 = sum(1 for k, v in i3.items() if isinstance(k, int) and v.get("kind") == "particle")
        pr3 = i3["_prior"]
        ok = n3 >= 16 and pr3["n"] >= 10 and pr3["corners"] == 4
        res["v210"] = (f"ok adaptive {n3}/18 squares (incl. light grey), prior n={pr3['n']} corners={pr3['corners']} "
                       f"w={pr3['weight']:.2f}" if ok else f"FAIL n={n3} prior={pr3.get('n')} corners={pr3.get('corners')}")
    except Exception:
        res["v210"] = "FAIL " + traceback.format_exc()[-400:]
    try:                                       # v2.11: any-node defect calculator, project files, update check, global style
        import io as _io4, zipfile as _zf4, urllib.request, ssl  # noqa: F401  (the update check needs HTTPS in the frozen app)
        from mfs import defects as D4, nodes as N4, chem as CH4, project as PJ4, update as UP4, style as ST4
        from mfs.chem import molar_mass as mm4
        d4 = CH4.FRAMEWORKS["HKUST-1"]; lk = d4["linker"]; f4, z4, _r4, nH4, _ = CH4.SPECIES[lk]
        inp4 = D4.DefectInput(sample="HKUST-1", framework="HKUST-1", metal=d4.get("metal", "Zr"),
                              species=[D4.Species(lk, f4, z4, "linker", nH4, 1.0 * nH4, False, True)],
                              tga_W_plateau=100 * mm4("Cu2(C9H3O6)1.3333333333") / mm4("Cu2O2"), tga_W_end=100.0)
        r4 = D4.method_a(inp4)
        rel = {"tag_name": "v9.9.0", "html_url": "x", "assets": [{"name": "Reticula_Setup_9.9.0.exe", "browser_download_url": "y"}]}
        a4 = UP4.pick_asset(rel, "Windows")
        pj = PJ4.export_bytes("Home", "selftest", {"_ps_pages": {}, "_ps_wp": {}, "_ps_up": {}, "gs_font": 7.0, "gs_reset": False})
        nm4 = _zf4.ZipFile(_io4.BytesIO(pj)).namelist()
        ok = (abs(r4.L - 4 / 3) < 1e-3 and len(N4.NODES) >= 9 and a4 and a4["name"].endswith(".exe") and UP4.is_newer("2.11.1", "2.11.0")
              and "manifest.json" in nm4 and len(ST4.JOURNAL_PRESETS) >= 3)
        res["v211"] = (f"ok HKUST-1 L={r4.L:.3f}/node, {len(N4.NODES)} nodes, update asset {a4['name']}, project {len(nm4)} files, "
                       f"{len(ST4.JOURNAL_PRESETS)} journal presets" if ok else f"FAIL L={r4.L} asset={a4} files={nm4}")
    except Exception:
        res["v211"] = "FAIL " + traceback.format_exc()[-400:]
    try:                                       # v2.12: sphere model (touching pair split by circle fits), label mover, macOS asset pick
        import numpy as _np5
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as _plt5
        from scipy import ndimage as _ndi5
        from mfs import particles as PA5, spheres as SP5, labelmove as LM5, update as UP5  # noqa: F401
        yy, xx = _np5.mgrid[:220, :320]; t5 = _np5.zeros((220, 320))
        for cx, cy in [(70, 70), (180, 110), (234, 110), (90, 170)]:
            t5 += 2 * _np5.sqrt(_np5.clip(28.0 ** 2 - (xx - cx) ** 2 - (yy - cy) ** 2, 0, None))
        im5 = _np5.clip(_ndi5.gaussian_filter(220 * _np5.exp(-t5 / 25.0), 1.2), 0, 255).astype(_np5.uint8)
        _l5, _f5, _x5, i5 = PA5.segment_adaptive(im5, sphere=True)
        nc5 = sum(1 for k, v in i5.items() if isinstance(v, dict) and v.get("circle") and v.get("kind") == "particle")
        fig5, ax5 = _plt5.subplots(); ax5.plot([0, 1], [0, 1]); ax5.text(0.5, 0.5, "S1")
        LM5.tag_originals(fig5); it5 = LM5.items(fig5)
        st5 = LM5.handle(fig5, {}, {"op": "move", "id": it5[0][0], "fx": 0.2, "fy": 0.8})
        _plt5.close(fig5)
        rel5 = {"tag_name": "v9.9.0", "html_url": "x", "assets": [{"name": "Reticula_9.9.0_macOS_x86_64.dmg", "browser_download_url": "a"},
                                                                {"name": "Reticula_9.9.0_macOS_arm64.dmg", "browser_download_url": "b"}]}
        a5 = UP5.pick_asset(rel5, "Darwin", "arm64")
        ok = nc5 == 4 and i5["_sphere"]["n_split"] >= 2 and len(st5) == 1 and a5 and "arm64" in a5["name"]
        res["v212"] = (f"ok sphere model {nc5}/4 circles ({i5['_sphere']['n_split']} split), label move stored, mac asset {a5['name']}"
                       if ok else f"FAIL circles={nc5} labels={st5} asset={a5}")
    except Exception:
        res["v212"] = "FAIL " + traceback.format_exc()[-400:]
    try:                                       # v2.13: page modules bind, BET window scan, 1H NMR prediction + assignment (RDKit), generic ratios
        import numpy as _np6
        from mfs import pages as PG6, analysis as A6, ratios as RT6, nmrpredict as HP6, nmrassign as NA6
        mods6 = PG6.modules()                                   # every page module imports (also inside the frozen app)
        ns6 = {}
        for m6 in mods6:
            for n6 in m6.INJECTED:
                ns6.setdefault(n6, object())
        PG6.bind(ns6)
        p6 = _np6.concatenate([_np6.logspace(-6, -2, 30), _np6.linspace(0.011, 0.3, 40)])
        q6 = 100 * 800 * p6 / ((1 - p6) * (1 + 799 * p6))
        b6 = A6.bet_rouquerol(p6, q6); sp6 = A6.bet_range_spread(b6)
        r6 = RT6.molar_ratios({"Zr": 20.0, "Hf": 10.0, "S": 2.0}, {}, ["Zr", "Hf"], per=6)
        pred6 = {x["label"]: x["shift"] for x in HP6.predict(HP6.mol_from_any("CCOC(C)=O"))}     # ethyl acetate 1.26 / 4.12 / 2.05
        g6 = NA6.groups(HP6.mol_from_any("Cc1ccc(C)cc1"))
        ok = (len(mods6) == len(PG6.MODULES) and abs(b6["n_m"] - 100) < 0.5 and sp6["spread_pct"] < 1 and "S per 6 (Zr+Hf)" in r6
              and abs(pred6[1] - 1.26) < 0.2 and abs(pred6[2] - 4.12) < 0.2 and sorted(g["nH"] for g in g6) == [4, 6])
        res["v213"] = (f"ok {len(mods6)} page modules, BET n_m={b6['n_m']:.1f} spread {sp6['spread_pct']:.2f} %, 1H prediction EtOAc "
                       f"{pred6[1]:.2f}/{pred6[2]:.2f}/{pred6[5]:.2f} ppm") if ok else f"FAIL pages={len(mods6)} bet={b6['n_m']} pred={pred6}"
    except Exception:
        res["v213"] = "FAIL " + traceback.format_exc()[-400:]
    try:                                       # v2.14: parallel particle segmentation inside the frozen app
        import numpy as _np7
        from mfs import segpool as SP7, segjobs as SJ7
        rng7 = _np7.random.default_rng(3)
        yy7, xx7 = _np7.mgrid[0:300, 0:300]
        imgs7 = []
        for k7 in range(4):
            im7 = _np7.full((300, 300), 200.0)
            for cx7, cy7 in rng7.integers(40, 260, (5, 2)):
                im7[(xx7 - cx7) ** 2 + (yy7 - cy7) ** 2 < 18 ** 2] = 70.0
            imgs7.append((im7 + rng7.normal(0, 2, im7.shape)).astype(_np7.uint8))
        a7 = (True, 1.0, 4.0, 0.0, 0.25, 0.8, 2.0, True, 1536, None, 0.7, 10.0, "auto", 0.5)
        jobs7 = [("adaptive", im7, a7) for im7 in imgs7]
        ser7 = [SJ7.run(*j) for j in jobs7]
        SP7.clear_cache(); last7 = dict(SP7.prefetch(jobs7))
        par7 = [SP7.get(*j) for j in jobs7]
        same7 = all(_np7.array_equal(x[0], y[0]) for x, y in zip(ser7, par7))
        res["v214_parallel"] = (f"ok {last7['mode']} ({last7['workers']} workers, {last7['seconds']:.1f} s), identical to serial" if same7
                                else "FAIL parallel result differs from serial") + (f"; fallback reason: {last7['error']}" if last7.get("error") else "")
    except Exception:
        res["v214_parallel"] = "FAIL " + traceback.format_exc()[-400:]
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()      # v2.14: the particle sizer starts worker processes (a frozen exe would otherwise start the whole app again)
    if len(sys.argv) >= 3 and sys.argv[1] == "--selftest":
        selftest(sys.argv[2])
    elif len(sys.argv) >= 3 and sys.argv[1] == "--serve":
        serve(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 0)
    else:
        main()
