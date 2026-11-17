"""v2.11 headless-Chrome checks: loading screen (launcher splash -> app), unified sample list (toggle / reorder),
project save/open, update banner against a local mock server."""
import os, re, sys, json, time, socket, subprocess, threading, urllib.request
import pytest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from test_sem_e2e_browser import _browser
from test_batchA_e2e_browser import _shot, _log

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAUNCHER = os.environ.get("MFS_LAUNCHER", "/workspace/exe_build/launcher/mfs_launcher.py")


@pytest.mark.skipif(not os.path.exists(LAUNCHER) or bool(os.environ.get("MFS_E2E_URL")), reason="launcher not available")
def test_v211_loading_screen_redirects_to_app(tmp_path):
    from playwright.sync_api import sync_playwright
    sf = tmp_path / "splash.txt"
    env = dict(os.environ, LOCALAPPDATA=str(tmp_path / "lad"), MFS_APP_DIR=ROOT, MFS_NO_GUI="1", MFS_NO_BROWSER="1",
               MFS_SPLASH_FILE=str(sf), MFS_EXPORT_DIR=str(tmp_path), MFS_STORE=str(tmp_path / "store.json"),
               MFS_DB=str(tmp_path / "r.db"), MFS_SLOW_START="4")
    proc = subprocess.Popen([sys.executable, LAUNCHER], env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                            start_new_session=True)
    try:
        for _ in range(100):
            if sf.exists() and sf.read_text().count("\n") >= 2:
                break
            time.sleep(0.1)
        splash, app = sf.read_text().split()[:2]
        with sync_playwright() as p:
            b = _browser(p); pg = b.new_page(viewport={"width": 1280, "height": 800})
            t0 = time.time(); pg.goto(splash)
            pg.get_by_text("Starting Reticula").first.wait_for(timeout=5000)
            assert pg.locator("img").count() >= 1
            pg.wait_for_timeout(1200); _shot(pg, "v211_loading_screen.png")
            stage = pg.locator("#stage").inner_text()
            pg.wait_for_url(app.rstrip("/") + "/**", timeout=120000)
            t_red = time.time() - t0
            pg.locator("[data-testid='stSidebar']").first.wait_for(timeout=120000)
            pg.get_by_text("Reticula").first.wait_for(timeout=60000)
            _log(f"loading screen: stage '{stage}', redirected to the app after {t_red:.1f} s")
            b.close()
    finally:
        import signal
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except Exception:
            proc.terminate()


@pytest.fixture(scope="module")
def gurl(tmp_path_factory):
    if os.environ.get("MFS_E2E_URL"):
        yield os.environ["MFS_E2E_URL"]; return
    s = socket.socket(); s.bind(("127.0.0.1", 0)); port = s.getsockname()[1]; s.close()
    tmp = tmp_path_factory.mktemp("v211srv")
    env = dict(os.environ, MFS_NO_BROWSER="1", MFS_STORE=str(tmp / "store.json"), MFS_DB=str(tmp / "r.db"), MFS_EXPORT_DIR=str(tmp),
               MFS_SETTINGS=str(tmp / "settings.json"))
    proc = subprocess.Popen([sys.executable, "-m", "streamlit", "run", os.path.join(ROOT, "app.py"), "--server.port", str(port),
                             "--server.headless", "true"], cwd=ROOT, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(120):
        try:
            urllib.request.urlopen(f"http://127.0.0.1:{port}/_stcore/health", timeout=2); break
        except Exception:
            time.sleep(1)
    yield f"http://127.0.0.1:{port}"
    proc.terminate(); proc.wait(10)


def _order_caption(pg):
    loc = pg.locator("[data-testid='stMain']").get_by_text("Plot order (").first
    loc.wait_for(timeout=60000)
    return loc.inner_text()


def _wait_caption_change(pg, old, timeout=30):
    t = time.time()
    while time.time() - t < timeout:
        c = _order_caption(pg)
        if c != old:
            return c
        pg.wait_for_timeout(300)
    return _order_caption(pg)


def test_v211_sample_list_toggle_reorder_rename(gurl):
    from playwright.sync_api import sync_playwright
    from _e2e_util import nav
    with sync_playwright() as p:
        b = _browser(p); pg = b.new_page(viewport={"width": 1400, "height": 1000})
        pg.goto(gurl); pg.locator("[data-testid='stSidebar']").first.wait_for(timeout=120000)
        nav(pg, "PXRD")
        c0 = _order_caption(pg)
        n_all = int(c0.split(" of ")[1].split(" ")[0]); n0 = int(c0.split("(")[1].split(" ")[0])
        assert n_all >= 3 and n0 == n_all, c0
        fr = pg.frame_locator("iframe[title*='sample_list']").first
        fr.locator(".row").first.wait_for(timeout=30000)
        from mfs.samplelist import plain; first_name = plain(fr.locator("input.nm").nth(0).input_value())
        # hide the 2nd sample
        fr.locator("input[type=checkbox]").nth(1).click()
        c1 = _wait_caption_change(pg, c0)
        assert f"{n_all - 1} of {n_all} shown" in c1, c1
        # shift-click range: hide 1st..3rd (click 1st unchecks it, shift-click 3rd unchecks range)
        fr.locator("input[type=checkbox]").nth(0).click(); pg.wait_for_timeout(150)
        fr.locator("input[type=checkbox]").nth(2).click(modifiers=["Shift"])
        for _ in range(60):                                   # wait for the FINAL state (the first click already changes the caption)
            c2 = _order_caption(pg)
            if f"{n_all - 3} of {n_all} shown" in c2:
                break
            pg.wait_for_timeout(250)
        assert f"{n_all - 3} of {n_all} shown" in c2, c2
        fr.get_by_role("button", name="All").click()
        c3 = _wait_caption_change(pg, c2)
        assert f"{n_all} of {n_all} shown" in c3
        # reorder: drag the first row below the third
        fr.locator(".row .h").nth(0).drag_to(fr.locator(".row").nth(2))
        c4 = _wait_caption_change(pg, c3)
        order = c4.split(": ", 1)[1].split(" → ")
        assert order[0] != first_name and first_name in order, (first_name, order)
        _shot(pg, "v211_sample_list.png")
        # rename
        nm = fr.locator("input.nm").nth(0); nm.fill("Renamed sample X"); nm.press("Enter")
        c5 = _wait_caption_change(pg, c4)
        assert c5.split(": ", 1)[1].startswith("Renamed sample X"), c5
        # stack evenly button exists and works
        pg.locator("[data-testid='stMain']").get_by_role("button", name="⇅ Stack evenly").first.click()
        pg.wait_for_timeout(2500)
        assert not pg.locator("[data-testid='stException']").count()
        # survives navigation (page state)
        nav(pg, "Home", answer="save"); pg.wait_for_timeout(1500); nav(pg, "PXRD", answer=None)
        c6 = _order_caption(pg)
        assert c6.split(": ", 1)[1].startswith("Renamed sample X"), c6
        _log(f"sample list: {n_all} samples; toggle, shift-range, All, drag reorder, rename, stack evenly, kept after navigation")
        b.close()


def _expand_sidebar(pg, title):
    ex = pg.locator("[data-testid='stSidebar'] [data-testid='stExpander']").filter(has_text=title).first
    ex.scroll_into_view_if_needed()
    if ex.locator("details").first.get_attribute("open") is None:
        ex.locator("summary").first.click(); pg.wait_for_timeout(800)
    return ex


def test_v211_node_library_undo_style_project(gurl, tmp_path):
    from playwright.sync_api import sync_playwright
    from _e2e_util import nav, bright_particles
    from test_sem_e2e_browser import _select
    import re
    with sync_playwright() as p:
        b = _browser(p); ctx = b.new_context(accept_downloads=True, viewport={"width": 1500, "height": 1100}); pg = ctx.new_page()
        pg.goto(gurl); pg.locator("[data-testid='stSidebar']").first.wait_for(timeout=120000)
        main = lambda: pg.locator("[data-testid='stMain']").inner_text()
        # --- node library (defect calculator, any metal)
        nav(pg, "Defect calculator"); pg.get_by_text("Framework, node and metal").first.wait_for(timeout=120000)
        _select(pg, "Framework", "HKUST-1")
        for _ in range(40):
            if "W residue (CuO)" in main():
                break
            pg.wait_for_timeout(500)
        assert "W residue (CuO)" in main()
        assert pg.locator("[data-testid='stSelectbox']").filter(has_text="Node (SBU)").first.locator("input").first.input_value() == "Cu2 paddlewheel"
        ex = pg.locator("[data-testid='stExpander']").filter(has_text="Node library").first
        ex.locator("summary").first.click(); pg.wait_for_timeout(2500)
        ex.scroll_into_view_if_needed(); pg.wait_for_timeout(600)
        _shot(pg, "v211_node_library.png")
        assert not pg.locator("[data-testid='stException']").count()
        # --- undo / redo in the sample list (PXRD)
        nav(pg, "PXRD", answer="discard")
        try:
            c0 = _order_caption(pg)
        except Exception:
            pg.screenshot(path="/tmp/p11/pxrd_fail.png"); raise
        fr = pg.frame_locator("iframe[title*='sample_list']").first
        fr.locator(".row").first.wait_for(timeout=30000)
        fr.locator("input[type=checkbox]").nth(1).click()
        c1 = _wait_caption_change(pg, c0); assert c1 != c0
        und = pg.locator("[data-testid='stSidebar'] button", has_text="Undo").first
        und.wait_for(timeout=30000); pg.wait_for_timeout(1500)
        und.click()
        c2 = _wait_caption_change(pg, c1); assert c2 == c0, (c0, c2)
        assert fr.locator("input[type=checkbox]").nth(1).is_checked()
        pg.locator("[data-testid='stSidebar'] button", has_text="Redo").first.click()
        c3 = _wait_caption_change(pg, c2); assert c3 == c1, (c1, c3)
        pg.locator(".st-key-ps_undo").first.scroll_into_view_if_needed()
        _shot(pg, "v211_undo_redo.png")
        # --- global style preset applies everywhere
        _expand_sidebar(pg, "Figure style")
        sb = pg.locator("[data-testid='stSidebar']")
        inp = sb.locator("[data-testid='stSelectbox']").filter(has_text="Journal preset").first.locator("input").first
        inp.click(); inp.fill("Nature"); pg.wait_for_timeout(600)
        pg.locator("[role='option']").filter(has_text="Nature / Springer Nature").first.click(); pg.wait_for_timeout(4000)
        assert sb.locator("[data-testid='stNumberInput']").filter(has_text="Font size").first.locator("input").input_value() == "7.00"
        _shot(pg, "v211_global_style.png")
        assert not pg.locator("[data-testid='stException']").count()
        # --- project save -> change -> open
        pex = _expand_sidebar(pg, "Project (.reticula)")
        with pg.expect_download(timeout=120000) as dl:
            pex.get_by_role("button", name=re.compile("Save project")).first.click()
        prj = tmp_path / "e2e.reticula"; dl.value.save_as(prj)
        import zipfile
        names = zipfile.ZipFile(prj).namelist(); assert "manifest.json" in names and any(n.startswith("workspace/") for n in names)
        fr.locator("input[type=checkbox]").nth(0).click()                 # change after saving
        c4 = _wait_caption_change(pg, c3); assert c4 != c3
        pex = _expand_sidebar(pg, "Project (.reticula)")
        pex.locator("[data-testid='stFileUploaderDropzoneInput']").first.set_input_files(str(prj)); pg.wait_for_timeout(2500)
        pex.get_by_role("button", name=re.compile("Open e2e.reticula")).first.click(); pg.wait_for_timeout(5000)
        try:
            c5 = _order_caption(pg)
        except Exception:
            pg.screenshot(path="/tmp/p11/open_fail.png", full_page=True); open("/tmp/p11/open_fail.txt","w", encoding="utf-8").write(pg.locator("body").inner_text()); raise
        assert c5 == c3, (c3, c5)
        assert "Opened" in pg.locator("[data-testid='stSidebar']").inner_text()
        _shot(pg, "v211_project_open.png")
        _log(f"v2.11 e2e: node library (HKUST-1 CuO), undo/redo sample list, Nature preset, project save/open ({len(names)} files)")
        b.close()


def test_v211_update_banner_mock(tmp_path):
    """Own server (update check on, mock feed): asked once -> Yes -> banner with Download / Later / Skip."""
    import http.server, socketserver
    from playwright.sync_api import sync_playwright
    rel = {"tag_name": "v9.9.0", "html_url": "https://github.com/Spru1tje/reticula-releases/releases/tag/v9.9.0", "body": "- test",
           "assets": [{"name": "Reticula_Setup_9.9.0.exe", "browser_download_url": "https://example.invalid/Reticula_Setup_9.9.0.exe", "size": 182000000}]}

    class H(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            bb = json.dumps(rel).encode(); self.send_response(200); self.send_header("Content-Length", str(len(bb))); self.end_headers(); self.wfile.write(bb)
        def log_message(self, *a):
            pass
    srv = socketserver.TCPServer(("127.0.0.1", 0), H); threading.Thread(target=srv.serve_forever, daemon=True).start()
    s = socket.socket(); s.bind(("127.0.0.1", 0)); port = s.getsockname()[1]; s.close()
    env = {k: v for k, v in os.environ.items() if k != "MFS_NO_UPDATE_CHECK"}
    env.update(MFS_UPDATE_URL=f"http://127.0.0.1:{srv.server_address[1]}/latest", MFS_STORE=str(tmp_path / "store.json"),
               MFS_DB=str(tmp_path / "r.db"), MFS_EXPORT_DIR=str(tmp_path / "ex"))
    proc = subprocess.Popen([sys.executable, "-m", "streamlit", "run", os.path.join(ROOT, "app.py"), "--server.port", str(port),
                             "--server.headless", "true"], cwd=ROOT, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        for _ in range(120):
            try:
                urllib.request.urlopen(f"http://127.0.0.1:{port}/_stcore/health", timeout=2); break
            except Exception:
                time.sleep(1)
        with sync_playwright() as p:
            b = _browser(p); pg = b.new_page(viewport={"width": 1400, "height": 900})
            pg.goto(f"http://127.0.0.1:{port}"); pg.get_by_text("Check for Reticula updates").first.wait_for(timeout=120000)
            _shot(pg, "v211_update_ask.png")
            pg.get_by_role("button", name="Yes, check").click()
            pg.get_by_text("Reticula 9.9.0 is available").first.wait_for(timeout=30000)
            href = pg.get_by_role("link", name=re.compile("Download")).first.get_attribute("href")
            import platform
            assert href.endswith("Reticula_Setup_9.9.0.exe") if platform.system() == "Windows" else href == rel["html_url"], href
            _shot(pg, "v211_update_banner.png")
            pg.get_by_role("button", name="Skip this version").click(); pg.wait_for_timeout(2500)
            assert not pg.get_by_text("is available").count()
            assert json.load(open(tmp_path / "update.json", encoding="utf-8"))["skip"] == "9.9.0"
            b.close()
    finally:
        proc.terminate(); srv.shutdown()
