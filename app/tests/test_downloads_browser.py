"""Real-browser download test (Chromium via Playwright): every download button on the main pages delivers a
file of the right type and size, also when clicked in quick succession; 'Save to folder' writes the files.
Skipped automatically if Playwright or a Chromium/Chrome browser is not installed."""
import os, socket, subprocess, sys, time, shutil
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAGIC = {"png": b"\x89PNG", "tif": b"II*\x00", "pdf": b"%PDF", "zip": b"PK\x03\x04", "json": b"{"}
PAGES = ["PXRD", "N2 sorption (BET/PSD)", "TGA", "Defect calculator", "SEM / TEM / EDX", "Kinetics & adsorption"]


def _port():
    s = socket.socket(); s.bind(("127.0.0.1", 0)); p = s.getsockname()[1]; s.close(); return p


@pytest.fixture(scope="module")
def server():
    pw = pytest.importorskip("playwright.sync_api")
    port = _port()
    env = dict(os.environ, MFS_NO_BROWSER="1")
    proc = subprocess.Popen([sys.executable, "-m", "streamlit", "run", os.path.join(ROOT, "app.py"), "--server.port", str(port),
                             "--server.headless", "true"], cwd=ROOT, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(120):
        time.sleep(0.5)
        with socket.socket() as s:
            if s.connect_ex(("127.0.0.1", port)) == 0:
                break
    yield f"http://localhost:{port}"
    proc.terminate(); proc.wait(10)


def _nav(pg, name):
    """Click a technique in the sidebar navigation (labels carry an icon prefix since the Reticula UI)."""
    from _e2e_util import nav
    nav(pg, name)


def _browser(p):
    for kw in ({"channel": "chrome"}, {}):
        try:
            return p.chromium.launch(headless=True, **kw)
        except Exception:
            continue
    pytest.skip("no Chromium/Chrome available for Playwright")


def _check(path):
    data = open(path, "rb").read(); ext = path.rsplit(".", 1)[-1].lower()
    if ext == "svg":
        return b"<svg" in data[:2000] and len(data) > 1000
    if ext == "csv":
        return len(data) > 50 and (b"," in data[:3000])
    return data.startswith(MAGIC[ext]) and len(data) > 1000


def test_every_download_arrives(server, tmp_path):
    from playwright.sync_api import sync_playwright
    results, retries = [], []
    with sync_playwright() as p:
        b = _browser(p)
        pg = b.new_context(accept_downloads=True).new_page()
        pg.goto(server, wait_until="networkidle")
        for name in PAGES:
            _nav(pg, name)
            pg.wait_for_timeout(9000)
            btns = pg.locator("[data-testid='stDownloadButton'] button:visible")
            for i in range(btns.count()):
                try:
                    with pg.expect_download(timeout=30000) as dl:
                        btns.nth(i).click()
                    d = dl.value
                except Exception:          # one slow rerun under load: retry once (a lost file fails twice)
                    retries.append((name, i)); pg.wait_for_timeout(2000)
                    with pg.expect_download(timeout=60000) as dl:
                        btns.nth(i).click()
                    d = dl.value
                fp = tmp_path / d.suggested_filename; d.save_as(fp)
                results.append((name, d.suggested_filename, _check(str(fp))))
                pg.wait_for_timeout(50)          # deliberately fast clicking (v1 lost files here)
        b.close()
    bad = [r for r in results if not r[2]]
    assert len(results) >= 60 and not bad, bad
    assert len(retries) <= 2, retries


def test_composer_and_save_to_folder(server, tmp_path):
    from playwright.sync_api import sync_playwright
    exp = os.path.join(ROOT, "exports")
    shutil.rmtree(exp, ignore_errors=True)
    with sync_playwright() as p:
        b = _browser(p)
        pg = b.new_context(accept_downloads=True).new_page()
        pg.goto(server, wait_until="networkidle")
        _nav(pg, "PXRD"); pg.wait_for_timeout(8000)
        pg.get_by_role("button", name="💾 Save to folder").first.click(); pg.wait_for_timeout(8000)
        assert os.path.exists(os.path.join(exp, "PXRD.png")) and os.path.exists(os.path.join(exp, "PXRD.pdf"))
        assert pg.get_by_text("Saved:", exact=False).count() >= 1
        pg.get_by_role("button", name="➕ Add this panel to the figure composer").first.click(); pg.wait_for_timeout(5000)
        _nav(pg, "TGA"); pg.wait_for_timeout(8000)
        pg.get_by_role("button", name="➕ Add this panel to the figure composer").first.click(); pg.wait_for_timeout(5000)
        _nav(pg, "Figure composer"); pg.wait_for_timeout(9000)
        btns = pg.locator("[data-testid='stDownloadButton'] button:visible")
        assert btns.count() >= 5
        for i in range(btns.count()):
            with pg.expect_download(timeout=60000) as dl:
                btns.nth(i).click()
            d = dl.value; fp = tmp_path / d.suggested_filename; d.save_as(fp)
            assert _check(str(fp)), d.suggested_filename
        b.close()
    shutil.rmtree(exp, ignore_errors=True)
