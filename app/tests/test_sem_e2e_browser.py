"""End-to-end UI test in headless Chrome on the user's real TEM images (Olympus SIS / JEOL JEM-1010):
upload -> ⚡ Auto-calibrate -> 200 nm scale bar -> download the PNG figure -> measure the drawn bar in the PNG and
compare with value / nm_per_px x (rendered image width / original width). Also checks the 'Applied to figure'
indicator, the crop of the old burned-in bar and the magnification label.
Set MFS_E2E_URL to test an already running app (e.g. the Windows exe under Wine); otherwise a server is started."""
import os, re, sys, time, hashlib, socket, subprocess
import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "tests", "data", "olympus_sis")
PIX = {"60K_MOF801SGC6.tif": 1.55075160758898, "100K_MOF801SGC6.tif": 0.8733188801503962,
       "120K_MOFSGDC.tif": 0.7167629607361451, "250K_MOF801SGDC.tif": 0.33103385433308367}
W = 1376


def _isolated():
    """Own store / export / workspace folder per test server (v2.9 saved page work would otherwise leak between runs)."""
    import tempfile
    d = tempfile.mkdtemp(prefix="mfs_e2e_")
    return {"MFS_STORE": os.path.join(d, "sample_store.json"), "MFS_DB": os.path.join(d, "reticula.db"), "MFS_EXPORT_DIR": os.path.join(d, "exports")}


@pytest.fixture(scope="module")
def url():
    pytest.importorskip("playwright.sync_api")
    if os.environ.get("MFS_E2E_URL"):
        yield os.environ["MFS_E2E_URL"]; return
    s = socket.socket(); s.bind(("127.0.0.1", 0)); port = s.getsockname()[1]; s.close()
    proc = subprocess.Popen([sys.executable, "-m", "streamlit", "run", os.path.join(ROOT, "app.py"), "--server.port", str(port),
                             "--server.headless", "true"], cwd=ROOT, env=dict(os.environ, MFS_NO_BROWSER="1", **_isolated()),
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(120):
        time.sleep(0.5)
        with socket.socket() as c:
            if c.connect_ex(("127.0.0.1", port)) == 0:
                break
    yield f"http://localhost:{port}"
    proc.terminate(); proc.wait(10)


def _browser(p):
    for kw in ({"channel": "chrome"}, {}):
        try:
            return p.chromium.launch(headless=True, **kw)
        except Exception:
            continue
    pytest.skip("no Chrome/Chromium for Playwright")


def _select(pg, label, option):
    inp = pg.locator("[data-testid='stSelectbox']").filter(has_text=label).first.locator("input").first
    inp.scroll_into_view_if_needed(); inp.click(); pg.wait_for_timeout(500)
    inp.fill(option); pg.wait_for_timeout(800)
    pg.locator("[role='option']").filter(has_text=option).first.click()
    pg.wait_for_timeout(3000)
    assert inp.input_value() == option, (label, inp.input_value())


def measure_bar(png_path):
    """(rendered image width px, longest black horizontal bar run px, row) in the exported PNG."""
    from PIL import Image
    a = np.asarray(Image.open(png_path).convert("L")).astype(int)
    h, w = a.shape
    mid = a[h // 2]
    cols = np.where(mid < 250)[0]
    img_w = cols.max() - cols.min() + 1
    best = (0, None)
    for y in range(int(h * 0.75), h):
        dark = np.concatenate([[0], (a[y] < 20).astype(int), [0]])
        d = np.diff(dark); st_, en = np.where(d == 1)[0], np.where(d == -1)[0]
        if len(st_):
            ln = (en - st_)
            ln = ln[ln < 0.8 * img_w]
            if not len(ln):
                continue
            L = ln.max()
            if L > best[0]:
                best = (int(L), y)
    return int(img_w), best[0], best[1]


@pytest.mark.parametrize("fname", list(PIX))
def test_auto_calibrate_to_exported_figure(url, fname, tmp_path):
    from playwright.sync_api import sync_playwright
    hk = "sem_" + hashlib.md5(fname.encode()).hexdigest()[:8]
    with sync_playwright() as p:
        b = _browser(p)
        pg = b.new_context(accept_downloads=True, viewport={"width": 1400, "height": 1000}).new_page()
        pg.goto(url); pg.wait_for_selector("[data-testid='stSidebar'] [data-testid='stRadioOption']", timeout=180000)
        pg.locator("[data-testid='stSidebar'] [data-testid='stRadioOption']", has_text="SEM / TEM").first.click()
        pg.wait_for_selector("input[type=file][accept*='.tif']", state="attached", timeout=120000); pg.wait_for_timeout(3000)
        pg.locator("input[type=file][accept*='.tif']").first.set_input_files(os.path.join(D, fname))
        auto = pg.get_by_role("button", name="⚡ Auto-calibrate").first
        auto.wait_for(timeout=180000); pg.wait_for_timeout(2000)
        _select(pg, "Colour", "black")
        auto.click(); pg.wait_for_timeout(5000)
        _select(pg, "Scale value", "200 nm")
        pg.get_by_text("Magnification label on figure").first.click(); pg.wait_for_timeout(5000)
        main = pg.locator("[data-testid='stMain']").inner_text()
        assert "Auto-calibrated" in main
        m = re.search(r"Applied to figure: ([\d.]+) nm/px.*?'200 nm' = ([\d.]+) px .*?cropped at row (\d+); magnification label '([^']+)'", main, re.S)
        assert m, main[-3000:]
        nmpx, px_orig, crop, magt = float(m.group(1)), float(m.group(2)), int(m.group(3)), m.group(4)
        assert nmpx == pytest.approx(PIX[fname], rel=1e-3)
        assert px_orig == pytest.approx(200 / PIX[fname], abs=0.11)
        assert crop == 976                                   # old analySIS bar (row 980) + label removed
        assert magt == "×" + f"{int(fname.split('K')[0]) * 1000:,}".replace(",", " ")
        assert not pg.locator("[data-testid='stException']").count()
        exp = pg.locator("[class*='st-key-exportbar_SEM_panel'] button", has_text="PNG").first
        with pg.expect_download(timeout=120000) as dl:
            exp.click()
        path = tmp_path / "fig.png"; dl.value.save_as(path)
        b.close()
    img_w, bar_px, row = measure_bar(str(path))
    expected = 200 / PIX[fname] * img_w / W
    if os.environ.get("MFS_E2E_LOG"):
        with open(os.environ["MFS_E2E_LOG"], "a") as fh:
            fh.write(f"{fname}\tnm/px {nmpx}\tbar(orig) {px_orig:.1f}px\tcrop {crop}\tmag {magt}\tPNG image width {img_w}\t"
                     f"bar in PNG {bar_px}px\texpected {expected:.1f}px\n")
    assert abs(bar_px - expected) <= max(2.0, 0.015 * expected), (fname, img_w, bar_px, expected, row)
