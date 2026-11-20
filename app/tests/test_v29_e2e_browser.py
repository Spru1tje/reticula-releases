"""v2.9 headless-Chrome checks (MFS_E2E_URL = running app, e.g. the exe under Wine; otherwise a local server with the
unsaved-work guard ON is started): particle picker click, unsaved-work dialog (Cancel / Save / Don't save), upload
restored after navigation, Rigaku .rasx in PXRD, UV/Vis page."""
import os, re, sys, socket, subprocess, time, urllib.request
import pytest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))      # repo root (mfs) when run from tests/
from test_sem_e2e_browser import _browser
from test_batchA_e2e_browser import _shot, _log
from _e2e_util import nav, SIDEBAR_OPT, DIALOG
from test_rigaku import _rasx

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture(scope="module")
def gurl(tmp_path_factory):
    if os.environ.get("MFS_E2E_URL"):
        yield os.environ["MFS_E2E_URL"]; return
    s = socket.socket(); s.bind(("127.0.0.1", 0)); port = s.getsockname()[1]; s.close()
    tmp = tmp_path_factory.mktemp("v29srv")
    env = dict(os.environ, MFS_NO_BROWSER="1", MFS_NAV_GUARD="1", MFS_STORE=str(tmp / "store.json"),
               MFS_DB=str(tmp / "r.db"), MFS_EXPORT_DIR=str(tmp))
    proc = subprocess.Popen([sys.executable, "-m", "streamlit", "run", os.path.join(ROOT, "app.py"), "--server.port", str(port),
                             "--server.headless", "true"], cwd=ROOT, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(120):
        try:
            urllib.request.urlopen(f"http://127.0.0.1:{port}/_stcore/health", timeout=2); break
        except Exception:
            time.sleep(1)
    yield f"http://127.0.0.1:{port}"
    proc.terminate(); proc.wait(10)


def _shot_el(loc, name):
    d = os.environ.get("MFS_SHOT_DIR")
    if d:
        os.makedirs(d, exist_ok=True); loc.screenshot(path=os.path.join(d, name))


def _main(pg):
    return pg.locator("[data-testid='stMain']").inner_text()


def _noexc(pg, what):
    m = _main(pg)
    assert not pg.locator("[data-testid='stException']").count(), (what, m[-1500:])
    return m


def _chip(pg):
    return pg.locator(".rt-header .chip").first.inner_text()


def _wait_text(pg, text, t=120000):
    pg.get_by_text(text).first.wait_for(timeout=t)


def _open(p, url):
    b = _browser(p)
    pg = b.new_context(viewport={"width": 1500, "height": 1100}).new_page()
    pg.goto(url); pg.wait_for_selector(SIDEBAR_OPT, timeout=180000); pg.wait_for_timeout(1500)
    return b, pg


def test_particle_picker_click(gurl):
    from playwright.sync_api import sync_playwright
    with sync_playwright() as p:
        b, pg = _open(p, gurl)
        nav(pg, "SEM / TEM"); _wait_text(pg, "Particle sizing", 180000); pg.wait_for_timeout(6000)
        from _e2e_util import bright_particles
        assert "dark (on a light background, TEM)" in pg.locator("[data-testid='stMain']").inner_text()
        bright_particles(pg)
        m0 = re.search(r"n = (\d+)\)", _noexc(pg, "SEM")); assert m0
        n0 = int(m0.group(1))
        fr = pg.frame_locator("iframe[title*='particle_picker']").first
        stats = fr.locator("#stats"); stats.wait_for(timeout=60000)
        s0 = int(re.search(r"N = (\d+)", stats.inner_text()).group(1))
        assert s0 == n0, (s0, n0)
        lab = fr.locator("svg polygon.p:not(.ax):not(.x)").nth(3)          # v2.10: numbers are hidden by default
        lab.scroll_into_view_if_needed(); lab.click(force=True)
        pg.wait_for_timeout(300)
        s1 = int(re.search(r"N = (\d+)", stats.inner_text()).group(1))
        assert s1 == s0 - 1                                     # instant, in the browser
        assert fr.locator("#hg rect").count() >= 5               # live histogram
        n1 = None
        for _ in range(20):
            pg.wait_for_timeout(1000)
            m = re.search(r"n = (\d+)\)", _main(pg))
            if m and int(m.group(1)) == n0 - 1:
                n1 = int(m.group(1)); break
        assert n1 == n0 - 1, _main(pg)[-800:]
        pg.locator("iframe[title*='particle_picker']").first.scroll_into_view_if_needed(); pg.wait_for_timeout(800)
        _shot(pg, "v29_picker_after_click.png")
        fr.locator("#undo").click(); pg.wait_for_timeout(300)
        assert int(re.search(r"N = (\d+)", stats.inner_text()).group(1)) == s0
        _log(f"picker: N {s0} -> {s1} (instant), app n {n0} -> {n1}")
        b.close()


def test_unsaved_dialog_and_restore(gurl, tmp_path):
    from playwright.sync_api import sync_playwright
    rx = tmp_path / "SmartLab_SampleA.rasx"; rx.write_bytes(_rasx([("Omega", 0, 10, 4.3), ("TwoThetaTheta", 4, 40, 8.6)]))
    with sync_playwright() as p:
        b, pg = _open(p, gurl)
        nav(pg, "PXRD"); _wait_text(pg, "Powder X-ray diffraction")
        dz = pg.locator("[data-testid='stMain'] [data-testid='stFileUploaderDropzoneInput']").first
        dz.wait_for(state="attached", timeout=120000); pg.wait_for_timeout(1500)
        for _attempt in range(3):          # the uploader can drop a file set before it is ready: set again until it is parsed
            dz.set_input_files(str(rx))
            try:
                _wait_text(pg, "Rigaku SmartLab", 60000); break
            except Exception:
                if _attempt == 2:
                    raise
                pg.wait_for_timeout(2000)
        pg.wait_for_timeout(4000)
        m = _noexc(pg, "PXRD rasx")
        assert "not 2θ/θ" in m and "Kα1 1.540593" in m
        _shot(pg, "v29_pxrd_rasx.png")
        # leave → dialog → Cancel
        assert not nav(pg, "TGA", answer=None)
        dlg = pg.locator(DIALOG).filter(has_text="Do you want to save")
        dlg.first.wait_for(timeout=15000); pg.wait_for_timeout(800)
        _shot(pg, "v29_unsaved_dialog.png")
        dlg.first.get_by_role("button", name="Cancel", exact=True).click(); pg.wait_for_timeout(2500)
        assert "PXRD" in _chip(pg) and "Rigaku SmartLab" in _main(pg)
        # leave → Save → TGA; come back → the upload is restored
        assert nav(pg, "TGA", answer="save"); pg.wait_for_timeout(3000)
        assert "TGA" in _chip(pg)
        nav(pg, "PXRD"); _wait_text(pg, "Rigaku SmartLab", 120000); pg.wait_for_timeout(3000)   # the uploader keeps the file; the "kept from earlier" caption shows only when it is empty
        m = _noexc(pg, "PXRD restored")
        assert "SmartLab_SampleA.rasx" in m and "Rigaku SmartLab" in m
        _shot(pg, "v29_pxrd_restored.png")
        # change a setting, leave → Don't save → back: reverted to the saved state (file still there)
        pg.get_by_text("Show reflection sticks").first.click(); pg.wait_for_timeout(2500)
        cb = pg.get_by_label("Show reflection sticks").first; cb.wait_for(state="attached", timeout=60000)
        assert not cb.is_checked()
        assert nav(pg, "XPS", answer="discard"); pg.wait_for_timeout(2500)
        assert "XPS" in _chip(pg)
        nav(pg, "PXRD"); _wait_text(pg, "Show reflection sticks"); pg.wait_for_timeout(3000)
        cb = pg.get_by_label("Show reflection sticks").first; cb.wait_for(state="attached", timeout=60000)
        assert cb.is_checked() and "SmartLab_SampleA.rasx" in _main(pg)
        # the defect calculator still sees data entered on other pages
        nav(pg, "TGA"); pg.wait_for_timeout(8000)
        nav(pg, "Defect calculator"); _wait_text(pg, "Pull data for sample"); pg.wait_for_timeout(3000)
        _noexc(pg, "defects")
        from test_sem_e2e_browser import _select
        _select(pg, "Pull data for sample", "Sample A")
        _wait_text(pg, "Pulled for"); pg.wait_for_timeout(3000)
        m = _noexc(pg, "defects pulled"); assert "Linkers / node" in m
        _log("dialog: cancel/save/don't-save ok; upload restored; defect calculator sees the store")
        b.close()


def test_uvvis_page(gurl):
    from playwright.sync_api import sync_playwright
    with sync_playwright() as p:
        b, pg = _open(p, gurl)
        nav(pg, "UV/Vis"); _wait_text(pg, "UV/Vis spectroscopy"); pg.wait_for_timeout(8000)
        m = _noexc(pg, "UV/Vis")
        assert "Put λmax into the sample store" in m
        _shot(pg, "v29_uvvis_spectra.png")
        pg.get_by_role("tab", name=re.compile("Calibration")).click(); pg.wait_for_timeout(4000)
        m = _noexc(pg, "UV cal"); assert "LOQ" in m
        _shot(pg, "v29_uvvis_calibration.png")
        pg.get_by_role("tab", name=re.compile("Tauc")).click(); pg.wait_for_timeout(4000)
        m = _noexc(pg, "UV tauc"); assert "Put the band gaps" in m
        _shot(pg, "v29_uvvis_tauc.png")
        b.close()


def test_v210_real_tem_numbers_and_shape_prior(gurl):
    """v2.10 on the user's TEM image: clean overlay (numbers hidden, toggle), shape prior off → on (dotted splits, panel)."""
    from playwright.sync_api import sync_playwright
    tif = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "TEM_SampleA_user5_sample1.tif")
    with sync_playwright() as p:
        b, pg = _open(p, gurl)
        nav(pg, "SEM / TEM"); _wait_text(pg, "Particle sizing", 180000); pg.wait_for_timeout(4000)
        pg.locator("[data-testid='stMain'] [data-testid='stFileUploaderDropzoneInput']").first.set_input_files(tif)
        for _ in range(240):                                      # the TEM image replaces the bundled examples
            pg.wait_for_timeout(1000)
            t_ = _main(pg)
            if "TEM_SampleA_user5_sample1" in t_ and "Particle sizing" in t_ and "Running" not in t_:
                break
        pg.wait_for_timeout(5000)
        fr = pg.frame_locator("iframe[title*='particle_picker']").first
        stats = fr.locator("#stats"); stats.wait_for(timeout=180000)
        n_off = int(re.search(r"N = (\d+)", stats.inner_text()).group(1))
        assert n_off >= 100, stats.inner_text()                          # adaptive default finds the cubes
        assert fr.locator("svg text.n").first.is_hidden()                 # numbers off by default
        fr.locator("#nums").click(); pg.wait_for_timeout(300)
        assert fr.locator("svg text.n").first.is_visible()
        fr.locator("#nums").click(); pg.wait_for_timeout(300)
        assert fr.locator("svg polygon.pr").count() == 0
        pg.locator("iframe[title*='particle_picker']").first.scroll_into_view_if_needed(); pg.wait_for_timeout(800)
        _shot(pg, "v210_picker_real_tem_prior_off.png")
        box = pg.locator(".st-key-ps_shape_box").first                           # v2.11: prominent toggle + hint
        box.scroll_into_view_if_needed()
        assert "similar particles found" in box.inner_text(), box.inner_text()
        _shot_el(box, "v211_shape_model_hint.png")
        cb = box.locator("[data-testid='stCheckbox'], [data-testid='stToggle']").filter(has_text="Shape model").first.locator("input")
        cb.click(force=True)
        _wait_text(pg, "split by the shape prior", 240000); pg.wait_for_timeout(4000)
        box = pg.locator(".st-key-ps_shape_box").first
        assert "Learned from" in box.inner_text() and box.locator("img").count() >= 1, box.inner_text()
        _shot_el(box, "v211_shape_model_toggle_on.png")
        assert "Learned shape prior" in _main(pg)
        fr = pg.frame_locator("iframe[title*='particle_picker']").first
        fr.locator("svg polygon.pr").first.wait_for(timeout=120000)
        npr = fr.locator("svg polygon.pr").count()
        assert npr >= 10, npr
        s_on = int(re.search(r"N = (\d+)", fr.locator("#stats").inner_text()).group(1))
        assert s_on > n_off, (n_off, s_on)
        fr.locator("svg polygon.pr").first.click(force=True); pg.wait_for_timeout(400)     # a prior split can be excluded
        s2 = int(re.search(r"N = (\d+)", fr.locator("#stats").inner_text()).group(1))
        assert s2 == s_on - 1
        _noexc(pg, "v2.10 prior")
        pg.locator("iframe[title*='particle_picker']").first.scroll_into_view_if_needed(); pg.wait_for_timeout(800)
        _shot(pg, "v210_picker_real_tem_prior_on.png")
        exp = pg.locator("[data-testid='stExpander']").filter(has_text="Learned shape prior").first
        exp.scroll_into_view_if_needed(); pg.wait_for_timeout(800)
        _shot(pg, "v210_prior_panel.png")
        _log(f"v2.10 real TEM: picker N {n_off} (prior off) -> {s_on} (prior on, {npr} dotted splits); exclude prior split -> {s2}")
        b.close()
