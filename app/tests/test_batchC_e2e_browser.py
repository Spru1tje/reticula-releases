"""Headless-Chrome check of Batch C in a running app (MFS_E2E_URL): SEM particle sizing (auto + click-to-measure),
matched scale bars / one calibration, EDX overlay + line scan, ICP/XRF ratios ± σ."""
import re
from test_sem_e2e_browser import url, _browser  # noqa: F401
from test_batchA_e2e_browser import _nav, _shot, _log
from test_batchB_e2e_browser import _noexc


def _scroll(pg, text):
    pg.get_by_text(text).first.scroll_into_view_if_needed(); pg.wait_for_timeout(1500)


def test_batchC_pages(url):
    from playwright.sync_api import sync_playwright
    with sync_playwright() as p:
        b = _browser(p)
        pg = b.new_context(viewport={"width": 1500, "height": 1100}).new_page()
        pg.goto(url); pg.wait_for_selector("[data-testid='stSidebar'] [data-testid='stRadioOption']", timeout=180000)
        _nav(pg, "SEM / TEM"); pg.get_by_text("Particle sizing").first.wait_for(timeout=180000); pg.wait_for_timeout(6000)
        from _e2e_util import bright_particles
        bright_particles(pg)
        pg.get_by_text("mean ± SD, n =").first.wait_for(timeout=120000); pg.wait_for_timeout(2000)   # the sizing result is computed after the click
        main = _noexc(pg, "SEM")
        m = re.search(r"ECD = (\d+) ± (\d+) nm \(mean ± SD, n = (\d+)\)", main); assert m, main[-2000:]
        auto = tuple(int(g) for g in m.groups())
        _scroll(pg, "Size reported"); _shot(pg, "exe_sem_particles_auto.png")
        _scroll(pg, "Element colours"); pg.wait_for_timeout(1000)
        main = _noexc(pg, "EDX"); assert "Overlay base image" in main
        pg.get_by_text("Line scan across the maps").first.click(); pg.wait_for_timeout(5000)
        _noexc(pg, "EDX line scan"); _scroll(pg, "Band width (px)"); _shot(pg, "exe_edx_overlay_linescan.png")
        # matched bars + one calibration
        pg.get_by_text("apply the first image's calibration to all").first.click(); pg.wait_for_timeout(6000)
        main = _noexc(pg, "SEM onecal"); assert "Batch calibration:" in main and "applied to all" in main
        # click-to-measure: three particles = six clicks
        pg.get_by_text("click-to-measure (two clicks per particle)").first.click(); pg.wait_for_timeout(5000)
        fr = pg.frame_locator("iframe[title*='streamlit_image_coordinates']").last
        img = fr.locator("img").first; img.wait_for(timeout=60000)
        for (x0, y0, x1, y1) in ((300, 300, 340, 300), (500, 400, 560, 400), (700, 200, 750, 200)):
            img.click(position={"x": x0, "y": y0}); pg.wait_for_timeout(3500)
            img.click(position={"x": x1, "y": y1}); pg.wait_for_timeout(3500)
        for _ in range(20):
            main = _noexc(pg, "SEM click")
            m = re.search(r"ECD = (\d+) ± (\d+) nm \(mean ± SD, n = (\d+)\)", main)
            if m and m.group(3) == "3":
                break
            pg.wait_for_timeout(2000)
        assert m and m.group(3) == "3", main[-1500:]
        _scroll(pg, "Clear clicks"); _shot(pg, "exe_sem_click_measure.png")
        _nav(pg, "ICP-OES"); pg.get_by_text("Zr : Hf : S molar ratios").first.wait_for(timeout=120000); pg.wait_for_timeout(2000)
        pg.get_by_text("XRF / any table").first.click(); pg.wait_for_timeout(5000)
        main = _noexc(pg, "XRF"); assert "Zr : Hf : S molar ratios" in main
        _scroll(pg, "Zr : Hf : S molar ratios"); _shot(pg, "exe_xrf_ratios_sigma.png")
        b.close()
    _log(f"Batch C: auto sizing ECD {auto[0]} ± {auto[1]} nm (n = {auto[2]}); click-to-measure n = 3; EDX overlay + line scan OK; "
         f"one calibration OK; XRF ratios ± σ OK")
