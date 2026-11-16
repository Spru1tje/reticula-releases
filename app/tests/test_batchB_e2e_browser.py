"""Headless-Chrome check of Batch B in a running app (MFS_E2E_URL): FTIR labels/inset, XPS multi-region fit + at.%,
TGA steps + TGA–MS, DSC cycle comparison + cp (sapphire)."""
import os, re
from test_sem_e2e_browser import url, _browser  # noqa: F401
from test_batchA_e2e_browser import _nav, _shot, _log


def _tick(pg, text):
    pg.get_by_text(text).first.wait_for(timeout=120000); pg.wait_for_timeout(1000)
    pg.get_by_text(text).first.click(); pg.wait_for_timeout(3500)


def _noexc(pg, where):
    main = pg.locator("[data-testid='stMain']").inner_text()
    assert not pg.locator("[data-testid='stException']").count(), (where, main[-1500:])
    return main


def test_batchB_pages(url):
    from playwright.sync_api import sync_playwright
    with sync_playwright() as p:
        b = _browser(p)
        pg = b.new_context(viewport={"width": 1500, "height": 1100}).new_page()
        pg.goto(url); pg.wait_for_selector("[data-testid='stSidebar'] [data-testid='stRadioOption']", timeout=180000)
        _nav(pg, "FTIR"); pg.get_by_text("Reference band, peak labels").first.wait_for(timeout=120000)
        pg.get_by_text("Reference band, peak labels").first.click(); pg.wait_for_timeout(1500)
        for t in ("Normalise to a reference band", "Label peaks", "Zoomed inset"):
            _tick(pg, t)
        main = _noexc(pg, "FTIR"); assert "Difference spectra" in main
        pg.locator("[data-testid='stImage'], [data-testid='stImageContainer'] img").first.scroll_into_view_if_needed(); pg.wait_for_timeout(1500)
        _shot(pg, "exe_ftir_labels_inset.png")
        _nav(pg, "XPS"); _tick(pg, "Use the multi-region EXAMPLE")
        pg.get_by_text("Atomic % (RSF-corrected)").first.wait_for(timeout=120000); pg.wait_for_timeout(3000)
        main = _noexc(pg, "XPS"); assert "Atomic %" in main and "Peak fitting" in main and "Region / doublet constraint" in main
        pg.get_by_text("Peak fitting").first.scroll_into_view_if_needed(); pg.wait_for_timeout(1500)
        _shot(pg, "exe_xps_fit.png")
        _nav(pg, "TGA"); pg.get_by_text("TGA–MS overlay").first.wait_for(timeout=120000); pg.wait_for_timeout(4000)
        main = _noexc(pg, "TGA"); assert "Automatic step labels" in main
        pg.get_by_text("TGA–MS overlay").first.scroll_into_view_if_needed(); pg.wait_for_timeout(1500)
        _shot(pg, "exe_tga_ms.png")
        _nav(pg, "DSC"); _tick(pg, "Use the heat–cool–heat EXAMPLE")
        pg.get_by_text("Cycle-to-cycle comparison").first.wait_for(timeout=120000); pg.wait_for_timeout(3000)
        main = _noexc(pg, "DSC"); assert "Cycle-to-cycle comparison" in main and "Heat 2" in main
        pg.get_by_text("Cycle-to-cycle comparison").first.scroll_into_view_if_needed(); pg.wait_for_timeout(1500)
        _shot(pg, "exe_dsc_cycle_compare.png")
        pg.get_by_text("cp from blank + sapphire").first.click(); pg.wait_for_timeout(5000)
        main = _noexc(pg, "DSC cp"); assert "sapphire" in main
        pg.get_by_text("cp from blank + sapphire").first.scroll_into_view_if_needed(); pg.wait_for_timeout(1500)
        _shot(pg, "exe_dsc_cp.png")
        b.close()
    _log(f"Batch B: FTIR labels/inset OK; XPS fit + at.% OK; TGA–MS OK; DSC cycle compare + cp OK")
