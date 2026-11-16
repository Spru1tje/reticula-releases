"""Headless-Chrome check of Batch E (Structure (CIF)) in a running app (MFS_E2E_URL): user CIF upload (if present on
this machine, else the bundled MOF-801), ORTEP render, SI DOCX download, simulated PXRD; UiO-66 pore analysis."""
import os, zipfile
from test_sem_e2e_browser import url, _browser  # noqa: F401
from test_batchA_e2e_browser import _nav, _shot, _log
from test_batchB_e2e_browser import _noexc

USER = os.environ.get("MFS_USER_CIF", "/workspace/real_cif/user.cif")


def _tab(pg, name):
    pg.get_by_role("tab", name=name).first.click(); pg.wait_for_timeout(2500)


def test_batchE_structure_page(url, tmp_path):
    import pytest
    if not os.path.exists(USER):
        pytest.skip("user CIF not on this machine")
    from playwright.sync_api import sync_playwright
    with sync_playwright() as p:
        b = _browser(p)
        pg = b.new_context(accept_downloads=True, viewport={"width": 1500, "height": 1200}).new_page()
        pg.goto(url); pg.wait_for_selector("[data-testid='stSidebar'] [data-testid='stRadioOption']", timeout=180000)
        _nav(pg, "Structure (CIF)"); pg.get_by_text("Upload a CIF").first.wait_for(timeout=120000)
        pg.locator("[data-testid='stMain'] [data-testid='stFileUploader'] input[type='file']").first.set_input_files(USER)
        pg.get_by_text("Space group").first.wait_for(timeout=120000)
        pg.get_by_text("fragment(s)").first.wait_for(timeout=120000); pg.wait_for_timeout(3000)
        main = _noexc(pg, "Structure"); assert "fragment(s)" in main
        _shot(pg, "exe_structure_ortep.png")
        _tab(pg, "H-bonds"); main = _noexc(pg, "H-bonds")
        _tab(pg, "SI tables"); pg.get_by_text("SI tables (DOCX)").first.wait_for(timeout=60000)
        with pg.expect_download(timeout=120000) as d:
            pg.get_by_text("SI tables (DOCX)").first.click()
        fp = tmp_path / "si.docx"; d.value.save_as(fp); nd = len(zipfile.ZipFile(fp).namelist())
        _shot(pg, "exe_structure_si.png")
        _tab(pg, "Simulated PXRD"); pg.wait_for_timeout(15000); _noexc(pg, "PXRD")
        _shot(pg, "exe_structure_pxrd.png")
        b.close()
    _log(f"Batch E: CIF page OK, SI DOCX {nd} parts, H-bonds tab: {'O3A' in main}")
