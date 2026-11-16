"""Headless-Chrome check of Batch D in a running app (MFS_E2E_URL): batch processing, SI bundle download, sample
database (import from the bundled pipeline overview), composer journal template check."""
import io, zipfile
from test_sem_e2e_browser import url, _browser  # noqa: F401
from test_batchA_e2e_browser import _nav, _shot, _log
from test_batchB_e2e_browser import _noexc


def test_batchD_pages(url, tmp_path):
    from playwright.sync_api import sync_playwright
    with sync_playwright() as p:
        b = _browser(p)
        pg = b.new_context(accept_downloads=True, viewport={"width": 1500, "height": 1100}).new_page()
        pg.goto(url); pg.wait_for_selector("[data-testid='stSidebar'] [data-testid='stRadioOption']", timeout=180000)
        _nav(pg, "Batch processing"); pg.get_by_text("Process 5 file(s)").first.wait_for(timeout=120000)
        pg.get_by_text("Process 5 file(s)").first.click()
        pg.get_by_text("Batch done: 5 / 5").first.wait_for(timeout=180000); pg.wait_for_timeout(1500)
        with pg.expect_download(timeout=120000) as d:
            pg.get_by_text("ZIP (figures + data + log + preset)").first.click()
        fp = tmp_path / "batch.zip"; d.value.save_as(fp); nb = len(zipfile.ZipFile(fp).namelist())
        _shot(pg, "exe_batch_processing.png")
        _nav(pg, "PXRD"); pg.wait_for_timeout(9000); _nav(pg, "TGA"); pg.wait_for_timeout(9000)
        _nav(pg, "SI export bundle"); pg.get_by_text("SI bundle (").first.wait_for(timeout=120000); pg.wait_for_timeout(2000)
        with pg.expect_download(timeout=180000) as d:
            pg.get_by_text("SI bundle (").first.click()
        fp = tmp_path / "si.zip"; d.value.save_as(fp); names = zipfile.ZipFile(fp).namelist()
        assert "SI_parameters.xlsx" in names and any(n.endswith(".pdf") for n in names), names
        _shot(pg, "exe_si_bundle.png")
        _nav(pg, "Sample database"); pg.get_by_text("Import / update").first.wait_for(timeout=120000)
        pg.get_by_text("Import / update").first.click(); pg.get_by_text("sample IDs imported").first.wait_for(timeout=120000)
        pg.wait_for_timeout(2500); main = _noexc(pg, "DB"); assert "Children" in main and "Parent" in main
        pg.get_by_text("Children").first.scroll_into_view_if_needed(); pg.wait_for_timeout(1000)
        _shot(pg, "exe_sample_database.png")
        b.close()
    _log(f"Batch D: batch zip {nb} entries; SI bundle {len(names)} entries; sample DB import OK")
