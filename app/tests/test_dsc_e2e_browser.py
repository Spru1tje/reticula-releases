"""Headless-Chrome check of the DSC page with the user's real Mettler STARe file: upload -> results table and the
STARe comparison box -> PNG export. MFS_E2E_URL tests a running app (e.g. the Windows exe under Wine)."""
import os, re
import pytest
from test_sem_e2e_browser import url, _browser  # noqa: F401  (fixture reuse)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REAL = os.path.join(ROOT, "tests", "data", "Mettler_STARe_SG_DSC_1.1_1K_REAL.txt")


def test_dsc_real_mettler_upload(url, tmp_path):
    if not os.path.exists(REAL):
        pytest.skip("private test data not in the public repository: Mettler_STARe_SG_DSC_1.1_1K_REAL.txt")
    from playwright.sync_api import sync_playwright
    with sync_playwright() as p:
        b = _browser(p)
        pg = b.new_context(accept_downloads=True, viewport={"width": 1400, "height": 1000}).new_page()
        pg.goto(url); pg.wait_for_selector("[data-testid='stSidebar'] [data-testid='stRadioOption']", timeout=180000)
        pg.locator("[data-testid='stSidebar'] [data-testid='stRadioOption']", has_text="DSC").first.click()
        sel = "[class*='st-key-data_dsc'] input[type=file]"
        pg.wait_for_selector(sel, state="attached", timeout=120000); pg.wait_for_timeout(3000)
        pg.locator(sel).first.set_input_files(REAL)
        pg.get_by_text("values reported by Mettler").first.wait_for(timeout=180000); pg.wait_for_timeout(3000)
        main = pg.locator("[data-testid='stMain']").inner_text()
        assert not pg.locator("[data-testid='stException']").count(), main[-2000:]
        m = re.search(r"Reticula: onset ([\d.]+) °C, peak ([\d.]+) °C, ΔH (-?[\d.]+) J/g", main)
        assert m, main[-3000:]
        on, pk, dh = map(float, m.groups())
        assert abs(on - 48.29) < 0.1 and abs(pk - 48.85) < 0.05 and abs(dh + 105.38) < 0.53
        assert "normalized" in main or "dH (J/g) -105.38" in main
        exp = pg.locator("[class*='st-key-exportbar_DSC'] button", has_text="PNG").first
        with pg.expect_download(timeout=120000) as dl:
            exp.click()
        path = tmp_path / "dsc.png"; dl.value.save_as(path)
        b.close()
    from PIL import Image
    assert Image.open(path).size[0] > 500
    if os.environ.get("MFS_E2E_LOG"):
        with open(os.environ["MFS_E2E_LOG"], "a", encoding="utf-8") as fh:
            fh.write(f"DSC real Mettler: onset {on} peak {pk} dH {dh} J/g; PNG {Image.open(path).size}\n")


def test_dsc_full_cycle_example(url, tmp_path):
    """Tick the full-cycle example (REAL heating + SIMULATED cooling, two files) -> combined cycle, summary, 3 views."""
    from playwright.sync_api import sync_playwright
    from test_sem_e2e_browser import _select
    with sync_playwright() as p:
        b = _browser(p)
        pg = b.new_context(accept_downloads=True, viewport={"width": 1400, "height": 1000}).new_page()
        pg.goto(url); pg.wait_for_selector("[data-testid='stSidebar'] [data-testid='stRadioOption']", timeout=180000)
        pg.locator("[data-testid='stSidebar'] [data-testid='stRadioOption']", has_text="DSC").first.click()
        pg.get_by_text("Use the full-cycle example").first.wait_for(timeout=120000); pg.wait_for_timeout(2000)
        pg.get_by_text("Use the full-cycle example").first.click()
        pg.get_by_text("Cycle summary").first.wait_for(timeout=120000); pg.wait_for_timeout(3000)
        main = pg.locator("[data-testid='stMain']").inner_text()
        assert "Combine files of one sample" in main and "Heat 1" in main and "Cool 1" in main
        sizes = {}
        for view in ["Full cycle vs temperature", "Full cycle vs time", "Stacked: heating above cooling"]:
            _select(pg, "View", view); pg.wait_for_timeout(2000)
            assert not pg.locator("[data-testid='stException']").count(), view
            exp = pg.locator("[class*='st-key-exportbar_DSC'] button", has_text="PNG").first
            with pg.expect_download(timeout=120000) as dl:
                exp.click()
            path = tmp_path / f"{len(sizes)}.png"; dl.value.save_as(path)
            from PIL import Image
            sizes[view] = Image.open(path).size
        b.close()
    assert all(w > 500 for w, h in sizes.values())
    assert sizes["Stacked: heating above cooling"][1] > sizes["Full cycle vs temperature"][1]    # two panels
    if os.environ.get("MFS_E2E_LOG"):
        with open(os.environ["MFS_E2E_LOG"], "a", encoding="utf-8") as fh:
            fh.write(f"DSC full cycle example: views OK {sizes}\n")
