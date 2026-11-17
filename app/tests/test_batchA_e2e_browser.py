"""Headless-Chrome check of Batch A in a running app (MFS_E2E_URL, e.g. the exe under Wine): kinetics temperature
series -> Arrhenius/Eyring/pD; N2 computed PSD; sample store -> defect calculator pull + SI table XLSX/DOCX."""
import os, re
import pytest
from test_sem_e2e_browser import url, _browser  # noqa: F401

SHOTS = os.environ.get("MFS_SHOT_DIR")


def _nav(pg, name):
    from _e2e_util import nav
    nav(pg, name)


def _shot(pg, name):
    if SHOTS:
        os.makedirs(SHOTS, exist_ok=True); pg.screenshot(path=os.path.join(SHOTS, name), full_page=False)


def _log(s):
    if os.environ.get("MFS_E2E_LOG"):
        with open(os.environ["MFS_E2E_LOG"], "a", encoding="utf-8") as fh:
            fh.write(s + "\n")


def test_batchA_kinetics_n2_store_si(url, tmp_path):
    from playwright.sync_api import sync_playwright
    with sync_playwright() as p:
        b = _browser(p)
        pg = b.new_context(accept_downloads=True, viewport={"width": 1500, "height": 1100}).new_page()
        pg.goto(url); pg.wait_for_selector("[data-testid='stSidebar'] [data-testid='stRadioOption']", timeout=180000)
        # kinetics
        _nav(pg, "Kinetics")
        pg.get_by_text("Use the temperature-series EXAMPLE").first.wait_for(timeout=120000); pg.wait_for_timeout(2000)
        pg.get_by_text("Use the temperature-series EXAMPLE").first.click()
        pg.get_by_text("Arrhenius").nth(1).wait_for(timeout=180000); pg.wait_for_timeout(4000)
        main = pg.locator("[data-testid='stMain']").inner_text()
        assert not pg.locator("[data-testid='stException']").count(), main[-1500:]
        m = re.search(r"Ea = ([\d.]+) ± ([\d.]+) kJ", main); assert m, main[-2500:]
        Ea = float(m.group(1)); assert abs(Ea - 72) < 3
        pk = re.search(r"pKa = ([\d.]+)", main); assert pk and abs(float(pk.group(1)) - 6.2) < 0.2
        pg.get_by_text("Arrhenius and Eyring analysis").first.scroll_into_view_if_needed(); pg.wait_for_timeout(1500)
        _shot(pg, "exe_kinetics_arrhenius.png")
        # N2: computed PSD
        _nav(pg, "N2 sorption")
        pg.get_by_text("Pore-size distribution computed from the isotherm").first.wait_for(timeout=180000); pg.wait_for_timeout(4000)
        assert not pg.locator("[data-testid='stException']").count()
        pg.get_by_text("Pore-size distribution computed from the isotherm").first.scroll_into_view_if_needed(); pg.wait_for_timeout(1500)
        _shot(pg, "exe_n2_computed_psd.png")
        main = pg.locator("[data-testid='stMain']").inner_text()
        assert "Rouquerol" in main and "t-plot" in main
        # TGA to fill the store, then the defect calculator
        _nav(pg, "TGA"); pg.wait_for_timeout(8000)
        _nav(pg, "Defect calculator")
        pg.get_by_text("Pull data for sample").first.wait_for(timeout=120000); pg.wait_for_timeout(3000)
        main = pg.locator("[data-testid='stMain']").inner_text()
        assert "write_stream" not in main      # no st.help output from a stray magic expression
        assert "Sample store" in main and not pg.locator("[data-testid='stException']").count()
        from test_sem_e2e_browser import _select
        _select(pg, "Pull data for sample", "MOF-801(Hf)")
        pg.get_by_text("Pulled for").first.wait_for(timeout=120000); pg.wait_for_timeout(4000)
        assert not pg.locator("[data-testid='stException']").count()
        _shot(pg, "exe_defects_store.png")
        pg.get_by_text("SI table – one consistent").first.click(); pg.wait_for_timeout(2500)
        got = {}
        for lab in ("XLSX", "DOCX"):
            bt = pg.locator("[data-testid='stDownloadButton'] button", has_text=f"SI table ({lab})")
            if bt.count():
                with pg.expect_download(timeout=120000) as dl:
                    bt.first.click()
                pth = tmp_path / f"si.{lab.lower()}"; dl.value.save_as(pth); got[lab] = pth.stat().st_size
                assert open(pth, "rb").read(2) == b"PK"
        _shot(pg, "exe_defects_si_table.png")
        b.close()
    assert got.get("XLSX") and got.get("DOCX"), got
    _log(f"Batch A: kinetics Ea {Ea} kJ/mol, pKa {pk.group(1)}; N2 computed PSD OK; SI table downloads {got}")
