"""v2.13 headless-Chrome checks: ¹H NMR characterisation page (Ketcher editor offline, drag-to-assign, ACS string, figure)."""
import os, re, sys
import pytest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from test_sem_e2e_browser import _browser
from test_v211_e2e_browser import gurl  # noqa: F401  (module fixture: a fresh app server)

SHOTS = os.environ.get("MFS_SHOTS") or os.path.join(__import__("tempfile").gettempdir(), "reticula_shots")


def _zones(fr):
    out = []
    for i in range(fr.locator(".zone").count()):
        e = fr.locator(".zone").nth(i)
        out.append((e.get_attribute("title"), e.bounding_box()))
    return out


def test_v213_nmr_assign_drag(gurl):
    from playwright.sync_api import sync_playwright
    from _e2e_util import nav
    os.makedirs(SHOTS, exist_ok=True)
    with sync_playwright() as p:
        b = _browser(p); pg = b.new_page(viewport={"width": 1400, "height": 1100})
        reqs = []
        pg.on("request", lambda r: reqs.append(r.url))
        pg.goto(gurl); pg.locator("[data-testid='stSidebar']").first.wait_for(timeout=120000)
        nav(pg, "¹H NMR"); pg.wait_for_timeout(3000)
        main = pg.locator("[data-testid='stMain']").first
        main.get_by_text("Characterisation string").first.wait_for(timeout=120000)
        txt = main.inner_text()
        assert "¹H NMR (300 MHz, CDCl3) δ" in txt and "(td, J = 7.8, 1.7 Hz, 1H)" in txt, txt[:3000]
        assert "Fulmer et al." in txt or "Every multiplet is assigned" in txt
        k = pg.locator("iframe[title*='ketcher']").first
        k.scroll_into_view_if_needed(); pg.wait_for_timeout(8000)
        k.screenshot(path=os.path.join(SHOTS, "e2e_ketcher.png"))
        fr = pg.frame_locator("iframe[title*='nmr_assign']").first
        fr.locator(".chip").first.wait_for(timeout=60000)
        a = pg.locator("iframe[title*='nmr_assign']").first; a.scroll_into_view_if_needed(); pg.wait_for_timeout(1500)
        a.screenshot(path=os.path.join(SHOTS, "e2e_assign_before.png"))
        z = _zones(fr)
        tz = [bb for t, bb in z if t.startswith("δ 7.27")]
        assert tz, z
        chip = fr.locator(".chip", has_text=re.compile(r"^6$")).first
        cb = chip.bounding_box()
        pg.mouse.move(cb["x"] + cb["width"] / 2, cb["y"] + cb["height"] / 2); pg.mouse.down()
        pg.mouse.move(tz[0]["x"] + tz[0]["width"] / 2, tz[0]["y"] + tz[0]["height"] / 2, steps=12); pg.mouse.up()
        ok = False
        for _ in range(120):
            pg.wait_for_timeout(500)
            try:
                if any(t.startswith("δ 7.27") and re.search(r"H-6\b", t) for t, _ in _zones(fr)):
                    ok = True; break
            except Exception:
                pass
        assert ok, _zones(fr)
        assert fr.locator(".chip.fix").count() == 1
        a.scroll_into_view_if_needed(); pg.wait_for_timeout(800)
        a.screenshot(path=os.path.join(SHOTS, "e2e_assign_after_drag.png"))
        fr.get_by_text("Automatic assignment").click()
        for _ in range(120):
            pg.wait_for_timeout(500)
            try:
                if fr.locator(".chip.fix").count() == 0:
                    break
            except Exception:
                pass
        assert fr.locator(".chip.fix").count() == 0
        h = main.get_by_text("6 · Labelled figure").first; h.scroll_into_view_if_needed(); pg.wait_for_timeout(2500)
        img = main.locator("[data-testid='stImage'] img, [data-testid='stImage']").first
        pg.screenshot(path=os.path.join(SHOTS, "e2e_figure.png"))
        ext = [u for u in reqs if u.startswith("http") and "localhost" not in u and "127.0.0.1" not in u]
        assert not ext, ext[:10]
        b.close()
