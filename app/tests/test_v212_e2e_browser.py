"""v2.12 headless-Chrome checks: movable / editable labels (PXRD) and the sphere model on the user's SiO2-NH2 TEM image."""
import os, re, sys
import pytest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from test_sem_e2e_browser import _browser
from test_batchA_e2e_browser import _shot, _log, SHOTS
from test_v211_e2e_browser import gurl  # noqa: F401  (module fixture: a fresh app server)

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _main(pg):
    return pg.locator("[data-testid='stMain']").first.inner_text()


def _boxes(fr):
    out = {}
    for i in range(fr.locator(".it").count()):
        e = fr.locator(".it").nth(i)
        out[e.get_attribute("data-id")] = (e.bounding_box(), e.get_attribute("title"))
    return out


def _wait_change(pg, fr, old, key, t=40):
    for _ in range(t * 4):
        pg.wait_for_timeout(250)
        try:
            b = _boxes(fr)
        except Exception:
            continue
        if key in b and old.get(key) and (abs(b[key][0]["x"] - old[key][0]["x"]) > 3 or abs(b[key][0]["y"] - old[key][0]["y"]) > 3
                                          or b[key][1] != old[key][1]):
            return b
    return _boxes(fr)


def test_v212_label_mover_pxrd(gurl):
    from playwright.sync_api import sync_playwright
    from _e2e_util import nav
    with sync_playwright() as p:
        b = _browser(p); pg = b.new_page(viewport={"width": 1400, "height": 1100})
        pg.goto(gurl); pg.locator("[data-testid='stSidebar']").first.wait_for(timeout=120000)
        nav(pg, "PXRD"); pg.wait_for_timeout(3000)
        tg = pg.locator("[data-testid='stMain'] [data-testid='stToggle'], [data-testid='stMain'] [data-testid='stCheckbox']").filter(
            has_text="Move / edit labels").first
        tg.scroll_into_view_if_needed(); tg.locator("input").click(force=True)
        fr = pg.frame_locator("iframe[title*='label_mover']").first
        fr.locator(".it").first.wait_for(timeout=60000); pg.wait_for_timeout(1500)
        b0 = _boxes(fr); keys = [k for k in b0 if not k.endswith("__legend__")]
        assert len(keys) >= 2, b0
        k0 = keys[0]; bb = b0[k0][0]
        pg.locator("iframe[title*='label_mover']").first.scroll_into_view_if_needed(); pg.wait_for_timeout(500)
        bb = _boxes(fr)[k0][0]; b0 = _boxes(fr)
        pg.mouse.move(bb["x"] + bb["width"] / 2, bb["y"] + bb["height"] / 2); pg.mouse.down()
        pg.mouse.move(bb["x"] + bb["width"] / 2 - 60, bb["y"] + bb["height"] / 2 - 25, steps=8); pg.mouse.up()
        b1 = _wait_change(pg, fr, b0, k0)
        dx1, dy1 = b1[k0][0]["x"] - b0[k0][0]["x"], b1[k0][0]["y"] - b0[k0][0]["y"]
        assert dx1 < -30 and dy1 < -10, (dx1, dy1)
        # the moved label is selected → apply the same offset to every sample
        e = fr.locator(f".it[data-id='{k0}']"); bb = e.bounding_box()
        pg.mouse.click(bb["x"] + 3, bb["y"] + 3)
        fr.get_by_role("button", name=re.compile("Save this offset")).click()
        b2 = _wait_change(pg, fr, b1, keys[1])
        dxa = b2[keys[1]][0]["x"] - b0[keys[1]][0]["x"]
        assert dxa < -30, (dxa, b0[keys[1]], b2[keys[1]])
        _shot(pg, "v212_label_editor.png")
        # double-click → edit text
        e = fr.locator(f".it[data-id='{keys[1]}']"); e.dblclick()
        ed = fr.locator("#ed"); ed.wait_for(timeout=5000); ed.fill("Renamed label Z"); ed.press("Enter")
        b3 = _wait_change(pg, fr, b2, keys[1])
        assert "Renamed label Z" in (b3[keys[1]][1] or ""), b3[keys[1]]
        assert not pg.locator("[data-testid='stException']").count()
        # switch off → the normal preview keeps the moved labels
        tg.locator("input").click(force=True); pg.wait_for_timeout(2500)
        cap = pg.locator("[data-testid='stMain']").get_by_text(re.compile(r"label\(s\) moved or renamed")).first
        cap.wait_for(timeout=30000); n_mv = int(re.search(r"(\d+) label", cap.inner_text()).group(1))
        assert n_mv >= len(keys)
        _shot(pg, "v212_labels_moved_preview.png")
        # reset
        tg.locator("input").click(force=True); fr = pg.frame_locator("iframe[title*='label_mover']").first
        fr.locator(".it").first.wait_for(timeout=60000)
        fr.get_by_role("button", name=re.compile("Reset labels")).click()
        b4 = _wait_change(pg, fr, b3, keys[1])
        assert abs(b4[keys[1]][0]["x"] - b0[keys[1]][0]["x"]) < 4 and "Renamed label Z" not in (b4[keys[1]][1] or "")
        _log(f"label mover: {len(keys)} sample labels; drag {dx1:.0f},{dy1:.0f} px; offset for all -> 2nd label {dxa:.0f} px; "
             f"text edit; {n_mv} moved in preview; reset")
        b.close()


def test_v212_spheres_real_tem(gurl):
    from playwright.sync_api import sync_playwright
    from _e2e_util import nav
    tif = os.path.join(D, "TEM_SiO2NH2_user6_20k.tif")
    with sync_playwright() as p:
        b = _browser(p); pg = b.new_page(viewport={"width": 1400, "height": 1100})
        pg.goto(gurl); pg.locator("[data-testid='stSidebar']").first.wait_for(timeout=120000)
        nav(pg, "SEM / TEM"); pg.get_by_text("Particle sizing").first.wait_for(timeout=180000); pg.wait_for_timeout(4000)
        pg.locator("[data-testid='stMain'] [data-testid='stFileUploaderDropzoneInput']").first.set_input_files(tif)
        for _ in range(240):
            pg.wait_for_timeout(1000)
            t_ = _main(pg)
            if "TEM_SiO2NH2_user6_20k" in t_ and "Spheres detected" in t_ and "Running" not in t_:
                break
        pg.wait_for_timeout(4000)
        box = pg.locator(".st-key-ps_shape_box").first; box.scroll_into_view_if_needed()
        tx = box.inner_text()
        m = re.search(r"(\d+) circle fits · (\d+) in touching / overlapping groups · (\d+) occluded · (\d+) low-confidence", tx)
        assert m, tx
        nfit, nspl, nocc, nlc = map(int, m.groups())
        assert nspl >= 30 and nlc >= 3, tx
        if SHOTS:
            box.screenshot(path=os.path.join(SHOTS, "v212_sphere_hint.png"))
        fr = pg.frame_locator("iframe[title*='particle_picker']").first
        fr.locator("#stats").wait_for(timeout=180000)
        n0 = int(re.search(r"N = (\d+)", fr.locator("#stats").inner_text()).group(1))
        lc0, lci0 = fr.locator("svg .p.lc").count(), fr.locator("svg .p.lci").count()
        assert lc0 + lci0 == nlc and lc0 >= 1, (lc0, lci0, nlc)   # lci0 > 0 if an earlier session put one back (kept in the sample store)
        pg.locator("iframe[title*='particle_picker']").first.scroll_into_view_if_needed(); pg.wait_for_timeout(800)
        _shot(pg, "v212_picker_spheres_orange.png")
        fr.locator("svg .p.lc").first.click(force=True); pg.wait_for_timeout(500)       # put a low-confidence fit back
        n1 = int(re.search(r"N = (\d+)", fr.locator("#stats").inner_text()).group(1))
        assert n1 == n0 + 1 and fr.locator("svg .p.lci").count() == lci0 + 1
        mean = re.search(r"([\d.]+) ± ([\d.]+) nm", _main(pg))
        assert not pg.locator("[data-testid='stException']").count()
        _log(f"v2.12 spheres (O10 20k): {nfit} fits, {nspl} split, {nocc} occluded, {nlc} low-conf; picker N {n0} -> {n1} "
             f"after including one orange fit; page stats {mean.group(0) if mean else '?'}")
        b.close()
