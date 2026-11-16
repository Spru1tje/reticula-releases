"""Shared Playwright helpers. Since v2.9 leaving a page with unsaved changes opens a Save / Don't save / Cancel dialog;
older e2e tests expect the pre-2.9 behaviour (a page starts fresh), so nav() answers "Don't save" by default."""
SIDEBAR_OPT = "[data-testid='stSidebar'] [data-testid='stRadioOption']"
DIALOG = "[data-testid='stDialog']"


def nav(pg, name, answer="discard", wait_ms=2500):
    pg.locator(SIDEBAR_OPT, has_text=name).first.click()
    if answer is None:
        return False
    dlg = pg.locator(DIALOG).filter(has_text="Do you want to save")
    try:
        dlg.first.wait_for(state="visible", timeout=wait_ms)
    except Exception:
        return False
    label = {"save": "💾 Save", "discard": "Don't save", "cancel": "Cancel"}[answer]
    dlg.first.get_by_role("button", name=label, exact=True).click()
    pg.wait_for_timeout(1200)
    return True


def bright_particles(pg, wait_ms=6000):
    """v2.9: the particle sizer defaults to dark particles (TEM); the bundled SEM example has bright particles."""
    opt = pg.locator("[data-testid='stMain']").get_by_text("bright (on a dark background, SEM)").first
    opt.wait_for(timeout=120000); opt.scroll_into_view_if_needed(); opt.click()
    pg.wait_for_timeout(wait_ms)
