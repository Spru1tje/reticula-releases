"""Shared Playwright helpers. Since v2.9 leaving a page with unsaved changes opens a Save / Don't save / Cancel dialog;
older e2e tests expect the pre-2.9 behaviour (a page starts fresh), so nav() answers "Don't save" by default.
v2.18.6: nav() waits for the sidebar option, polls for the save dialog (up to 15 s) and answers it when it appears; a dialog left
open before a click is answered first. The old fixed 2.5 s wait missed a dialog that appeared later and blocked the next click."""
import time

# v2.19.2: the sidebar pages are buttons in categorised expanders (sidebar_nav.py), no longer radio options.
SIDEBAR_OPT = "[data-testid='stSidebar'] button"
CURRENT_OPT = "[data-testid='stSidebar'] button[data-testid='stBaseButton-primary']"
DIALOG = "[data-testid='stDialog']"


def _answer(pg, answer):
    """Click the button of the open save dialog; returns True when a dialog was answered."""
    dlg = pg.locator(DIALOG).filter(has_text="Do you want to save")
    if not (dlg.count() and dlg.first.is_visible()):
        return False
    label = {"save": "💾 Save", "discard": "Don't save", "cancel": "Cancel"}[answer]
    dlg.first.get_by_role("button", name=label, exact=True).click()
    pg.wait_for_timeout(1500)
    return True


def nav(pg, name, answer="discard", wait_ms=15000):
    """Switch the sidebar page. answer=None: click only and return False (the caller checks the dialog itself).
    Otherwise poll until either the save dialog appears (answered) or the page has switched (no dialog)."""
    opt = pg.locator(SIDEBAR_OPT, has_text=name).first
    opt.wait_for(state="attached", timeout=60000)                      # the sidebar renders after the first load
    if not opt.is_visible():
        # the page sits in a category that is collapsed: open the category that holds this page, as a user would
        cat = pg.locator("[data-testid='stSidebar'] details", has=pg.locator("button", has_text=name)).first
        if cat.count() and not cat.evaluate("d => d.open"):
            cat.locator("summary").first.click()
            pg.wait_for_timeout(600)
    opt.wait_for(state="visible", timeout=60000)
    if answer is not None:
        _answer(pg, answer)           # a dialog left open by an earlier step would block the click
    opt.click()
    if answer is None:
        return False
    deadline = time.time() + wait_ms / 1000
    while time.time() < deadline:
        if _answer(pg, answer):
            return True
        try:
            if pg.locator(CURRENT_OPT, has_text=name).count():           # the current page is the primary button
                pg.wait_for_timeout(800)
                if _answer(pg, answer):          # the dialog can appear just after the switch
                    return True
                return False
        except Exception:
            pass
        pg.wait_for_timeout(300)
    return False


def bright_particles(pg, wait_ms=6000):
    """v2.9: the particle sizer defaults to dark particles (TEM); the bundled SEM example has bright particles."""
    opt = pg.locator("[data-testid='stMain']").get_by_text("bright (on a dark background, SEM)").first
    opt.wait_for(timeout=120000); opt.scroll_into_view_if_needed(); opt.click()
    pg.wait_for_timeout(wait_ms)