"""Every Streamlit page renders without exception on the bundled example data; the composer collects panels."""
import pytest
from apptest_all_pages import PAGES, run, composer_check


@pytest.mark.parametrize("page", PAGES)
def test_page_renders(page):
    at = run(page)
    assert not at.exception, [e.value for e in at.exception]
    assert not at.error, [e.value for e in at.error]


def test_composer_collects_panels():
    n, ndl, errs = composer_check()
    assert n == 2 and not errs and ndl >= 4


def _goto(page, attempts=3):
    """Open a page in AppTest. On Windows, asyncio's socketpair can fail with WinError 10055 (loopback buffer/port exhaustion
    after many test runs in a row); that is an environment error, so the session is retried a few times. Other errors are raised."""
    import time
    from streamlit.testing.v1 import AppTest
    import os
    for i in range(attempts):
        try:
            at = AppTest.from_file(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app.py"), default_timeout=300)
            at.run(); at.button(key="nv_" + page).click().run()
            return at
        except OSError as e:
            if getattr(e, "winerror", None) != 10055 or i == attempts - 1:
                raise
            time.sleep(2.0)


def _ok(at):
    assert not at.exception, [e.value for e in at.exception]
    assert not [w.value for w in at.warning if "Customisation not applied" in w.value]


def test_figure_controls_in_app():
    at = _goto("N2 sorption (BET/PSD)")
    at.checkbox(key="fc_N2_isotherms_leg_custom").check().run(); _ok(at)
    at.selectbox(key="fc_N2_isotherms_leg_loc").set_value("outside right").run(); _ok(at)
    at.number_input(key="fc_N2_isotherms_leg_ncol").set_value(2).run(); _ok(at)
    at.selectbox(key="fc_N2_isotherms_preset").set_value("2 sides (bottom + left)").run(); _ok(at)
    at.radio(key="fc_N2_isotherms_direction").set_value("out").run(); _ok(at)
    at.number_input(key="fc_N2_isotherms_p_0_ymajor_n").set_value(25.0).run(); _ok(at)
    assert at.slider(key="fc_N2_isotherms_p_0_ymajor_s").value == 25.0      # slider follows the number box
    at.selectbox(key="fc_N2_isotherms_p_0_y_unit").set_value("mmol g⁻¹").run(); _ok(at)
    at.text_input(key="fc_N2_isotherms_title").set_value("N_2 at 77 K").run(); _ok(at)


def test_scalebar_controls_in_app():
    import hashlib
    at = _goto("SEM / TEM / EDX")
    hk = "sem_" + hashlib.md5("SampleA_SEM.tif".encode()).hexdigest()[:8]
    at.selectbox(key=f"{hk}_pre").set_value("1 µm").run(); _ok(at)
    assert any("label '1 µm'" in c.value for c in at.caption)
    at.radio(key=f"{hk}_mode").set_value("length (bar size)").run(); _ok(at)
    at.number_input(key=f"{hk}_lpx").set_value(400.0).run(); _ok(at)
    assert any("400 px = 1040 nm" in c.value for c in at.caption)
    at.radio(key=f"{hk}_src").set_value("manual (nm/px)").run()
    at.number_input(key=f"{hk}_pix").set_value(5.0).run(); _ok(at)
    assert any("400 px = 2000 nm" in c.value and "'2 µm'" in c.value for c in at.caption)
    at.selectbox(key="sem_all_pos").set_value("upper left").run(); _ok(at)
    at.checkbox(key="sem_all_box").check().run(); _ok(at)


def test_goto_retries_a_winerror_10055_and_then_succeeds(monkeypatch):
    """The retry in _goto is proven with a simulated WinError 10055: the first two AppTest sessions fail, the third works."""
    from streamlit.testing.v1 import AppTest
    import test_app_smoke as me
    real = AppTest.from_file
    calls = {"n": 0}

    def flaky(*a, **k):
        calls["n"] += 1
        if calls["n"] <= 2:
            err = OSError("simulated socket buffer full")
            err.winerror = 10055
            raise err
        return real(*a, **k)

    monkeypatch.setattr(AppTest, "from_file", staticmethod(flaky))
    # the 2 s pause in _goto is kept real: patching time.sleep globally would stall Streamlit's own run loop
    at = me._goto("SEM / TEM / EDX")
    assert calls["n"] == 3
    assert not at.exception


def test_goto_does_not_retry_other_oserrors(monkeypatch):
    from streamlit.testing.v1 import AppTest
    import test_app_smoke as me
    calls = {"n": 0}

    def broken(*a, **k):
        calls["n"] += 1
        raise OSError("disk on fire")

    monkeypatch.setattr(AppTest, "from_file", staticmethod(broken))
    with pytest.raises(OSError):
        me._goto("SEM / TEM / EDX")
    assert calls["n"] == 1
