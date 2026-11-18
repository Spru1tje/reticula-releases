"""v2.16: the guided Home page."""
import os, sys, warnings
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _app():
    warnings.filterwarnings("ignore")
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=300)
    at.run()
    return at


def test_home_has_a_start_for_every_tool_and_navigates():
    from mfs import ui
    at = _app()
    assert not at.exception, at.exception
    html = " ".join(m.value for m in at.markdown)
    assert "hm-hero" in html and "What do you want to do" in html and "How it works" in html and "All tools" in html
    opens = [b for b in at.button if b.key.startswith("hm_open_")]
    assert len(opens) == len(ui.PAGES) - 1                                   # one tile per tool
    [b for b in at.button if b.key == "hm_i_1_0"][0].click().run()          # 'I want several techniques in one figure' -> Combined figure
    assert at.session_state["nav"] == "Figure composer" and not at.exception


def test_home_level_cards_switch_the_level(monkeypatch):
    monkeypatch.setenv("MFS_DEFAULT_LEVEL", "Standard")
    at = _app()
    assert [b for b in at.button if b.key == "hm_setlvl_0"][0].disabled          # Standard is the current level
    [b for b in at.button if b.key == "hm_setlvl_2"][0].click().run()
    assert not at.exception and at.session_state["_lv"]["global"] == "Developer"
    assert [b for b in at.button if b.key == "hm_setlvl_2"][0].disabled
