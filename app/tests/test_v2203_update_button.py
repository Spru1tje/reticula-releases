"""v2.20.3: the 'Check now' button works without the automatic-check opt-in, and shows the result at once (no network: fetch is faked)."""
import textwrap
import os

from streamlit.testing.v1 import AppTest


def _script(tmp_path):
    here = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    p = tmp_path / "upd_script.py"
    p.write_text(textwrap.dedent(f"""
        import sys, streamlit as st
        sys.path.insert(0, {here!r})
        from mfs import update as U
        U.settings_box({str(tmp_path)!r}, "2.20.3")
    """))
    return str(p)


def test_check_now_works_with_opt_in_off(tmp_path, monkeypatch):
    monkeypatch.delenv("MFS_NO_UPDATE_CHECK", raising=False)
    from mfs import update as U
    # fake the network for this test only; monkeypatch restores it afterwards (no leak into other tests)
    monkeypatch.setattr(U, "fetch", lambda url, timeout=0: {"tag_name": "v9.9.9", "html_url": "x", "body": "", "assets": [
        {"name": "Reticula_Setup_9.9.9.exe", "browser_download_url": "u"}]})
    at = AppTest.from_file(_script(tmp_path), default_timeout=60)
    at.run()
    btn = [b for b in at.button if b.key == "_upd_now"][0]
    assert not btn.disabled
    btn.click()
    at.run()
    texts = " ".join(x.value for x in at.success) + " ".join(x.value for x in at.warning)
    assert "9.9.9" in texts and "newer" in texts
