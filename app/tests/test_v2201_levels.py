"""v2.20.x: the detail level hides Advanced-only expanders on the structure and microscopy pages (Standard) and shows them from Advanced."""
import os
import textwrap

import pytest


def _run(level, tmp_path):
    from streamlit.testing.v1 import AppTest
    script = tmp_path / "gate_script.py"
    script.write_text(textwrap.dedent("""
        import sys, streamlit as st
        sys.path.insert(0, %r)
        from mfs import levels as LV
        st.session_state["nav"] = "Probe"
        with LV.gate("Advanced", globals()), st.expander("Advanced details"):
            st.write("inside")
        st.write("always")
    """ % os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))))
    old = os.environ.get("MFS_DEFAULT_LEVEL")
    os.environ["MFS_DEFAULT_LEVEL"] = level
    try:
        at = AppTest.from_file(str(script), default_timeout=60)
        at.run()
    finally:
        if old is None:
            os.environ.pop("MFS_DEFAULT_LEVEL", None)
        else:
            os.environ["MFS_DEFAULT_LEVEL"] = old
    return [m.value for m in at.markdown] + [e.label for e in at.expander]


def test_standard_hides_advanced_expander(tmp_path):
    shown = _run("Standard", tmp_path)
    assert "Advanced details" not in shown and "inside" not in shown


def test_developer_shows_advanced_expander(tmp_path):
    shown = _run("Developer", tmp_path)
    assert "Advanced details" in shown
