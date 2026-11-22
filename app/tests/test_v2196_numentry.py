"""v2.19.6: every slider also takes a typed number (mfs/numentry.py). Known answers: a typed number moves the slider, a dragged slider sets
the number box, typed values are clamped to the range and ranges are kept in order, the select_slider index box picks the option."""
import os
import sys

import pytest

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, APP)

from streamlit.testing.v1 import AppTest  # noqa: E402


def _script():
    import streamlit as st
    from mfs import numentry as NE
    a = NE.slider("Single", 0.0, 10.0, 4.0, 0.5, key="s1")
    b = NE.slider("Range", 0.0, 10.0, (2.0, 6.0), 0.5, key="s2")
    c = NE.select_slider("Pick", options=[293.0, 303.0, 313.0, 323.0], value=313.0, format_func=lambda v: f"{v:.1f} K", key="s3")
    st.session_state["out"] = (a, b, c)


def _app():
    at = AppTest.from_function(_script, default_timeout=30)
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    return at


def _num(at, label_part):
    for n in at.number_input:
        if label_part in (n.label or ""):
            return n
    raise AssertionError(f"no number box with {label_part!r}: {[n.label for n in at.number_input]}")


def test_initial_values_are_the_defaults():
    at = _app()
    assert at.slider[0].value == 4.0
    assert tuple(at.slider[1].value) == (2.0, 6.0)
    assert at.session_state["out"][2] == 313.0
    assert _num(at, "Single").value == 4.0


def test_typed_number_moves_the_slider():
    at = _app()
    _num(at, "Single").set_value(7.5).run()
    assert at.slider[0].value == 7.5


def test_dragged_slider_sets_the_number_box():
    at = _app()
    at.slider[0].set_value(2.5).run()
    assert _num(at, "Single").value == 2.5


def test_typed_value_outside_the_range_keeps_the_slider_inside_it():
    """Streamlit itself returns a number box value above its maximum as the box's minimum, and a browser refuses it, so the
    clamp in numentry is a second guard: the slider must always stay inside its range and no error may appear."""
    at = _app()
    for typed in (99.0, -3.0):
        _num(at, "Single").set_value(typed).run()
        assert not at.exception
        assert 0.0 <= at.slider[0].value <= 10.0


def test_range_boxes_are_kept_in_order():
    at = _app()
    # 'from' typed as 8 while 'to' is 6: the pair is put in order -> (6, 8), and both boxes show it
    _num(at, "Range from").set_value(8.0).run()
    assert tuple(at.slider[1].value) == (6.0, 8.0)
    assert _num(at, "Range from").value == 6.0 and _num(at, "Range to").value == 8.0
    # 'to' typed as 1 with the boxes at (6, 8): the pair becomes (1, 6)
    _num(at, "Range to").set_value(1.0).run()
    assert tuple(at.slider[1].value) == (1.0, 6.0)


def test_select_slider_index_box_picks_the_option():
    at = _app()
    idx = [n for n in at.number_input if "index" in (n.label or "")][0]
    idx.set_value(3).run()
    assert at.session_state["out"][2] == 323.0
    idx.set_value(0).run()
    assert at.session_state["out"][2] == 293.0


def test_no_exception_for_any_of_it():
    at = _app()
    _num(at, "Single").set_value(3.0).run()
    at.slider[1].set_value((1.0, 9.0)).run()
    assert not at.exception


def test_value_set_by_other_code_moves_the_number_box():
    """A value changed by other code (a click in a chart, a reset) must reach the number box as well."""
    at = _app()
    at.session_state["s1"] = 8.0
    at.run()
    assert _num(at, "Single").value == 8.0
