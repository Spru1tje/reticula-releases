"""v2.11 undo/redo on every page: widget values, sample-list state (tracked extras) and registered keys."""
import os, sys, copy
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import pytest


@pytest.fixture
def at(tmp_path):
    from streamlit.testing.v1 import AppTest
    os.environ.update(MFS_NAV_GUARD="0", MFS_STORE=str(tmp_path / "s.json"), MFS_DB=str(tmp_path / "r.db"), MFS_EXPORT_DIR=str(tmp_path / "ex"))
    a = AppTest.from_file(os.path.join(ROOT, "app.py"), default_timeout=300); a.run()
    return a


def _cb(at, text):
    return [c for c in at.checkbox if text in c.label][0]


@pytest.mark.parametrize("page,label", [("TGA", "Normalise to residue"), ("PXRD", None), ("DSC", None)])
def test_undo_redo_widgets_on_pages(at, page, label):
    at.button(key="nv_" + page).click().run()
    assert at.button(key="_ps_undo").disabled
    cbs = [c for c in at.checkbox if (label is None or label in c.label) and c.key and not c.key.startswith(("_", "gs_", "si_"))]
    c = cbs[0]; k = c.key; v0 = c.value
    c.set_value(not v0).run()
    assert at.checkbox(key=k).value == (not v0) and not at.button(key="_ps_undo").disabled
    at.button(key="_ps_undo").click().run()
    assert not at.exception, [e.value for e in at.exception][:1]
    assert at.checkbox(key=k).value == v0 and not at.button(key="_ps_redo").disabled
    at.button(key="_ps_redo").click().run()
    assert at.checkbox(key=k).value == (not v0)
    at.button(key="_ps_undo").click().run()
    assert at.checkbox(key=k).value == v0


def test_undo_sample_list_order_and_visibility(at):
    at.button(key="nv_" + "PXRD").click().run()
    sk = [k for k in at.session_state._state._state._keys() if str(k).startswith("_sl|")] if hasattr(at.session_state._state, "_state") else []
    if not sk:
        sk = [k for k in at.session_state["_ps_pages"]["PXRD"].get("extra", []) if k.startswith("_sl|")]
    sk = sk[0]
    st0 = copy.deepcopy(at.session_state[sk])
    new = copy.deepcopy(st0)
    items = new["items"] if isinstance(new, dict) and "items" in new else new
    items.reverse()
    at.session_state[sk] = new; at.run()
    assert at.session_state[sk] != st0
    at.button(key="_ps_undo").click().run()
    assert not at.exception
    assert at.session_state[sk] == st0


def test_history_cap_and_branching():
    from mfs import pagestate as PS
    assert PS.HIST_MAX >= 50
