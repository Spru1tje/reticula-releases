"""v2.9: page state survives navigation; unsaved-work dialog; defect calculator unchanged."""
import os, json, sys, tempfile
import pytest
from streamlit.testing.v1 import AppTest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
APP = os.path.join(HERE, "app.py")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def _app(tmp, guard="1"):
    os.environ.update(MFS_STORE=os.path.join(tmp, "store.json"), MFS_DB=os.path.join(tmp, "r.db"), MFS_EXPORT_DIR=tmp,
                      MFS_NAV_GUARD=guard)
    at = AppTest.from_file(APP, default_timeout=400); at.run()
    return at


def _nav(at, p):
    at.sidebar.radio[0].set_value(p).run()
    assert not at.exception, [e.value for e in at.exception][:2]
    return at


def _page(at):
    return at.session_state["nav"]


def _slider(at, label):
    return next(s for s in at.slider if s.label.startswith(label))


def test_settings_survive_navigation(tmp_path):
    at = _app(str(tmp_path), guard="0")
    _nav(at, "PXRD")
    _slider(at, "2θ range").set_value((6.0, 40.0)).run()
    sb = next(s for s in at.selectbox if s.label == "Metal in simulation"); sb.set_value("Zr (as in CIF)").run()
    _nav(at, "TGA"); _nav(at, "XPS"); _nav(at, "PXRD")
    assert _slider(at, "2θ range").value == (6.0, 40.0)
    assert next(s for s in at.selectbox if s.label == "Metal in simulation").value == "Zr (as in CIF)"


def test_dialog_blocks_navigation_and_cancel(tmp_path):
    at = _app(str(tmp_path))
    _nav(at, "PXRD")
    assert _page(at) == "PXRD"
    _slider(at, "2θ range").set_value((6.0, 40.0)).run()
    _nav(at, "TGA")
    assert _page(at) == "PXRD"                              # switch held back
    labels = [b.label for b in at.button]
    assert "💾 Save" in labels and "Don't save" in labels and "Cancel" in labels
    at.button(key="_psd_cancel").click().run()
    assert _page(at) == "PXRD" and _slider(at, "2θ range").value == (6.0, 40.0)


def test_dialog_save_then_return_restores(tmp_path):
    at = _app(str(tmp_path))
    _nav(at, "PXRD")
    _slider(at, "2θ range").set_value((6.0, 40.0)).run()
    _nav(at, "TGA")
    at.button(key="_psd_save").click().run()
    assert _page(at) == "TGA"
    assert os.path.exists(os.path.join(str(tmp_path), "workspace", "pxrd.json"))
    _nav(at, "PXRD")                                       # TGA was not changed → no dialog
    assert _page(at) == "PXRD" and _slider(at, "2θ range").value == (6.0, 40.0)
    _nav(at, "TGA"); assert _page(at) == "TGA"              # saved → leaving is silent
    # a fresh session restores the saved work from the workspace folder
    at2 = _app(str(tmp_path)); _nav(at2, "PXRD")
    assert _slider(at2, "2θ range").value == (6.0, 40.0)


def test_dialog_dont_save_reverts(tmp_path):
    at = _app(str(tmp_path))
    _nav(at, "PXRD")
    _slider(at, "2θ range").set_value((6.0, 40.0)).run()
    _nav(at, "TGA")
    at.button(key="_psd_discard").click().run()
    assert _page(at) == "TGA"
    _nav(at, "PXRD")
    assert _slider(at, "2θ range").value == (4.0, 50.0)     # back to the default
    assert not os.path.exists(os.path.join(str(tmp_path), "workspace", "pxrd.json"))


def test_home_tiles_route_through_guard(tmp_path):
    at = _app(str(tmp_path))
    btn = [b for b in at.button if b.label.endswith("Open PXRD")]
    if not btn:
        pytest.skip("no gallery tile on Home")
    btn[0].click().run()
    assert _page(at) == "PXRD"


def _cmp(a, b):
    assert set(a) == set(b), set(a) ^ set(b)
    for k in a:
        assert a[k] == b[k], k


@pytest.mark.parametrize("guard", ["0", "1"])
def test_defect_calculator_identical_to_v28(tmp_path, guard):
    import _defect_capture as C
    base = json.load(open(os.path.join(HERE, "tests", "data", "defect_baseline_v28.json"), encoding="utf-8"))
    now = json.loads(json.dumps(C.capture(str(tmp_path), nav_guard=guard), ensure_ascii=False))
    _cmp(now, base)
    assert len(now["_store_ids"]) >= 4 and now["pulled::consensus"]       # sees TGA/ICP/NMR/N2 data after navigating


def test_defect_methods_pure_regression():
    """methods a–h + consensus straight from mfs.defects on the stored example inputs (no UI)."""
    from mfs import defects as D, report as REP
    import glob, numpy as np
    base = json.load(open(os.path.join(HERE, "tests", "data", "defect_baseline_v28.json"), encoding="utf-8"))
    for p in sorted(glob.glob(os.path.join(HERE, "example_data", "Defects", "defect_input_*.json"))):
        t = D.results_table(D.run_all(REP.input_from_dict(json.load(open(p)))))
        ref = base[os.path.basename(p)]
        assert len(t) == len(ref)
        for (_, r), rr in zip(t.iterrows(), ref):
            for col in ("linkers/node L", "missing linkers/node x", "caps/node m", "open sites/node"):
                v = r[col]
                if rr.get(col) is None:
                    assert v is None or (isinstance(v, float) and not np.isfinite(v)) or str(v) in ("", "None", "nan")
                elif isinstance(rr[col], float):
                    assert abs(float(v) - rr[col]) < 1e-6, (p, r["method"], col)
