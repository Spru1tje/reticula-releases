"""Shared helper: run the defect calculator page in AppTest on (a) each bundled example input, (b) data pulled from the
sample store after visiting TGA / ICP / NMR / N2 pages with their example data. Returns {case: [table rows]}."""
import os, json, glob
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _table(at):
    for d in at.dataframe:
        v = d.value
        if hasattr(v, "columns") and "linkers/node L" in list(map(str, v.columns)):
            return v
    return None


def _rows(df):
    out = []
    for _, r in df.iterrows():
        out.append({str(k): (None if (isinstance(v, float) and not np.isfinite(v)) else (round(float(v), 6) if isinstance(v, (int, float, np.floating, np.integer)) and not isinstance(v, bool) else str(v)))
                    for k, v in r.items()})
    return out


def capture(tmp, nav_guard="0"):
    from streamlit.testing.v1 import AppTest
    os.environ.update(MFS_STORE=os.path.join(tmp, "store.json"), MFS_DB=os.path.join(tmp, "r.db"), MFS_EXPORT_DIR=tmp, MFS_NAV_GUARD=nav_guard)
    res = {}
    for p in sorted(glob.glob(os.path.join(HERE, "example_data", "Defects", "defect_input_*.json"))):
        at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=400); at.run()
        at.session_state["def_in"] = json.load(open(p))
        at.sidebar.radio[0].set_value("Defect calculator").run()
        assert not at.exception, [e.value for e in at.exception][:2]
        res[os.path.basename(p)] = _rows(_table(at))
        res[os.path.basename(p) + "::consensus"] = [m.label + "=" + str(m.value) for m in at.metric]
    # cross-page: visit the data pages (example data), then pull the stored results in the defect calculator
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=400); at.run()
    for pg in ["TGA", "ICP-OES / XRF", "1H NMR (digest)", "N2 sorption (BET/PSD)"]:
        if pg == "1H NMR (digest)":
            at.session_state["nmr_mode"] = "⚗️ Digestion / composition (ratios, qNMR)"   # v2.13
        at.sidebar.radio[0].set_value(pg).run(); assert not at.exception, (pg, [e.value for e in at.exception][:2])
    at.sidebar.radio[0].set_value("Defect calculator").run()
    opts = [o for o in at.selectbox(key="def_pull").options if not o.startswith("(none")]
    res["_store_ids"] = sorted(opts)
    if opts:
        at.selectbox(key="def_pull").set_value(at.selectbox(key="def_pull").options[1]).run()
        assert not at.exception, [e.value for e in at.exception][:2]
        res["pulled::" + str(opts[0])] = _rows(_table(at)) if _table(at) is not None else []
        res["pulled::consensus"] = [m.label + "=" + str(m.value) for m in at.metric]
    return res
