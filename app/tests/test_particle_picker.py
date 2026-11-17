"""v2.9 click-to-exclude particle picker: logic (mirrors the JS), exports respect exclusions, persistence in the store."""
import os, re, json
import numpy as np
from mfs import particle_picker as PK

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P = [{"id": 1, "cx": 10, "cy": 10, "r": 5, "kind": "particle"}, {"id": 2, "cx": 50, "cy": 10, "r": 5, "kind": "aggregate"},
     {"id": 3, "cx": 50, "cy": 60, "r": 5, "kind": "particle"}]


def test_toggle_nearest_box():
    assert PK.toggle([], 3) == [3] and PK.toggle([3], 3) == []
    assert PK.nearest(P, 12, 11) == 1 and PK.nearest(P, 200, 200, snap=30) is None
    assert PK.box(P, [], 0, 0, 60, 20) == [1, 2]
    assert PK.box(P, [1, 2, 3], 0, 0, 60, 20, include=True) == [3]
    rows = [{"id": 1, "kind": "particle"}, {"id": 2, "kind": "aggregate"}, {"id": 3, "kind": "particle"}]
    assert PK.auto_excluded(rows) == [2]
    assert [r["id"] for r in PK.apply(rows, [2, 3])] == [1]


def test_payload_outlines():
    lab = np.zeros((40, 40), int); lab[5:15, 5:15] = 1; lab[20:35, 20:30] = 2
    rows = [{"id": 1, "ECD (nm)": 10.0, "kind": "particle", "x (px)": 9.5, "y (px)": 9.5, "area (nm²)": 100.0},
            {"id": 2, "ECD (nm)": 14.0, "kind": "aggregate", "x (px)": 24.5, "y (px)": 27.0, "area (nm²)": 150.0}]
    pl = PK.particles_payload(rows, lab, 1, 1.0, 1.0, "ECD (nm)")
    assert [p["id"] for p in pl] == [1, 2] and all(len(p["poly"]) >= 4 for p in pl)
    assert abs(pl[0]["cx"] - 9.5) < 1.5 and pl[1]["kind"] == "aggregate"


def _n(at):
    for m in at.markdown:
        g = re.search(r"n = (\d+)\)", m.value)
        if g:
            return int(g.group(1))


def test_exclusions_drive_stats_store_and_survive_navigation(tmp_path, monkeypatch):
    monkeypatch.setenv("MFS_STORE", str(tmp_path / "store.json")); monkeypatch.setenv("MFS_EXPORT_DIR", str(tmp_path))
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=400); at.run()
    at.sidebar.radio[0].set_value("SEM / TEM / EDX").run(); assert not at.exception
    at.radio(key="ps_dark").set_value("bright").run(); assert not at.exception   # bundled SEM example = bright particles
    n0 = _n(at); assert n0
    sig = list(at.session_state["_ps_excl"].values())[0][0]          # the current segmentation's exclusion list
    xk = "psx_" + sig
    assert xk in at.session_state, "exclusion key"
    ex0 = list(at.session_state[xk])
    # exclude two more particles (what two clicks in the picker send back)
    sem = json.load(open(tmp_path / "store.json", encoding="utf-8"))
    ids = list(at.session_state["_ps_ids"].values())[0]
    inc = [i for i in ids if i not in ex0]
    at.session_state[xk] = sorted(set(ex0) | set(inc[:2]))
    at.run()
    n1 = _n(at)
    assert n1 == n0 - 2, (n0, n1)
    df = next(d.value for d in at.dataframe if "n" in d.value.columns)
    assert int(df["n"].iloc[0]) == n1                       # the summary table (and every export built from df) uses only included
    store = json.load(open(tmp_path / "store.json", encoding="utf-8"))
    sem = next(v for v in store.values() if "sem" in v.get("sections", v))
    sec = sem.get("sections", sem)["sem"]
    vals = sec.get("values", sec)
    exk = [k for k in vals if k.startswith("excluded::")]
    assert exk and len(vals[exk[0]]) == len(ex0) + 2
    # navigation away and back keeps the exclusions
    at.sidebar.radio[0].set_value("TGA").run(); at.sidebar.radio[0].set_value("SEM / TEM / EDX").run()
    assert _n(at) == n1
