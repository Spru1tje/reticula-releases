"""v2.8 particle watershed: support holes, splitting of touching faceted particles, aggregates, manual edits."""
import sys, os, time
import numpy as np
from mfs import particles as PA

HERE = os.path.dirname(os.path.abspath(__file__))


def _synth():
    sys.path.insert(0, os.path.join(HERE, "data"))
    import importlib.util
    p = os.path.join(HERE, "data", "synth_support.py")
    spec = importlib.util.spec_from_file_location("synth_support", p); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m.make(1)


def test_watershed_beats_legacy_on_holey_carbon():
    img, lab, ecd, nm = _synth(); crop = img[:-70]; L = lab[:-70]
    edge = set(np.unique(np.concatenate([L[0], L[-1], L[:, 0], L[:, -1]]))) - {0}
    truth = [i for i in np.unique(L) if i and i not in edge]
    t0 = time.time(); l1, f, t, info = PA.segment_ws(crop, True); dt = time.time() - t0
    n1 = len(PA.measure(l1, nm, f))
    l0, f0, _ = PA.segment(crop, dark_particles=True, threshold=t); n0 = len(PA.measure(l0, nm, f0))
    assert info["_holes"] >= 6                       # the bright round holes are found and excluded
    assert n1 > n0 and n1 >= 0.7 * len(truth), (n0, n1, len(truth))
    e = np.array([r["ECD (nm)"] for r in PA.measure(l1, nm, f)])
    assert abs(np.median(e) - 2 * np.sqrt(np.median([(L == i).sum() for i in truth]) / np.pi) * nm) < 0.25 * np.median(e)
    assert dt < 20


def test_touching_discs_split_and_aggregate():
    yy, xx = np.mgrid[:300, :400]
    img = np.full((300, 400), 200, np.uint8)
    for cx, cy in [(100, 150), (160, 150), (220, 150)]:              # three touching discs (necks)
        img[(xx - cx) ** 2 + (yy - cy) ** 2 < 34 ** 2] = 40
    lab, f, t, info = PA.segment_ws(img, True, exclude_holes=False)
    assert len([i for i in np.unique(lab) if i]) == 3
    img2 = np.full((300, 400), 200, np.uint8)                          # a blob that cannot be split: a thick 'L'
    img2[60:240, 60:120] = 40; img2[180:240, 60:330] = 40
    for x0 in range(70, 330, 40):
        img2[180:240, x0:x0 + 2] = 40
    lab2, f2, t2, info2 = PA.segment_ws(img2, True, exclude_holes=False, min_dist_factor=1.5)
    kinds = [v["kind"] for k, v in info2.items() if isinstance(v, dict)]
    assert kinds and ("aggregate" in kinds or len(kinds) >= 2)


def test_manual_edits():
    lab = np.zeros((50, 100), int); lab[10:40, 5:45] = 1; lab[10:40, 55:95] = 2
    m = PA.apply_edits(lab, [("merge", (20, 20), (70, 20))]); assert set(np.unique(m)) == {0, 1}
    s = PA.apply_edits(lab, [("split", (25, 5), (25, 45))]); assert len(np.unique(s)) == 4
    d = PA.apply_edits(lab, [("delete", (70, 20))]); assert set(np.unique(d)) == {0, 1}
    a = PA.apply_edits(np.zeros((50, 50), int), [("add", (25, 25), 6)]); assert (a > 0).sum() > 100


def test_app_watershed_options_and_editor(tmp_path, monkeypatch):
    monkeypatch.setenv("MFS_STORE", str(tmp_path / "store.json"))
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(os.path.join(os.path.dirname(HERE), "app.py"), default_timeout=400); at.run()
    at.sidebar.radio[0].set_value("SEM / TEM / EDX").run(); assert not at.exception
    assert at.radio(key="ps_dark").value == "dark"                 # v2.9 default: dark particles (TEM)
    at.radio(key="ps_dark").set_value("bright").run(); assert not at.exception     # bundled SEM example: bright particles
    assert at.radio(key="ps_algo").value.startswith("adaptive")          # v2.10 default
    assert any("typical radius" in c.value for c in at.caption)
    assert at.toggle(key="ps_prior").value is False                     # shape prior: off by default
    at.toggle(key="ps_prior").set_value(True).run(); assert not at.exception, [e.value for e in at.exception][:2]
    assert any("split by the shape prior" in c.value for c in at.caption)
    at.radio(key="ps_algo").set_value("watershed v2.8 (h-maxima markers, aggregates, support holes)").run(); assert not at.exception
    at.selectbox(key="ps_wsm").set_value("combined").run(); assert not at.exception, [e.value for e in at.exception][:2]
    assert any("excluded by click/auto" in c.value for c in at.caption)
    at.radio(key="ps_pmode").set_value("edit outlines (merge / split / delete / add)").run()
    assert not at.exception, [e.value for e in at.exception][:2]
    assert any("manual edit" in c.value for c in at.caption)
    at.radio(key="ps_algo").set_value("legacy v2.6 (distance markers + nearest-marker split)").run(); assert not at.exception
    assert any("± " in m.value and "n = " in m.value for m in at.markdown)


def test_default_polarity_is_dark():
    import inspect, re
    from mfs import particles as PA
    assert inspect.signature(PA.segment).parameters["dark_particles"].default is True
    assert inspect.signature(PA.segment_ws).parameters["dark_particles"].default is True
    src = open(os.path.join(os.path.dirname(HERE), "mfs", "pages", "microscopy.py"), encoding="utf-8").read()     # v2.13: page code lives in mfs/pages
    assert 'PARTICLE_POLARITY = ["dark", "bright"]' in src and "auto" not in re.search(r'PARTICLE_POLARITY = \[.*\]', src).group(0)
    fix = "/workspace/upgrades/batchC/particles_user_fix.py"
    if os.path.exists(fix):
        assert 'DARK = "--bright" not in sys.argv' in open(fix).read()
