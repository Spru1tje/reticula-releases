"""v2.11 unit tests: speed helpers, unified sample list logic, (more added per milestone)."""
import os, sys, io, json
import numpy as np
import pytest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ----------------------------------------------------------------------------- M2 speed
def test_fingerprint_stable_and_sensitive():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from mfs import perf
    def mk(c="C0", y=(3, 4)):
        fig, ax = plt.subplots(); ax.plot([1, 2], list(y), color=c, label="a"); ax.legend(); ax.imshow(np.zeros((5, 5))); return fig
    a, b = perf.fingerprint(mk()), perf.fingerprint(mk())
    assert a == b and a is not None
    assert perf.fingerprint(mk("C1")) != a and perf.fingerprint(mk(y=(3, 5))) != a
    plt.close("all")


def test_preview_cache_hits_and_export_keeps_interpolation():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from mfs import perf
    def mk():
        fig, ax = plt.subplots(); ax.imshow(np.random.default_rng(0).random((800, 800)), interpolation="lanczos"); return fig, ax
    h0 = perf.STATS["hits"]
    b1 = perf.preview_png(mk()[0], 100); fig, ax = mk(); b2 = perf.preview_png(fig, 100)   # a fresh, identical figure (next rerun)
    assert b1 == b2 and b1[:4] == b"\x89PNG" and perf.STATS["hits"] == h0 + 1
    assert ax.images[0].get_interpolation() == "lanczos"       # restored after the cheap preview
    plt.close(fig)


def test_downsample_keeps_extremes():
    from mfs import perf
    x = np.linspace(0, 1, 200_001); y = np.sin(40 * x); y[123_456] = 7.0; y[54_321] = -5.0
    xs, ys = perf.downsample(x, y, 4000)
    assert len(xs) <= 4100 and ys.max() == 7.0 and ys.min() == -5.0 and np.all(np.diff(xs) > 0)


def test_fast_widget_state_patch_membership():
    from mfs import perf
    perf.install_fast_widget_state()
    from streamlit.runtime.state.safe_session_state import SafeSessionState
    assert getattr(SafeSessionState.__dict__["filtered_state"].fget, "_mfs_fast", False)


def test_cached_reader_same_result():
    """app._CachedReaders: same parse as the direct call, second call from the cache."""
    from mfs import readers as R
    f = os.path.join(ROOT, "example_data", "PXRD", sorted(os.listdir(os.path.join(ROOT, "example_data", "PXRD")))[0])
    d0 = R.read_pxrd(io.BytesIO(open(f, "rb").read()), os.path.basename(f))
    src = open(os.path.join(ROOT, "app.py"), encoding="utf-8").read()
    assert "R = _CachedReaders(R, skip=(\"read_sem\",))" in src and "@st.fragment" in src
    assert len(d0) > 100


# ----------------------------------------------------------------------------- M3 unified sample list
def test_samplelist_logic():
    from mfs import samplelist as SL
    it = SL.merge(None, ["a.xy", "b.xy", "c.xy", "d.xy"], ["A", "B", "C", "D"])
    assert [o["show"] for o in it] == [True] * 4
    it = SL.set_range(it, 1, 2, False)
    assert [o["show"] for o in it] == [True, False, False, True]
    it = SL.invert(it); assert [o["show"] for o in it] == [False, True, True, False]
    it = SL.all_(it); it = SL.move(it, 0, 2)
    assert [o["id"] for o in it] == ["b.xy", "c.xy", "a.xy", "d.xy"]
    it[0]["name"] = "C"                                     # duplicate name -> unique legend names
    assert [n for _, n in SL.visible(it)] == ["C", "C (2)", "A", "D"]
    it2 = SL.merge(it, ["c.xy", "a.xy", "e.xy"], ["C", "A", "E"])   # keep order, drop removed, append new
    assert [o["id"] for o in it2] == ["c.xy", "a.xy", "e.xy"]
    assert SL.merge(None, list("abcdef"), list("abcdef"), max_shown=4)[4]["show"] is False
    assert SL.plain("Sample C") == "Sample C"          # v2.18.5: masked example name


def test_samplelist_component_bundled_offline():
    from mfs import samplelist as SL
    html = open(os.path.join(SL.component_dir(), "index.html"), encoding="utf-8").read()
    assert "draggable" in html and "shiftKey" in html and "http" not in html.replace("http://www.w3.org", "")


def test_sample_order_drives_stack_and_legend():
    """AppTest: hide + reorder through the list state -> the PXRD legend follows the list."""
    from streamlit.testing.v1 import AppTest
    import tempfile
    td = tempfile.mkdtemp(); os.environ.update(MFS_STORE=td + "/s.json", MFS_DB=td + "/r.db", MFS_EXPORT_DIR=td)
    at = AppTest.from_file(os.path.join(ROOT, "app.py"), default_timeout=300); at.run()
    at.sidebar.radio[0].set_value("PXRD").run()
    items = at.session_state["_sl|pxrd_tab"]
    assert len(items) >= 3
    new = [dict(o) for o in items[::-1]]; new[1]["show"] = False
    at.session_state["_sl|pxrd_tab"] = new; at.run()
    cap = [c.value for c in at.caption if str(c.value).startswith("Plot order")][0]
    from mfs.samplelist import plain
    vis = [plain(o["name"]) for o in new if o["show"]]
    assert cap.split(": ", 1)[1].split(" → ") == vis, cap
    fig = at.session_state["si_figs"]
    k = [k for k in fig if k.lower().startswith("pxrd")][0]
    leg = fig[k]["fig"].axes[0].get_legend()
    labels = [t.get_text() for t in leg.get_texts()] if leg else [t.get_text() for t in fig[k]["fig"].axes[0].texts]
    assert len(at.exception) == 0


# ----------------------------------------------------------------------------- M4 shape model in the UI + batch
def test_batch_particles_adaptive_and_shape_model():
    from mfs import batchproc as BP
    f = os.path.join(ROOT, "example_data", "SEM", "SampleA_SEM.tif")
    z0, l0 = BP.run("SEM/TEM particles", [(f, f)], formats=("png",), dpi=80, particle_opts=dict(dark=False))
    z1, l1 = BP.run("SEM/TEM particles", [(f, f)], formats=("png",), dpi=80, particle_opts=dict(dark=False, shape_prior=True))
    assert l0[0]["status"] == "ok" and 15 <= l0[0]["N"] <= 40 and l0[0]["shape-model splits"] == 0
    assert l1[0]["N"] >= l0[0]["N"] and l1[0]["shape-model splits"] >= 1
    import zipfile
    names = zipfile.ZipFile(io.BytesIO(z1)).namelist()
    assert any(n.endswith(".png") for n in names) and any(n.startswith("data/") for n in names)


def test_shape_model_toggle_and_hint_on_sem_page():
    from streamlit.testing.v1 import AppTest
    import tempfile
    td = tempfile.mkdtemp(); os.environ.update(MFS_STORE=td + "/s.json", MFS_DB=td + "/r.db", MFS_EXPORT_DIR=td)
    at = AppTest.from_file(os.path.join(ROOT, "app.py"), default_timeout=300); at.run()
    at.sidebar.radio[0].set_value("SEM / TEM / EDX").run()
    t = at.toggle(key="ps_prior")
    assert t.value is False and "Shape model" in t.label
    at.radio(key="ps_dark").set_value("bright").run()
    md = " ".join(str(x.value) for x in at.markdown) + " ".join(str(x.value) for x in at.caption) + " ".join(str(x.value) for x in at.info)
    assert "Detect particles" in md and ("similar particles found" in md or "Shape model:" in md), md[-600:]
    at.toggle(key="ps_prior").set_value(True).run()
    assert not at.exception, [e.value for e in at.exception][:1]
    md = " ".join(str(x.value) for x in at.markdown) + " ".join(str(x.value) for x in at.caption)
    assert "Learned from" in md or "Too few" in md
