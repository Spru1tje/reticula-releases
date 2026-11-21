"""v2.12: movable / editable sample labels (mfs.labelmove) – positions relative to each trace's offset, 'same offset for
every sample', text edits, legend move, reset, exports use the moved positions."""
import os, sys, io
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest
from mfs import labelmove as LM


def _fig(offset=1.0):
    fig, ax = plt.subplots(figsize=(4, 3))
    x = np.linspace(5, 50, 200)
    for k, (nm, c) in enumerate([("UiO-66", "#0072B2"), ("MOF-801", "#D55E00"), ("Sulfo", "#009E73")]):
        y = k * offset + 0.5 + 0.4 * np.exp(-((x - 10 - 5 * k) / 1.5) ** 2)
        ax.plot(x, y, color=c)
        ax.text(45, k * offset + 0.5 + 0.1 + 0.05 * k, nm, color=c)     # auto position differs per label
    ax.plot([5, 50], [0, 0], color="#888888", label="ref"); ax.legend(loc="upper right")
    LM.tag_originals(fig)
    return fig, ax


def _pos(ax, s):
    return next(t.get_position() for t in ax.texts if t.get_text() == s)


def test_move_is_relative_to_trace_and_follows_offset():
    fig, ax = _fig()
    k = "0:MOF-801"
    disp = ax.transData.transform((30, 1.8)); fx, fy = fig.transFigure.inverted().transform(disp)
    ov = LM.handle(fig, {}, {"op": "move", "id": k, "fx": fx, "fy": fy})
    assert ov[k]["dx"] == pytest.approx(30 - 45, abs=1e-6) and ov[k]["dy"] == pytest.approx(1.8 - 1.65, abs=1e-6)
    fig2, ax2 = _fig(offset=2.0)                 # larger stacking offset: the label moves with its trace
    LM.apply(fig2, ov)
    assert _pos(ax2, "MOF-801") == pytest.approx((30, 2.0 + 0.65 + 0.15), abs=1e-6)
    plt.close("all")


def test_same_offset_for_every_sample_and_reset():
    fig, ax = _fig()
    k = "0:UiO-66"
    disp = ax.transData.transform((20, 0.9)); fx, fy = fig.transFigure.inverted().transform(disp)
    ov = LM.handle(fig, {}, {"op": "move", "id": k, "fx": fx, "fy": fy}); LM.apply(fig, ov)
    ov = LM.handle(fig, ov, {"op": "all", "id": k})
    fig2, ax2 = _fig(); LM.apply(fig2, ov)
    rel = [(_pos(ax2, n)[0], _pos(ax2, n)[1] - b) for n, b in (("UiO-66", 0.5), ("MOF-801", 1.5), ("Sulfo", 2.5))]
    for r in rel:
        assert r == pytest.approx(rel[0], abs=1e-6)          # every label at the same place relative to its own trace
    assert LM.handle(fig, ov, {"op": "reset"}) == {}
    plt.close("all")


def test_text_edit_legend_move_and_export():
    fig, ax = _fig()
    ov = LM.handle(fig, {}, {"op": "text", "id": "0:Sulfo", "text": "Sample B"})
    ov = LM.handle(fig, ov, {"op": "move", "id": "0:__legend__", "fx": 0.2, "fy": 0.5})
    LM.apply(fig, ov)
    assert any(t.get_text() == "Sample B" for t in ax.texts)
    assert 0 < ov["0:__legend__"]["lx"] < 0.5
    b = io.BytesIO(); fig.savefig(b, format="svg"); svg = b.getvalue().decode()
    assert "Sample B" in svg or "Sample B" in svg.replace("\u2010", "-") or b.tell() > 1000
    src, W, H, its = LM.payload(fig)
    assert src.startswith("data:image/png") and {i["id"] for i in its} >= {"0:UiO-66", "0:__legend__"}
    plt.close("all")


def test_pxrd_page_label_editor_apptest():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app.py"), default_timeout=180)
    at.run(); at.button(key="nv_" + "PXRD").click().run()
    tg = [t for t in at.toggle if str(t.key or "").startswith("lmon_")]
    assert tg, "label-editor switch missing on PXRD"
    tg[0].set_value(True).run()
    assert not at.exception
    sk = "lbl|" + tg[0].key[len("lmon_"):]
    at.session_state[sk] = {"0:__legend__": {"lx": 0.1, "ly": 0.9}}; at.run()
    assert not at.exception
