"""v2.13 BET/N2 graph redesign (t94u): plot functions + sorption page options."""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from mfs import analysis as A, plots as P, readers as R  # noqa: E402

EX = os.path.join(HERE, "example_data", "Sorption")


def _load():
    out = []
    for f in sorted(glob.glob(os.path.join(EX, "*_3Flex_*.xls*")))[:3]:
        d = R.read_sorption(f, os.path.basename(f))
        out.append((os.path.basename(f).split("_")[0], d))
    return out


def test_bet_panel_shows_whole_fit_range():
    nm, d = _load()[0]
    a = d["iso"][d["iso"].branch == "ads"]
    b = A.bet_rouquerol(a.p_rel.values, a.q.values)
    fig, ax = plt.subplots()
    P.draw_bet(ax, b, nm)
    y0, y1 = ax.get_ylim(); x0, x1 = ax.get_xlim()
    pts = [l for l in ax.lines if len(l.get_xdata()) == b["n"]][0]
    assert np.all(pts.get_ydata() <= y1) and np.all(pts.get_xdata() <= x1)
    assert x1 < 0.35 and x1 >= b["p_hi"]
    txt = " ".join(t.get_text() for t in ax.texts)
    assert "S_\\mathrm{BET}" in txt or "BET" in txt
    assert any(p.get_label() != "" for p in ax.patches) or ax.patches or ax.collections   # shaded range
    plt.close(fig)


def test_overlays_and_rouquerol_tplot():
    items = []
    for nm, d in _load():
        a = d["iso"][d["iso"].branch == "ads"]
        b = A.bet_rouquerol(a.p_rel.values, a.q.values)
        if b:
            items.append((nm, b, "C%d" % len(items)))
    fig, (a1, a2, a3) = plt.subplots(1, 3)
    P.draw_bet(a1, items[0][1], overlay=items)
    P.draw_rouquerol(a2, items[0][1], overlay=items)
    a = _load()[0][1]["iso"]; a = a[a.branch == "ads"]
    tp = A.t_plot(a.p_rel.values, a.q.values, (3.5, 5.0))
    P.draw_tplot(a3, tp, (3.5, 5.0), name="x")
    assert len([l for l in a1.lines if l.get_linestyle() == "-"]) >= len(items)
    assert a2.get_xscale() == "log"
    plt.close(fig)


def test_isotherm_options():
    data = _load()
    isos = [(nm, d["iso"]) for nm, d in data]
    psds = [(nm, d["psd"]) for nm, d in data if d["psd"] is not None]
    fig, ax = plt.subplots()
    P.draw_isotherms(ax, isos, hysteresis=True, direct_labels=True, thin=True, inset_psd=psds or None)
    names = {t.get_text() for t in ax.texts}
    assert {nm for nm, _ in isos} <= names                     # labels at the curve ends
    assert ax.get_legend() is not None                          # ads/des key only
    assert {t.get_text() for t in ax.get_legend().get_texts()} == {"adsorption", "desorption"}
    if psds:
        assert len(fig.axes[0].child_axes) == 1
    fig2, ax2 = plt.subplots()
    if psds:
        P.draw_psd(ax2, psds, cumulative=True, legend=True)
        assert len(fig2.axes) == 2                              # twin axis
    plt.close("all")


@pytest.fixture
def at(tmp_path, monkeypatch):
    monkeypatch.setenv("MFS_STORE", str(tmp_path / "store.json"))
    from streamlit.testing.v1 import AppTest
    a = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=300); a.run()
    return a


def test_sorption_page_redesign_options(at):
    at.sidebar.radio[0].set_value("N2 sorption (BET/PSD)").run()
    assert not at.exception, [e.value for e in at.exception][:2]
    keys = [c.key for c in at.checkbox]
    assert "sorp_inset" not in keys          # v2.13: the PSD inset was removed (no graphs inside graphs)
    for k in ("sorp_hyst", "sorp_dlab", "sorp_thin", "sorp_cum", "sorp_bet_all"):
        at.checkbox(key=k).check()
    at.run()
    assert not at.exception, [e.value for e in at.exception][:2]


def test_bet_range_spread_reports_sensitivity():
    f = os.path.join(EX, "SampleA_3Flex_MicroActive.xlsx")
    d = R.read_sorption(f, os.path.basename(f))
    a = d["iso"][d["iso"].branch == "ads"]
    b = A.bet_rouquerol(a.p_rel.values, a.q.values)
    assert b["windows"] and all(w["area"] > 0 for w in b["windows"])
    sp = A.bet_range_spread(b)
    assert sp["min"] <= sp["median"] <= sp["max"] and sp["n_windows"] == len(b["windows"]) and sp["spread_pct"] >= 0
    assert A.bet_range_spread(None) is None
    # the headline result is unchanged by the extra bookkeeping (regression against the stored v2.12 behaviour)
    b2 = A._bet_rouquerol(a.p_rel.values, a.q.values)
    assert abs(b["area"] - b2["area"]) < 1e-9 and b["i"] == b2["i"] and b["j"] == b2["j"]


def test_bet_synthetic_langmuir_like_type_ii_has_tight_spread():
    p = np.array([0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.15, 0.18, 0.22, 0.26, 0.30])
    nm_, C = 100.0, 120.0
    q = nm_ * C * p / ((1 - p) * (1 + (C - 1) * p))         # exact BET isotherm -> every window gives the same n_m
    b = A.bet_rouquerol(p, q)
    sp = A.bet_range_spread(b)
    assert abs(b["n_m"] - nm_) / nm_ < 1e-6 and sp["spread_pct"] < 0.01
