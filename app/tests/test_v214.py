"""v2.14: Zetasizer Excel import, per-figure font family, legend title size, DLS annotations, progress helper."""
import io
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from mfs import figctl as F, readers as R, style as S  # noqa: E402

CSV = os.path.join(HERE, "example_data", "Zeta-DLS", "Zetasizer_export_pD7.4_D2O.csv")


def _xlsx(preamble_rows=3, sheet="Results"):
    df = pd.read_csv(CSV)
    pre = pd.DataFrame([["Zetasizer export (test)"] + [None] * (df.shape[1] - 1)] * preamble_rows)
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as w:
        pre.to_excel(w, sheet_name=sheet, header=False, index=False, startrow=0)
        df.to_excel(w, sheet_name=sheet, index=False, startrow=preamble_rows)
    return buf.getvalue()


def test_zetasizer_excel_equals_csv():
    a = R.read_zetasizer(CSV, "x.csv")
    b = R.read_zetasizer(_xlsx(), "x.xlsx")
    pd.testing.assert_frame_equal(a["summary"].reset_index(drop=True), b["summary"].reset_index(drop=True), check_dtype=False)
    assert b["dist"] is not None and len(b["dist"]) == len(a["dist"])


def test_zetasizer_excel_without_results_table_explains():
    buf = io.BytesIO()
    pd.DataFrame({"a": [1, 2], "b": [3, 4]}).to_excel(buf, index=False)
    with pytest.raises(ValueError, match="No Zetasizer results table"):
        R.read_zetasizer(buf.getvalue(), "bad.xlsx")


def _fig_with_legend():
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], label="a"); ax.plot([0, 1], [1, 0], label="b")
    ax.legend(title="Samples")
    return fig, ax


def test_font_family_is_applied_to_everything():
    fam = next(f for f in S.available_fonts() if f not in ("auto", "DejaVu Sans"))
    fig, ax = _fig_with_legend(); ax.set_title("t"); ax.set_xlabel("x"); ax.text(0.5, 0.5, "note", transform=ax.transAxes)
    cfg = F.default_settings(); cfg["font"]["family"] = fam
    F.apply(fig, cfg)
    want = S._font_family(fam)
    got = {ax.xaxis.label.get_fontfamily()[0], ax.title.get_fontfamily()[0], ax.texts[0].get_fontfamily()[0],
           ax.get_legend().get_texts()[0].get_fontfamily()[0], ax.get_xticklabels()[0].get_fontfamily()[0]}
    assert got == {want}, got
    plt.close(fig)


def test_legend_title_size_and_family():
    fig, ax = _fig_with_legend()
    cfg = F.default_settings(); cfg["legend"]["title_size"] = 13.0; cfg["legend"]["family"] = "DejaVu Serif"
    F.apply(fig, cfg)
    lg = ax.get_legend()
    assert lg.get_title().get_fontsize() == 13.0
    assert lg.get_title().get_fontfamily()[0] == "DejaVu Serif" and lg.get_texts()[0].get_fontfamily()[0] == "DejaVu Serif"
    plt.close(fig)


def test_house_style_keeps_legend_title_size():
    cfg = F.default_settings(); cfg["legend"]["title_size"] = 11.0; cfg["legend"]["family"] = "Arial"
    h = F.house_style(cfg)
    assert h["legend"]["title_size"] == 11.0 and h["legend"]["family"] == "Arial"


def test_available_fonts_contains_the_classics():
    f = S.available_fonts()
    assert f[0] == "auto" and "Arial" in f and "Times New Roman" in f and len(f) == len(set(f))


@pytest.fixture
def at(tmp_path, monkeypatch):
    monkeypatch.setenv("MFS_STORE", str(tmp_path / "store.json"))
    from streamlit.testing.v1 import AppTest
    a = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=300); a.run()
    return a


def test_zeta_page_annotations_and_legend(at):
    at.sidebar.radio[0].set_value("Zeta / DLS").run()
    assert not at.exception, [e.value for e in at.exception][:2]
    for k in ("zeta_show_z", "zeta_show_pdi", "zeta_show_zp", "zeta_bar_leg"):
        at.checkbox(key=k).check()
    at.run()
    assert not at.exception, [e.value for e in at.exception][:2]


def test_progress_helper_learns_speed_and_returns_result():
    import numpy as np
    from mfs import progress as PG

    class _Box:
        def progress(self, *a, **k): pass
        def empty(self): pass
    p = PG.SegProgress.__new__(PG.SegProgress)
    p.n, p.label, p.box, p.done, p.rates, p._pending = 2, "x", _Box(), 0, {}, 0.0
    import time
    r = p.call(np.zeros((1000, 1000), np.uint8), "auto", lambda img, k: (time.sleep(0.6), k)[1], 7)
    assert r == 7 and p.done == 1 and "auto" in p.rates and p.rates["auto"] > 0.3
    # a fast (cache-like) call does not change the learned rate
    before = dict(p.rates)
    p.call(np.zeros((1000, 1000), np.uint8), "auto", lambda img: None)
    assert p.rates == before
