"""v2.15 tests."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest
from mfs import style as S


@pytest.mark.parametrize("fmt", ["png", "pdf", "svg", "tiff"])
def test_provenance_roundtrip(fmt):
    fig, ax = plt.subplots(); ax.plot([1, 2, 3])
    prov = {"version": "2.15.0", "page": "FTIR", "sources": ["a (1).csv", "b.csv"], "customisation": {"x": [1, 2], "t": r"(a) \ b"}}
    data = S.export(fig, fmt, 100, prov)
    assert S.read_provenance(data, fmt) == prov
    assert S.export(fig, fmt, 100) == S.export(fig, fmt, 100)       # without provenance: reproducible bytes
    plt.close(fig)


def test_journal_check_and_fix():
    from mfs import journal as J
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot([0, 1], [0, 1], color="#d62728", lw=0.3, label="a")
    ax.plot([0, 1], [1, 0], color="#2ca02c", lw=0.3, label="b")
    ax.set_xlabel("x", fontsize=8); ax.legend(fontsize=7)
    iss = {i["key"]: i for i in J.check(fig, 3.25, 5.0, 0.5, png_dpi=300)}
    assert iss["font"]["level"] == "error" and iss["line"]["level"] == "error" and iss["colour"]["level"] == "warn" and iss["dpi"]["level"] == "warn"
    done = J.fix(fig, 3.25, 5.0, 0.5)
    assert abs(fig.get_figwidth() - 3.25) < 1e-6 and done
    iss = {i["key"]: i for i in J.check(fig, 3.25, 5.0, 0.5)}
    assert all(iss[k]["level"] == "ok" for k in ("width", "font", "line", "colour")), iss
    plt.close(fig)


def test_figure_data_xlsx():
    import io, pandas as pd
    from mfs import journal as J
    figs = {"fig/one": {"data": pd.DataFrame({"x": [1, 2], "y": [3, 4]}), "page": "A"}, "fig:one": {"data": pd.DataFrame({"z": [1]}), "page": "B"},
            "nodata": {"data": None, "page": "C"}}
    x = pd.read_excel(io.BytesIO(J.figure_data_xlsx(figs, "2.15.0")), sheet_name=None)
    assert list(x)[0] == "Index" and len(x) == 3 and len(x["Index"]) == 3
