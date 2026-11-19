"""v2.17: simple mode (all black, one symbol per sample)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mfs import figctl as F


def test_simple_black_recolours_and_gives_each_series_its_own_symbol():
    fig, ax = plt.subplots()
    for k in range(4):
        ax.plot([0, 1, 2, 3, 4, 5], [k, k + 1, k + 0.5, k + 2, k + 1, k], color=f"C{k}", label=f"s{k}")
    ax.scatter([1, 2], [1, 2], color="red")
    cfg = F.default_settings()
    cfg["simple_black"] = True
    F.apply(fig, cfg)
    lines = [ln for ln in ax.get_lines() if ln.get_label().startswith("s")]
    assert all(ln.get_color() in ("black", "#000000") for ln in lines)
    # v2.18.2: continuous curves get one line style per sample (no markers, so the graph stays sleek)
    assert len({ln.get_linestyle() for ln in lines}) == 4
    assert all(ln.get_marker() in ("None", None, "") for ln in lines)
    assert all(np.allclose(c.get_facecolor()[0][:3], 0) for c in ax.collections)
    plt.close(fig)


def test_default_keeps_colours():
    fig, ax = plt.subplots()
    ax.plot([0, 1, 2], [0, 1, 0], color="C3", label="x")
    F.apply(fig, F.default_settings())
    assert ax.get_lines()[0].get_color() == "C3"
    plt.close(fig)
