"""v2.17: the structure picture is centred on the molecule and its labels have a halo."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mfs import cif as CX, cifplot as CP

CIF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "references", "MOF-801_Zr_COD4121459.cif")


def test_structure_is_centred_and_labels_have_halo():
    S = CX.read(CIF, "MOF-801")[0]
    sites, bonds = S.grow(mode="molecules")
    fig, ax = plt.subplots(figsize=(4, 4))
    V = CP.view_matrix(S, "best", sites=sites)
    P = CP.draw(ax, S, sites, bonds, view=V, show_cell=True, labels="non-H")
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    cx, cy = (P[:, 0].min() + P[:, 0].max()) / 2, (P[:, 1].min() + P[:, 1].max()) / 2
    assert abs(cx - (x0 + x1) / 2) < 0.08 * (x1 - x0) and abs(cy - (y0 + y1) / 2) < 0.08 * (y1 - y0)
    assert any(t.get_bbox_patch() is not None for t in ax.texts)
    plt.close(fig)
