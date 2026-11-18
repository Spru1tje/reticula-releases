"""NMR 'simple mode': with simple=True every ink of the 1H figure is black (spectrum, labels, leaders, integrals, marks)."""
from types import SimpleNamespace as N

import matplotlib
import matplotlib.image
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mfs import nmrfig as F


def _data():
    ppm = np.linspace(0, 8, 4000)
    y = np.exp(-((ppm - 3.5) / 0.02) ** 2) + 0.5 * np.exp(-((ppm - 7.2) / 0.03) ** 2)
    sp = N(ppm=ppm, y=y, nucleus="1H")
    mults = [N(centre=3.5, lo=3.4, hi=3.6, lines=[3.55, 3.45], nH=3.0, label="CH3", note=None, mult="s"),
             N(centre=7.2, lo=7.1, hi=7.3, lines=[7.25, 7.15], nH=2.0, label=None, note=None, mult="d")]
    return sp, mults


def _colours(ax):
    return {str(t.get_color()) for t in ax.texts} | {str(ln.get_color()) for ln in ax.lines}


def test_colour_mode_keeps_palette():
    sp, mults = _data()
    fig, ax = plt.subplots()
    F.draw(ax, sp, mults, simple=False)
    plt.close(fig)
    assert _colours(ax) - {"black"}, "the coloured figure should use its palette"


def test_simple_mode_is_all_black():
    sp, mults = _data()
    fig, ax = plt.subplots()
    F.draw(ax, sp, mults, simple=True)
    plt.close(fig)
    assert _colours(ax) == {"black"}


def _max_saturation(png_bytes):
    """Largest channel spread (max - min) over all pixels of a PNG: 0 means pure greys/black/white."""
    import io
    img = matplotlib.image.imread(io.BytesIO(png_bytes), format="png")[..., :3]
    return float((img.max(axis=2) - img.min(axis=2)).max())


def test_simple_mode_figure_with_structure_has_no_colour():
    import io
    import pytest
    rdkit = pytest.importorskip("rdkit")
    from rdkit import Chem
    sp, mults = _data()
    mol = Chem.MolFromSmiles("OC(=O)c1ccc(cc1)C(=O)O")       # O and carbonyl: RDKit colours heteroatoms by default
    buf = io.BytesIO()
    fig, ax = plt.subplots(figsize=(6, 4))
    F.draw(ax, sp, mults, mol=mol, simple=True, inset_box=(0.55, 0.45, 0.4, 0.5))
    fig.savefig(buf, format="png", dpi=100)
    plt.close(fig)
    assert _max_saturation(buf.getvalue()) < 0.02
    buf2 = io.BytesIO()
    fig, ax = plt.subplots(figsize=(6, 4))
    F.draw(ax, sp, mults, mol=mol, simple=False, inset_box=(0.55, 0.45, 0.4, 0.5))
    fig.savefig(buf2, format="png", dpi=100)
    plt.close(fig)
    assert _max_saturation(buf2.getvalue()) > 0.1       # the coloured figure is really coloured
