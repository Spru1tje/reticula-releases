"""v2.16: automatic indexing + Le Bail of the indexed cell."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest
from mfs import indexing as IX, phaseid as PID, rietveld as RV


def lines(name, tt=(10, 90), noise=0.01, seed=2, min_I=0.02):
    e = [e for e in PID.library() if e["name"].startswith(name)][0]
    t, I = PID.sticks(e["path"], 1.540593, tt)
    t = t[I > min_I]
    return t + np.random.default_rng(seed).normal(0, noise, len(t)), e


@pytest.mark.parametrize("name,system,cell", [("Zincite", "hexagonal", (3.2494, 5.2054)), ("Corundum", "hexagonal", (4.7606, 12.994)), ("Rutile", "tetragonal", (4.5937, 2.9587)),
                                                ("Ceria", "cubic", (5.4113,)), ("Silicon", "cubic", (5.4310,)), ("Fluorapatite", "hexagonal", (9.3973, 6.88))])
def test_known_structures_are_indexed(name, system, cell):
    pk, e = lines(name)
    r = IX.index(pk)
    assert r, name
    b = r[0]
    assert b["system"] == system and b["quality"] == "reliable", (name, b)
    assert abs(b["cell"][0] - cell[0]) < 0.02 and (len(cell) == 1 or abs(b["cell"][2] - cell[1]) < 0.02), (name, b["cell"], cell)


def test_tiny_cells_and_few_peaks_are_rejected():
    assert IX.index([20.0, 30.0, 40.0]) == []
    pk, _ = lines("Zincite")
    assert all(min(s["cell"][:3]) >= IX.MIN_EDGE and s["volume"] >= IX.MIN_VOLUME for s in IX.index(pk))


def test_le_bail_of_the_indexed_cell_fits_the_pattern():
    pk, e = lines("Zincite")
    cell = IX.index(pk)[0]["cell"]
    # the 'measured' pattern: the real structure; the Le Bail phase knows only the indexed cell
    x = np.arange(20, 80, 0.02)
    t, I = PID.sticks(e["path"], 1.540593, (15, 85))
    y = np.zeros_like(x)
    for pos, inten in zip(t, I):
        y += inten * np.exp(-4 * np.log(2) * ((x - pos) / 0.12) ** 2)
    y = np.random.default_rng(3).poisson(y / y.max() * 8000 + 100).astype(float)
    ph = RV.Phase.from_cell(cell, "unknown")
    r = RV.quick_refine(x, y, [ph], mode="lebail", offset="zero")
    assert r["stats"]["Rwp"] < 0.12 and abs(r["phases"][0]["cell"][0] - 3.2494) < 0.01 and abs(r["phases"][0]["cell"][2] - 5.2054) < 0.02
