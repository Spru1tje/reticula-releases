"""v2.24.1: end members from a pattern: indexing + Le Bail (cell-only route) on a SYNTHETIC cubic pattern with a known
cell (fluorite-type a = 5.4113 Å, the repo reference CeO2 cell). The forward model is the app's own Le Bail model, so this
checks the chain (peak picking, indexing, fit, symmetrisation), not a real sample. pymatgen and gemmi are required."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest
from mfs import vegard as VG
from mfs import rietveld as RV

A_TRUE = 5.4113


def test_symmetric_cell_removes_the_asymmetry_of_a_fit():
    assert VG.symmetric_cell("cubic", 5.41) == (5.41, 5.41, 5.41, 90.0, 90.0, 90.0)
    assert VG.symmetric_cell("hexagonal", 4.76, 12.99) == (4.76, 4.76, 12.99, 90.0, 90.0, 120.0)
    with pytest.raises(ValueError):
        VG.symmetric_cell("orthorhombic", 4.0)


def test_cell_from_a_synthetic_fluorite_pattern_recovers_the_cell():
    cell = (A_TRUE, A_TRUE, A_TRUE, 90.0, 90.0, 90.0)
    x = np.arange(10.0, 90.0, 0.02)
    ph = RV.Phase.from_cell(cell, "true")
    R = RV.Refinement(x, np.ones_like(x) * 1000.0, [ph], "Cu Kα1", 6, "lebail", counts=True)
    D = R._design(R._values())
    y = D @ np.ones(D.shape[1])
    y = y / y.max() * 20000 + 300
    y = np.random.default_rng(5).poisson(y).astype(float)
    r = VG.cell_from_pattern(x, y, 10.0, 90.0)
    assert r["status"] == "ok", r
    assert r["system"] == "cubic"
    assert r["cell"][0] == pytest.approx(A_TRUE, abs=0.02)
    assert r["cell"][0] == r["cell"][1] == r["cell"][2]
