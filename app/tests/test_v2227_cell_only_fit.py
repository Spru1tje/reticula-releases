"""v2.22.7: positive control for the cell-only Le Bail route. A SYNTHETIC I432 pattern with a known cell is generated
with the same forward model as the Le Bail refinement, then refined with quick_refine(mode='lebail') from a CIF that
has a cell and a space group but no atoms (cell-only reserve). The fit must recover the cell.

Limits (stated on purpose): the pattern is synthetic and uses the app's own forward model, so this checks the fit
chain (reflection list, absences, profile, refinement), not the physics of a real sample. No real XRD data is used."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest
from mfs import rietveld as RV

A_TRUE = 20.00
CELL_CIF = f"""data_cell_only_positive_control
_cell_length_a {A_TRUE}
_cell_length_b {A_TRUE}
_cell_length_c {A_TRUE}
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'I 4 3 2'
_symmetry_Int_Tables_number 211
"""


@pytest.fixture
def cell_cif(tmp_path):
    p = tmp_path / "cell_only_i432.cif"
    p.write_text(CELL_CIF, encoding="utf-8")
    return str(p)


def test_cell_only_lebail_fit_recovers_the_cell_from_a_synthetic_pattern(cell_cif):
    pytest.importorskip("gemmi")
    x = np.arange(5.0, 60.0, 0.02)
    ph = RV.Phase.from_cif(cell_cif)
    R = RV.Refinement(x, np.ones_like(x) * 1000.0, [ph], "Cu Kα1", 6, "lebail", counts=True)
    # true parameters: the cell of the phase, a sharp profile, a flat background
    for k, v in {"c0.a": A_TRUE, "W": 0.004, "X": 0.05, "B0": 1.0}.items():
        if k in R.P:
            R.P[k].value = v
    y_model = R._design(R._values()) @ np.ones(R._design(R._values()).shape[1])
    y_model = y_model / y_model.max() * 20000 + 300
    rng = np.random.default_rng(11)
    y = rng.poisson(y_model).astype(float)
    fit = RV.quick_refine(x, y, [cell_cif], radiation="Cu Kα1", mode="lebail", tt_range=(5, 60))
    a_fit = fit["phases"][0]["cell"][0]
    assert abs(a_fit - A_TRUE) < 0.02
    assert fit["stats"]["Rwp"] < 0.15
