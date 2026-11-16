import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tools.make_example_data import composition_masses


def design_mass(tr):
    """Molar mass (RT, with water) of the designed example composition, incl. extra-framework ZrO2 share."""
    Mrt, *_ = composition_masses(tr["L"], tr["formate"], tr["sulfosuccinate"], tr["Cl"], tr["water"])
    fz = tr.get("extra_ZrO2_wt_frac", 0) or 0
    return Mrt / (1 - fz)
