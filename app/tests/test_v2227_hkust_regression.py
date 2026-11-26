"""v2.22.7: regression test of the HKUST-1 net check with atomic (Cromer–Mann) form factors, on the repo reference
HKUST-1.cif (COD 2300380). Synthetic: the 'observed' intensities are the full-structure intensities of that CIF.
The test records what the check gives; it does NOT prove that tbo is the net of HKUST-1.

Stated limits: no measured data, no Debye–Waller factor, and the match uses the node-only net model. The tbo node set
comes from the same CIF (topology.py), so a tbo match is not independent. The release claim about tbo
(helpdocs.py r.207-208) is NOT covered by this test and is decided by the coordinator and the user."""
import os, sys, warnings
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest
from mfs import netmatch as NM

CIF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "references", "library", "HKUST-1.cif")

# values measured when this test was written (2026-11-26, same CIF, R up to |h| <= 12): recorded to detect drift
EXPECTED_R = {"pcu": 0.6582, "fcu": 0.6582, "dia": 0.7703, "tbo": 0.8078}


def _structure_factor_sq():
    pymatgen = pytest.importorskip("pymatgen.core")
    from pymatgen.analysis.diffraction.xrd import ATOMIC_SCATTERING_PARAMS as P
    warnings.filterwarnings("ignore")
    s = pymatgen.Structure.from_file(CIF)
    a_rec = s.lattice.reciprocal_lattice_crystallographic.matrix
    R = 12
    g = np.arange(-R, R + 1)
    H = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T
    H = H[np.abs(H).sum(1) > 0]
    Hf = H.astype(float)
    s2 = (np.linalg.norm(Hf @ a_rec, axis=1) / 2.0) ** 2
    F = np.zeros(len(H), complex)
    for site in s:
        fj = np.zeros(len(H))
        for sp, occ in site.species.items():
            co = np.array(P[sp.symbol], float)
            fj += occ * (sp.Z - 41.78214 * s2 * np.sum(co[:, 0] * np.exp(-co[:, 1] * s2[:, None]), axis=1))
        F += fj * np.exp(2j * np.pi * (Hf @ site.frac_coords))
    Fsq = np.abs(F) ** 2
    keep = Fsq > 1e-9
    return H[keep].astype(float), Fsq[keep]


def test_hkust_ranking_with_atomic_form_factors_matches_the_recorded_values():
    hkl, Fsq = _structure_factor_sq()
    res = NM.rank_nets(hkl, Fsq)
    got = {r["net"]: r["R"] for r in res["ranking"]}
    for net, r_expected in EXPECTED_R.items():
        assert got[net] == pytest.approx(r_expected, abs=2e-4), net


def test_tbo_is_not_the_best_net_under_atomic_form_factors():
    hkl, Fsq = _structure_factor_sq()
    res = NM.rank_nets(hkl, Fsq)
    assert res["best"] != "tbo"
    assert res["confident"] is False
