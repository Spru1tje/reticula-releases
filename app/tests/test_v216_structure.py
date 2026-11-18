"""v2.16: symmetry-aware structure refinement (atom positions, B per element, soft bond restraints)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest
from mfs import rietveld as RV, structmodel as SM

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
FAP_CIF = os.path.join(DATA, "FAP_REAL.cif")
FAP_XRA = os.path.join(DATA, "FAP_REAL.xra")
needs_data = pytest.mark.skipif(not os.path.exists(FAP_XRA), reason="real data not present (kept out of the repository)")
RAD = [(1.5405, 1.0), (1.5443, 0.5)]


def _fap_pattern():
    L = open(FAP_XRA).read().splitlines()
    vals = []
    for ln in L[2:]:
        vals += [float(ln[8 * k:8 * k + 8]) for k in range(10) if ln[8 * k:8 * k + 8].strip()]
    y = np.array(vals[:5753])
    return 15.0 + 0.02 * np.arange(len(y)), y


@needs_data
def test_model_finds_the_space_group_special_positions_and_consistent_structure_factors():
    from pymatgen.core import Structure
    m = SM.StructureModel(Structure.from_file(FAP_CIF))
    assert m.sg_number == 176 and m.spacegroup.replace(" ", "") == "P6_3/m" and m.n_ops == 12 and len(m.reps) == 7
    assert sorted(m.basis[r].shape[1] for r in m.reps) == [0, 1, 2, 2, 2, 2, 3]            # F on the 6_3 axis: no freedom; Ca1: z only; O3 general
    assert np.abs(m.site_positions() - m.site_frac).max() < 1e-9
    ph = RV.Phase.from_cif(FAP_CIF).build(1.0)
    inv, reps, cnt = m.classify(ph.m_hkl.astype(int))
    Gs = np.linalg.inv(RV._cell_metric(ph.cell0))
    s2 = np.einsum("ki,ij,kj->k", reps, Gs, reps) / 4
    F2c = m.f2(reps, s2, None, np.zeros(m.n_atoms))
    ref = ph.m_F2 / 2
    assert np.max(np.abs(F2c[inv] - ref) / (ref + 1e-9 * ref.max())) < 1e-9               # classes reproduce the full P1 sum
    # a symmetry-allowed random shift keeps the space group: members of one class keep the same |F|²
    rng = np.random.default_rng(0)
    X = m.atom_positions(rng.normal(0, 0.01, m.n_pos))
    Z = np.array([RV.Phase.from_cif(FAP_CIF).elems.count(e) for e in m.atom_el], float) * 0 + 1.0
    members = np.where(inv == inv[5])[0][:5]
    vals = [abs((Z * np.exp(2j * np.pi * (X @ ph.m_hkl[i]))).sum()) ** 2 for i in members]
    assert np.allclose(vals, vals[0], rtol=1e-9)


@needs_data
def test_real_fluorapatite_structure_refinement_agrees_with_gsas():
    """Real laboratory pattern (GSAS-II tutorial): fully automatic refinement of atom positions + B from the COD structure. GSAS reference (same data):
    Ca2 (0.241976, 0.992603, 1/4), P (0.397416, 0.367704, 1/4), O1 (0.325053, 0.484763, 1/4), O2 (0.591494, 0.469954, 1/4), O3 (0.339510, 0.258126, 0.070641);
    the COD setting has x and y exchanged."""
    x, y = _fap_pattern()
    r0 = RV.quick_refine(x, y, [FAP_CIF], radiation=RAD)
    r = RV.quick_refine(x, y, [FAP_CIF], radiation=RAD, refine_structure=True)
    assert r["stats"]["Rwp"] < r0["stats"]["Rwp"] - 0.003 and r["phases"][0]["Rb"] < r0["phases"][0]["Rb"] and np.sqrt(r["stats"]["chi2"]) < 1.7
    st = r["phases"][0]["structure"]
    assert st["bond_rms_deviation"] < 0.04 and st["bond_max_deviation"] < 0.08                  # chemically sensible geometry
    at = {a["label"]: a for a in st["atoms"]}
    gsas = {"Ca2": (0.241976, 0.992603, 0.25), "P1": (0.397416, 0.367704, 0.25), "O1": (0.325053, 0.484763, 0.25), "O2": (0.591494, 0.469954, 0.25),
            "O3": (0.339510, 0.258126, 0.070641)}
    for lab, (gx, gy, gz) in gsas.items():
        a = at[lab]
        assert abs(a["x"] - gy) < 0.002 and abs(a["y"] - gx) < 0.002 and abs(a["z"] - gz) < 0.002, (lab, a["x"], a["y"], a["z"])
        assert a["x_esd"] > 0 or a["y_esd"] > 0
    assert 0.1 < at["Ca1"]["B"] < 1.5 and at["F1"]["x_esd"] == 0 and at["F1"]["free_coordinates"] == 0


def test_synthetic_positions_are_recovered_and_restraints_keep_geometry():
    cif = FAP_CIF if os.path.exists(FAP_CIF) else None
    if cif is None:
        pytest.skip("needs the fluorapatite CIF")
    x = np.arange(15, 90, 0.02)
    ph = RV.Phase.from_cif(cif).enable_structure()
    R = RV.Refinement(x, np.ones_like(x) * 1000.0, [ph], RAD, 6, "rietveld", counts=True)
    for k, v in {"c0.a": 9.372, "c0.c": 6.886, "W": 0.004, "X": 0.06}.items():
        R.P[k].value = v
    rng = np.random.default_rng(3)
    true = rng.normal(0, 0.004, ph.model.n_pos)                                        # the "real" structure differs from the CIF
    for k, t in enumerate(true):
        R.P[f"sx0.{k}"].value = float(t)
    for e in ph.model.elements:
        R.P[f"sB0.{e}"].value = 0.6
    yc = R._design(R._values())[:, 0]
    y = np.random.default_rng(4).poisson(yc / yc.max() * 30000 + 200).astype(float)
    r = RV.quick_refine(x, y, [cif], radiation=RAD, refine_structure=True, offset="zero")
    a = r["phases"][0]["structure"]["atoms"]
    # compare through the model: recovered shifts of the free coordinates
    ph2 = RV.Phase.from_cif(cif).enable_structure()
    m2 = ph2.model
    got = {a_["label"]: np.array([a_["x"], a_["y"], a_["z"]]) for a_ in a}
    want = m2.site_frac.copy()
    for (rr, c), v in zip(m2.pos_params, true):
        want[rr] = want[rr] + m2.basis[rr][:, c] * v
    labs = m2.labels()
    for rr in m2.reps:
        assert np.abs(got[labs[rr]] - want[rr]).max() < 0.002, (labs[rr], got[labs[rr]], want[rr])
    assert r["stats"]["chi2"] < 1.6
    assert r["phases"][0]["structure"]["bond_rms_deviation"] < 0.08
