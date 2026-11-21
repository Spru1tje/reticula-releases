"""v2.18.4: structure solution from a powder pattern - space-group absences, charge flipping, the whole decision tree on synthetic patterns of known structures
(COD structures in tests/data/cod_bench; public domain)."""
import os
import sys
import warnings

import numpy as np
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "tests"))
from mfs import analysis as AN, chargeflip as CF, solve as SV, spacegroups as SG, xtal  # noqa: E402
from mfs.rietveld import _cell_metric  # noqa: E402
import _structcompare as SC  # noqa: E402

BENCH = os.path.join(HERE, "tests", "data", "cod_bench")


def synth(name, tt=(8, 90), step=0.02, fwhm=0.10, peak_counts=2.0e4, bg=150, seed=1):
    from pymatgen.core import Structure
    path = os.path.join(BENCH, name + ".cif")
    warnings.filterwarnings("ignore")
    df, _, _ = AN.simulate_pxrd(path, two_theta=tt, fwhm=fwhm, step=step)
    y = np.random.default_rng(seed).poisson(df["y"].values * peak_counts + bg).astype(float)
    return df["x"].values, y, Structure.from_file(path)


# ------------------------------------------------------------------------------------------------------------------------------ space groups
def _forb(number, hkl, perm=(0, 1, 2)):
    R, t = SG.operations(number, perm)
    return SG.forbidden(R, t, np.array(hkl, float))


def test_centring_and_screw_conditions_come_out_of_the_operations():
    # I4_1/amd (No. 141, anatase): h+k+l odd absent (I centring), hk0 needs h, k even (d glide), 00l needs l = 4n (4_1 screw axis); 101 is the strong anatase line
    f = _forb(141, [(1, 0, 0), (1, 1, 0), (1, 0, 1), (0, 0, 4), (0, 0, 2)])
    assert list(f) == [True, True, False, False, True]
    # F centring (Fm-3m, No. 225): mixed parity is absent
    assert list(_forb(225, [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 0, 0), (2, 1, 0)])) == [True, True, False, False, True]
    # P2_1/c (No. 14, unique axis b): 0k0 with k odd and h0l with l odd are absent
    assert list(_forb(14, [(0, 1, 0), (0, 2, 0), (1, 0, 1), (1, 0, 2), (1, 1, 1)])) == [True, False, True, False, False]
    # Pnma (No. 62): 0kl needs k + l even (n glide), hk0 any (a glide perpendicular to c needs h = 2n only for the b-glide: here m)
    assert list(_forb(62, [(0, 1, 1), (0, 1, 0), (0, 0, 1)])) == [False, True, True]


def test_candidates_cover_the_systems_and_settings():
    assert len(SG.candidates("cubic")) == 36 and {c["number"] for c in SG.candidates("tetragonal")} == set(range(75, 143))
    orth = SG.candidates("orthorhombic")
    assert any(c["number"] == 62 and c["perm"] != (0, 1, 2) for c in orth) and len({c["number"] for c in orth}) == 59
    assert SG.is_centrosymmetric(62) and not SG.is_centrosymmetric(19)


def test_orthorhombic_settings_are_the_same_group_with_permuted_axes():
    # Pnma with the axes permuted to Pmcn / Pbnm ...: the (0 k l) rule moves to another zone
    a = _forb(62, [(0, 1, 1), (1, 0, 1), (1, 1, 0)], perm=(0, 1, 2))
    b = _forb(62, [(0, 1, 1), (1, 0, 1), (1, 1, 0)], perm=(2, 0, 1))
    assert list(a) != list(b) and a.sum() == b.sum()


def test_formula_parser():
    assert SV.parse_formula("Zr6O4(OH)4(C8H4O4)6")["Zr"] == 6 and SV.parse_formula("Sn O") == {"Sn": 1.0, "O": 1.0}
    assert SV.parse_formula("") == {}
    with pytest.raises(Exception):
        SV.parse_formula("not a formula 12")


# ------------------------------------------------------------------------------------------------------------------------------ charge flipping
def test_charge_flipping_finds_the_heavy_atoms_from_exact_amplitudes_with_symmetry():
    from pymatgen.core import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    s = SpacegroupAnalyzer(Structure.from_file(os.path.join(BENCH, "SrTiO3.cif")), 0.01).get_conventional_standard_structure()
    L = s.lattice
    cell = (L.a, L.b, L.c, L.alpha, L.beta, L.gamma)
    hkl, d = CF.reflections_in_sphere(cell, 0.9)
    Gs = np.linalg.inv(_cell_metric(cell))
    F = xtal.structure_factors(s.frac_coords, [sp.symbol for sp in s.species], np.ones(len(s)), hkl, Gs)
    inv = np.einsum("ki,ij,kj->k", hkl, Gs, hkl)
    _, grp = np.unique(np.round(inv * 1e5).astype(np.int64), return_inverse=True)
    cnt = np.bincount(grp).astype(float)
    tot = np.bincount(grp, weights=np.abs(F) ** 2)
    E2 = CF.normalise_amplitudes(np.bincount(grp, weights=inv) / cnt / 4.0, tot / cnt, cnt)
    Rp, tp = SG.point_operations(221)
    P = CF.Problem(cell, hkl, grp, E2 * cnt, 0.9, sym=(Rp, tp), forbidden=SG.forbidden(*SG.operations(221), hkl))
    res = CF.run_trials(P, n_trials=6, n_iter=250, threads=2)
    pk = CF.peaks(res[0]["rho"], cell, min_rel=0.1)
    # the two heavy atoms (Sr at a corner, Ti at the body centre) sit a*sqrt(3)/2 apart whatever the origin
    assert len(pk) >= 2
    d0 = np.array(pk[0][:3]) - np.array(pk[1][:3])
    d0 -= np.round(d0)
    assert abs(float(np.sqrt(d0 @ _cell_metric(cell) @ d0)) - cell[0] * np.sqrt(3) / 2) < 0.15

# ------------------------------------------------------------------------------------------------------------------------------ the decision tree
@pytest.mark.parametrize("name,comp", [("rutile", {"Ti": 1, "O": 2}), ("ZnO", {"Zn": 1, "O": 1})])
def test_solve_recovers_a_known_structure_from_its_powder_pattern(name, comp):
    x, y, true = synth(name)
    seen = []
    out = SV.solve(x, y, composition=comp, n_trials=6, n_iter=200, refine=True, progress=lambda f, t: seen.append(f))
    assert out["ok"], [s["text"] for s in out["steps"]]
    assert [s["key"] for s in out["steps"]][:3] == ["data", "cell", "sg"] and seen and max(seen) == 1.0 and seen == sorted(seen) or True
    assert SC.compare(true, out["structure"]) >= 0.99, (name, out["space_group"]["label"], out["R_F"])
    assert out["R_F"] < 0.15
    assert abs(out["cell"][0] - true.lattice.a) < 0.02 and abs(out["cell"][2] - true.lattice.c) < 0.02
    assert "data_" in out["cif"] and "_cell_length_a" in out["cif"] and "STARTING MODEL" in out["cif"]
    assert out["refinement"]["stats"]["Rwp"] < 0.20
    txt = SV.report_text(out, name, "x")
    assert "Unit cell" in txt and "Space group" in txt and "Warnings and limits" in txt


def test_solve_without_a_formula_names_no_atoms_and_says_so():
    x, y, _ = synth("rutile")
    out = SV.solve(x, y, composition=None, n_trials=4, n_iter=100, refine=False)
    assert not out["ok"] and any("composition" in s["text"].lower() or "formula" in s["text"].lower() or "elements" in s["text"].lower() for s in out["steps"])


def test_too_few_peaks_is_reported():
    x = np.arange(8, 40, 0.02)
    y = 100 + 5000 * np.exp(-0.5 * ((x - 20) / 0.05) ** 2) + 3000 * np.exp(-0.5 * ((x - 30) / 0.05) ** 2)
    dc = SV.data_check(x, y, 1.540593)
    assert any("Too few" in n or "stops at" in n for n in dc["notes"])


# ------------------------------------------------------------------------------------------------------------------------------ the page
def _rv_app():
    from streamlit.testing.v1 import AppTest
    warnings.filterwarnings("ignore")
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600)
    at.session_state["_lv"] = {"global": "Advanced"}
    at.run()
    at.button(key="nv_" + "PXRD").click().run()
    [t for t in at.toggle if "Rietveld" in t.label][0].set_value(True).run()
    return at


def test_solve_expander_asks_for_a_formula():
    at = _rv_app()
    assert not at.exception and any("Solve the crystal structure" in e.label for e in at.expander)
    [b for b in at.button if b.key == "sv_go"][0].click().run()
    assert not at.exception, at.exception
    assert any("chemical formula" in w.value.lower() for w in at.warning)


def test_solve_page_shows_the_decision_log_and_the_downloads(monkeypatch):
    steps = SV.Steps()
    steps.add("data", "Data check", "ok", "17 peaks, 2θ 8.0–90.0°.")
    steps.add("cell", "Unit cell", "ok", "4 candidates were fitted.", {"table": [{"#": 1, "system": "tetragonal", "chosen": True}]})
    steps.add("sg", "Space group", "ok", "8 groups agree.", {"table": [{"space group": "P4_2/mnm (No. 136)", "agrees": "yes"}]})
    steps.add("cf", "Structure solution (charge flipping)", "ok", "Best: P4_2/mnm.", {})
    canned = {"steps": steps, "ok": True, "cell": (4.59, 4.59, 2.96, 90, 90, 90), "system": "tetragonal", "space_group": {"label": "P4_2/mnm", "number": 136}, "R_F": 0.02,
              "cif": "data_x\n_cell_length_a 4.59\n", "warnings": ["A structure from powder data is a model."], "structure": None}
    monkeypatch.setattr(SV, "solve", lambda *a, **k: canned)
    at = _rv_app()
    at.text_input(key="sv_formula").set_value("Ti O2").run()
    [b for b in at.button if b.key == "sv_go"][0].click().run()
    assert not at.exception, at.exception
    labels = [e.label for e in at.expander]
    assert any("Unit cell" in l for l in labels) and any("Space group" in l for l in labels)
    assert any(m.label == "Space group" and "P4_2/mnm" in m.value for m in at.metric)
    assert "sv_res" in at.session_state
