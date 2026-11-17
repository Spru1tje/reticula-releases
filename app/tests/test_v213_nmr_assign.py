"""v2.13 M3/M4: ¹H prediction (nmrshiftdb2 HOSE table), assignment, impurity check, qNMR, labelled figure, SI page."""
import sys, pathlib
import numpy as np
import pytest
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from mfs import nmrproc as N, nmrpeaks as P, nmrpredict as HP, nmrassign as A

D = pathlib.Path(__file__).parent / "data" / "nmr"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


def _aspirin():
    fn = "aspirin_300_XWIN35.zip"
    sp = N.process(N.load((D / fn).read_bytes(), fn)[0])
    ms = P.analyse(sp); P.normalise(ms)
    return sp, ms


def test_prediction_heldout_accuracy():
    r = HP.validate(limit=400)
    assert r["MAE"] < 0.3 and r["within_0.3"] > 0.7, r


def test_prediction_simple_molecules():
    p = {x["label"]: x["shift"] for x in HP.predict(HP.mol_from_any("CCOC(C)=O"))}   # ethyl acetate (CDCl3: 1.26, 4.12, 2.05)
    assert abs(p[1] - 1.26) < 0.15 and abs(p[2] - 4.12) < 0.15 and abs(p[5] - 2.05) < 0.15
    p = {x["label"]: x["shift"] for x in HP.predict(HP.mol_from_any("c1ccccc1"))}
    assert abs(p[1] - 7.34) < 0.15


def test_groups_symmetry():
    G = A.groups(HP.mol_from_any("Cc1ccc(C)cc1"))       # p-xylene: CH3 (6H) and ArH (4H)
    assert sorted(g["nH"] for g in G) == [4, 6]


def test_aspirin_auto_assignment():
    sp, ms = _aspirin()
    mol = HP.mol_from_any(ASPIRIN); G = A.groups(mol); a = A.assign(G, ms)
    at = {tuple(g["labels"]): (None if a[g["gid"]] is None else round(ms[a[g["gid"]]].centre, 1)) for g in G}
    assert at[(1,)] == 2.3 and at[(9,)] == 8.0 and at[(7,)] == 7.5 and at[(6,)] == 7.0 and at[(8,)] == 7.3, at


def test_fixed_assignment_respected():
    sp, ms = _aspirin()
    mol = HP.mol_from_any(ASPIRIN); G = A.groups(mol)
    g6 = next(g["gid"] for g in G if g["labels"] == [6])
    m73 = next(i for i, m in enumerate(ms) if abs(m.centre - 7.27) < 0.03)
    a = A.assign(G, ms, fixed={g6: m73})
    assert a[g6] == m73
    a = A.assign(G, ms, fixed={g6: None})
    assert a[g6] is None


def test_impurity_check_synthetic():
    mk = lambda c, mult, nH: P.Multiplet(lo=c - 0.02, hi=c + 0.02, centre=c, integral=nH, nH=nH, mult=mult, J=[7.1] if mult in "tq" else [])
    ms = [mk(4.12, "q", 0.2), mk(2.05, "s", 0.3), mk(1.26, "t", 0.3), mk(5.30, "s", 0.1), mk(3.10, "s", 0.05)]
    r = A.impurity_check(ms, list(range(len(ms))), "CDCl3")
    names = {b["name"]: b for b in r["best"]}
    assert "ethyl acetate" in names and set(names["ethyl acetate"]["peaks"]) == {0, 1, 2}
    assert abs(names["ethyl acetate"]["mol_ratio"] - 0.1) < 0.01
    assert "dichloromethane" in names
    assert 4 in r["unexplained"]


def test_qnmr_formula():
    from mfs import nmr_page as NPG
    # 10 mg analyte (M 180.16, 3H signal, I 3.0) vs 5 mg maleic acid (M 116.07, 2H, I 2.0, 99.9 %)
    p = NPG.qnmr_purity(3.0, 3, 180.16, 10.0, 2.0, 2, 116.07, 5.0, 99.9)
    assert abs(p - (1.0 * 180.16 / 116.07 * 0.5 * 99.9)) < 1e-9


def test_figure_and_si_pdf():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from mfs import nmrfig as F, nmr_page as NPG
    sp, ms = _aspirin()
    mol = HP.mol_from_any(ASPIRIN); G = A.groups(mol); A.assign(G, ms)
    fig, ax = plt.subplots(); F.draw(ax, sp, ms, mol, G, {}); plt.close(fig)
    acs = P.acs_string(ms, mhz=300, solvent="CDCl3")
    assert "(br s, 1H, Ha)" in acs                    # COOH (most downfield group = Ha): nH from the assignment, not the broad integral
    pdf = NPG.si_pdf(sp, ms, mol, G, {}, acs, "aspirin", "CDCl3", 300, 98.7)
    assert pdf[:4] == b"%PDF" and len(pdf) > 10000


def test_component_payload_zones():
    from mfs import nmrui as UI
    sp, ms = _aspirin()
    mol = HP.mol_from_any(ASPIRIN); G = A.groups(mol); a = A.assign(G, ms)
    p = UI.payload(sp, ms, mol, G, a)
    assert len(p["groups"]) == len(G) and len(p["zones"]) >= 6
    assert all(0 <= z["x0"] < z["x1"] <= p["sw"] for z in p["zones"])
    assert UI.handle({}, {"op": "assign", "gid": 2, "mi": 4}) == {2: 4}
    assert UI.handle({2: 4}, {"op": "unassign", "gid": 3}) == {2: 4, 3: None}
    assert UI.handle({2: 4}, {"op": "reset"}) == {}
