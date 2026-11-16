"""Batch E: Structure (CIF) — user's molecular CIF (simvastatin, P2₁) and a MOF CIF (UiO-66)."""
import os, io, zipfile
import numpy as np
import pytest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mfs import cif as CX, cifplot as CP, cifsi as CS, cifpxrd as CPX, poretools as PT

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
USER = os.environ.get("MFS_USER_CIF", "/workspace/real_cif/user.cif")
UIO = os.path.join(HERE, "references", "UiO-66_Zr_COD4512074.cif")
need_user = pytest.mark.skipif(not os.path.exists(USER), reason="user CIF not on this machine")


@need_user
def test_user_cif_numbers_and_tables():
    S = CX.read(USER)[0]
    t = dict(CS.crystal_table(S).values)
    assert t["Space group"].replace(" ", "") in ("P1211", "P2_1", "P21")
    assert t["Z"] in ("2", 2)
    assert "0.0440" in t["Final R indexes [I > 2σ(I)]"] and "0.1218" in t["Final R indexes [I > 2σ(I)]"]
    assert "0.0458" in t["Final R indexes [all data]"]
    assert "-0.03(5)" in t["Flack parameter"].replace("−", "-")
    assert "1.079" in str(t["Goodness-of-fit on F²"])
    assert abs(S.V - 2256.96) < 0.1 and abs(S.density() - 1.232) < 0.002
    g = CX.geom_tables(S)
    assert len(g["bonds"]) == 64 and len(g["angles"]) == 94 and len(g["torsions"]) == 118
    assert any("(" in str(v) for v in g["bonds"].iloc[:, 1])          # ESDs present
    hb = CX.hbonds(S)
    assert len(hb) == 2 and all(abs(r["D···A (Å)"] - 2.905) < 0.01 for r in hb)
    assert {r["symmetry (A)"] for r in hb} and all(r["symmetry (A)"] for r in hb)
    secs = CS.si_sections(S, g, hb)
    d = CS.docx(secs); x = CS.xlsx(secs)
    assert zipfile.ZipFile(io.BytesIO(d)).namelist() and len(d) > 3000
    import openpyxl
    wb = openpyxl.load_workbook(io.BytesIO(x)); assert len(wb.sheetnames) == 5


@need_user
def test_user_cif_render_and_pxrd():
    S = CX.read(USER)[0]
    for mode, style in [("molecules", "ortep"), ("packing", "ball"), ("asu", "space"), ("supercell", "stick")]:
        sites, bonds = S.grow(mode, supercell=(2, 1, 1))
        assert sites and (bonds or style == "space")
        fig, ax = plt.subplots()
        CP.draw(ax, S, sites, bonds, style=style, prob=0.5, hbonds=CX.hbond_pairs(S, sites), view=CP.view_matrix(S, "a"))
        b = io.BytesIO(); fig.savefig(b, format="png"); plt.close(fig); assert b.tell() > 5000
    mol, mb = S.grow("molecules"); assert len(CX._components(mol, mb)) == 2   # Z' = 2
    pk, _ = S.grow("packing"); assert len(CX.hbond_pairs(S, pk)) >= 1
    assert abs(CP.prob_scale(0.5) - 1.5382) < 1e-3
    sim, st_ = CPX.simulate(open(USER, "rb").read(), (5, 30))
    top = st_.nlargest(1, "I").iloc[0]
    assert abs(top.two_theta - 9.32) < 0.05
    fr = CPX.fit_to(sim.x, 3 * sim.y + 10, sim); assert fr["Rwp"] < 1e-6 and abs(fr["scale"] - 3) < 1e-6


def test_mof_cif_render_and_pores():
    S = CX.read(UIO)[0]
    sites, bonds = S.grow("cell", part="major")
    fig, ax = plt.subplots()
    CP.draw(ax, S, sites, bonds, style="ball", polyhedra=True, view=CP.view_matrix(S, "uvw", uvw=[1, 1, 0]))
    plt.close(fig)
    r = PT.analyse(S, 0.5, n_mc=60)
    assert 6.5 < r["Di (Å)"] < 9.5 and 2.0 < r["Df (Å)"] < 5.0
    assert 0.45 < r["void fraction (probe 0)"] < 0.75 and 600 < r["ASA (m² g⁻¹)"] < 2000
    assert "Monte Carlo" in r["method"]


def test_percolation_simple_channel():
    m = np.zeros((10, 10, 10), bool); m[4:6, 4:6, :] = True
    assert PT._percolates(m)
    m[:, :, 5] = False
    assert not PT._percolates(m)
