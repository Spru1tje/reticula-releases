"""v2.18.4: the Lázaro 2020 multifaceted TGA method (mfs.lazaro). The expected values are the tables of the paper's Supporting Information
(Eur. J. Inorg. Chem. 2020, 4284-4294): S9 / S10 (pristine and defective UiO-66 with formate), S14 / S15 (multivariate with a partly decomposing linker),
S22 / S23 (inert TiO2 composite), S24 / S25 (a polymer that burns away before the linker)."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mfs import lazaro as LZ  # noqa: E402
from mfs.defects import DefectInput, Species, run_all  # noqa: E402


def _uio(W200=None, Wplat=None, Wend=100.0, fa=None, extra_species=(), **kw):
    sp = [Species("BDC (terephthalate)", "C8H4O4", 2, "linker", 4, 4.0, True, True)]
    if fa is not None:
        sp.append(Species("Formate", "CHO2", 1, "cap", 1, fa, False, False))
    sp += list(extra_species)
    return DefectInput(sample="t", framework="UiO-66", metal="Zr", species=sp, tga_W_200=W200, tga_W_plateau=Wplat, tga_W_end=Wend, **kw)


def _i(inp):
    return {m.key: m for m in LZ.method_i(inp)}


# Table S9: %MOF (start structure, residue = 100 %) and %DH for Zr6O4(OH)4(L)x(FA)y with x = 6, 5, 4, 3 and the NMR ratio FA/L
S9 = [(6, 0.0, 225.10, 220.22), (5, 0.4, 215.08, 200.19), (4, 1.0, 205.06, 180.15), (3, 2.0, 195.04, 160.11)]


@pytest.mark.parametrize("L,nmr,pct_mof,pct_dh", S9)
def test_pristine_and_defective_uio66_start_structure_table_S9(L, nmr, pct_mof, pct_dh):
    r = _i(_uio(W200=pct_mof, Wplat=pct_dh, fa=nmr))
    assert r["i1"].L == pytest.approx(L, abs=0.02)
    assert (r["i1"].m or 0) == pytest.approx(L * nmr, abs=0.04)
    assert r["i1"].extra["R_calc"] == pytest.approx(pct_mof / 100, abs=1e-9)         # the solution reproduces the measured ratio
    assert abs(r["i1"].extra["charge_left_per_node"]) < 1e-9                          # balanced by the balancing species


@pytest.mark.parametrize("L,nmr,pct_mof,pct_dh", S9)
def test_dehydroxylated_structure_table_S9(L, nmr, pct_mof, pct_dh):
    r = _i(_uio(W200=pct_mof, Wplat=pct_dh, fa=nmr))
    assert r["i2"].L == pytest.approx(L, abs=0.02)                                    # exact, also for the strongly defective samples (older methods: up to 21 % error)


@pytest.mark.parametrize("L,nmr,pct_mof,pct_dh", [(6, 0.0, 225.10, 220.22), (5, 0.4, 215.08, 198.02), (4, 1.0, 205.06, 175.82), (3, 2.0, 195.04, 153.62), (2, 4.0, 185.0, 131.42)])
def test_dehydroxylated_structure_without_oxygen_compensation_table_S10(L, nmr, pct_mof, pct_dh):
    """ZrO(L)x DH structures (Table S10): exact with this method, while the older plateau method (a) is up to 21.6 % off (Table S7: L = 1.57 for L = 2)."""
    from mfs.defects import method_a
    inp = _uio(W200=pct_mof, Wplat=pct_dh, fa=nmr, lazaro_dh_compensation="none")
    assert _i(inp)["i2"].L == pytest.approx(L, abs=0.02)
    if L == 2:
        assert method_a(inp).L == pytest.approx(1.568, abs=0.01)          # the error the paper describes for method 2 / this convention


def test_round_trip_is_exact_with_chloride_and_hydroxide():
    """Forward model -> solve: the amounts come back exactly (charge balance closed by OH-, Cl- measured)."""
    from mfs.chem import M_CL, molar_mass
    nM = 6
    inp0 = _uio(fa=0.5)
    Mres = inp0.MMO2(); Msbu = inp0.M_clust() / nM; q = inp0.Q() / nM
    lk = [LZ.Part("BDC (terephthalate)", molar_mass("C8H4O4"), 2, 1.0)]
    md = [LZ.Part("Formate", molar_mass("CHO2"), 1, 0.5)]
    cl = [(LZ.Part("Cl-", M_CL, 1.0), 0.15)]
    bal = LZ.Part("OH-", 17.007, 1.0)
    x = 0.61
    R = LZ.forward_R(x, Mres, Msbu, q, lk, md, cl, bal)
    s = LZ.solve_per_metal(R, Mres, Msbu, q, lk, md, cl, bal)
    assert s["x"] == pytest.approx(x, abs=1e-12) and s["R_calc"] == pytest.approx(R, abs=1e-12)
    # through the method (Cl per node = 6 * 0.15)
    inp = _uio(W200=100 * R, Wplat=None, fa=0.5, cl_per_node=0.9)
    r = _i(inp)["i1"]
    assert r.L == pytest.approx(6 * x, abs=1e-6) and r.Cl == pytest.approx(0.9, abs=1e-9)


def test_multivariate_mof_with_a_partly_decomposing_linker_table_S14():
    """L1 = BDC, L2 = NH2-BDC (50 % stays in the residue), L3 = NO2-BDC; DH structures Zr6O6(L1)4(L2)1(L3)1 (%DH 203.66) and (L1)3(L2)1.5(L3)1.5 (196.66)."""
    for (L1, a2, pct) in ((4, 0.25, 203.66), (3, 0.5, 196.66), (2, 1.0, 190.34), (1, 2.5, 184.61)):
        extra = [Species("NH2-BDC", "C8H5NO4", 2, "linker", 3, 3 * a2, True, False, 0.5), Species("NO2-BDC", "C8H3NO6", 2, "linker", 3, 3 * a2, True, False)]
        r = _i(_uio(Wplat=pct, extra_species=extra))["i2"]
        assert r.L == pytest.approx(L1, abs=0.03), (L1, r.L)
        assert r.extra["R_calc"] == pytest.approx(pct / 100, abs=1e-9)


def test_inert_composite_tio2_at_uio66_table_S23():
    """TiO2@UiO-66: Zr6O4(OH)4(L)5(FA)2 with 17.8 % TiO2 / 82.2 % ZrO2 in the residue; %DH of the composite 182.4."""
    r = _i(_uio(Wplat=182.4, fa=0.4, mof_residue_fraction=0.822, composite_kind="inert"))["i2"]
    assert r.L == pytest.approx(5.0, abs=0.03)
    r3 = _i(_uio(Wplat=100 * (0.5 * 2.0019 + 0.5), fa=0.4, mof_residue_fraction=0.5, composite_kind="inert"))["i2"]     # 50 % of the residue is TiO2
    assert r3.L == pytest.approx(5.0, abs=0.03)


def test_polymer_that_burns_before_the_linker_table_S24_S25():
    """PVAc@UiO-66 Zr6O4(OH)4(L)5(FA)2 (molar fraction 0.75 MOF): the DH mass is the MOF alone; the start ratio gives the polymer weight fraction."""
    M_mof, M_pvac, S = 1590.092, 86.09, 0.75
    R_start = (M_mof * S + M_pvac * (1 - S)) / (739.34 * S)                          # = 2.190 (219.0 %)
    r = _i(_uio(W200=100 * R_start, Wplat=200.19, fa=0.4, composite_kind="before_dh"))
    assert r["i1"].status == "n/a" and "before the linker" in r["i1"].message
    assert r["i2"].L == pytest.approx(5.0, abs=0.03)
    w_true = (1 - S) * M_pvac / (M_mof * S + M_pvac * (1 - S))                        # mass fraction of the polymer
    assert r["i2"].extra["composite_wt_pct"] == pytest.approx(100 * w_true, abs=0.15)


def test_composite_burning_with_the_linker_uses_its_known_weight_fraction():
    # a material (30 wt %) that burns together with the linker: the MOF mass at the plateau is 70 % of the measured mass
    inp0 = _i(_uio(Wplat=200.19, fa=0.4))["i2"]
    r = _i(_uio(Wplat=200.19 / 0.7, fa=0.4, composite_kind="with_dh", composite_wt_pct=30.0))["i2"]
    assert r.L == pytest.approx(inp0.L, abs=1e-9) and r.L == pytest.approx(5.0, abs=0.03)


def test_methods_are_reported_as_not_available_without_data_and_run_all_keeps_old_results():
    r = _i(_uio(fa=0.4))
    assert r["i1"].status == "n/a" and r["i2"].status == "n/a" and r["i1"].message
    res = run_all(_uio(W200=215.08, Wplat=200.19, fa=0.4))
    assert "lazaro" not in res                                        # default: the earlier methods and the consensus are unchanged
    res2 = run_all(_uio(W200=215.08, Wplat=200.19, fa=0.4), lazaro=True)
    assert [m.key for m in res2["lazaro"]] == ["i1", "i2"]
    assert res2["consensus"] == res["consensus"]


def test_other_nodes_work_through_the_node_library():
    """MOF-5 (Zn4O): charge per metal 1.5, SBU ZnO1/4; exact round trip."""
    from mfs.chem import molar_mass
    inp = DefectInput(sample="m5", framework="MOF-5", metal="Zn", species=[Species("BDC (terephthalate)", "C8H4O4", 2, "linker", 4, 4.0, True, True)],
                      tga_W_end=100.0)
    nM = inp.nM(); Mres = inp.MMO2(); Msbu = inp.M_clust() / nM; q = inp.Q() / nM
    assert q == pytest.approx(1.5) and nM == 4
    lk = [LZ.Part("BDC (terephthalate)", molar_mass("C8H4O4"), 2, 1.0)]
    x = 0.7
    R = LZ.forward_R(x, Mres, Msbu, q, lk, (), (), LZ.Part("OH-", 17.007, 1.0))
    inp.tga_W_200 = 100 * R
    r = _i(inp)["i1"]
    assert r.L == pytest.approx(nM * x, abs=1e-6)
