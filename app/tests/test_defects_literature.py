"""Defect calculator vs. worked examples from the literature (and exact round-trip checks).

Each test docstring states the source and whether the published number is reproduced."""
import math
import pytest
from mfs.chem import molar_mass, oxide_mass
from mfs.defects import DefectInput, Species, method_a, method_b, method_c, method_d, method_e, method_g, method_f, run_all, open_sites
from mfs.chem import SPECIES


def sp(name, integral=None, at_plateau=False, main=False, **kw):
    f, z, role, nH, _ = SPECIES[name]
    d = dict(name=name, formula=f, z=z, role=role, nH=nH, integral=integral, at_plateau=at_plateau, main=main)
    d.update(kw)
    return Species(**d)


# ---------------------------------------------------------------- Valenzano 2011 / Shearer 2014-2016
def test_ideal_uio66_plateau_220():
    """Shearer et al. 2014/2016: ideal dehydroxylated UiO-66, Zr6O6(BDC)6, plateau = 220.2 % of the ZrO2 end mass.
    PUBLISHED: 220.2 % -> reproduced."""
    inp = DefectInput(framework="UiO-66", metal="Zr", species=[sp("BDC (terephthalate)", main=True)],
                      tga_W_plateau=220.2, tga_W_end=100)
    r = method_a(inp)
    assert r.extra["W_plat_ideal_%"] == pytest.approx(220.2, abs=0.05)
    assert r.L == pytest.approx(6.0, abs=0.01)


def test_valenzano_dehydroxylated_mass():
    """Valenzano 2011: Zr6O6(BDC)6 (dehydroxylated) M = 1628 g/mol, Zr6O4(OH)4(BDC)6 M = 1664 g/mol.
    Reproduced to 0.1 g/mol (IUPAC atomic weights)."""
    assert molar_mass("Zr6O6(C8H4O4)6") == pytest.approx(1628.0, abs=0.2)
    assert molar_mass("Zr6O4(OH)4(C8H4O4)6") == pytest.approx(1664.1, abs=0.1)


def test_defective_uio66_example():
    """Shearer-type calculation: plateau at 200.0 % (ZrO2 = 100 %) -> L = (2.000-1)*739.33/148.12 = 4.99."""
    inp = DefectInput(framework="UiO-66", metal="Zr", species=[sp("BDC (terephthalate)", main=True)],
                      tga_W_plateau=200.0, tga_W_end=100)
    r = method_a(inp)
    assert r.L == pytest.approx(739.332 / 148.121, rel=1e-3)
    assert r.x == pytest.approx(6 - 4.991, abs=0.01)


def test_hf_and_mof801_ideal_plateaus():
    """Internal consistency for Hf-UiO-66 (HfO2 end product) and MOF-801 (fumarate)."""
    hf = DefectInput(framework="UiO-66", metal="Hf", species=[sp("BDC (terephthalate)", main=True)], tga_W_plateau=100, tga_W_end=100)
    assert method_a(hf).extra["W_plat_ideal_%"] == pytest.approx(100 * (1 + 6 * 148.121 / (6 * 210.488)), rel=1e-4)
    m801 = DefectInput(framework="MOF-801", metal="Zr", species=[sp("Fumarate", main=True)], tga_W_plateau=100, tga_W_end=100)
    assert method_a(m801).extra["W_plat_ideal_%"] == pytest.approx(179.6, abs=0.1)


# ---------------------------------------------------------------- De Roo 2025 (Adv. Sci.) MOF-808 worked example
def _deroo_mof808():
    # [std] = 2 mM TMSP-d4 (integral 9), [MOF] = 1.5 mM in 0.5 mL -> 0.75 µmol digested (as in the paper)
    Mexp = 100 / 45.6 * 6 * oxide_mass(0)
    return DefectInput(sample="MOF-808 (De Roo 2025)", framework="MOF-808", metal="Zr",
                       species=[sp("BTC (trimesate)", integral=4.5, main=True), sp("Acetate", integral=4.5),
                                sp("Methanol (in pores)", integral=0.6)],
                       tga_W_initial=100, tga_W_end=45.6, nmr_std_conc_mM=2.0, nmr_std_integral=9, nmr_std_nH=9,
                       nmr_sample_mg=0.75e-3 * Mexp, nmr_volume_mL=0.5, cl_conc_mM=4.2)


def test_deroo_mof808_experimental_molar_mass():
    """De Roo 2025, eq. 4: W_ZrO2 = 45.6 % -> M_exp = 1619.5 g/mol (paper). We get 1621.3 (the paper's 1619.5
    corresponds to W = 45.65 %, i.e. rounding of the residue). Reproduced within 0.15 %."""
    r = method_d(_deroo_mof808())
    assert r.extra["M_exp"] == pytest.approx(1619.5, rel=0.0015)


def test_deroo_mof808_minimal_formula():
    """De Roo 2025, eqs. 5-9: x(BTC) = 2, y(OAc) = 2, z(MeOH) = 0.3, i(Cl) = 2.8, w(OH) = 1.2,
    M_anhy = 1341 g/mol, n(H2O) = 15.5. All reproduced."""
    r = method_d(_deroo_mof808())
    assert r.L == pytest.approx(2.0, abs=1e-6)
    assert r.m == pytest.approx(2.0, abs=1e-6)
    assert r.Cl == pytest.approx(2.8, abs=1e-6)
    assert r.extra["w_OH"] == pytest.approx(1.2, abs=1e-6)
    assert r.extra["M_anhy"] == pytest.approx(1341, abs=1.5)
    assert r.extra["n_H2O_total"] == pytest.approx(15.5, abs=0.15)


def test_deroo_table3_molar_masses():
    """De Roo 2025, Table 3: ideal MOF-808 1447.9, ideal UiO-66(Zr) 1664.1, experimental MOF-808 1620.2 (33.8 % Zr),
    experimental UiO-66(Zr) 1852.6 (29.6 % Zr). Masses reproduced within 0.5 g/mol; %Zr of UiO-66 computes to 29.55 % (paper rounds to 29.6)."""
    assert molar_mass("Zr6O4(OH)4(C9H3O6)2(C2H3O2)6") == pytest.approx(1447.9, abs=0.1)
    assert molar_mass("Zr6O4(OH)4(C8H4O4)6") == pytest.approx(1664.1, abs=0.1)
    m = molar_mass("Zr6O4(OH)4(OH)1.2(C9H3O6)2(C2H3O2)2(Cl)2.8(CH4O)0.3(H2O)15.5")
    assert m == pytest.approx(1620.2, abs=0.5)
    assert 100 * 6 * 91.224 / m == pytest.approx(33.8, abs=0.05)
    u = molar_mass("Zr6O4(OH)4(C8H4O4)5.1(Cl)1.8(C2H8NCl)0.9(C3H7NO)0.8(H2O)7.8")
    assert u == pytest.approx(1852.6, abs=0.5)
    # paper: 29.6 %; exact value from its own formula is 29.55 % (rounding in the paper)
    assert 100 * 6 * 91.224 / u == pytest.approx(29.6, abs=0.06)


# ---------------------------------------------------------------- Sannes 2023 M-200: exact round trip
def test_sannes_m200_round_trip():
    """Sannes 2023 eqs. 2, 6, 7: build Zr6O4(OH)4(BDC)5(HCOO)1(OH)1(H2O)1 + 0.4 Cl, compute m200/m800,
    feed NMR ratios -> recover L = 5, m = 1, beta. (No raw numbers published in the main text -> exact round-trip.)"""
    L, mf, cl = 5.0, 1.0, 0.4
    w = 12 - 2 * L - mf - cl
    M = molar_mass("Zr6O4(OH)4") + L * molar_mass("C8H4O4") + mf * molar_mass("CHO2") + cl * 35.45 + w * (17.007 + 18.015)
    W200 = 100.0
    Wend = W200 * 6 * oxide_mass(0) / M
    inp = DefectInput(framework="UiO-66", metal="Zr", species=[sp("BDC (terephthalate)", integral=4 * L * 0.37, main=True),
                      sp("Formate", integral=1 * mf * 0.37)], tga_W_200=W200, tga_W_end=Wend, cl_per_node=cl)
    r = method_b(inp)[0]
    assert r.L == pytest.approx(L, abs=1e-6)
    assert r.m == pytest.approx(mf, abs=1e-6)
    assert r.extra["beta"] == pytest.approx(w, abs=1e-6)
    assert open_sites(inp, r.L, r.m, r.Cl) == pytest.approx(2 * 1 - 1 - 0.4)


def test_nmr_corrected_plateau_benzoate():
    """Method (c): benzoate surviving 400 °C. Zr6O(6+...)(BDC)5(BA)1 built exactly; (a) overestimates L, (c) recovers 5."""
    L, nb = 5.0, 1.0
    Mpl = 6 * oxide_mass(0) + L * (molar_mass("C8H4O4") - 15.999) + nb * (molar_mass("C7H5O2") - 15.999 / 2)
    Wpl = 100 * Mpl / (6 * oxide_mass(0))
    inp = DefectInput(framework="UiO-66", metal="Zr", species=[sp("BDC (terephthalate)", integral=4 * L, main=True),
                      sp("Benzoate", integral=5 * nb, at_plateau=True)], tga_W_plateau=Wpl, tga_W_end=100)
    assert method_c(inp).L == pytest.approx(5.0, abs=1e-6)
    assert method_a(inp).L > 5.5


# ---------------------------------------------------------------- Klet 2016
def test_klet_convention_and_pka_reference():
    """Klet 2016: '1.60 missing linkers (of ideally 12)' for Zr-UiO-66 = missing connections per node.
    With p = 1.60 OH/OH2 pairs per node the calculator returns L = 5.20 and x = 0.80 per formula unit (ideal 6).
    The raw titrant volumes are only in the ESI, so this is a convention check, not a raw-data reproduction."""
    M = molar_mass("Zr6O4(OH)4(C8H4O4)5.2(OH)1.6(H2O)1.6")
    n_node = 50.0 / M  # mmol in 50 mg
    c = 0.1
    v1 = 1.0
    v2 = v1 + 1.6 * n_node / c
    v3 = v2 + 1.6 * n_node / c
    inp = DefectInput(framework="UiO-66", metal="Zr", species=[sp("BDC (terephthalate)", main=True)],
                      tit_mass_mg=50.0, tit_c_NaOH=c, tit_V_eq=[v1, v2, v3], tit_M_node=M)
    r = method_g(inp)
    assert r.extra["missing_connections_Klet"] == pytest.approx(1.60, abs=1e-6)
    assert r.L == pytest.approx(5.20, abs=1e-6)
    assert r.x == pytest.approx(0.80, abs=1e-6)


def test_composite_extraframework_zro2_correction():
    """ZrO2@Hf-MOF-801: 10 wt% extra-framework ZrO2 must not lower the Hf-MOF linker count."""
    L = 5.5
    M = 6 * oxide_mass(1) + L * (molar_mass("C4H2O4") - 15.999)
    hfo2 = 6 * oxide_mass(1)
    # per 100 mg plateau framework + ZrO2
    m_fw_pl = 90.0; m_fw_end = m_fw_pl * hfo2 / M; m_zro2 = 10.0
    w_zr = m_zro2 / oxide_mass(0) * 91.224
    # Hf wt% relative to the same initial mass (initial = plateau here for simplicity)
    w_hf = m_fw_end / oxide_mass(1) * 178.49
    inp = DefectInput(framework="MOF-801", metal="Hf", extraframework_Zr=True, species=[sp("Fumarate", main=True)],
                      tga_W_plateau=m_fw_pl + m_zro2, tga_W_end=m_fw_end + m_zro2, icp_Zr_wt=w_zr, icp_Hf_wt=w_hf)
    assert method_a(inp).L == pytest.approx(L, abs=1e-6)
    inp.extraframework_Zr = False
    inp.metal = "Hf"
    assert method_a(inp).L < L - 0.5   # uncorrected -> strongly underestimated


def test_run_all_consensus():
    inp = _deroo_mof808()
    res = run_all(inp)
    assert res["consensus"]["L"] == pytest.approx(2.0, abs=0.05)
