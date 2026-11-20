"""Regression test on the user's REAL Micromeritics TriStar II 3020 'Full Report Set' export (MOF-801, SGA2).
All reports sit side by side on one sheet, in column bands separated by '|' columns, with decimal commas in the
header/summary values. v1 merged every table (154 points, zig-zag); the reader must return exactly the
Isotherm Tabular Report (33 adsorption + 26 desorption points, identical to the instrument's own Isotherm Linear Plot)."""
import os
import numpy as np
import pytest
import xlrd
from mfs import readers as R, analysis as A

F = os.path.join(os.path.dirname(__file__), "data", "SampleASGA2_TriStar_full_report_REAL.xls")


@pytest.fixture(scope="module")
def d():
    return R.read_sorption(open(F, "rb"), os.path.basename(F))


def test_isotherm_is_tabular_report_split_like_instrument(d):
    iso = d["iso"]
    ads, des = iso[iso.branch == "ads"], iso[iso.branch == "des"]
    assert len(ads) == 33 and len(des) == 26
    assert d["meta"]["isotherm_table"].startswith("Isotherm Tabular Report")
    # instrument's own split (its desorption column starts with the highest-pressure point again): Isotherm Linear Plot columns L/M (ads) and N/O (des), rows 32+
    sh = xlrd.open_workbook(F).sheet_by_index(0)
    ref_ads = [(sh.cell_value(r, 11), sh.cell_value(r, 12)) for r in range(31, sh.nrows) if isinstance(sh.cell_value(r, 11), float)]
    ref_des = [(sh.cell_value(r, 13), sh.cell_value(r, 14)) for r in range(31, sh.nrows) if isinstance(sh.cell_value(r, 13), float)]
    assert np.allclose(ads[["p_rel", "q"]].values, ref_ads) and np.allclose(des[["p_rel", "q"]].values, ref_des[1:])   # instrument repeats the turning point
    assert np.all(np.diff(ads.p_rel.values) > 0) and np.all(np.diff(des.p_rel.values) < 0)   # no zig-zag
    assert iso.p_rel.min() > 0.01 and iso.p_rel.max() < 0.99


def test_meta_per_report_band_and_decimal_commas(d):
    m = d["meta"]
    assert m["sample"] == "SampleASGA2" and m["instrument"] == "TriStar"
    assert float(m["mass_g"]) == pytest.approx(0.0425) and float(m["bath_K"]) == pytest.approx(77.35)
    assert float(m["bet_reported"]) == pytest.approx(335.0806)
    assert float(m["bet_C_reported"]) == pytest.approx(-82.007886)          # BET band, not another report's value
    assert float(m["bet_R_reported"]) == pytest.approx(0.998715)            # not Freundlich's 0.112715
    assert m["bet_reported_range"] == pytest.approx((0.05666, 0.19946), abs=1e-4) and m["bet_reported_npts"] == 6
    kinds = {t["kind"] for t in d["tables"]}
    assert {"iso_tabular", "iso_plot", "bet", "tplot", "langmuir"} <= kinds
    assert d["psd"] is None                                                  # only BJH/D-H summary values in this file


def test_bet_app_method_and_file_reference(d):
    a = d["iso"][d["iso"].branch == "ads"]
    b = A.bet_rouquerol(a.p_rel.values, a.q.values)
    assert b is not None and b["few_points"] and b["n"] == 3
    assert b["area"] == pytest.approx(400.9, abs=0.5) and b["C"] == pytest.approx(1127.5, rel=0.01)
    assert b["R2"] > 0.9999 and b["all_criteria"]
    assert (b["p_lo"], b["p_hi"]) == pytest.approx((0.010726, 0.044454), abs=1e-5)
    # same data on the file's window reproduces the instrument's number (and its negative C)
    f = A.bet_fixed_range(a.p_rel.values, a.q.values, *d["meta"]["bet_reported_range"])
    assert f["area"] == pytest.approx(335.08, abs=0.05) and f["C"] < 0


def test_tplot_matches_report_and_quality_flags(d):
    a = d["iso"][d["iso"].branch == "ads"]
    tp = A.t_plot(a.p_rel.values, a.q.values)
    assert tp["V_micro"] == pytest.approx(0.138845, rel=0.01) and tp["S_ext"] == pytest.approx(33.85, rel=0.02)
    w = A.isotherm_quality(a.p_rel.values, a.q.values)
    assert any("DECREASES" in x for x in w) and any("not resolved" in x for x in w)


def test_synthetic_examples_unaffected():
    ex = os.path.join(os.path.dirname(os.path.dirname(__file__)), "example_data", "Sorption", "SampleA_3Flex_full_report_EXAMPLE.xls")
    dd = R.read_sorption(open(ex, "rb"), os.path.basename(ex))
    assert len(dd["iso"]) > 20 and dd["meta"].get("bet_reported")
