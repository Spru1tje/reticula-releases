"""v2.16: plausibility checks between techniques."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mfs import crosscheck as X, samplestore as SS


def rec(**sections):
    st = {}
    for sec, vals in sections.items():
        SS.put(st, "S1", sec, vals, "test")
    return st["sample_store"]["s1"]


def levels(r):
    return {f["rule"]: f["level"] for f in X.check_sample(r)}


def test_size_rules():
    assert levels(rec(pxrd={"crystallite_size_WH_nm": 40.0}, sem={"mean (nm)": 20.0}))["size"] == "warn"          # particle smaller than crystal
    assert levels(rec(pxrd={"crystallite_size_WH_nm": 40.0}, sem={"mean (nm)": 45.0}))["size"] == "ok"
    assert levels(rec(pxrd={"crystallite_size_scherrer_nm": 20.0}, sem={"mean (nm)": 50.0}))["size"] == "info"
    f = [f for f in X.check_sample(rec(pxrd={"crystallite_size_scherrer_nm": 10.0}, sem={"mean (nm)": 100.0}))][0]
    assert f["level"] == "info" and "polycrystalline" in f["title"] and f["values"]["ratio"] == 10.0
    assert X.check_sample(rec(pxrd={"crystallite_size_scherrer_nm": 10.0})) == []                   # one technique only: nothing to compare


def test_bet_and_dls_rules():
    # spheres of 100 nm, density 1.8 g/cm3 -> external area 33 m2/g
    r = rec(sem={"mean (nm)": 100.0}, pxrd={"density_g_cm3": 1.8}, n2={"S_BET": 10.0})
    f = {x["rule"]: x for x in X.check_sample(r)}
    assert f["bet"]["level"] == "warn" and abs(f["bet"]["values"]["S_geometric"] - 33.3) < 0.2
    assert levels(rec(sem={"mean (nm)": 100.0}, pxrd={"density_g_cm3": 1.8}, n2={"S_BET": 900.0}))["bet"] == "info"
    assert levels(rec(sem={"mean (nm)": 100.0}, pxrd={"density_g_cm3": 1.8}, n2={"S_BET": 50.0}))["bet"] == "ok"
    assert levels(rec(sem={"mean (nm)": 100.0}, dls={"z_avg_nm": 40.0}))["dls"] == "warn"
    assert levels(rec(sem={"mean (nm)": 100.0}, dls={"z_avg_nm": 160.0}))["dls"] == "ok"
    assert levels(rec(sem={"mean (nm)": 100.0}, dls={"z_avg_nm": 600.0, "pdi": 0.5}))["dls"] == "info"


def test_defect_and_surface_rules_and_ordering():
    assert levels(rec(defects={"L_range": "4.20–5.80"}))["defects"] == "warn"
    assert levels(rec(defects={"L_range": "5.10–5.60"}))["defects"] == "ok"
    assert levels(rec(icp={"Hf/(Zr+Hf)": 0.50}, xps={"Hf/(Zr+Hf)": 0.80}))["surface"] == "info"
    assert levels(rec(icp={"Hf/(Zr+Hf)": 0.50}, xps={"Hf/(Zr+Hf)": 0.52}))["surface"] == "ok"
    st = {}
    SS.put(st, "A", "pxrd", {"crystallite_size_WH_nm": 10.0}, "")
    SS.put(st, "A", "sem", {"mean (nm)": 11.0}, "")
    SS.put(st, "B", "pxrd", {"crystallite_size_WH_nm": 50.0}, "")
    SS.put(st, "B", "sem", {"mean (nm)": 10.0}, "")
    res = X.check(st["sample_store"])
    assert [f["level"] for f in res] == ["warn", "ok"] and res[0]["sample"] == "B"
    md = X.as_markdown(res)
    assert "⚠️" in md and "B – Particles smaller" in md and "No cross" in X.as_markdown([])
