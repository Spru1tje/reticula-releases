"""v2.15: NMR library matching."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mfs import nmrlib as L


def test_digest_in_naod_bdc_formate_acetate_ratios():
    # 1 BDC (4 H, integral 4.0), 0.5 formate (1 H, 0.5), 0.3 acetate (3 H, 0.9)
    mults = [(7.76, 4.0), (8.44, 0.5), (1.90, 0.9), (4.79, 50.0)]
    m = L.match(mults, "NaOD/D2O")
    names = [x["compound"] for x in m]
    assert "terephthalic acid (BDC)" in names[:3] and "formic acid / formate" in names and "acetic acid / acetate" in names
    r = L.ratios(m, "terephthalic acid (BDC)")
    assert abs(r["formic acid / formate"] - 0.5) < 0.02 and abs(r["acetic acid / acetate"] - 0.3) < 0.02


def test_incomplete_matches_rank_lower_and_own_rows_win():
    m = L.match([(7.95, 1.0), (2.89, 3.0)], "DMSO-d6")
    dmf = next(x for x in m if x["compound"] == "DMF")
    assert dmf["found"] == "2/3" and not dmf["complete"]
    own = [{"compound": "my linker", "solvent": "DMSO", "δ (ppm)": 7.30, "mult": "s", "nH": 2}]
    m2 = L.match([(7.31, 2.0)], "DMSO", own)
    assert m2[0]["compound"] == "my linker" and m2[0]["source"] == "own" and m2[0]["complete"]
    assert L.solvent_key("DMSO-d6") == "DMSO" and L.solvent_key(None) == "D2O"
