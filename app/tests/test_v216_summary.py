"""v2.16: summary with sources."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mfs import samplestore as SS, summary as SUM


def test_sources_travel_with_every_value():
    state = {}
    SS.put(state, "MOF-801", "tga", {"residue": 41.2, "T": 600}, "TGA page", files=["SampleA_TGA.txt"], version="2.16.0")
    SS.put(state, "MOF-801", "n2", {"S_BET": 1050.0}, "N2 page", files=["SampleA_N2.xls"], version="2.16.0")
    SS.put(state, "Other", "tga", {"residue": 30.0}, "TGA page", files=["o.txt"], version="2.16.0")
    long = SUM.summary_long(state["sample_store"])
    row = long[(long.Sample == "MOF-801") & (long.Quantity == "residue")].iloc[0]
    assert row["Value"] == 41.2 and row["Source"] == "TGA page" and row["Data files"] == "SampleA_TGA.txt" and row["Version"] == "2.16.0" and row["Technique"] == "TGA"
    notes, ref = SUM.summary_notes(state["sample_store"])
    assert len(notes) == 3 and ref[("MOF-801", "tga")] == "a" and "SampleA_TGA.txt" in notes[0][1]
    tex = SUM.summary_latex_sourced(state["sample_store"])
    assert r"\toprule" in tex and "$^{a}$" in tex and r"SampleA\_N2.xls" in tex and "footnotesize" in tex and tex.count("\\\\") >= 5
    assert SUM.summary_latex_sourced({}) == ""
