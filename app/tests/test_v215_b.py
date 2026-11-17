"""v2.15: methods text + sample summary."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mfs import summary as M, ui


def test_methods_and_summary():
    skip = ("Home", "Figure composer", "Batch processing", "SI export bundle", "Sample database", "Export guide")
    assert all(p in M.TECHNIQUE_METHODS for p in ui.PAGES if p not in skip), "a technique has no methods text"
    t = M.methods_text(["PXRD", "TGA", "PXRD", "Home"], "2.15.0")
    assert t.count("**PXRD.**") == 1 and "Reticula 2.15.0" in t and "Home" not in t
    store = {"a": {"name": "MOF_1", "sections": {"tga": {"values": {"residue": 41.25, "x": {"y": 2}}}}},
             "b": {"name": "MOF 2", "sections": {"n2": {"values": {"bet": 1050.0}}}}}
    df = M.sample_summary(store)
    assert list(df["Sample"]) == ["MOF_1", "MOF 2"] and "tga: x.y" in df and df.loc[1, "n2: bet"] == 1050.0
    tex = M.summary_latex(df)
    assert "MOF\\_1" in tex and "\\toprule" in tex and tex.count("\\\\") == 3
