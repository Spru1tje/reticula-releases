"""v2.15.1: dossier per material."""
import io, os, sys, zipfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mfs import dossier as D, style as S


def test_guess_group_and_build():
    assert D.guess_material(["MOF-801"]) == "MOF-801" and D.guess_material(["A", "B"]) == "" and D.guess_material([]) == ""
    figs = {"PXRD_a": {"page": "PXRD"}, "TGA_a": {"page": "TGA"}, "other": {"page": "FTIR"}}
    g = D.group_by_material(figs, {"PXRD_a": "MOF-801", "TGA_a": "MOF-801", "other": ""})
    assert list(g) == ["MOF-801"] and set(g["MOF-801"]) == {"PXRD_a", "TGA_a"}
    fig, ax = plt.subplots(); ax.plot([1, 2], [3, 4])
    rend = {"PXRD_a": {"page": "PXRD", "time": "10:00", "params": {"lo": 4}, "cfg": None, "title": "PXRD", "n_rows": 2,
                       "png": S.export(fig, "png", 100), "pdf": S.export(fig, "pdf"), "csv": b"x,y\n1,3\n", "data": pd.DataFrame({"x": [1, 2], "y": [3, 4]})},
            "TGA_a": {"page": "TGA", "time": "10:01", "params": {}, "cfg": None, "title": "TGA", "n_rows": 0, "png": S.export(fig, "png", 100), "data": None}}
    store = {"mof801": {"name": "MOF-801", "sections": {"tga": {"values": {"residue": 41.2}, "source": "t", "updated": "x"}}},
             "other": {"name": "Other", "sections": {"tga": {"values": {"residue": 9.9}}}}}
    z = zipfile.ZipFile(io.BytesIO(D.build_dossier("MOF-801", rend, store, "2.15.1")))
    names = z.namelist()
    assert "figure_data.xlsx" in names and "methods.md" in names and "sample_summary.csv" in names and any(n.startswith("figures/") for n in names)
    assert "41.2" in z.read("sample_summary.csv").decode() and "9.9" not in z.read("sample_summary.csv").decode()
    assert "PXRD" in z.read("methods.md").decode() and "TGA" in z.read("methods.md").decode()
    allz = zipfile.ZipFile(io.BytesIO(D.build_all({"MOF-801": rend, "B b": rend}, store, "2.15.1")))
    assert any(n.startswith("MOF-801/") for n in allz.namelist()) and any(n.startswith("B_b/") for n in allz.namelist())
    plt.close(fig)
