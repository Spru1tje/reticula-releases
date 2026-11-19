"""Batch D: journal templates, batch processing, SI bundle, SQLite sample database."""
import io, os, zipfile
import numpy as np
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX = os.path.join(HERE, "example_data")
MASTER = os.path.join(EX, "SampleDB", "Sample_overview_EXAMPLE.xlsx")     # synthetic sample list (bundled example)


def test_templates_check():
    import matplotlib.pyplot as plt
    from mfs import templates as TP, composer as C
    fig, axs = C.compose([lambda ax: ax.plot([1, 2], [1, 2])] * 2, ncols=2, width=TP.width("RSC", "double"), font_size=7, label_style="(a)")
    assert any(l == "ok" for l, _ in TP.check(fig, "RSC", "double"))
    axs[0].text(0.5, 0.5, "tiny", fontsize=4)
    msgs = TP.check(fig, "Wiley", "single")
    assert any("width" in m for _, m in msgs) and any("below 6 pt" in m for _, m in msgs)
    assert abs(TP.width("Wiley", "double") - 17.5 / 2.54) < 1e-9 and TP.width("ACS", "single") == 3.25
    plt.close(fig)


def test_batchproc_all_kinds():
    from mfs import batchproc as BP
    for kind, d in (("PXRD", "PXRD"), ("TGA", "TGA"), ("DSC", "DSC"), ("N2 isotherm", "Sorption"), ("FTIR", "FTIR")):
        fs = BP.list_folder(os.path.join(EX, d), kind)
        assert fs and not any("README" in f for f in fs)
        z, log = BP.run(kind, [(f, f) for f in fs], cfg={"ticks": {"preset": "Bottom + left"}}, layout="one figure per file",
                        formats=("png", "pdf"), dpi=100)
        names = zipfile.ZipFile(io.BytesIO(z)).namelist()
        assert all(l["status"] == "ok" for l in log), (kind, log)
        assert sum(n.endswith(".pdf") for n in names) >= len(fs) and "batch_log.csv" in names and "preset_used.json" in names
    z, log = BP.run("PXRD", [(f, f) for f in BP.list_folder(os.path.join(EX, "PXRD"), "PXRD")], layout="stacked", formats=("svg",))
    assert [n for n in zipfile.ZipFile(io.BytesIO(z)).namelist() if n.endswith(".svg")] == ["PXRD_all_stacked.svg"]


def test_si_bundle():
    import matplotlib.pyplot as plt
    from mfs import si_bundle as SB, samplestore as SS
    st = {}; SS.put(st, "MOF-801", "defects", {"L": 5.0, "x": 1.0, "OS": 1.6}); SS.put(st, "MOF-801", "tga", {"W_end": 40.0})
    z = SB.build({"PXRD": {"page": "PXRD", "pdf": b"%PDF-1.4", "png": b"\x89PNG", "csv": b"x,y\n1,2", "params": {"pxrd_off": 1.1},
                           "cfg": {"ticks": {"preset": "All 4 sides"}}}}, st["sample_store"], "2.7.0")
    zf = zipfile.ZipFile(io.BytesIO(z)); n = zf.namelist()
    assert {"figures/Figure_S1_PXRD.pdf", "data/Figure_S1_PXRD.csv", "SI_parameters.xlsx", "SI_Table_defects.docx", "README.txt"} <= set(n)
    import pandas as pd
    par = pd.read_excel(io.BytesIO(zf.read("SI_parameters.xlsx")), sheet_name=None)
    assert set(par) == {"Figures", "Parameters", "Results per sample"}
    assert "pxrd_off" in set(par["Parameters"].parameter) and "figure.ticks.preset" in set(par["Parameters"].parameter)
    assert len(par["Results per sample"]) == 4


def test_sampledb(tmp_path):
    from mfs import sampledb as DB, samplestore as SS
    con = DB.connect(str(tmp_path / "r.db"))
    assert DB.import_master(con, MASTER) == 8
    s = DB.sample(con, "EX-U66-SO3H")
    assert s["parent"] == "EX-U66" and "NP@EX-U66-SO3H" in s["children"] and s["protocol_pdf"].endswith("EX-U66-SO3H.pdf")
    assert os.path.join("Paper1", "Protocols") in s["protocol_pdf"]
    assert DB.match_id(con, "EX-M801-F25-B1_PXRD.xy") == "EX-M801-F25" and DB.match_id(con, "notes.txt") is None
    d = tmp_path / "data" / "TGA"; d.mkdir(parents=True)
    (d / "NP@EX-U66-SO3H-B2_TGA.txt").write_text("x"); (d / "unknown.txt").write_text("x")
    n, un = DB.link_folder(con, str(tmp_path / "data"))
    assert n == 1 and len(un) == 1
    f = DB.sample(con, "NP@EX-U66-SO3H")["files"][0]; assert f["technique"] == "TGA" and f["batch"] == "2"
    st = {}; SS.put(st, "Ex-U66-SO3H", "xrf", {"S_per_M6": 0.6}); SS.put(st, "New sample X", "tga", {"W_end": 41})
    assert DB.sync_store(con, st["sample_store"]) == 2
    assert DB.sample(con, "EX-U66-SO3H")["results"][0]["value"] == 0.6 and DB.sample(con, "New sample X")
    ov = DB.overview(con); assert "protocol" in ov.columns and len(ov) >= 9


def test_sampledb_generic_columns(tmp_path):
    """A user's own list: only 'Sample ID' is required; 'Project' is accepted instead of 'Paper'."""
    import openpyxl
    from mfs import sampledb as DB
    wb = openpyxl.Workbook(); ws = wb.active
    ws.append(["Sample ID", "Name", "Project"]); ws.append(["A-1", "first", "Thesis 2"]); ws.append(["A-2", None, None])
    p = tmp_path / "mine.xlsx"; wb.save(p)
    con = DB.connect(str(tmp_path / "g.db"))
    assert DB.import_master(con, str(p)) == 2
    assert DB.sample(con, "A-1")["paper"] == "Thesis 2" and DB.sample(con, "A-2")["name"] is None


def test_no_private_pipeline_bundled():
    """The author's synthesis pipeline must never ship again."""
    assert not os.path.exists(os.path.join(HERE, "references", "Synthesis_Pipeline"))


# ---------------- app ----------------
def _goto(at, page):
    at.sidebar.radio[0].set_value(page).run()
    assert not at.exception, [e.value for e in at.exception][:2]
    return at


@pytest.fixture
def at(tmp_path, monkeypatch):
    monkeypatch.setenv("MFS_STORE", str(tmp_path / "store.json")); monkeypatch.setenv("MFS_DB", str(tmp_path / "r.db"))
    monkeypatch.setenv("MFS_EXPORT_DIR", str(tmp_path / "exports"))
    from streamlit.testing.v1 import AppTest
    a = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=400); a.run()
    return a


def test_app_batch_si_db_composer(at):
    _goto(at, "Batch processing")
    at.selectbox(key="bp_kind").set_value("TGA").run()
    at.button(key="bp_run").click().run(); assert not at.exception, [e.value for e in at.exception][:2]
    from mfs import batchproc as BP
    n = len(BP.list_folder(os.path.join(EX, "TGA"), "TGA"))   # the folder decides the count; do not hard-code it
    assert any(f"Batch done: {n} / {n}" in s.value for s in at.success)
    _goto(at, "PXRD"); _goto(at, "TGA")
    figs = at.session_state["si_figs"]; assert len(figs) >= 2 and {f["page"] for f in figs.values()} >= {"PXRD", "TGA"}
    _goto(at, "SI export bundle"); assert not at.info
    _goto(at, "Sample database")
    at.button(key="db_imp").click().run(); assert not at.exception, [e.value for e in at.exception][:2]
    assert any("sample IDs imported" in s.value for s in at.success)
    assert at.selectbox(key="db_pick").options
    _goto(at, "Figure composer")
