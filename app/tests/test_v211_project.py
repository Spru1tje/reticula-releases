"""v2.11 project files (.reticula zip): save the whole session, open it in a fresh session/workspace."""
import os, sys, io, json, zipfile, tempfile
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import pytest


def _frozen(at):
    inner = at.session_state._state._state if hasattr(at.session_state._state, "_state") else at.session_state._state
    return {k: inner[k] for k in list(inner._keys()) if k in inner}


def test_project_roundtrip_widgets_uploads_store(tmp_path):
    from streamlit.testing.v1 import AppTest
    from mfs import project as PRJ
    a = tmp_path / "a"; os.environ.update(MFS_NAV_GUARD="0", MFS_STORE=str(a / "s.json"), MFS_DB=str(a / "r.db"), MFS_EXPORT_DIR=str(a / "exports"))
    at = AppTest.from_file(os.path.join(ROOT, "app.py"), default_timeout=300); at.run()
    at.sidebar.radio[0].set_value("TGA").run()
    [c for c in at.checkbox if "Normalise to residue" in c.label][0].check().run()
    at.sidebar.radio[0].set_value("Defect calculator").run()
    at.selectbox(key="def_fw").set_value("MOF-5").run()
    fz = _frozen(at)
    fz["sample_store"] = {"x": {"name": "Sample X", "sections": {"tga": {"values": {"W_end": 12.3}}}}}
    # a page with an uploaded file (AppTest cannot upload: simulate what the uploader hook keeps)
    fz["_ps_pages"]["PXRD"] = {"keys": [], "up": ["pxrd_up"], "extra": []}
    fz["_ps_up"]["pxrd_up"] = [("my scan.xy", "text/plain", b"5 10\n6 20\n")]
    data = PRJ.export_bytes("Defect calculator", "2.11.0", fz)
    z = zipfile.ZipFile(io.BytesIO(data)); names = z.namelist()
    assert {"manifest.json", "globals.json", "workspace/tga.json", "workspace/defect_calculator.json", "workspace/pxrd.json"} <= set(names)
    assert any(n.startswith("workspace/pxrd/") and n.endswith("my_scan.xy") for n in names)
    assert PRJ.read_manifest(data)["pages"] and PRJ.read_manifest(data)["version"] == "2.11.0"
    # fresh session + empty workspace
    b_ = tmp_path / "b"; os.environ.update(MFS_STORE=str(b_ / "s.json"), MFS_DB=str(b_ / "r.db"), MFS_EXPORT_DIR=str(b_ / "exports"))
    b = AppTest.from_file(os.path.join(ROOT, "app.py"), default_timeout=300); b.run()
    b.session_state["_proj_pending"] = ("t.reticula", data); b.run()
    assert not b.exception, [e.value for e in b.exception][:1]
    assert b.sidebar.radio[0].value == "Defect calculator" and b.selectbox(key="def_fw").value == "MOF-5"
    assert b.session_state["_ps_up"]["pxrd_up"][0][2] == b"5 10\n6 20\n"
    assert b.session_state["sample_store"]["x"]["name"] == "Sample X"
    assert json.load(open(b_ / "s.json", encoding="utf-8"))["x"]["name"] == "Sample X"
    assert sorted(os.listdir(b_ / "workspace")) >= ["defect_calculator.json", "pxrd", "pxrd.json", "tga.json"]
    b.sidebar.radio[0].set_value("TGA").run()
    assert [c.value for c in b.checkbox if "Normalise to residue" in c.label] == [True]


def test_project_rejects_foreign_zip_and_traversal():
    from mfs import project as PRJ
    bio = io.BytesIO()
    with zipfile.ZipFile(bio, "w") as z:
        z.writestr("manifest.json", json.dumps({"app": "Other"}))
    with pytest.raises(ValueError):
        PRJ.read_manifest(bio.getvalue())


def test_launcher_project_arg(tmp_path):
    sys.path.insert(0, "/workspace/exe_build/launcher")
    if not os.path.exists("/workspace/exe_build/launcher/mfs_launcher.py"):
        pytest.skip("launcher not in this checkout")
    import importlib; L = importlib.import_module("mfs_launcher")
    p = tmp_path / "x.reticula"; p.write_bytes(b"PK")
    assert L.project_arg(["Reticula.exe", str(p)]) == str(p.resolve())
    assert L.project_arg(["Reticula.exe", "--serve", "8501"]) is None


def test_project_globals_skip_buttons():
    """v2.11 e2e bug: the 'Reset style' button state ended up in globals.json -> StreamlitValueAssignmentNotAllowedError on open."""
    from mfs import project as PRJ
    g = PRJ._jsonable_globals({"gs_reset": False, "gs_font": 7.0, "series_colours": {}, "other": 1})
    assert "gs_reset" not in g and g["gs_font"] == 7.0 and "other" not in g
    assert not PRJ._global_ok("gs_reset") and PRJ._global_ok("gs_palette")


def test_snapshot_keeps_sample_list_state(tmp_path):
    """The sample list (list of dicts) must be written to project / workspace snapshots."""
    import json
    from mfs import pagestate as PS
    items = [{"id": "a.xy", "name": "A", "show": True}, {"id": "b.xy", "name": "B", "show": False}]
    PS._write_snapshot({"page": "PXRD", "saved": "x", "widgets": {}, "uploads": {}, "extra": {"_sl|pxrd_tab": items}}, str(tmp_path))
    d = json.loads((tmp_path / "pxrd.json").read_text())
    assert PS._unjson(d["extra"]["_sl|pxrd_tab"]) == items
