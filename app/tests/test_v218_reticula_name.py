"""v2.18: the program is called Reticula; project files are .reticula (projects saved as .reticula still open)."""
import io, json, os, sys, zipfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pytest
from mfs import project as PRJ


def _zip(manifest: dict) -> bytes:
    bio = io.BytesIO()
    with zipfile.ZipFile(bio, "w") as z:
        z.writestr("manifest.json", json.dumps(manifest))
        z.writestr("globals.json", "{}")
    return bio.getvalue()


def test_project_extension_is_reticula_and_old_reticula_files_still_open():
    assert PRJ.EXT == ".reticula"
    assert ".reticula" in PRJ.OPEN_EXTS and ".reticula" in PRJ.OPEN_EXTS
    new = _zip({"app": "Reticula", "format": 1, "version": "2.18.0", "pages": []})
    old = _zip({"app": "Reticula", "format": 1, "version": "2.16.0", "pages": []})
    assert PRJ.read_manifest(new)["app"] == "Reticula"
    assert PRJ.read_manifest(old)["app"] == "Reticula"


def test_foreign_zip_is_rejected_with_the_program_name():
    with pytest.raises(ValueError, match="Reticula"):
        PRJ.read_manifest(_zip({"app": "Something else"}))


def test_saved_project_is_a_reticula_manifest(tmp_path):
    from streamlit.testing.v1 import AppTest
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.environ.update(MFS_STORE=str(tmp_path / "s.json"), MFS_DB=str(tmp_path / "r.db"), MFS_EXPORT_DIR=str(tmp_path / "exports"))
    at = AppTest.from_file(os.path.join(root, "app.py"), default_timeout=300).run()
    assert not at.exception
    fz = {k: at.session_state[k] for k in at.session_state.keys()}
    man = PRJ.read_manifest(PRJ.export_bytes("Home", "2.18.0", fz))
    assert man["app"] == "Reticula" and man["version"] == "2.18.0"


def test_home_page_has_the_project_box():
    from streamlit.testing.v1 import AppTest
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    at = AppTest.from_file(os.path.join(root, "app.py"), default_timeout=300).run()
    assert not at.exception, [e.value for e in at.exception][:1]
    text = " ".join(m.value for m in at.markdown)
    assert "Continue a project" in text and "Save your work" in text
    assert any("Reticula" in h.value for h in at.markdown)
