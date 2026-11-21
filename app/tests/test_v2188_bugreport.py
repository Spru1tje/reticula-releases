"""v2.18.8: an error on a page must end up in the problem-report ZIP.

Before: the report panel was drawn after the page call, so when the page raised, the panel was not drawn in that run.
Now: the error is recorded at the page dispatcher, and the ZIP is built from the session when the user clicks."""
import io
import json
import zipfile

from mfs import problem_report as PR


def _zip_text(data: bytes):
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        return json.loads(z.read("report.json")), z.read("report.md").decode("utf-8")


def test_error_from_a_page_is_in_the_zip():
    session = {"nav": "Photoluminescence"}
    try:
        raise ValueError("Singular matrix at power fit (test)")
    except ValueError as exc:
        PR.record_error(session, "Photoluminescence", exc)
    rep, md = _zip_text(PR.build(session, "Photoluminescence", "2.18.8", [], "test"))
    assert rep["last_error"]["type"] == "ValueError"
    assert "Singular matrix at power fit (test)" in rep["last_error"]["message"]
    assert "Traceback" in rep["last_error"]["traceback"] or "ValueError" in rep["last_error"]["traceback"]
    assert "Singular matrix at power fit (test)" in md


def test_no_error_says_so():
    rep, md = _zip_text(PR.build({}, "Home", "2.18.8", [], ""))
    assert rep["last_error"] is None
    assert "No error was recorded" in md


def test_home_paths_are_anonymised_in_the_error():
    import os
    home = os.path.expanduser("~")
    session = {}
    try:
        raise OSError(f"cannot open {home}{os.sep}data.csv")
    except OSError as exc:
        PR.record_error(session, "PXRD", exc)
    rep, _ = _zip_text(PR.build(session, "PXRD", "2.18.8", [], ""))
    assert home not in rep["last_error"]["message"]
