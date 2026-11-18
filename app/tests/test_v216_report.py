"""v2.16: problem report."""
import io, json, os, sys, zipfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mfs import problem_report as PR


def test_report_has_no_data_and_anonymises_paths():
    home = os.path.expanduser("~")
    ss = {"nav": "PXRD", "pxrd_lo": 4.0, "name_x": "MOF-801", "big_table": "x" * 500, "_lv": {"global": "Advanced"}, "ak|pxrd|x|rv_off": "zero",
          "fc_PXRD_title": "secret title", "si_figs": {"PXRD": {"page": "PXRD", "title": "PXRD", "n_rows": 10, "data": "DATA"}},
          "sample_store": {"mof801": {"name": "MOF-801", "sections": {"tga": {"values": {"r": 1}}}}}}
    try:
        open(os.path.join(home, "nope_nope", "x.csv"))
    except Exception as e:
        PR.record_error(ss, "PXRD", e)
    z = zipfile.ZipFile(io.BytesIO(PR.build(ss, "PXRD", "2.16.0", [os.path.join(home, "data", "a.xy")], "it crashed when I clicked Refine")))
    rep = json.loads(z.read("report.json"))
    txt = z.read("report.json").decode() + z.read("report.md").decode()
    assert rep["version"] == "2.16.0" and rep["page"] == "PXRD" and rep["detail_level"] == "Advanced"
    assert rep["last_error"]["type"] == "FileNotFoundError" and "it crashed" in z.read("report.md").decode()
    assert home not in txt and "DATA" not in txt and "secret title" not in txt and "x" * 500 not in txt
    assert rep["settings_on_this_page"].get("pxrd_lo") == 4.0 and rep["settings_on_this_page"].get("rv_off") == "zero"
    assert rep["figures_made"][0]["figure"] == "PXRD" and rep["sample_store_sections"] == {"mof801": ["tga"]}
    assert rep["loaded_files"][0].startswith("~") and rep["loaded_files"][0].endswith("a.xy") and home not in rep["loaded_files"][0]


def test_report_button_in_the_app():
    import warnings
    warnings.filterwarnings("ignore")
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "app.py"), default_timeout=300)
    at.run()
    assert not at.exception
    assert any("Report a problem" in e.label for e in at.sidebar.expander)
