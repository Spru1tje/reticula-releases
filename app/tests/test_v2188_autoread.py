"""v2.18.8: automatic reading of instrument exports with unit recognition, each on a synthetic file with a known answer."""
import numpy as np

from mfs import autoread as AR


def test_units_are_recognised_and_converted():
    assert AR.detect_unit("Temperature (°C)")["status"] == "converted"
    t, info = AR.convert_column([25.0, 100.0], "Temperature (°C)")
    assert np.allclose(t, [298.15, 373.15])
    p, _ = AR.convert_column([250.0], "Power / mW")
    assert np.allclose(p, [0.25])
    m, _ = AR.convert_column([2.0], "Mass [mg]")
    assert np.allclose(m, [0.002])
    tm, _ = AR.convert_column([5.0], "Time (min)")
    assert np.allclose(tm, [300.0])


def test_unknown_unit_is_kept_and_reported_not_guessed():
    v, info = AR.convert_column([1.0, 2.0], "Signal (arb)")
    assert info["status"] == "unknown unit"
    assert np.allclose(v, [1.0, 2.0])


TGA_LIKE = """Instrument: TGA-synthetic
Operator: test
Sample: MOF-A
Temperature (°C)\tMass (mg)\tTime (min)
25.0\t10.00\t0
100.0\t9.80\t15
200.0\t9.10\t30
300.0\t8.00\t45
"""

SEMICOLON = """Export v2
Temperature [K];Power / mW
298,15;250
373,15;260
448,15;270
"""


def test_tab_separated_export_with_metadata_and_units():
    res = AR.read_text(TGA_LIKE)
    assert res["skipped_lines"] == 3
    assert res["delimiter"] == "\t"
    assert np.allclose(res["data"]["Temperature (°C)"], [298.15, 373.15, 473.15, 573.15])
    assert np.allclose(res["data"]["Mass (mg)"], [0.01, 0.0098, 0.0091, 0.008])
    assert np.allclose(res["data"]["Time (min)"], [0, 900, 1800, 2700])


def test_semicolon_export_with_decimal_comma():
    res = AR.read_text(SEMICOLON)
    assert res["delimiter"] == ";"
    assert np.allclose(res["data"]["Temperature [K]"], [298.15, 373.15, 448.15])
    assert np.allclose(res["data"]["Power / mW"], [0.25, 0.26, 0.27])


def test_file_without_a_table_raises_a_clear_error():
    try:
        AR.read_text("only text here\nno numbers\n")
    except ValueError as exc:
        assert "no table" in str(exc)
    else:
        raise AssertionError("expected ValueError")


def test_read_file_from_disk(tmp_path):
    p = tmp_path / "export.txt"
    p.write_bytes(TGA_LIKE.encode("utf-8"))
    res = AR.read_file(str(p))
    assert res["units"][0]["quantity"] == "temperature"
    assert res["units"][0]["status"] == "converted"
