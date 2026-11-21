"""SEM/TEM calibration from the original image: metadata sidecars, scale-bar detection, OCR of the bar label and
the magnification, manual two-point calibration, magnification overlay and the calibration panel in the app."""
import os, sys, glob, hashlib
import numpy as np
import pytest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from semcal_synth import make, STYLES
from mfs import semcal, readers as R, plots as P

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CALEX = os.path.join(ROOT, "example_data", "SEM_calibration", "JEOL-style_noMetadata_x25000.tif")
needs_ocr = pytest.mark.skipif(not semcal.ocr_available(), reason="Tesseract not installed")


@pytest.mark.parametrize("txt,val", [("x10,000", 1e4), ("10.0kx", 1e4), ("50 000 x", 5e4), ("Mag = 10.00 K X", 1e4),
                                     ("×5.00k", 5e3), ("10000x", 1e4), ("MAG: 25.0 kx", 2.5e4), ("Magn 40000x", 4e4),
                                     ("X 1,500", 1500), ("x25,000", 2.5e4), ("mag 50 000 x", 5e4), ("Mag = 2.50 KX", 2500)])
def test_parse_magnification(txt, val):
    assert semcal.parse_magnification(txt) == pytest.approx(val)


@pytest.mark.parametrize("txt", ["WD 10.1 mm", "HV 5.00 kV", "EXAMPLE DATA", ""])
def test_parse_magnification_none(txt):
    assert semcal.parse_magnification(txt) is None


@pytest.mark.parametrize("txt,nm", [("1 µm", 1000), ("500nm", 500), ("1um", 1000), ("2 μm", 2000), ("0.5 µm", 500),
                                    ("5.00um", 5000), ("200 nm", 200), ("10 Å", 1.0)])
def test_parse_scale_label(txt, nm):
    assert semcal.parse_scale_label(txt)[2] == pytest.approx(nm)


def test_ocr_fixes_and_exclusions():
    assert semcal.parse_scale_label(semcal._fix_ocr("lum"))[2] == 1000
    assert semcal._excluded("WD 8.0mm", None) and not semcal._excluded("500 nm", None)


@pytest.mark.parametrize("style", STYLES)
@pytest.mark.parametrize("seed", [0, 1])
def test_bar_detected_exactly(style, seed):
    a, t = make(style, seed=seed)
    strip = semcal.find_info_strip(a)
    assert (strip == 884) if style != "overlay" else (strip is None)
    b = semcal.detect_bar(a, strip)
    assert b and b["length_px"] == t["bar_px"]


def test_bundled_examples_bar_and_strip():
    for f in sorted(glob.glob(os.path.join(ROOT, "example_data", "SEM", "*.tif"))):
        d = R.read_sem(open(f, "rb"), f)
        b = semcal.detect_bar(d["image"], semcal.find_info_strip(d["image"]))
        assert semcal.find_info_strip(d["image"]) == 968
        L_nm = 5000 if "lowmag" in f else 1000
        assert d["pixel_nm"] * b["length_px"] == pytest.approx(L_nm, rel=0.01)


@needs_ocr
@pytest.mark.parametrize("style", STYLES)
def test_calibrate_synthetic_with_ocr(style):
    a, t = make(style)
    c = semcal.calibrate(a)
    sb = [x for x in c["candidates"] if x["method"].startswith("scale bar")][0]
    assert sb["nm_per_px"] == pytest.approx(t["nm_per_px"], rel=1e-3)
    assert (c["mag_ocr"] or {}).get("value") == t["mag"]


@needs_ocr
def test_calibrate_examples_cross_check_metadata():
    for f in sorted(glob.glob(os.path.join(ROOT, "example_data", "SEM", "*.tif"))):
        d = R.read_sem(open(f, "rb"), f)
        c = semcal.calibrate(d["image"], d["pixel_nm"], d["meta"].get("magnification"), d["meta"].get("calibration"))
        assert c["method"].startswith("metadata") and c["nm_per_px"] == d["pixel_nm"]
        sb = [x for x in c["candidates"] if x["method"].startswith("scale bar")][0]
        assert sb["nm_per_px"] == pytest.approx(d["pixel_nm"], rel=0.01)
        assert c["mag_ocr"]["value"] == d["meta"]["magnification"]
        assert any("agree" in n for _, n in c["notes"])


@needs_ocr
def test_builtin_calibration_example_no_metadata():
    d = R.read_sem(open(CALEX, "rb"), CALEX)
    assert d["pixel_nm"] is None
    c = semcal.calibrate(d["image"])
    assert c["method"] == "scale bar (detected) + OCR label"
    assert c["nm_per_px"] == pytest.approx(1000 / 252, rel=1e-6)
    assert c["bar"]["length_px"] == 252 and c["label"]["nm"] == 1000 and c["mag_ocr"]["value"] == 25000
    assert c["strip_row"] == 896


def test_metadata_mag_vs_ocr_disagreement_flagged(monkeypatch):
    a, t = make("zeiss")
    monkeypatch.setattr(semcal, "read_strip", lambda *a_, **k: {"ocr": True, "text": "", "label": None,
                                                              "mag": {"value": 10000.0, "text": "Mag = 10.00 K X", "box": (0, 0, 1, 1)}})
    c = semcal.calibrate(a, 3.3333, 20000.0, "test")
    assert any(l == "warn" and "differ" in n for l, n in c["notes"])


def test_manual_two_points():
    assert semcal.nm_per_px_from_points(400, 652, 1, "µm") == pytest.approx(1000 / 252)
    assert semcal.nm_per_px_from_points(10, 10, 1, "µm") is None


def test_sidecars():
    jeol = "$CM_FORMAT JEOL/SEM\n$CM_MAG 25000\n$CM_FULL_SIZE 2560 1920\n$$SM_MICRON_BAR 504\n$$SM_MICRON_MARKER 1um\n"
    s = R.parse_sem_sidecar(jeol, 1280)
    assert s["source"] == "JEOL sidecar" and s["pixel_nm"] == pytest.approx(1000 / 252) and s["magnification"] == 25000
    hit = "[SemImageFile]\nMagnification=10000\nPixelSize=12.4023\n"
    s = R.parse_sem_sidecar(hit)
    assert s["source"] == "Hitachi sidecar" and s["pixel_nm"] == pytest.approx(12.4023) and s["magnification"] == 10000
    tes = "[MAIN]\nDevice=TESCAN MIRA\nMagnification=50000\nPixelSizeX=2.2727e-09\n"
    s = R.parse_sem_sidecar(tes)
    assert s["source"] == "TESCAN .hdr" and s["pixel_nm"] == pytest.approx(2.2727) and s["magnification"] == 50000
    z = "AP_MAG\nMag = 10.00 K X\nAP_IMAGE_PIXEL_SIZE\nImage Pixel Size = 3.333 nm\n"
    s = R.parse_sem_sidecar(z)
    assert s["pixel_nm"] == pytest.approx(3.333) and s["magnification"] == 10000


def test_read_sem_with_sidecar():
    d = R.read_sem(open(CALEX, "rb"), CALEX, sidecar="$CM_MAG 25000\n$$SM_MICRON_BAR 252\n$$SM_MICRON_MARKER 1um\n")
    assert d["pixel_nm"] == pytest.approx(1000 / 252) and d["meta"]["calibration"] == "JEOL sidecar"
    assert d["meta"]["magnification"] == 25000


def test_mag_overlay_on_figure():
    a, _ = make("fei_box")
    fig = P.fig_sem_panel([{"image": a, "pixel_nm": 2.27, "crop_row": 884, "bar": {"mode": "value", "value": 500, "unit": "nm"},
                            "mag_text": "×50 000", "mag_pos": "lower left"}], 1, 3.3, ["(a)"])
    texts = [t.get_text() for t in fig.axes[0].texts]
    assert "×50 000" in texts and "500 nm" in texts


def _goto_sem():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(os.path.join(ROOT, "app.py"), default_timeout=300)
    at.run(); at.button(key="nv_" + "SEM / TEM / EDX").click().run()
    return at


@needs_ocr
def test_calibration_panel_in_app():
    at = _goto_sem()
    assert not at.exception
    met = {m.label: m.value for m in at.metric}
    assert met["Scale bar"] == "385 px" and met["Label (OCR)"] == "1 µm" and met["Metadata (nm/px)"] == "2.6"
    at.checkbox(key="sem_calex").check().run(); assert not at.exception
    hk = "sem_" + hashlib.md5("JEOL-style_noMetadata_x25000.tif".encode()).hexdigest()[:8]
    assert at.radio(key=f"{hk}_src").value == "scale bar (detected)"
    assert any("252 px = 3.968 nm/px" in m.value for m in at.markdown)
    at.number_input(key=f"{hk}_blen").set_value(250.0).run()               # override the detected length
    assert any("250 px = 4 nm/px" in m.value for m in at.markdown)
    at.radio(key=f"{hk}_src").set_value("scale bar (mark ends)").run(); assert not at.exception
    assert any("Select the scale bar length" in i.value for i in at.info)
    at.number_input(key=f"{hk}_m0").set_value(400.0).run(); at.number_input(key=f"{hk}_m1").set_value(652.0).run()
    assert any("252.0 px = 3.968 nm/px" in m.value for m in at.markdown)
    at.select_slider(key=f"{hk}_zoom").set_value(4).run(); assert not at.exception
    at.selectbox(key=f"{hk}_pre").set_value("500 nm").run()
    assert any("126 px = 500 nm" in c.value for c in at.caption)           # new bar uses this calibration
    at.checkbox(key=f"{hk}_magon").check().run(); assert not at.exception
    assert at.text_input(key=f"{hk}_magtxt").value == "×25 000"
