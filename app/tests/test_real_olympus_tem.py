"""The user's 4 real images: JEOL JEM-1010 TEM ('TEM1010'), Olympus SIS MegaView III camera, analySIS 5.0 palette TIFFs
with a burned-in black scale bar + label (no info strip, no printed magnification). Truth read visually/pixel-exactly."""
import os, pytest
from mfs import readers as R, semcal

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "olympus_sis")
#            file                   bar px  label nm  metadata nm/px          magnification
TRUTH = {"60K_MOF801SGC6.tif":  (354, 200.0, 1.55075160758898, 60000),
         "100K_MOF801SGC6.tif": (322, 500.0, 0.8733188801503962, 100000),
         "120K_MOFSGDC.tif":    (248, 2000.0, 0.7167629607361451, 120000),
         "250K_MOF801SGDC.tif": (279, 200.0, 0.33103385433308367, 250000)}
needs_ocr = pytest.mark.skipif(not semcal.ocr_available(), reason="Tesseract not installed")


@pytest.mark.parametrize("f", TRUTH)
def test_olympus_sis_metadata_and_palette(f):
    L, nm, pix, mag = TRUTH[f]
    d = R.read_sem(open(os.path.join(D, f), "rb"), f)
    assert d["pixel_nm"] == pytest.approx(pix) and d["meta"]["magnification"] == mag
    assert d["meta"]["calibration"].startswith("Olympus SIS") and d["meta"]["instrument"] == "TEM1010"
    a = d["image"]
    assert a[980, 1200] == 0 and a.mean() > 100          # palette mapped to grey: burned bar black, image bright


@pytest.mark.parametrize("f", TRUTH)
def test_bar_detected_pixel_exact(f):
    d = R.read_sem(open(os.path.join(D, f), "rb"), f)
    assert semcal.find_info_strip(d["image"]) is None    # overlay bar, no strip to crop
    b = semcal.detect_bar(d["image"], None)
    assert b["length_px"] == TRUTH[f][0] and (b["x0"], b["y0"], b["y1"]) == (1365 - TRUTH[f][0], 980, 981)


@needs_ocr
@pytest.mark.parametrize("f", TRUTH)
def test_calibrate_label_and_crosscheck(f):
    L, nm, pix, mag = TRUTH[f]
    d = R.read_sem(open(os.path.join(D, f), "rb"), f)
    c = semcal.calibrate(d["image"], d["pixel_nm"], d["meta"]["magnification"], d["meta"]["calibration"])
    assert c["method"].startswith("metadata") and c["nm_per_px"] == pytest.approx(pix)
    assert c["label"]["nm"] == nm
    sb = [x for x in c["candidates"] if x["method"].startswith("scale bar")][0]
    assert sb["nm_per_px"] == pytest.approx(nm / L)
    assert c["mag_ocr"] is None                            # analySIS printed no magnification
    assert any(l == "warn" for l, _ in c["notes"])          # burned bar disagrees with the metadata


def test_stale_bar_detected_across_images():
    info = {f: {"meta": pix, "bar": nm / L, "label": f"{nm:g} nm", "bar_px": L} for f, (L, nm, pix, mag) in TRUTH.items()}
    notes = semcal.stale_bar_notes(info)
    assert len(notes) == 4
    assert any("100K_MOF801SGC6.tif" in n.split(":")[0] and "60K_MOF801SGC6.tif" in n.split("matches")[1] for n in notes if "matches" in n)
    assert any(n.startswith("**250K") and "120K_MOFSGDC.tif" in n for n in notes)
