"""OCR scale label versus metadata pixel size: a warning with the numbers (metadata value, measured bar length) when they
disagree by more than 2 %. Nothing is changed automatically. Synthetic images and a fake OCR only; no private data."""
import numpy as np
import pytest

from mfs import semcal as SC


X0 = 500          # bar starts here: it must not touch the picture edge (detect_bar rejects edge-touching bars)


def _strip(bar_px, extra_px=0):
    """Dark info strip with a 2 px thick white bar of bar_px at row 980; optional touching light stroke of extra_px."""
    img = np.full((1032, 1376), 30, np.uint8)
    img[980:982, X0:X0 + bar_px] = 255
    if extra_px:
        img[980:982, X0 + bar_px:X0 + bar_px + extra_px] = 255
    return img


def _fake_label_ocr(monkeypatch, text, bar_px=240):
    cx = X0 + bar_px / 2                                  # label box centred on the bar
    def fake(region, x_off=0, y_off=0, psm=11, scale=2, whitelist=None, ink=None):
        if scale == 3 and psm == 7 and whitelist and ink is not None:
            return [{"text": text, "x0": cx - 35, "y0": 985.0, "x1": cx + 35, "y1": 1000.0, "conf": 90.0, "line": (1, 1, 1)}]
        return []
    monkeypatch.setattr(SC, "ocr_words", fake)
    monkeypatch.setattr(SC, "ocr_available", lambda: True)
    monkeypatch.setattr(SC, "find_info_strip", lambda img8: None)


def _mismatch(cal):
    return [t for lvl, t in cal["notes"] if lvl == "warn" and "does not match the metadata pixel size" in t]


def test_agreeing_label_gives_no_mismatch_warning(monkeypatch):
    _fake_label_ocr(monkeypatch, "100 nm")
    cal = SC.calibrate(_strip(240), meta_pixel_nm=100 / 240)       # 240 px bar is exactly 100 nm
    assert cal["label"]["value"] == 100.0
    assert _mismatch(cal) == []


def test_outlier_100nm_with_touching_stroke_is_flagged_with_numbers(monkeypatch):
    _fake_label_ocr(monkeypatch, "100 nm", bar_px=240)
    cal = SC.calibrate(_strip(240, extra_px=35), meta_pixel_nm=0.4174)   # metadata expects ~240 px; bar measures 275 px
    msg = _mismatch(cal)
    assert len(msg) == 1
    assert "0.4174 nm/px" in msg[0] and "275 px" in msg[0]
    assert cal["bar"]["length_px"] == pytest.approx(275.0)
    assert cal["label"]["value"] == 100.0                  # label kept as read
    assert cal["nm_per_px"] == pytest.approx(0.4174)        # metadata stays the value used; nothing changed automatically


def test_outlier_1um_read_as_7um_is_flagged_and_not_changed(monkeypatch):
    _fake_label_ocr(monkeypatch, "7 µm", bar_px=418)
    cal = SC.calibrate(_strip(418), meta_pixel_nm=2.3941)                # 418 px x 2.3941 nm/px = 1 µm
    msg = _mismatch(cal)
    assert len(msg) == 1 and "418 px" in msg[0] and "2.394 nm/px" in msg[0]
    assert cal["label"]["value"] == 7.0                    # label kept as read
    assert cal["nm_per_px"] == pytest.approx(2.3941)        # metadata stays the value used; nothing changed automatically


def test_correct_1um_label_gives_no_warning(monkeypatch):
    _fake_label_ocr(monkeypatch, "1 µm", bar_px=418)
    cal = SC.calibrate(_strip(418), meta_pixel_nm=2.3941)
    assert cal["label"]["value"] == 1.0 and _mismatch(cal) == []


def test_no_duplicate_generic_cross_check_for_the_scale_bar(monkeypatch):
    _fake_label_ocr(monkeypatch, "100 nm", bar_px=240)
    cal = SC.calibrate(_strip(240, extra_px=35), meta_pixel_nm=0.4174)
    generic = [t for _, t in cal["notes"] if t.startswith("Cross-check") and "scale bar" in t]
    assert generic == []
