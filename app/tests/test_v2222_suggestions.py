"""v2.22.2 follow-up: (1) scale-label check against standard bar lengths, (2) precise sphere mode checked on synthetic
discs with a known radius (own model, not validated against the instrument). Synthetic images only; no private measurement data is used or named here.

Test environment (documentation only; no test depends on it beyond the skips already in the suite):
  - Tesseract is needed for the real OCR tests. These tests replace OCR with a fake, so they run without it.
  - threadpoolctl is optional for this module; without it test_v2195_plfit is skipped (seen with Python 3.12 venv, 2026-11-26).
  - The CI totals differ between machines when these optional packages or MFS_TEST_DTS* files differ; record them with the result."""
import numpy as np
import pytest
from scipy import special

from mfs import semcal as SC
from mfs import particles as PA

_BAR = {"x0": 993, "x1": 1365, "y0": 980, "y1": 981, "length_px": 372.0, "thickness_px": 2, "polarity": "black", "region_top": 774}


def _fake_label_ocr(monkeypatch, text):
    """Synthetic OCR: the single-line label passes return `text`, everything else returns nothing."""
    def fake(region, x_off=0, y_off=0, psm=11, scale=2, whitelist=None, ink=None):
        if scale == 3 and psm == 7 and whitelist and ink is not None:
            return [{"text": text, "x0": 1130.0, "y0": 985.0, "x1": 1200.0, "y1": 1000.0, "conf": 90.0, "line": (1, 1, 1)}]
        return []
    monkeypatch.setattr(SC, "ocr_words", fake)
    monkeypatch.setattr(SC, "ocr_available", lambda: True)
    monkeypatch.setattr(SC, "find_info_strip", lambda img8: None)
    monkeypatch.setattr(SC, "detect_bar", lambda img8, strip=None: dict(_BAR))


def _warnings(cal):
    return [t for lvl, t in cal["notes"] if lvl == "warn"]


# ------------------------------------------------------------------ (1) scale-label check
@pytest.mark.parametrize("text, standard", [("100 nm", True), ("200 nm", True), ("500 nm", True), ("1 µm", True)])
def test_standard_label_gets_no_scale_warning(monkeypatch, text, standard):
    _fake_label_ocr(monkeypatch, text)
    cal = SC.calibrate(np.zeros((1032, 1376), np.uint8))
    assert not any("not a standard bar length" in w for w in _warnings(cal))


def test_non_standard_label_is_flagged_and_value_is_unchanged(monkeypatch):
    _fake_label_ocr(monkeypatch, "400 nm")
    cal = SC.calibrate(np.zeros((1032, 1376), np.uint8))
    assert cal["label"]["value"] == 400.0                                 # the value is used as read ...
    assert cal["nm_per_px"] == pytest.approx(400.0 / 372.0)
    flagged = [w for w in _warnings(cal) if "not a standard bar length" in w]
    assert len(flagged) == 1 and "400 nm" in flagged[0]                  # ... and the user is asked to check it


def test_non_standard_label_with_metadata_is_not_flagged_here(monkeypatch):
    _fake_label_ocr(monkeypatch, "400 nm")
    cal = SC.calibrate(np.zeros((1032, 1376), np.uint8), meta_pixel_nm=1.0, meta_source="test")
    assert not any("not a standard bar length" in w for w in _warnings(cal))


# ------------------------------------------------------------------ (2) precise sphere mode, synthetic discs
def _discs(R, s_blur=1.0, H=200, W=200, centers=((50, 50), (150, 50), (50, 150), (150, 150)), noise=2.0, seed=0):
    """Blurred absorbing discs of known radius R (px): absorption = 0.6 x (50 % point at R, edge width s_blur)."""
    yy, xx = np.mgrid[:H, :W]
    absorb = np.zeros((H, W))
    for cy, cx in centers:
        d = np.hypot(yy - cy, xx - cx)
        absorb += 0.6 * 0.5 * special.erfc((d - R) / (np.sqrt(2) * s_blur))
    bg = 200 - 20 * xx / W
    img = bg * np.exp(-absorb) + noise * np.random.default_rng(seed).standard_normal((H, W))
    return np.clip(img, 0, 255).astype(np.uint8)


@pytest.mark.parametrize("R", [10.0, 14.0])
def test_precise_mode_within_tolerance_on_synthetic_discs(R):
    # Tolerance check only: each radius must lie within 1.0 px of the known R. This is not a validation of the
    # instrument radius; the discs come from our own absorption model.
    img = _discs(R)
    lab, f, hi, info = PA.segment_adaptive(img, sphere="precise", clear_border=True)
    radii = [v["d_px"] / 2 for k, v in info.items() if isinstance(k, int) and isinstance(v, dict) and v.get("d_px")]
    assert len(radii) == 4, radii
    assert all(abs(r - R) < 1.0 for r in radii), (R, radii)                # tolerance 1.0 px


def test_precise_mode_returns_no_radius_on_empty_film():
    rng = np.random.default_rng(5)
    img = np.clip(200 + 2.0 * rng.standard_normal((200, 200)), 0, 255).astype(np.uint8)
    lab, f, hi, info = PA.segment_adaptive(img, sphere="precise", clear_border=True)
    assert not any(isinstance(v, dict) and v.get("d_px") for k, v in info.items() if isinstance(k, int))
