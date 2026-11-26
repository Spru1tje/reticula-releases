"""Calibration must read the unstretched array 'raw', never the contrast-stretched 'image' (a stretch can push a thin
white scale bar below the white threshold). Behaviour tests: the readers return 'raw', and calibration_panel calibrates
on it and refuses a reader result without it. Synthetic images and synthetic reader output only; no private data."""
import io

import numpy as np
import pytest
from PIL import Image

from mfs import readers as R


class _Stop(Exception):
    """Raised by the fake calibration to end calibration_panel right after the array was handed over."""


def _raw_with_bar():
    img = np.full((1032, 1376), 30, np.uint8)
    img[980:982, 500:872] = 255                    # 2 px thick white bar, 372 px long
    return img


def test_bar_is_found_in_raw_but_lost_after_a_halving_stretch():
    from mfs import semcal as SC
    raw = _raw_with_bar()
    assert SC.detect_bar(raw, None)["length_px"] == 372.0
    assert SC.detect_bar((raw * 0.5).astype(np.uint8), None) is None   # why the page must use 'raw'


def test_read_sem_returns_raw_for_a_picture_file():
    buf = io.BytesIO()
    Image.fromarray(_raw_with_bar()).save(buf, format="PNG")
    buf.seek(0)
    d = R.read_sem(buf, "bar.png")
    assert "raw" in d and d["raw"].dtype == np.uint8 and d["raw"].shape == d["image"].shape


def test_read_em_returns_raw_for_an_ncempy_file(monkeypatch):
    import ncempy.io as nio
    data = _raw_with_bar().astype(np.float32)

    def fake_read(path):                           # stands in for the ncempy reader; no real file is opened
        return {"data": data, "pixelSize": [0.5], "pixelUnit": ["nm"]}
    monkeypatch.setattr(nio, "read", fake_read)
    d = R.read_em(io.BytesIO(b"not-a-real-file"), "fake.dm4")
    assert "raw" in d and d["raw"].dtype == np.uint8 and d["raw"].shape == d["image"].shape


def test_calibration_panel_calibrates_on_raw_not_on_image(monkeypatch):
    import mfs.pages.microscopy as MS
    seen = {}

    def fake_cal(img8, pixel_nm, meta_mag, meta_src):
        seen["arr"] = img8
        raise _Stop
    monkeypatch.setattr(MS, "_sem_calibrate", fake_cal)
    raw = np.zeros((4, 4), np.uint8)
    img = np.full((4, 4), 9, np.uint8)
    m = {"image": img, "raw": raw, "pixel_nm": 1.0, "meta": {}, "databar_row": None}
    with pytest.raises(_Stop):
        MS.calibration_panel("hk", m, "bar.png")
    assert seen["arr"] is raw and seen["arr"] is not img


def test_calibration_panel_refuses_a_reader_result_without_raw(monkeypatch):
    import mfs.pages.microscopy as MS
    called = []
    monkeypatch.setattr(MS, "_sem_calibrate", lambda *a, **k: called.append(a) or None)
    m = {"image": np.full((4, 4), 9, np.uint8), "pixel_nm": 1.0, "meta": {}, "databar_row": None}
    with pytest.raises(KeyError, match="raw"):
        MS.calibration_panel("hk", m, "bar.png")
    assert called == []                            # no silent fallback to the stretched image
