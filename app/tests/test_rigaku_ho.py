"""Rigaku .HO spectrum reader (mfs.rigaku_ho): spectra are located by the length marker, not by fixed offsets.

Synthetic tests (always run): buffers with the same layout, optionally shifted by padding. They check the mechanics:
schema refusal, search, output columns and flags, the second-spectrum match, ambiguous energy setting, no-spectrum error.
Real-file test: only if MFS_TEST_HO is set to a path (private measurement, not in the repository). It checks structure
only; no values or positions from that file are written here.
The reference spectrum is chosen by the user per measurement; the reader does not decide it.
The energy calibration is UNDETERMINED; the reader derives no keV axis.
"""
import os
import struct

import numpy as np
import pytest

from mfs import rigaku_ho as RH


def _buf(shift=0, second=True, schema=RH.SCHEMA_TAG, peak=300, height=90000):
    n = RH.NCH
    head = bytearray(256)
    head[:2 * len(schema)] = schema.encode("utf-16-le")
    body = bytearray(b"\x00" * shift + bytes(head))
    body += b"\x00" * (17452 + shift - len(body))
    body += struct.pack("<i", n) + struct.pack("<d", 50.0)           # numCh, maxenergy
    a = np.full(n, 100, dtype="<i4"); a[peak] = height; a[811] = 0
    for arr in ([a, a.copy()] if second else [a]):
        body += b"\x00" * 4 + struct.pack("<i", 0) + struct.pack("<i", RH.LEN_WORD)   # length word
        body += b"\x00" * (RH.DATA_AFTER_LEN - 4)                                       # 20 bytes to the data
        body += arr.tobytes()
    return bytes(body), a


def test_synthetic_finds_both_spectra_and_flags():
    raw, a = _buf()
    r = RH.read_rigaku_ho(raw)
    assert r["numCh"] == 4096 and r["maxenergy_keV"] == 50.0
    df = r["spectra"]
    assert list(df.columns) == ["channel", "spectrum_1", "spectrum_2"] and len(df) == 4096
    assert r["flags"] == {"spectrum_2_found": True} and r["peak_channel_1"] == 300
    assert int(df["spectrum_1"].argmax()) == 300 and int(df["spectrum_2"].argmax()) == 300


def test_no_default_copy_label_in_output():
    raw, _ = _buf()
    r = RH.read_rigaku_ho(raw)
    assert "copy_A" not in r["spectra"].columns and "reference" not in r


def test_shifted_offsets_are_still_found():
    raw, _ = _buf(shift=1000)
    r = RH.read_rigaku_ho(raw)
    assert int(r["spectra"]["spectrum_1"].argmax()) == 300 and r["flags"]["spectrum_2_found"]


def test_single_spectrum_has_no_second_column():
    raw, _ = _buf(second=False)
    r = RH.read_rigaku_ho(raw)
    assert "spectrum_2" not in r["spectra"].columns and r["flags"]["spectrum_2_found"] is False


def test_peak_channel_is_reported_as_int():
    raw, _ = _buf(peak=777, height=80000)
    r = RH.read_rigaku_ho(raw)
    assert r["peak_channel_1"] == 777 and isinstance(r["peak_channel_1"], int)


def test_same_peak_but_other_sum_is_not_spectrum_2():
    raw, _ = _buf(second=False)
    other = np.full(4096, 100, dtype="<i4"); other[300] = 90000; other[2000:2100] = 9000   # same peak, much larger sum
    raw = raw + b"\x00" * 20 + struct.pack("<i", RH.LEN_WORD) + b"\x00" * 20 + other.tobytes()
    r = RH.read_rigaku_ho(raw)
    assert r["flags"]["spectrum_2_found"] is False


@pytest.mark.parametrize("tag", ["Module SxSchema 30", "Module SxSchema 3.1", "Module SxSchema 3a", "Module SxSchema 4"])
def test_other_schema_spellings_are_refused(tag):
    raw, _ = _buf(schema=tag)
    with pytest.raises(RH.RigakuHOError, match="SxSchema 3"):
        RH.read_rigaku_ho(raw)


def test_ambiguous_energy_pairs_are_refused():
    raw, _ = _buf(second=False)
    raw = raw + struct.pack("<i", 4096) + struct.pack("<d", 30.0)
    with pytest.raises(RH.RigakuHOError, match="ambiguous"):
        RH.read_rigaku_ho(raw)


def test_no_spectrum_is_a_clear_error():
    raw, _ = _buf(second=False)
    raw = raw.replace(struct.pack("<I", RH.LEN_WORD), struct.pack("<I", 1234))
    with pytest.raises(RH.RigakuHOError, match="no spectrum"):
        RH.read_rigaku_ho(raw)


@pytest.mark.skipif(not os.environ.get("MFS_TEST_HO"), reason="no private .HO file (MFS_TEST_HO); not part of the repository")
def test_real_file_structure():
    with open(os.environ["MFS_TEST_HO"], "rb") as fh:
        r = RH.read_rigaku_ho(fh.read())
    assert r["numCh"] > 0 and r["maxenergy_keV"] > 0
    df = r["spectra"]
    assert "spectrum_1" in df.columns and len(df) == r["numCh"]
    assert set(r["flags"]) == {"spectrum_2_found"} and isinstance(r["peak_channel_1"], int)
