"""Rigaku .HO (XRF measurement, object store 'Module SxSchema 3'): spectrum reader.

Layout used (found by search, not by fixed offsets):
  * the schema text 'Module SxSchema 3' (UTF-16) must occur in the first 4096 bytes; other versions are refused.
  * each spectrum object (RowIntensity, int count[4096]) has the length word 16396 (= 12 + 4096*4) 24 bytes before its
    data. Every such place is a candidate; a candidate is accepted only if its 4096 int32 values are non-negative,
    below 1e7 each, sum below 1e9, with at least 100 non-zero channels.
  * numCh (int, 4096) and maxenergy (double, keV) are found as a pair: int 4096 followed by a double in (0, 200].
    The whole file is searched; pairs with different maxenergy values raise an error.
    Only 4-byte-aligned positions are searched.

Output (no default reference):
  * spectrum_1 = the first accepted candidate in file order.
  * spectrum_2 = a second accepted candidate with the same peak channel, height within 10 % and sum within 10 % of
    spectrum_1; present only if found. This matching is a heuristic; shape is not checked.
  * peak_channel_1: the peak channel of spectrum_1 (int); flags: spectrum_2_found (True/False).
  * extra_candidates: other accepted candidates, unlabelled.
The choice of which spectrum is the REFERENCE is left to the user for each measurement. The reader does not decide it.

Not done: the channel-to-energy mapping. The energy calibration is UNDETERMINED (no reference line was verified),
so no keV axis is derived.
"""
from __future__ import annotations

import re
import struct

import numpy as np
import pandas as pd

SCHEMA_TAG = "Module SxSchema 3"
LEN_WORD = 16396                 # 12 + 4096 * 4 bytes of a RowIntensity object
DATA_AFTER_LEN = 24              # data starts 24 bytes after the length word
NCH = 4096


class RigakuHOError(ValueError):
    pass


def _schema_ok(b: bytes) -> bool:
    head = b[:4096].decode("utf-16-le", errors="ignore")
    return re.search(re.escape(SCHEMA_TAG) + r"(?![\d.\w])", head) is not None


def _candidates(b: bytes):
    """All accepted spectrum arrays (offset, int64 array) in file order."""
    out = []
    for p in range(0, len(b) - 4, 4):
        if struct.unpack_from("<I", b, p)[0] != LEN_WORD:
            continue
        d = p + DATA_AFTER_LEN
        if d + 4 * NCH > len(b):
            continue
        a = np.frombuffer(b, "<i4", NCH, d)
        if a.min() < 0 or a.max() >= 1e7 or a.sum() >= 1e9 or (a > 0).sum() < 100:
            continue
        out.append((d, a.astype(np.int64)))
    return out


def _num_ch_and_maxenergy(b: bytes):
    vals = set()
    for p in range(0, len(b) - 11, 4):
        if struct.unpack_from("<i", b, p)[0] == NCH:
            e = struct.unpack_from("<d", b, p + 4)[0]
            if 0 < e <= 200:
                vals.add(float(e))
    if not vals:
        raise RigakuHOError("numCh / maxenergy pair not found in this .HO file")
    if len(vals) > 1:
        raise RigakuHOError(f"ambiguous maxenergy: {len(vals)} different (numCh, maxenergy) pairs; not read")
    return NCH, vals.pop()


def read_rigaku_ho(src) -> dict:
    """Read the spectra of a .HO file. `src` is a path or bytes. Raises RigakuHOError when the schema differs, when no
    spectrum is found, or when the energy setting is ambiguous."""
    if isinstance(src, (bytes, bytearray)):
        b = bytes(src)
    else:
        with open(src, "rb") as fh:
            b = fh.read()
    if not _schema_ok(b):
        raise RigakuHOError(f"not a Rigaku .HO file of schema '{SCHEMA_TAG}' (other version or other format); not read")
    n, max_e = _num_ch_and_maxenergy(b)
    cands = _candidates(b)
    if not cands:
        raise RigakuHOError("no spectrum (RowIntensity int count[4096]) found in this .HO file")
    off_1, s1 = cands[0]
    peak_1, max_1, sum_1 = int(s1.argmax()), int(s1.max()), int(s1.sum())
    off_2, s2 = None, None
    for off, arr in cands[1:]:
        if (int(arr.argmax()) == peak_1 and abs(int(arr.max()) - max_1) <= 0.10 * max_1
                and abs(int(arr.sum()) - sum_1) <= 0.10 * sum_1):
            off_2, s2 = off, arr
            break
    cols = {"spectrum_1": s1}
    if s2 is not None:
        cols["spectrum_2"] = s2
    df = pd.DataFrame({"channel": np.arange(n), **{k: v[:n] for k, v in cols.items()}})
    extra = [{"offset": off, "peak_channel": int(arr.argmax()), "max_counts": int(arr.max())}
             for off, arr in cands[1:] if off != off_2]
    return {"numCh": n, "maxenergy_keV": max_e, "spectra": df,
            "peak_channel_1": peak_1, "flags": {"spectrum_2_found": s2 is not None},
            "spectrum_1_offset": off_1, "spectrum_2_offset": off_2, "extra_candidates": extra,
            "note": "reference choice left to the user per measurement; energy calibration UNDETERMINED"}
