"""v2.24.0: (1) the grid estimate of x works for hexagonal end members (round trip and status cases);
(2) MIL-53 (exclude_matching in index.json) is not offered by the phase-identification matching.
Synthetic numbers for (1); repo reference cells only. pymatgen and gemmi are required: no importorskip."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest
from mfs import vegard as VG

WL = 1.540593
CORUNDUM = (4.7606, 4.7606, 12.994, 90.0, 90.0, 120.0)     # repo reference, hexagonal
B_HEX = (4.9, 4.9, 13.5, 90.0, 90.0, 120.0)                # synthetic hexagonal end member B
HKLS = [(1, 0, 4), (1, 1, 0), (0, 1, 2), (1, 1, 3)]


def _peaks_for(x, cell_A=CORUNDUM, cell_B=B_HEX):
    cell = VG.cell_from_composition(x, cell_A, cell_B)
    return [(VG.two_theta_for_cell(cell, h, WL), h) for h in HKLS]


def test_hexagonal_round_trip_recovers_x():
    peaks = _peaks_for(0.37)
    r = VG.estimate_x_grid(peaks, WL, CORUNDUM, B_HEX)
    assert r["status"] == "ok"
    assert r["x"] == pytest.approx(0.37, abs=0.002)
    assert r["rms_deg"] < 1e-3


def test_pure_end_member_is_reported_as_end_member():
    r = VG.estimate_x_grid(_peaks_for(1.0), WL, CORUNDUM, B_HEX)
    assert r["status"] == "end member"
    assert r["x"] == pytest.approx(1.0, abs=0.002)


def test_peaks_not_explained_by_either_end_member_are_outside_range():
    # reflections of a cell much larger than both end members: no x can fit them
    big = (6.0, 6.0, 20.0, 90.0, 90.0, 120.0)
    peaks = [(VG.two_theta_for_cell(big, h, WL), h) for h in HKLS]
    r = VG.estimate_x_grid(peaks, WL, CORUNDUM, B_HEX)
    assert r["status"] == "outside range"


def test_one_peak_is_too_few():
    r = VG.estimate_x_grid(_peaks_for(0.5)[:1], WL, CORUNDUM, B_HEX)
    assert r["status"] == "too few peaks"


def test_mil53_is_excluded_from_phase_identification_matching():
    """A library entry flagged exclude_matching must not appear in phaseid.match results."""
    from mfs import phaseid as PID
    lib = PID.library()
    mil = [e for e in lib if e["file"] == "MIL-53_Cr.cif"]
    assert mil and mil[0].get("exclude_matching") is True
    # a pattern built from CeO2 peaks: MIL-53 must not be among the ranked phases, even when its entry is passed in
    x = np.linspace(10, 80, 4000)
    y = np.full_like(x, 50.0)
    for t in (28.5, 33.1, 47.5, 56.3, 69.4, 76.7):
        y += 800 * np.exp(-((x - t) / 0.12) ** 2)
    res = PID.match(x, y, lib=lib)
    assert all(r["file"] != "MIL-53_Cr.cif" for r in res)
