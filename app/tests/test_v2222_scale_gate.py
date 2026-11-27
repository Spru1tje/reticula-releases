"""Scale check before a figure export: blocks when the nm/px in use disagrees with the file's metadata by more than 2 %,
and when there is no calibration or no metadata. Pure logic; no image or private data is used."""
import pytest

from mfs import semcal as SC


def test_agreeing_scale_passes():
    ok, msg = SC.scale_gate(0.2743, 0.2742656749)
    assert ok and "passed" in msg


def test_mismatch_above_two_percent_blocks_with_the_numbers():
    ok, msg = SC.scale_gate(0.2743 * 1.05, 0.2743)
    assert not ok and "+5.0 %" in msg and "Confirm the scale to export" in msg


def test_exactly_two_percent_passes():
    ok, _ = SC.scale_gate(1.02, 1.0)
    assert ok


def test_confirmed_mismatch_passes_and_says_so():
    ok, msg = SC.scale_gate(1.05, 1.0, confirmed=True)
    assert ok and "confirmed by the user" in msg


def test_mismatch_within_two_percent_passes():
    ok, _ = SC.scale_gate(0.2743 * 1.015, 0.2743)
    assert ok


def test_no_metadata_blocks():
    ok, msg = SC.scale_gate(0.5, None)
    assert not ok and "no pixel-size metadata" in msg


@pytest.mark.parametrize("used", [None, 0, -1.0])
def test_no_calibration_blocks(used):
    ok, msg = SC.scale_gate(used, 0.5)
    assert not ok and "No scale is set" in msg
