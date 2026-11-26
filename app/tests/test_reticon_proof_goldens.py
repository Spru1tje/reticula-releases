"""Reticon proof set: current SVGs must match the stored goldens byte for byte.

A deliberate change means updating the golden copy in the same commit, with the reason in the message.
"""
from pathlib import Path

PROOF = Path(__file__).resolve().parents[1] / "assets" / "icons_proof"
GOLDEN = Path(__file__).resolve().parent / "golden" / "reticon_proof"


def test_proof_svgs_match_goldens():
    current = sorted(p.name for p in PROOF.glob("*.svg"))
    golden = sorted(p.name for p in GOLDEN.glob("*.svg"))
    assert current == golden, f"file set differs: current={current} golden={golden}"
    for name in current:
        assert (PROOF / name).read_bytes() == (GOLDEN / name).read_bytes(), f"{name} differs from golden"
