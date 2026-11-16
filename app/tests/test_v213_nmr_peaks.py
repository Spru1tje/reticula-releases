"""v2.13 M2: 1H peak picking, multiplets, integrals, ACS string (public Bruker test data, see tests/data/nmr/README)."""
import sys, pathlib
import numpy as np
import pytest
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from mfs import nmrproc as N, nmrpeaks as P

D = pathlib.Path(__file__).parent / "data" / "nmr"


def _an(fn, **kw):
    e = N.load((D / fn).read_bytes(), fn)[0]
    sp = N.process(e, **kw)
    ms = P.analyse(sp); P.normalise(ms)
    return sp, ms


def test_multiplet_patterns_synthetic():
    assert P.analyse_multiplet([0.0], [1.0])[0] == "s"
    assert P.analyse_multiplet([0, 7.0], [1, 1])[:2] == ("d", [7.0])
    nm, J, fit = P.analyse_multiplet([0, 7, 14], [1, 2, 1]); assert nm == "t" and abs(J[0] - 7) < 0.2
    nm, J, _ = P.analyse_multiplet([0, 2, 8, 10], [1, 1, 1, 1]); assert nm == "dd" and J == [8.0, 2.0]
    nm, J, _ = P.analyse_multiplet([0, 7, 14, 21], [1, 3, 3, 1]); assert nm == "q"


def test_aspirin_acs():
    sp, ms = _an("aspirin_300_XWIN35.zip")
    s = P.acs_string(ms, mhz=sp.info["spectrometer (MHz)"], solvent=sp.info["solvent"])
    assert s.startswith("¹H NMR (300 MHz, CDCl3) δ")
    assert "(dd, J = 7.8, 1.7 Hz, 1H)" in s and "(td, J = 7.8, 1.7 Hz, 1H)" in s and "(s, 3H)" in s


def test_naphthoic_matches_topspin_multiplets():
    sp, ms = _an("naphtoic_500_TS35_dspfvs12.zip", reference="none")
    d = [m for m in ms if m.mult == "d" and (m.nH or 0) >= 0.9]
    assert len(d) == 4
    J = sorted(m.J[0] for m in d)
    assert np.allclose(J, [7.2, 8.2, 8.2, 8.6], atol=0.3)          # TopSpin multiplet.txt: 7.2, 8.2, 8.2, 8.6 Hz
    m3 = [m for m in ms if m.mult == "m" and 7.5 < m.centre < 8.0]
    assert m3 and round(m3[0].nH) == 3


def test_solvent_tagged():
    sp, ms = _an("aspirin_300_XWIN35.zip")
    notes = {m.note for m in ms}
    assert any("CDCl3" in n for n in notes)
    s = P.acs_string(ms, mhz=300, solvent="CDCl3")
    assert "7.26" not in s
    assert any(m.overlap for m in ms if m.note)          # aspirin H5 td hides under CHCl3 → hint


def test_normalise_reference_and_total():
    sp, ms = _an("aspirin_300_XWIN35.zip")
    i = [k for k, m in enumerate(ms) if m.mult == "s" and 2.0 < m.centre < 2.4][0]
    P.normalise(ms, ref_index=i, ref_n=3)
    assert abs(ms[i].nH - 3) < 1e-6
    P.normalise(ms, total_h=8)
    assert abs(sum(m.nH for m in ms if not m.note and m.nH and m.integral >= 0.06 * max(x.integral for x in ms if not x.note)) - 8) < 0.6
