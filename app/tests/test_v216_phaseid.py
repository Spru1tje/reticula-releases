"""v2.16: phase identification against the reference library."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest
from mfs import phaseid as PID

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
EX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "example_data")


def pattern_from_sticks(entries, weights, lam=1.540593, scale=None, zero=0.0, noise=0.01, seed=1, fwhm=0.12, rng_=(10, 80)):
    x = np.arange(rng_[0], rng_[1], 0.02)
    y = np.zeros_like(x)
    for e, w, sc in zip(entries, weights, scale or [1.0] * len(entries)):
        tt, I = PID.sticks(e["path"], lam, rng_)
        s = np.clip(np.sin(np.radians(tt / 2)) / sc, 0, 0.999)
        t2 = 2 * np.degrees(np.arcsin(s)) + zero
        for pos, inten in zip(t2, I):
            y += w * inten * np.exp(-4 * np.log(2) * ((x - pos) / fwhm) ** 2)
    y = y / y.max() * 10000 + 150
    return x, np.random.default_rng(seed).poisson(y).astype(float)


def test_library_is_complete_and_loadable():
    lib = PID.library()
    names = {e["name"] for e in lib}
    assert len(lib) >= 16 and any("UiO-66" in n for n in names) and any("HKUST-1" in n for n in names) and any("Silicon" in n for n in names)
    assert all(os.path.exists(e["path"]) and e["source"] for e in lib)


def test_single_phase_with_shifted_cell_and_a_mixture():
    lib = {e["name"].split(" ")[0]: e for e in PID.library()}
    cef = [e for e in PID.library() if e["name"].startswith("Ceria")][0]
    x, y = pattern_from_sticks([cef], [1.0], scale=[1.004], zero=0.05)                  # cell 0.4 % larger, zero shift 0.05°
    r = PID.match(x, y)
    assert r and r[0]["name"].startswith("Ceria") and r[0]["confidence"] in ("very likely", "likely") and abs(r[0]["cell_change_pct"] - 0.4) < 0.15
    corundum = [e for e in PID.library() if e["name"].startswith("Corundum")][0]
    zn = [e for e in PID.library() if e["name"].startswith("Zincite")][0]
    x, y = pattern_from_sticks([corundum, zn, cef], [1.0, 0.6, 0.4])
    r = PID.match(x, y, n_phases=4)
    found = [p["name"].split(" ")[0] for p in r]
    assert {"Corundum", "Zincite", "Ceria"} <= set(found), found
    assert len(r) <= 4 and r[0]["explained"] > 0.2


def test_silicon_and_its_cell_change():
    si = [e for e in PID.library() if e["name"].startswith("Silicon")][0]
    x, y = pattern_from_sticks([si], [1.0], scale=[0.9985], zero=-0.03)
    r = PID.match(x, y)
    assert r and r[0]["name"].startswith("Silicon") and abs(r[0]["cell_change_pct"] + 0.15) < 0.06 and abs(r[0]["zero"] + 0.03) < 0.04


def test_unknown_pattern_gives_no_confident_match():
    rng = np.random.default_rng(5)
    x = np.arange(10, 80, 0.02)
    y = np.zeros_like(x)
    for pos in rng.uniform(12, 78, 14):                                                   # peaks at random places: not in the library
        y += rng.uniform(0.3, 1.0) * np.exp(-4 * np.log(2) * ((x - pos) / 0.15) ** 2)
    y = np.random.default_rng(6).poisson(y / y.max() * 5000 + 100).astype(float)
    r = PID.match(x, y)
    assert not r or r[0]["confidence"] == "possible"


@pytest.mark.skipif(not os.path.exists(os.path.join(DATA, "FAP_REAL.xra")), reason="real data not present")
def test_real_laboratory_pattern_is_identified_as_fluorapatite():
    L = open(os.path.join(DATA, "FAP_REAL.xra")).read().splitlines()
    vals = []
    for ln in L[2:]:
        vals += [float(ln[8 * k:8 * k + 8]) for k in range(10) if ln[8 * k:8 * k + 8].strip()]
    y = np.array(vals[:5753])
    x = 15.0 + 0.02 * np.arange(len(y))
    r = PID.match(x, y)
    assert r and r[0]["name"].startswith("Fluorapatite") and r[0]["score"] > 0.9 and r[0]["confidence"] == "very likely"


def test_example_data_mof801():
    from mfs import readers as R
    f = os.path.join(EX, "PXRD", "MOF801_STOE.xy")
    d = R.read_pxrd(open(f, "rb"), f)
    r = PID.match(d.x.values, d.y.values)
    assert r and r[0]["name"].startswith("MOF-801") and r[0]["score"] > 0.8
