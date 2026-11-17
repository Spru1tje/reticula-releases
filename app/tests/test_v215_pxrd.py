"""v2.15: crystallite size on synthetic patterns with a known answer."""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from mfs import crystsize as CS, xrd as X

LAM = X.CU_KA1


def pattern(D_nm, eps, inst=0.05, seed=1, noise=0.004, K=0.9):
    """fcc-like peak positions; FWHM = Lorentzian sum of size + strain + instrument, pseudo-Voigt (eta 1 = Lorentzian)."""
    rng = np.random.default_rng(seed)
    x = np.arange(15, 85, 0.02)
    y = np.zeros_like(x)
    for tt, h in [(28.4, 1000), (32.9, 480), (47.3, 700), (56.1, 520), (58.8, 150), (69.1, 110), (76.4, 200)]:
        th = math.radians(tt / 2)
        beta = math.degrees(K * LAM / (D_nm * 10 * math.cos(th)) + 4 * eps * math.tan(th))
        w = beta + inst
        y += h * (1 / (1 + ((x - tt) / (w / 2)) ** 2))
    return x, y + 20 + rng.normal(0, noise * 1000, x.size)


def test_auto_peaks_finds_the_peaks():
    x, y = pattern(30, 0.0)
    pk = CS.auto_peaks(x, y, 20, 80)
    assert len(pk) >= 6 and min(abs(p - 28.4) for p in pk) < 0.1


def test_scherrer_and_wh_recover_known_values():
    for D, eps in [(25.0, 0.0), (18.0, 0.003)]:
        x, y = pattern(D, eps)
        r = CS.analyse(x, y, shape="lorentzian", lo=20, hi=80)
        s, w = r["scherrer"], r["wh"]
        assert w is not None and len(r["peaks"]) >= 6
        if eps == 0.0:
            assert abs(s["D_mean"] - D) / D < 0.06, s
        assert abs(w["D_nm"] - D) / D < 0.10, (w["D_nm"], D)
        assert abs(w["strain"] - eps) < 0.0012, (w["strain"], eps)
        assert w["D_err_nm"] > 0 and (eps == 0 or w["r2"] > 0.9)
        if eps > 0:
            assert r["hw"] is not None and abs(r["hw"]["D_nm"] - D) / D < 0.25


def test_notes_for_large_crystals_and_few_peaks():
    x, y = pattern(2000.0, 0.0, noise=0.002)                           # peaks as narrow as the instrument
    r = CS.analyse(x, y, shape="lorentzian", lo=20, hi=80)
    assert r["scherrer"] is None or r["scherrer"]["D_mean"] > 100
    assert r["notes"]
    x, y = pattern(20.0, 0.0)
    r = CS.analyse(x, y, peaks=[28.4, 47.3], shape="lorentzian")
    assert r["wh"] is None and any("≥ 3 peaks" in n for n in r["notes"])
    assert "Scherrer" in CS.summary_line(r)
