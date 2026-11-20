"""v2.15.1: PL extras."""
import io, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from mfs import plextras as PX


def test_correction_file_and_application():
    for txt in ("# wl factor\n400 1.0\n500 2.0\n600 4.0\n", "400;1,0\n500;2,0\n600;4,0\n", "400,1.0\n500,2.0\n600,4.0\n"):
        if ";1,0" in txt:
            txt = txt.replace(",", ".")
        wl, f = PX.read_correction(txt.encode())
        assert list(wl) == [400, 500, 600] and list(f) == [1, 2, 4]
    y = PX.apply_correction([450, 550, 700], [10, 10, 10], wl, f)
    assert np.allclose(y, [15, 30, 40])
    try:
        PX.read_correction(b"400 1\n500 2\n")
        raise AssertionError("read_correction accepted non-monotonic input")
    except ValueError:
        pass


def test_repeat_uncertainty():
    T = np.array([293, 293.1, 293, 313, 313.2, 323.0])
    D = np.array([0.30, 0.31, 0.29, 0.40, 0.41, 0.45])
    r = PX.repeat_uncertainty(T, D, lambda t: 100 * 0.092 / (8.617e-5 * t ** 2))
    assert [x["n"] for x in r] == [3, 2] and abs(r[0]["mean"] - 0.30) < 1e-9
    assert abs(r[0]["std"] - 0.01) < 1e-9 and r[0]["sigma_T"] > 0 and abs(r[0]["rel_pct"] - 100 * 0.01 / 0.30) < 1e-6


def test_power_parse_and_fit():
    assert PX.parse_power("sample_250mW.csv") == 0.25 and PX.parse_power("P = 1.5 W") == 1.5 and PX.parse_power("a 12 mW b") == 0.012
    assert PX.parse_power("no power here") is None
    P = np.array([5, 10, 20, 40, 80, 160, 320.0])
    I = 3.0 * P ** 1.9 * (1 + np.random.default_rng(1).normal(0, 0.01, P.size))
    f = PX.fit_power(P, I)
    assert abs(f["n"] - 1.9) < 0.05 and f["n_err"] > 0 and f["R2"] > 0.999 and not f["saturating"]
    Isat = I / (1 + (P / 150.0) ** 2)
    assert PX.fit_power(P, Isat)["saturating"]
    assert PX.fit_power([1, 2], [1, 2])["status"] != "ok"


def test_cie_known_colours():
    wl = np.arange(380.0, 781.0, 1.0)
    for lam, (x, y) in {520: (0.0743, 0.8338), 600: (0.6270, 0.3725), 460: (0.1440, 0.0297)}.items():
        I = np.exp(-0.5 * ((wl - lam) / 0.6) ** 2)
        c = PX.cie_xy(wl, I, photon_counts=False)
        assert abs(c["x"] - x) < 0.02 and abs(c["y"] - y) < 0.02, (lam, c)
    w = PX.cie_xy(wl, np.ones_like(wl), photon_counts=False)
    assert abs(w["x"] - 1 / 3) < 0.015 and abs(w["y"] - 1 / 3) < 0.015 and 4500 < w["CCT"] < 6500
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(); PX.draw_cie(ax, [(0.3, 0.3)]); plt.close(fig)
