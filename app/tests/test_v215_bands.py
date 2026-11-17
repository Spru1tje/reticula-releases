"""v2.15: band assignment + multi-band fit on synthetic spectra with known bands."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from mfs import bands as B


def spectrum(seed=2):
    x = np.arange(1300.0, 1800.0, 1.0)
    y = 0.02 + 1e-5 * (x - 1300)
    true = [(1710.0, 18.0, 0.8), (1655.0, 22.0, 0.5), (1585.0, 30.0, 1.0), (1395.0, 28.0, 0.9)]
    for c, w, h in true:
        y = y + h * (0.5 / (1 + 4 * ((x - c) / w) ** 2) + 0.5 * np.exp(-4 * np.log(2) * ((x - c) / w) ** 2))
    return x, y + np.random.default_rng(seed).normal(0, 0.003, x.size), true


def test_multiband_fit_recovers_overlapping_bands():
    x, y, true = spectrum()
    r = B.fit_bands(x, y, [1708, 1650, 1590, 1400], 1330, 1760)
    for b, (c, w, h) in zip(r["bands"], sorted(true)):
        assert abs(b["center"] - c) < 1.5 and abs(b["fwhm"] - w) / w < 0.1 and abs(b["height"] - h) / h < 0.08, (b, c, w, h)
        assert b["center_err"] > 0
    assert r["r2"] > 0.999
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    B.draw_fit(ax, r)
    plt.close(fig)


def test_assignment_own_rows_first_and_ranges():
    a = B.assign([1585.0, 1395.0, 2925.0], technique="ir")
    assert "ν_as(COO⁻)" in a[0]["candidates"][0][0] and "ν_s(COO⁻)" in a[1]["candidates"][0][0] and "aliphatic" in a[2]["candidates"][0][0]
    own = [{"from (cm⁻¹)": 1580, "to (cm⁻¹)": 1595, "assignment": "my sulfonate–Zr band"}]
    b = B.assign([1585.0], own, "ir")[0]
    assert b["candidates"][0][1] == "own" and b["ambiguous"]
    assert B.assign([1000.0], technique="raman")[0]["candidates"][0][0].startswith("ring breathing")
    assert B.assign([5000.0])[0]["candidates"] == []
