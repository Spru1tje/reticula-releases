"""v2.15: robust ratio of overlapping / dirty luminescence bands, on synthetic series with a known answer."""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from mfs import plratio as PR

KB = 8.617333262e-5


def series(overlap=1.0, shoulder=True, drift=True, seed=4, scale=30000.0, temps=None):
    rng = np.random.default_rng(seed)
    x = np.arange(500.0, 580.01, 0.4)
    temps = np.arange(293.0, 334.0, 5.0) if temps is None else temps
    out, true = [], []
    for T in temps:
        t = (T - temps.mean()) / np.ptp(temps)
        w1, w2 = 9.0 + (1.5 * t if drift else 0), 10.0 + (1.0 * t if drift else 0)
        c1, c2 = 524.0 + (0.3 * t if drift else 0), 524.0 + overlap * 0.5 * (w1 + w2) / 1.0 + 0.0
        D = 10.7 * math.exp(-(741.0 / 8065.54) / (KB * T))
        a2 = 1.0
        a1 = D * a2
        y = a1 * PR.pv_unit(x, c1, w1, 0.4) + a2 * PR.pv_unit(x, c2, w2, 0.4)
        if shoulder:
            y = y + 0.10 * PR.pv_unit(x, 562.0, 14.0, 0.3) + 0.04 * PR.pv_unit(x, 500.0, 12.0, 0.2)
        y = y / y.max() * scale + 200 + 3.0 * (x - 500)
        out.append({"T": float(T), "x": x, "y": rng.poisson(y).astype(float)})
        true.append(D)
    return out, np.array(true)


def test_global_fit_recovers_ratio_where_windows_are_biased():
    sp, true = series(overlap=0.9)
    mean = np.mean([s["y"] for s in sp], axis=0)
    bands = PR.auto_bands(sp[0]["x"], mean, 505, 575)
    assert bands is not None and len(bands) >= 2
    r = PR.auto_fit(sp, counts=True, lo=500, hi=580)
    ratio = np.array([row["ratio"] for row in r["rows"]])
    assert np.nanmax(np.abs(ratio / true - 1)) < 0.06, (ratio / true)
    assert all(row["flag"] != "poor" for row in r["rows"])
    # windows split at the valley: strongly biased because of the overlap
    c = [b[0] for b in bands[:2]]
    mid = 0.5 * (c[0] + c[1])
    I1, I2 = PR.window_integration(sp, (c[0] - 25, mid), (mid, c[1] + 25))
    assert np.nanmedian(np.abs(I1 / I2 / true - 1)) > 0.08
    assert np.nanmedian(np.abs(ratio / true - 1)) < 0.03


def test_extra_band_removes_residual_and_suggestion_finds_it():
    sp, true = series(overlap=1.4, shoulder=True)
    mean = np.mean([s["y"] for s in sp], axis=0)
    bands = PR.auto_bands(sp[0]["x"], mean, 505, 575)[:2]
    r0 = PR.fit_global(sp, bands, window=(500, 580), counts=True, drift="none")
    sug = PR.suggest_extra_band(r0)
    assert sug is not None and (abs(sug[0] - 560) < 8 or abs(sug[0] - 500) < 10)
    r1 = PR.fit_global(sp, bands + [sug], window=(500, 580), counts=True, drift="none")
    assert r1["rss"] < 0.6 * r0["rss"]


def test_uncertainty_is_realistic_and_robust_check_adds_systematic():
    sp, true = series(overlap=1.2, seed=7)
    bands = PR.auto_bands(sp[0]["x"], np.mean([s["y"] for s in sp], axis=0), 505, 575)[:2]
    r = PR.robust_ratio(sp, bands, counts=True, main=PR.auto_fit(sp, counts=True, lo=500, hi=580))
    ratio = np.array([row["ratio"] for row in r["rows"]])
    tot = r["ratio_total"]
    z = (ratio - true) / np.maximum(tot, 1e-9)
    assert np.sqrt(np.mean(z ** 2)) < 2.0                      # errors are neither wildly too small …
    assert np.median(tot / ratio) < 0.06                       # … nor useless
    assert r["n_alternatives"] >= 2 and np.all(r["ratio_sys"] >= 0)


def test_model_choice_and_quality_flags():
    sp, true = series(overlap=1.3, drift=False, shoulder=False)
    bands = PR.auto_bands(sp[0]["x"], np.mean([s["y"] for s in sp], axis=0), 505, 575)[:2]
    r = PR.fit_global(sp, bands, counts=True)
    assert r["model"] == "shared shape"                        # no drift in the data: the simpler model wins on BIC
    sp2, _ = series(overlap=1.3, drift=True, shoulder=False, scale=200000.0)
    r2 = PR.fit_global(sp2, bands, counts=True)
    assert r2["model"] == "drifting shape"
    bad = [dict(s) for s in sp]
    rng = np.random.default_rng(1)
    bad[3] = dict(bad[3], y=bad[3]["y"] + 4000 * np.sin(bad[3]["x"] / 2.0) + rng.normal(0, 300, len(bad[3]["x"])))   # wavy / dirty spectrum
    rb = PR.fit_global(bad, bands, counts=True)
    assert rb["rows"][3]["flag"] in ("check", "poor") and rb["rows"][0]["flag"] == "good"
