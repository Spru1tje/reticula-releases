"""v2.18.12: known-answer tests for the round's changes (Avrami standard errors, vectorised Tauc edge search, the
kinetics input checks, model selection for Sips data)."""
import numpy as np

from mfs import avrami as AV
from mfs import uvvis as UV
from mfs import isotherm_models as IM
from mfs import kinetics_nonisothermal as KN


# ---- Avrami: standard errors of n and ln k
def test_avrami_exact_data_has_zero_residual_and_recovers_n_and_k():
    t = np.linspace(1.0, 60.0, 30)
    X = 1 - np.exp(-0.02 * t ** 2.5)                      # exact: n = 2.5, k = 0.02
    r = AV.avrami_fit(t, X)
    assert abs(r["n"] - 2.5) < 1e-9
    assert abs(r["k"] - 0.02) < 1e-9
    assert r["n_se"] < 1e-8                               # no scatter, so no uncertainty


def test_avrami_standard_error_covers_the_true_exponent_on_noisy_data():
    rng = np.random.default_rng(12)
    t = np.linspace(1.0, 60.0, 40)
    X = np.clip(1 - np.exp(-0.02 * t ** 2.5) + rng.normal(0, 0.004, t.size), 0, 1)
    r = AV.avrami_fit(t, X)
    assert r["n_se"] > 0
    assert abs(r["n"] - 2.5) < 4 * r["n_se"]              # the true value lies within four standard errors


# ---- Tauc: the vectorised edge search matches the original per-window loop
def _tauc_edge_loop_reference(E, T, min_pts=6, min_frac=0.1):
    """The original (v2.18.11) window loop, kept here as the reference."""
    Tm = np.nanmax(T)
    best = None
    n = len(E)
    for w in sorted({min_pts, max(min_pts, n // 60), max(min_pts, n // 30), max(min_pts, n // 15)}):
        for i in range(0, n - w):
            seg = slice(i, i + w)
            ys = T[seg]
            if ys.min() < min_frac * Tm or ys.max() > 0.9 * Tm:
                continue
            s, b = np.polyfit(E[seg], ys, 1)
            if s <= 0:
                continue
            r2 = np.corrcoef(E[seg], ys)[0, 1] ** 2
            if r2 < 0.98:
                continue
            if best is None or s > best[0]:
                best = (s, b, r2, (float(E[seg].min()), float(E[seg].max())))
    return best


def test_tauc_noiseless_direct_edge_is_recovered_exactly():
    # (F * h nu)^2 = c (h nu - 2.5) above a direct gap of 2.5 eV (synthetic, no noise): the answer is 2.5 eV
    E = np.linspace(2.0, 4.5, 400)
    T = np.clip(E - 2.5, 0, None) * 3.0
    r = UV.band_gap(E, T)
    assert abs(r["Eg"] - 2.5) < 1e-6
    assert r["R2"] > 0.999999


def test_tauc_vectorised_search_gives_the_same_edge_as_the_reference_loop():
    rng = np.random.default_rng(7)
    E = np.linspace(1.8, 4.2, 240)
    T = np.clip(E - 2.9, 0, None) * 1.7 + rng.normal(0, 0.01, E.size)
    ref = _tauc_edge_loop_reference(E, T)
    new = UV.band_gap(E, T)
    assert ref is not None
    assert abs(new["Eg"] - (-ref[1] / ref[0])) < 1e-9
    assert abs(new["slope"] - ref[0]) < 1e-9 * abs(ref[0])
    assert abs(new["R2"] - ref[2]) < 1e-9
    assert new["window"] == ref[3]


# ---- Sips data is identified as Sips
def test_model_comparison_picks_sips_for_sips_data():
    p = np.geomspace(0.01, 10.0, 25)
    q = IM.sips(p, 5.0, 0.8, 0.6)                          # exact Sips, n = 0.6 (not Langmuir, not Toth)
    r = IM.compare_isotherms(p, q)
    assert r["status"] == "ok"
    assert r["best"] == "Sips"
    assert abs(r["fits"]["Sips"]["n"] - 0.6) < 1e-3
    assert abs(r["fits"]["Sips"]["qm"] - 5.0) < 1e-3


# ---- Kinetics input check: peak temperatures given in °C instead of K
def test_kissinger_warns_when_peak_temperatures_look_like_celsius():
    betas = [5.0, 10.0, 20.0]
    in_kelvin = KN.kissinger(betas, [600.0, 615.0, 630.0])
    assert in_kelvin["status"] == "ok"
    assert "warning" not in in_kelvin
    in_celsius = KN.kissinger(betas, [85.0, 100.0, 115.0])      # peaks entered in °C (85 to 115 °C), not K
    assert in_celsius["status"] == "ok"
    assert "273.15" in in_celsius.get("warning", "")


# ---- Avrami: the result carries a warning when the fit is not reliable
def test_avrami_warns_with_few_points_and_not_on_clean_data():
    t = np.linspace(0.5, 7.0, 30)                        # X runs from about 0.05 to 0.95 here (the fitting window)
    clean = AV.avrami_fit(t, 1 - np.exp(-0.02 * t ** 2.5))
    assert "warning" not in clean
    few = AV.avrami_fit([1.0, 2.0, 3.0, 4.0], [0.1, 0.2, 0.5, 0.8])
    assert "fewer" in few.get("warning", "") or "not reliable" in few.get("warning", "")
    rng = np.random.default_rng(3)
    noisy = AV.avrami_fit(t, np.clip(rng.uniform(0.05, 0.95, t.size), 0, 1))
    assert "warning" in noisy
