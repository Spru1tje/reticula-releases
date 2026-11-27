"""PL thermometry on user windows (mfs.plwindows): synthetic spectra with a KNOWN ratio. Only the data inside the windows are used.

The fit is a JOINT fit of both windows with both bands in both windows, so the neighbour band's tail is modelled. The synthetic data
use the same pseudo-Voigt shape as the model; a shape mismatch in real spectra is NOT tested here."""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest

from mfs import plwindows as PW
from mfs.plratio import pv_unit

W1, W2 = (512.0, 540.0), (552.0, 585.0)


def _spectrum(A1, A2, seed, c1=524.0, c2=566.0, w1=9.0, w2=10.0, scale=1.0):
    rng = np.random.default_rng(seed)
    x = np.arange(500.0, 600.0, 0.4)
    y = A1 * pv_unit(x, c1, w1, 0.5) + A2 * pv_unit(x, c2, w2, 0.5)
    y = scale * y + 40.0 + 0.2 * (x - 500.0)
    return x, rng.poisson(y).astype(float)


def _noise_free(A1=3000.0, A2=1000.0, scale=30.0):
    # 480-600 nm: both bands' tails (band 1 reaches ~497 nm) lie inside the measured range
    x = np.arange(480.0, 600.0, 0.4)
    y = scale * (A1 * pv_unit(x, 524.0, 9.0, 0.5) + A2 * pv_unit(x, 566.0, 10.0, 0.5)) + 40.0 + 0.2 * (x - 500.0)
    return x, y


def test_ratio_recovered_on_clean_synthetic_pair():
    # areas: I1 = 3 * I2 for a clean pair with plenty of counts
    x, y = _spectrum(3000.0, 1000.0, seed=1, scale=30.0)
    rows = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)
    r = rows[0]
    assert not any("poor fit" in f or "not significant" in f for f in r["flags"])
    assert abs(r["ratio"] - 3.0) / 3.0 < 0.03
    assert np.isfinite(r["ratio_err"]) and r["ratio_err"] > 0


def test_neighbour_tail_is_modelled_no_bias_on_noise_free_pair():
    # the previous separate-window fit was biased by +1.06 % here; the joint fit removes it (noise-free, so only the model can bias)
    x, y = _noise_free()
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
    assert abs(100 * (r["ratio"] - 3.0) / 3.0) < 0.05


def test_error_is_calibrated_over_seeds():
    # z = (ratio - truth) / ratio_err over 60 noise realisations should have spread ~1 (the error is neither too small nor too large)
    z = []
    for seed in range(60):
        x, y = _spectrum(3000.0, 1000.0, seed=100 + seed, scale=8.0)
        r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
        if np.isfinite(r["ratio_err"]):
            z.append((r["ratio"] - 3.0) / r["ratio_err"])
    z = np.asarray(z)
    assert len(z) >= 55
    assert 0.7 < np.std(z) < 1.3
    assert abs(np.mean(z)) < 0.5


def test_only_window_data_are_used():
    # changing the spectrum OUTSIDE both windows must not change the result
    x, y = _spectrum(3000.0, 1000.0, seed=7, scale=30.0)
    r1 = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
    y2 = y.copy()
    outside = (x < W1[0]) | ((x > W1[1]) & (x < W2[0])) | (x > W2[1])
    y2[outside] = 1e6 * np.arange(outside.sum()) % 97 + 5e5
    r2 = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y2}], W1, W2)[0]
    assert r1["ratio"] == pytest.approx(r2["ratio"], rel=1e-12)


def test_series_keeps_order_and_temperature():
    specs = []
    for k, T in enumerate([290.0, 310.0, 330.0]):
        x, y = _spectrum(3000.0 - 500 * k, 1000.0, seed=20 + k, scale=30.0)
        specs.append({"T": T, "x": x, "y": y})
    rows = PW.window_ratio_series(specs, W1, W2)
    assert [r["T"] for r in rows] == [290.0, 310.0, 330.0]
    assert rows[0]["ratio"] > rows[1]["ratio"] > rows[2]["ratio"]


def test_tail_inside_the_used_range_is_not_flagged():
    # windows wide enough that both bands' tails (band 1 to ~497 nm, band 2 to ~596 nm) lie inside the used range: no tail warning
    x, y = _noise_free()
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], (490.0, 540.0), (552.0, 600.0))[0]
    assert not any("tail" in f for f in r["flags"])


def test_band_tail_beyond_the_measured_range_is_flagged():
    # the measured range ends at 530 nm; band 1 (centre 524, FWHM 9) extends to ~551 nm, beyond the data: the area is an extrapolation
    x = np.arange(500.0, 530.0, 0.4)
    y = 30.0 * 3000.0 * pv_unit(x, 524.0, 9.0, 0.5) + 40.0 + 0.2 * (x - 500.0)
    rng = np.random.default_rng(3)
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": rng.poisson(y).astype(float)}], (512.0, 520.0), (521.0, 530.0))[0]
    assert any("I1: tail outside measured range" in f for f in r["flags"])


def test_too_few_points_is_reported_not_fitted():
    x, y = _spectrum(3000.0, 1000.0, seed=4, scale=30.0)
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], (524.0, 525.0), W2)[0]
    assert any("few points" in f for f in r["flags"])
    assert math.isnan(r["ratio"])
    assert any("ratio undefined" in f for f in r["flags"])


def test_no_signal_in_a_window_gives_no_ratio():
    x, y = _spectrum(3000.0, 0.0, seed=5, scale=30.0)     # band 2 absent
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
    assert math.isnan(r["ratio"])
    assert any("ratio undefined" in f for f in r["flags"])


def test_overlapping_windows_are_refused():
    x, y = _spectrum(3000.0, 1000.0, seed=6, scale=30.0)
    with pytest.raises(ValueError, match="overlap"):
        PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], (520.0, 560.0), (550.0, 580.0))


@pytest.mark.parametrize("bad", [None, (540.0, 512.0, 3.0), (np.nan, 540.0), (530.0, 530.0)])
def test_bad_window_input_is_refused(bad):
    x, y = _spectrum(3000.0, 1000.0, seed=8, scale=30.0)
    with pytest.raises(ValueError):
        PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], bad, W2)


def test_insignificant_band_gives_no_ratio_and_a_flag():
    x, y = _spectrum(3000.0, 0.0, seed=9, scale=30.0)     # band 2 absent: the joint fit has no area for band 2
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
    assert math.isnan(r["ratio"])
    assert any(f.startswith("I2:") and ("not significant" in f or "no area" in f) for f in r["flags"])


def test_window_systematic_is_reported_and_grows_with_the_error_budget():
    x, y = _spectrum(3000.0, 1000.0, seed=1, scale=30.0)
    r = PW.window_systematics([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
    assert r["n_variants"] == 8
    assert np.isfinite(r["ratio_sys"]) and r["ratio_sys"] > 0
    assert r["ratio_err_stat_sys"] == pytest.approx(math.hypot(r["ratio_err"], r["ratio_sys"]))
    assert r["ratio_err_stat_sys"] >= r["ratio_err"]
    assert "bias_note_pct" not in r                                  # no separate bias term any more: the neighbour tail is in the model


def test_window_systematic_stays_small_on_a_clean_pair():
    # the window choice changes the ratio by a few per cent on this pair, not by a factor
    x, y = _spectrum(3000.0, 1000.0, seed=1, scale=30.0)
    r = PW.window_systematics([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
    assert r["ratio_sys"] / r["ratio"] < 0.10


def test_window_systematic_skips_variants_that_overlap():
    # windows that touch each other: the +10 % shifts of window 1 overlap window 2 and are skipped, the rest are fitted
    x, y = _spectrum(3000.0, 1000.0, seed=2, scale=30.0)
    r = PW.window_systematics([{"T": 300.0, "x": x, "y": y}], (512.0, 552.0), (552.0, 585.0))[0]
    assert 0 < r["n_variants"] < 8


# CHARACTERISATION (not a correctness check): the bias that a band-SHAPE mismatch causes in the joint fit, noise-free, measured when the
# fit was built. The model fixes eta = 0.5; real bands with another eta or width give a bias of this size. A real spectrum is not tested.
SHAPE_MISMATCH = [  # (true eta, true FWHM 1, FWHM 2, bias % measured)
    (0.00, 9.0, 10.0, -1.125),
    (0.25, 9.0, 10.0, -0.453),
    (1.00, 9.0, 10.0, +1.135),
    (0.00, 6.0, 7.0, -3.345),
    (1.00, 14.0, 15.0, +2.817),
]


@pytest.mark.parametrize("eta,w1,w2,bias_pct", SHAPE_MISMATCH)
def test_characterised_shape_mismatch_bias(eta, w1, w2, bias_pct):
    # the FIXED-eta model (free_eta=False, the earlier default): this is the bias that the shape mismatch causes
    x = np.arange(500.0, 600.0, 0.4)
    y = 30.0 * (3000.0 * pv_unit(x, 524.0, w1, eta) + 1000.0 * pv_unit(x, 566.0, w2, eta)) + 40.0 + 0.2 * (x - 500.0)
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2, True, False)[0]
    assert abs(100 * (r["ratio"] - 3.0) / 3.0 - bias_pct) < 0.1


@pytest.mark.parametrize("eta,w1,w2,bias_pct", SHAPE_MISMATCH)
def test_free_eta_removes_the_shape_bias_on_noise_free_data(eta, w1, w2, bias_pct):
    # S-PL-1: eta is fitted. On data generated with the SAME pseudo-Voigt family (any eta) the bias is gone. A real band that is not a
    # pseudo-Voigt is NOT tested here, and the error grows (see the ratio_err of these cases)
    x = np.arange(500.0, 600.0, 0.4)
    y = 30.0 * (3000.0 * pv_unit(x, 524.0, w1, eta) + 1000.0 * pv_unit(x, 566.0, w2, eta)) + 40.0 + 0.2 * (x - 500.0)
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2, True, True)[0]  # opt-in free eta
    assert abs(100 * (r["ratio"] - 3.0) / 3.0) < 0.1


# REGRESSION GUARD ONLY: these values are the CURRENT output of the joint fit, pinned so that a change is noticed. They are NOT evidence
# that the ratio is correct: the synthetic data share the model's shape, and the residual error on real spectra is not known.
PINNED_SEEDS = {1: 0.009, 2: 0.464, 3: 0.916}


@pytest.mark.parametrize("seed,bias_pct", sorted(PINNED_SEEDS.items()))
def test_regression_pinned_relative_bias_on_known_overlap(seed, bias_pct):
    # pins the current output (regression guard, not a correctness check); a change of more than 0.5 percentage points is to be looked at
    x, y = _spectrum(3000.0, 1000.0, seed=seed, scale=30.0)
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
    assert abs(100 * (r["ratio"] - 3.0) / 3.0 - bias_pct) < 0.5


# S-PL-2: one joint fit of the series with SHARED band widths and eta. Assumption: the widths do not change with temperature.
def _series_const_width(n=5, w1=9.0, w2=10.0, scale=30.0, noise=False):
    x = np.arange(500.0, 600.0, 0.4)
    ratios = [3.0, 2.5, 2.0, 1.5, 1.0]
    out = []
    for k in range(n):
        A2 = 1000.0
        y = scale * (ratios[k] * A2 * pv_unit(x, 524.0, w1, 0.5) + A2 * pv_unit(x, 566.0, w2, 0.5)) + 40.0 + 0.2 * (x - 500.0)
        if noise:
            y = np.random.default_rng(10 + k).poisson(y).astype(float)
        out.append({"T": 290.0 + 10.0 * k, "x": x, "y": y})
    return out, ratios


def test_shared_fit_is_exact_on_noise_free_series_with_constant_widths():
    specs, ratios = _series_const_width()
    rows = PW.window_ratio_series_shared(specs, W1, W2)
    for r, truth in zip(rows, ratios):
        assert abs(100 * (r["ratio"] - truth) / truth) < 0.05


def test_shared_fit_error_is_smaller_than_per_spectrum_on_constant_width_series():
    specs, _ = _series_const_width()
    per = PW.window_ratio_series(specs, W1, W2)
    sh = PW.window_ratio_series_shared(specs, W1, W2)
    assert np.mean([r["ratio_err"] / r["ratio"] for r in sh]) < np.mean([r["ratio_err"] / r["ratio"] for r in per])


def test_shared_fit_biases_a_series_whose_width_changes_with_temperature():
    # characterisation of the stated assumption: a width that drifts 0.5 nm per step is NOT described by one shared width.
    x = np.arange(500.0, 600.0, 0.4)
    ratios = [3.0, 2.5, 2.0, 1.5, 1.0]
    specs = []
    for k in range(5):
        w1 = 8.0 + 0.5 * k
        y = 30.0 * (ratios[k] * 1000.0 * pv_unit(x, 524.0, w1, 0.5) + 1000.0 * pv_unit(x, 566.0, 10.0, 0.5)) + 40.0 + 0.2 * (x - 500.0)
        specs.append({"T": 290.0 + 10.0 * k, "x": x, "y": y})
    rows = PW.window_ratio_series_shared(specs, W1, W2)
    assert max(abs(100 * (r["ratio"] - t) / t) for r, t in zip(rows, ratios)) > 2.0


def test_shared_fit_returns_one_row_per_spectrum_in_order():
    specs, _ = _series_const_width(n=3)
    rows = PW.window_ratio_series_shared(specs, W1, W2)
    assert [r["T"] for r in rows] == [290.0, 300.0, 310.0]
    assert all("flags" in r and "ratio_err" in r for r in rows)


# S-Figuryte-1: a window with too few counts per point gives no ratio (no fit is trusted there)
def _dim_spectrum(seed=3):
    rng = np.random.default_rng(seed)
    x = np.arange(500.0, 600.0, 0.4)
    mean = 0.2 * (3000.0 * pv_unit(x, 524.0, 9.0, 0.5) + 1000.0 * pv_unit(x, 566.0, 10.0, 0.5)) + 0.5
    return x, rng.poisson(mean).astype(float)


def test_too_little_signal_gives_no_ratio_and_a_flag():
    x, y = _dim_spectrum()
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
    assert math.isnan(r["ratio"])
    assert any("too little signal" in f for f in r["flags"])


def test_too_little_signal_also_applies_to_the_shared_series_fit():
    x, y = _dim_spectrum()
    r = PW.window_ratio_series_shared([{"T": 300.0, "x": x, "y": y}, {"T": 310.0, "x": x, "y": y}], W1, W2)[0]
    assert math.isnan(r["ratio"])
    assert any("too little signal" in f for f in r["flags"])


def test_enough_signal_is_not_flagged():
    x, y = _spectrum(3000.0, 1000.0, seed=1, scale=30.0)
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
    assert not any("too little signal" in f for f in r["flags"])


def test_min_mean_counts_is_a_parameter_and_the_default_is_the_measured_value():
    # the default is the value measured on the synthetic series (see plwindows.py); a lower threshold lets the dim spectrum through
    assert PW.MIN_MEAN_COUNTS == 100.0
    x, y = _dim_spectrum()
    r = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2, True, False, 0.0)[0]
    assert not any("too little signal" in f for f in r["flags"])


# S-Figuryte-2: diagnosis of the neighbour band's share of each window's band signal
def test_neighbour_share_on_the_clean_pair_matches_the_known_overlap():
    x, y = _noise_free()
    d = PW.neighbour_share([{"T": 300.0, "x": x, "y": y}], (490.0, 540.0), (552.0, 600.0))[0]
    assert 0.03 < d["share2"] < 0.07        # band 1 reaches window 2 (about 4.9 % of that window's band signal)
    assert d["share1"] < 0.02                # band 2 barely reaches window 1
    assert not d["flags"]                    # below the 5 % flag


def test_neighbour_share_flags_a_window_dominated_by_the_neighbour():
    x, y = _noise_free()
    d = PW.neighbour_share([{"T": 300.0, "x": x, "y": y}], (490.0, 540.0), (552.0, 600.0), True, False, 0.03)[0]
    assert any(f.startswith("I2: neighbour band gives") for f in d["flags"])


# S-Figuryte-1 (opt-in): Poisson deviance weighting. Default stays sqrt; this checks the opt-in on fixed seeds.
def test_poisson_weighting_removes_the_low_count_bias_on_fixed_seeds():
    base = 3000.0 * pv_unit(np.arange(500.0, 600.0, 0.4), 524.0, 9.0, 0.5) + 1000.0 * pv_unit(np.arange(500.0, 600.0, 0.4), 566.0, 10.0, 0.5)
    x = np.arange(500.0, 600.0, 0.4)
    m1 = (x >= W1[0]) & (x <= W1[1])
    mu = (10.0 / base[m1].mean()) * base + 0.5
    rel = []
    for seed in range(1000, 1060):            # the 60 benchmark seeds; per-seed scatter is about 20 %, so the mean has an SE of about 2.6 %
        y = np.random.default_rng(seed).poisson(mu).astype(float)
        f = PW.fit_pair(x, y, W1, W2, None, True, False, "poisson")
        r, e, ok = PW._ratio_and_error(f["A1"], f["eA1"], f["A2"], f["eA2"], f["cov12"], True)
        if ok:
            rel.append(100.0 * (r - 3.0) / 3.0)
    assert len(rel) >= 55
    assert abs(np.mean(rel)) < 4.0           # measured +1.5 % (poisson) vs +20.6 % (sqrt) on these seeds, see tools/benchmark_plwindows.py


def test_default_weighting_is_still_sqrt():
    x, y = _noise_free()
    a = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2)[0]
    b = PW.window_ratio_series([{"T": 300.0, "x": x, "y": y}], W1, W2, True, False, 100.0, "sqrt")[0]
    assert a["ratio"] == b["ratio"]


# S-Figuryte-2: warm start from the previous temperature gives the same ratios (same minimum)
def test_warm_start_gives_the_same_ratios_as_the_cold_fit():
    x = np.arange(500.0, 600.0, 0.4)
    ratios = [3.0, 2.5, 2.0, 1.5, 1.0]
    specs = []
    for k, R in enumerate(ratios):
        base = 30.0 * (R * 1000.0 * pv_unit(x, 524.0, 9.0, 0.5) + 1000.0 * pv_unit(x, 566.0, 10.0, 0.5)) + 40.0 + 0.2 * (x - 500.0)
        specs.append({"T": 290.0 + 10.0 * k, "x": x, "y": np.random.default_rng(50 + k).poisson(base).astype(float)})
    cold = PW.window_ratio_series(specs, W1, W2)
    warm = PW.window_ratio_series(specs, W1, W2, warm_start=True)
    for c, w in zip(cold, warm):
        # measured largest relative difference 1.4e-7 on this series; the tolerance is 1e-6
        assert w["ratio"] == pytest.approx(c["ratio"], rel=1e-6)


# S-Figuryte-3: input check before the fit
def test_input_problems_flags_too_few_points_and_an_edge_peak():
    x = np.arange(500.0, 600.0, 0.4)
    y = 30.0 * (3000.0 * pv_unit(x, 524.0, 9.0, 0.5) + 1000.0 * pv_unit(x, 566.0, 10.0, 0.5)) + 40.0
    spec = [{"T": 300.0, "x": x, "y": y}]
    assert PW.input_problems(spec, W1, W2) == []
    msgs = PW.input_problems(spec, (524.0, 524.5), W2)
    assert any("window 1 holds" in m for m in msgs)
    msgs = PW.input_problems(spec, (500.0, 524.0), W2)   # the peak of band 1 (524 nm) is at the window edge
    assert any("lies on its edge" in m for m in msgs)


# S-Figuryte-2: the default weighting is a decision for the user; a change of it must be made on purpose (this test then fails)
def test_default_weighting_of_the_public_functions_is_sqrt():
    import inspect
    assert inspect.signature(PW.fit_pair).parameters["weighting"].default == "sqrt"
    assert inspect.signature(PW.window_ratio_series).parameters["weighting"].default == "sqrt"


# S-Figuryte-1: the shared-width fit flags a series whose band widths change; a constant-width series is not flagged
def test_shared_fit_flags_a_series_whose_width_changes():
    x = np.arange(500.0, 600.0, 0.4)
    specs = []
    for k in range(5):
        w1 = 8.0 + 0.5 * k
        y = 30.0 * (3000.0 * pv_unit(x, 524.0, w1, 0.5) + 1000.0 * pv_unit(x, 566.0, 10.0, 0.5)) + 40.0 + 0.2 * (x - 500.0)
        specs.append({"T": 290.0 + 10.0 * k, "x": x, "y": y})
    rows = PW.window_ratio_series_shared(specs, W1, W2)
    assert all(any("shared width not supported" in f for f in r["flags"]) for r in rows)


def test_shared_fit_does_not_flag_a_constant_width_series():
    specs, _ = _series_const_width()
    rows = PW.window_ratio_series_shared(specs, W1, W2)
    assert not any("shared width not supported" in f for r in rows for f in r["flags"])


# S-Figuryte-3: the fitted parameters per spectrum as text lines for the analysis record
def test_window_fit_record_gives_one_line_per_spectrum_with_the_fitted_shape():
    specs, _ = _series_const_width(n=3)
    lines = PW.window_fit_record(specs, W1, W2)
    assert len(lines) == 3
    assert lines[0].startswith("290 K: band 1 centre")
    assert "eta 0.50 (fixed)" in lines[0]


# S-Figuryte-1: the per-window fit is cached; the numbers are the same as the uncached fit
def test_fit_cache_gives_the_same_numbers_and_a_hit_on_the_second_call():
    x, y = _spectrum(3000.0, 1000.0, seed=11, scale=30.0)
    PW._FIT_CACHE.clear()
    PW._FIT_CACHE_STATS.update(hits=0, misses=0)
    a = PW.fit_pair_cached(x, y, W1, W2)
    b = PW.fit_pair_cached(x, y, W1, W2)
    c = PW.fit_pair(x, y, W1, W2)
    assert PW._FIT_CACHE_STATS["hits"] == 1 and PW._FIT_CACHE_STATS["misses"] == 1
    assert a is not b and a["flags1"] is not b["flags1"]   # copies: changing one result does not change the cache
    assert a["A1"] == pytest.approx(c["A1"], rel=1e-12) and a["A2"] == pytest.approx(c["A2"], rel=1e-12)


# S-Figuryte-3: the comparison (window vs global) is computed in one function and agrees with the cross-check
def test_comparison_table_agrees_with_the_crosscheck():
    specs = _series_const_width(n=3)[0]
    rows = PW.window_ratio_series(specs, W1, W2)
    T = [q["T"] for q in specs]
    cmp_ = PW.comparison_table(rows, T, [3.0, 2.5, 2.0], [0.5, 0.5, 0.5])
    cc = PW.crosscheck(rows, T, [3.0, 2.5, 2.0], [0.5, 0.5, 0.5])
    assert cmp_["n_compared"] == cc["n_compared"] and cmp_["n_disagree"] == cc["n_disagree"]
    assert len(cmp_["disagree"]) == 3 and len(cmp_["window ratio"]) == 3


# S-Figuryte-3: the cache key covers every input; a hit never returns a result computed with other options
def test_cache_key_covers_counts_weighting_eta_warm_start_and_windows():
    x, y = _spectrum(3000.0, 1000.0, seed=21, scale=30.0)
    PW._FIT_CACHE.clear()
    PW._FIT_CACHE_STATS.update(hits=0, misses=0)
    calls = [
        dict(),                                   # base
        dict(counts=False),                       # other counts option
        dict(weighting="poisson"),                # other weighting
        dict(free_eta=True),                      # other eta option
        dict(warm=[1000.0] * 10),                 # other warm start
        dict(w1=(510.0, 540.0)),                  # other window
    ]
    results = []
    for kw in calls:
        args = dict(x=x, y=y, w1=W1, w2=W2)
        args.update(kw)
        results.append(PW.fit_pair_cached(**args))
    assert PW._FIT_CACHE_STATS["misses"] == len(calls) and PW._FIT_CACHE_STATS["hits"] == 0
    # the weighting option changes the numbers, so the cache must not have mixed them up
    plain = PW.fit_pair(x, y, W1, W2)
    poisson = PW.fit_pair(x, y, W1, W2, None, True, False, "poisson")
    assert results[0]["A1"] == pytest.approx(plain["A1"], rel=1e-12)
    assert results[2]["A1"] == pytest.approx(poisson["A1"], rel=1e-12)
    assert results[0]["A1"] != pytest.approx(results[2]["A1"], rel=1e-6)


def test_cache_hits_only_for_identical_inputs():
    x, y = _spectrum(3000.0, 1000.0, seed=22, scale=30.0)
    PW._FIT_CACHE.clear()
    PW._FIT_CACHE_STATS.update(hits=0, misses=0)
    PW.fit_pair_cached(x, y, W1, W2)
    PW.fit_pair_cached(x, y.copy(), W1, W2)          # same values, different array object: still a hit
    PW.fit_pair_cached(x, y + 1.0, W1, W2)           # different data: a miss
    assert PW._FIT_CACHE_STATS == {"hits": 1, "misses": 2}


# the analysis record shows the weighting that was actually used
def test_window_fit_record_states_the_weighting_used():
    specs, _ = _series_const_width(n=2)
    assert all("weighting sqrt" in line for line in PW.window_fit_record(specs, W1, W2))
    assert all("weighting poisson" in line for line in PW.window_fit_record(specs, W1, W2, True, False, "poisson"))

