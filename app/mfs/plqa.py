r"""v2.18.4: quality, confidence and uncertainty of photoluminescence results, in particular of ratiometric luminescence thermometry.

Everything here answers one question: CAN THIS NUMBER BE TRUSTED, and how well? The page shows the answers as a traffic light with plain reasons.
  * data checks: detector saturation (clipped / flat-topped peaks, analogue-to-digital limits), different settings inside one series;
  * information: does the ratio change with temperature by more than its error (otherwise there is no thermometer), is it monotonic;
  * physics checks: Boltzmann residuals (reduced chi-squared, runs of residual signs), the fitted energy gap against the spectroscopic gap of the two
    bands and against the usual range of thermally coupled levels;
  * the uncertainty budget of the energy gap (Monte Carlo over the statistical ratio errors and the temperature-sensor uncertainty; refits of every
    alternative processing for the systematic parts; spread between independent runs), and the two different temperature uncertainties:
    resolution (precision of one reading) and calibration (accuracy of the T(Delta) curve);
  * experiment checks: heating against cooling branch (hysteresis), laser heating from a power series, repeats;
  * a known-answer self-check: series simulated from the user's own fit (same bands, temperatures, noise) with a known energy gap are analysed
    again - the recovered gap and the spread of its errors tell whether the analysis is unbiased and its error bars honest;
  * an analysis record (JSON) with settings, file fingerprints and results for the supporting information.
numpy / scipy only."""
from __future__ import annotations

import datetime as _dt
import json
import math

import numpy as np

KB_EV = 8.617333262e-5
EV_TO_CM = 8065.54429
ADC_LIMITS = (4095, 16383, 65535, 262143, 1048575, 16777215)       # 12, 14, 16, 18, 20, 24-bit converters (array detectors)


def check(name, status, text, detail=""):
    """One line of the quality report. status: 'ok' | 'warn' | 'fail' | 'info'."""
    return {"name": name, "status": status, "text": text, "detail": detail}


ICON = {"ok": "✅", "warn": "⚠️", "fail": "⛔", "info": "ℹ️"}


def verdict(checks):
    """Overall confidence from the checks: 'high' (no warning), 'moderate' (warnings only), 'low' (at least one failure)."""
    st = [c["status"] for c in checks]
    if "fail" in st:
        return "low"
    if "warn" in st:
        return "moderate"
    return "high"


# ------------------------------------------------------------------------------------------------------------------------------ data
def saturation(spectra, rate_limit_cps=1.0e6):
    """Saturation of every spectrum: a flat top (three or more points at the maximum), a maximum at an analogue-to-digital limit, or a photon-counting rate
    above `rate_limit_cps` when the dwell time is known. spectra: [{'name', 'I', 'meta'}]. Returns [{'name', 'max', 'status', 'reason'}]."""
    out = []
    for s in spectra:
        y = np.asarray(s.get("I_raw", s["I"]), float)
        mx = float(np.nanmax(y)) if len(y) else float("nan")
        status, reason = "ok", ""
        det = str((s.get("meta") or {}).get("detector", "")).upper()
        counting = "PMT" in det or "PHOTOMULT" in det or "SPAD" in det            # photon counting: no converter limit (dead time is checked below)
        if np.isfinite(mx) and mx > 500:
            at = np.abs(y - mx) <= max(1e-4 * mx, 0.5)
            run, best = 0, 0
            for a in at:
                run = run + 1 if a else 0
                best = max(best, run)
            lim = None if counting else next((L for L in ADC_LIMITS if 0.995 * L <= mx <= L * 1.0001), None)
            near = None if counting else next((L for L in ADC_LIMITS if 0.90 * L <= mx < 0.995 * L), None)
            if lim is not None:
                status, reason = "fail", f"the maximum ({mx:,.0f}) sits at the {lim + 1:,}-level converter limit: the peak is clipped"
            elif best >= 3:
                status, reason = "fail", f"{best} points share the maximum ({mx:,.0f}): a flat-topped (saturated) peak"
            elif near is not None:
                status, reason = "warn", f"the maximum ({mx:,.0f}) is above 90 % of the {near + 1:,}-level converter range: close to saturation"
        meta = s.get("meta") or {}
        dwell = meta.get("dwell_s")
        if dwell and np.isfinite(mx) and mx / float(dwell) > rate_limit_cps and "PMT" in str(meta.get("detector", "")).upper():
            status = "warn" if status == "ok" else status
            reason = (reason + "; " if reason else "") + f"peak count rate {mx / float(dwell):.2g} s⁻¹ is above {rate_limit_cps:.0e} s⁻¹ (dead-time losses of photon counting)"
        out.append({"name": s.get("name", ""), "max": mx, "status": status, "reason": reason})
    return out


SETTING_KEYS = [("detector", "detector"), ("excitation", "excitation"), ("excitation_nm", "excitation (nm)"), ("dwell_s", "dwell / integration (s)"),
                ("scans", "scans"), ("slit_ex_nm", "excitation slit (nm)"), ("slit_em_nm", "emission slit (nm)"), ("step_nm", "step (nm)")]


def settings_table(spectra):
    """Rows with the measurement settings that the files tell, and the list of settings that DIFFER inside the series."""
    rows, vals = [], {k: set() for k, _ in SETTING_KEYS}
    for s in spectra:
        m = dict(s.get("meta") or {})
        wl = np.asarray(s["wl"], float)
        if "step_nm" not in m or m.get("step_nm") in (None, ""):
            m["step_nm"] = round(float(np.median(np.diff(wl))), 4) if len(wl) > 2 else None
        r = {"spectrum": s.get("name", ""), "range (nm)": f"{wl.min():.1f}–{wl.max():.1f}" if len(wl) else ""}
        for k, lab in SETTING_KEYS:
            v = m.get(k)
            if v not in (None, ""):
                r[lab] = v
                vals[k].add(round(v, 6) if isinstance(v, float) else v)
        rows.append(r)
    differ = [lab for k, lab in SETTING_KEYS if len(vals[k]) > 1]
    return rows, differ


# ------------------------------------------------------------------------------------------------------------------------------ the ratio
def window_ratio(spectra, w1, w2, baseline="none"):
    """Integration-window ratio with its counting error. spectra: [{'x', 'y', 'var'}] ('var' = variance of every point, the raw counts for photon counting;
    without it the noise is estimated from the point-to-point scatter). Returns I1, I2, D, sD (arrays)."""
    I1, I2, sD = [], [], []
    for s in spectra:
        x = np.asarray(s["x"], float)
        y = np.asarray(s["y"], float)
        if s.get("var") is not None:
            var = np.maximum(np.asarray(s["var"], float), 1.0)
        else:
            d = np.diff(y)
            sg = float(np.median(np.abs(d - np.median(d))) / 0.6745 / math.sqrt(2)) if len(d) > 3 else 1.0
            var = np.full(len(y), max(sg, 1e-12) ** 2)
        res = []
        for lo, hi in (w1, w2):
            m = (x >= min(lo, hi)) & (x <= max(lo, hi))
            if m.sum() < 3:
                res.append((0.0, float("inf")))
                continue
            xs, v, vv = x[m], y[m], var[m]
            base = np.interp(xs, [xs[0], xs[-1]], [v[0], v[-1]]) if baseline == "linear" else 0.0
            dx = np.gradient(xs)
            res.append((float(np.trapezoid(v - base, xs)), float(np.sum(vv * dx ** 2))))
        (a, va), (b, vb) = res
        I1.append(a)
        I2.append(b)
        sD.append(abs(a / b) * math.sqrt(va / a ** 2 + vb / b ** 2) if a > 0 and b > 0 else float("nan"))
    I1, I2 = np.array(I1), np.array(I2)
    with np.errstate(divide="ignore", invalid="ignore"):
        D = I1 / I2
    return I1, I2, D, np.array(sD)


def information(T, D, sD):
    """Does the ratio carry temperature information? Weighted straight line of ln D against 1/T: the slope over its error (t), the total change of the ratio
    over the range in units of its typical error, and the rank correlation (monotonic)."""
    from scipy.stats import spearmanr
    T, D, sD = (np.asarray(v, float) for v in (T, D, sD))
    ok = np.isfinite(T) & np.isfinite(D) & (D > 0) & np.isfinite(sD) & (sD > 0)
    T, D, sD = T[ok], D[ok], sD[ok]
    if len(T) < 3:
        return {"t": float("nan"), "change_sigma": float("nan"), "rho": float("nan"), "n": int(len(T))}
    x, y = 1.0 / T, np.log(D)
    w = (D / sD) ** 2
    W, Wx, Wy = w.sum(), (w * x).sum(), (w * y).sum()
    Wxx, Wxy = (w * x * x).sum(), (w * x * y).sum()
    den = W * Wxx - Wx ** 2
    b = (W * Wxy - Wx * Wy) / den
    a = (Wy - b * Wx) / W
    chi2 = float(np.sum(w * (y - a - b * x) ** 2) / max(len(x) - 2, 1))
    se = math.sqrt(W / den * max(chi2, 1.0))
    total = abs(float(b) * (x.max() - x.min()))
    rel = float(np.median(sD / D))
    rho = float(spearmanr(T, D).correlation) if np.ptp(D) > 0 else 0.0
    return {"t": float(abs(b) / se) if se > 0 else float("inf"), "change_sigma": total / max(rel * math.sqrt(2), 1e-12), "rho": rho, "n": int(len(T)),
            "change_pct": 100 * (math.exp(total) - 1)}


def runs_test(res):
    """Wald-Wolfowitz runs of the residual signs: z < -2 means too few runs (a systematic trend: wrong model, background or heating)."""
    s = np.sign(np.asarray(res, float))
    s = s[s != 0]
    n1, n2 = int(np.sum(s > 0)), int(np.sum(s < 0))
    if n1 == 0 or n2 == 0 or len(s) < 6:
        return float("nan"), 0
    runs = 1 + int(np.sum(s[1:] != s[:-1]))
    mu = 2 * n1 * n2 / (n1 + n2) + 1
    var = 2 * n1 * n2 * (2 * n1 * n2 - n1 - n2) / ((n1 + n2) ** 2 * (n1 + n2 - 1))
    return (runs - mu) / math.sqrt(var) if var > 0 else float("nan"), runs


def spectroscopic_gap_cm1(c1_nm, c2_nm):
    """Energy difference of two emission bands (cm-1) from their centres (or barycentres) - for two thermally coupled levels emitting to the same
    final level this is close to the energy gap of the Boltzmann law."""
    return abs(1e7 / float(c1_nm) - 1e7 / float(c2_nm))


def boltzmann_checks(T, D, sD, bf, gap_spec_cm1=None, n_flag_bad=0, n_total=None):
    """The physics and statistics checks of a Boltzmann calibration. bf = result of pl.fit_boltzmann. Returns a list of check()."""
    out = []
    T, D, sD = (np.asarray(v, float) for v in (T, D, sD))
    n = len(T)
    if n < 5:
        out.append(check("Temperatures", "warn" if n >= 3 else "fail", f"{n} temperatures. A Boltzmann fit through fewer than 5 points cannot show whether the law holds.",
                         "General practice: 6 or more temperatures over at least 30–50 K."))
    else:
        out.append(check("Temperatures", "ok", f"{n} temperatures over {np.ptp(T):.0f} K."))
    inf = information(T, D, sD)
    if inf["t"] == inf["t"]:
        if inf["t"] < 3:
            out.append(check("Thermometric information", "fail", f"The ratio changes by only {inf['change_pct']:.1f} % over the range – {inf['t']:.1f} times its error. "
                             "There is no reliable temperature dependence: the energy gap and sensitivity are not determined.",
                             "Use a wider temperature range, other bands (thermally coupled levels), or reduce the ratio error."))
        elif inf["t"] < 10:
            out.append(check("Thermometric information", "warn", f"The ratio changes by {inf['change_pct']:.1f} % over the range ({inf['t']:.0f} times its error): weakly determined."))
        else:
            out.append(check("Thermometric information", "ok", f"The ratio changes by {inf['change_pct']:.1f} % ({inf['t']:.0f} times its error)."))
        if abs(inf["rho"]) < 0.8 and n >= 4:
            out.append(check("Monotonic", "warn", f"The ratio is not monotonic in temperature (rank correlation {inf['rho']:+.2f}). A thermometer must be monotonic: look for drift, "
                             "laser heating, a changing background, or a wrong temperature assignment."))
    if n_total:
        sh = n_flag_bad / n_total
        if sh > 0.3:
            out.append(check("Band fits", "fail" if sh > 0.6 else "warn", f"{n_flag_bad} of {n_total} spectra are not described well by the band model (wavy residual).",
                             "Crowded bands: use 'Manifold fit (crowded / Stark-split bands)' and choose the two ranges, or narrow the fit range."))
        else:
            out.append(check("Band fits", "ok", f"{n_total - n_flag_bad} of {n_total} spectra fitted well."))
    if bf.get("status") == "ok":
        res = (D - bf["f"](T)) / np.where(sD > 0, sD, np.nan)
        chi2 = float(np.nanmean(res ** 2) * n / max(n - 2, 1))
        z, _ = runs_test(res)
        if chi2 > 4:
            out.append(check("Boltzmann residuals", "warn", f"The points scatter {math.sqrt(chi2):.1f}× more around the Boltzmann curve than their error bars say (reduced χ² = {chi2:.1f}).",
                             "The ratio errors miss something (drift, temperature errors), or the law does not describe the data. The energy-gap error was scaled up accordingly."))
        elif chi2 > 2:
            out.append(check("Boltzmann residuals", "info", f"Reduced χ² = {chi2:.1f}: the points scatter somewhat more than their error bars; the energy-gap error was scaled up (Birge ratio)."))
        else:
            out.append(check("Boltzmann residuals", "ok", f"Reduced χ² = {chi2:.2f}: the scatter agrees with the error bars."))
        if z == z and z < -2:
            out.append(check("Systematic trend", "warn", "The residuals of the Boltzmann fit form long runs of the same sign: a curvature that the law does not have "
                             "(non-thermalised levels, heating, a background that changes with temperature)."))
        de = bf["dE_cm1"]
        if bf.get("_swapped"):
            out.append(check("Sign of the gap", "info", "The ratio I₁/I₂ fell with temperature: band 2 is the upper (thermally populated) level, so Δ = I₂/I₁ is used."))
        de = abs(de)
        if gap_spec_cm1:
            r = de / gap_spec_cm1
            if 0.7 <= r <= 1.3:
                out.append(check("Gap vs spectrum", "ok", f"Fitted ΔE = {de:.0f} cm⁻¹ agrees with the spacing of the two bands ({gap_spec_cm1:.0f} cm⁻¹, ratio {r:.2f})."))
            else:
                out.append(check("Gap vs spectrum", "warn", f"Fitted ΔE = {de:.0f} cm⁻¹ differs from the spacing of the two bands ({gap_spec_cm1:.0f} cm⁻¹, ratio {r:.2f}). "
                                 "For two thermally coupled levels they should be close.",
                                 "Possible reasons: the levels are not in thermal equilibrium (non-radiative or energy-transfer channels), other bands contribute, the temperature is not the sample "
                                 "temperature, or the background changes with temperature."))
        if de < 200:
            out.append(check("Coupling range", "info", f"ΔE = {de:.0f} cm⁻¹ is below 200 cm⁻¹ (Stark-level thermometry): a low sensitivity is expected."))
        elif de > 2000:
            out.append(check("Coupling range", "warn", f"ΔE = {de:.0f} cm⁻¹ is above 2000 cm⁻¹: the upper level is hardly populated and usually not in thermal equilibrium (General practice: 200–2000 cm⁻¹)."))
    return [c for c in out if c]


# ------------------------------------------------------------------------------------------------------------------------------ uncertainty
def _linfit_batch(Tm, Dm, w):
    """Weighted straight line ln D = ln a - b / T for many draws at once. Tm, Dm: (n_draw, n_pt); w: (n_pt,) or (n_draw, n_pt). Returns (b, ln a)."""
    x = 1.0 / Tm
    y = np.log(np.clip(Dm, 1e-300, None))
    w = np.broadcast_to(w, x.shape)
    W, Wx, Wy = w.sum(1), (w * x).sum(1), (w * y).sum(1)
    Wxx, Wxy = (w * x * x).sum(1), (w * x * y).sum(1)
    den = W * Wxx - Wx ** 2
    slope = (W * Wxy - Wx * Wy) / den
    icpt = (Wy - slope * Wx) / W
    return -slope, icpt


def predict_T(D, alpha, dE_eV, C=0.0):
    """Temperature from a measured ratio with a Boltzmann calibration: T = dE / (kB ln(alpha / (D - C)))."""
    D = np.asarray(D, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return dE_eV / (KB_EV * np.log(alpha / (D - C)))


def monte_carlo(T, D, sD, sT_rand=0.0, sT_sys=0.0, n=3000, seed=0, offset=False):
    """Monte Carlo of the Boltzmann calibration: every draw adds Gaussian noise to the ratios (sD), to every temperature (random sensor error sT_rand) and one common
    temperature offset (sensor calibration sT_sys) and refits. Returns {dE_eV (draws), alpha (draws), T_pred (draws x points)} - the T predicted for the measured ratios."""
    rng = np.random.default_rng(seed)
    T, D, sD = (np.asarray(v, float) for v in (T, D, sD))
    k = len(T)
    Dm = D[None, :] + rng.standard_normal((n, k)) * sD[None, :]
    Tm = T[None, :] + rng.standard_normal((n, k)) * sT_rand + rng.standard_normal((n, 1)) * sT_sys
    good = np.all(Dm > 0, axis=1)
    Dm, Tm = Dm[good], Tm[good]
    if not offset:
        w = (D / np.where(sD > 0, sD, np.inf)) ** 2
        b, la = _linfit_batch(Tm, Dm, w)
        dE = b * KB_EV
        al = np.exp(la)
        C = np.zeros_like(dE)
    else:
        from scipy.optimize import curve_fit
        mdl = lambda t, a, de, c: a * np.exp(-de / (KB_EV * t)) + c
        dE, al, C = [], [], []
        for i in range(min(len(Dm), 400)):
            try:
                p, _ = curve_fit(mdl, Tm[i], Dm[i], p0=[float(np.max(D)) * 2, 0.05, 0.0], maxfev=4000)
                al.append(p[0]); dE.append(p[1]); C.append(p[2])
            except Exception:
                continue
        dE, al, C = np.array(dE), np.array(al), np.array(C)
    Tp = predict_T(D[None, :], al[:, None], dE[:, None], C[:, None])
    return {"dE_eV": dE, "alpha": al, "C": C, "T_pred": Tp, "n_ok": int(len(dE))}


def _robust_sd(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if len(v) < 3:
        return float("nan")
    q = np.percentile(v, [15.87, 84.13])
    return float((q[1] - q[0]) / 2)


def budget(T, D, s_stat, bf, sT_rand=0.0, sT_sys=0.0, alt_sets=(), runs_sd_cm1=None, n_runs=0, offset=False, resp_rel=0.0, seed=0):
    """The uncertainty budget of the energy gap (cm-1) and of alpha. Components:
      statistical (ratio errors, Monte Carlo), temperature sensor random, temperature sensor calibration, every alternative processing (refit of the Boltzmann law
      on each alternative ratio set: fit range, peak model, background choice), between-run reproducibility, and the spectral-response correction (alpha only:
      a wavelength-independent error of the response multiplies every ratio by the same factor and so cancels in the gap and in Sr).
    Returns {rows: [{component, u_dE_cm1, u_alpha_rel_pct, note}], u_dE_cm1 (total), u_alpha_rel_pct, T_cal (calibration uncertainty per point, K), T_res (resolution, K)}."""
    from mfs import pl as PL
    T, D, s_stat = (np.asarray(v, float) for v in (T, D, s_stat))
    rows = []

    def mc(sD, a, b):
        r = monte_carlo(T, D, sD, a, b, 2500, seed, offset)
        return _robust_sd(r["dE_eV"]) * EV_TO_CM, _robust_sd(np.log(np.clip(r["alpha"], 1e-300, None))) * 100, r
    u_stat, a_stat, r_stat = mc(s_stat, 0.0, 0.0)
    rows.append({"component": "statistical (counting noise of the ratios)", "u_dE_cm1": u_stat, "u_alpha_rel_pct": a_stat, "note": "Monte Carlo over the ratio errors"})
    if sT_rand > 0:
        u, a, _ = mc(s_stat * 1e-9, sT_rand, 0.0)
        rows.append({"component": f"temperature sensor, random (± {sT_rand:g} K per point)", "u_dE_cm1": u, "u_alpha_rel_pct": a, "note": "Monte Carlo"})
    if sT_sys > 0:
        u, a, _ = mc(s_stat * 1e-9, 0.0, sT_sys)
        rows.append({"component": f"temperature sensor, calibration / sample offset (± {sT_sys:g} K, common)", "u_dE_cm1": u, "u_alpha_rel_pct": a, "note": "Monte Carlo"})
    base_de, base_la = bf["dE_cm1"], math.log(max(bf["alpha"], 1e-300))
    for lab, Dalt in alt_sets:
        Dalt = np.asarray(Dalt, float)
        if Dalt.shape != D.shape or not np.all(np.isfinite(Dalt)) or np.any(Dalt <= 0):
            continue
        b2 = PL.fit_boltzmann(T, Dalt, offset=offset)
        if b2.get("status") != "ok":
            continue
        rows.append({"component": f"processing: {lab}", "u_dE_cm1": abs(abs(b2["dE_cm1"]) - abs(base_de)), "u_alpha_rel_pct": 100 * abs(math.log(max(b2["alpha"], 1e-300)) - base_la),
                     "note": "refit of the Boltzmann law on the alternative ratios"})
    if runs_sd_cm1 is not None and n_runs >= 2:
        rows.append({"component": f"reproducibility between {n_runs} runs", "u_dE_cm1": runs_sd_cm1 / math.sqrt(n_runs), "u_alpha_rel_pct": float("nan"),
                     "note": f"standard deviation of ΔE between runs {runs_sd_cm1:.0f} cm⁻¹, divided by √{n_runs}"})
    if resp_rel > 0:
        rows.append({"component": "spectral-response correction", "u_dE_cm1": 0.0, "u_alpha_rel_pct": 100 * resp_rel,
                     "note": "a constant response error cancels in ΔE and Sr; it shifts α (and every absolute ratio)"})
    tot = math.sqrt(sum((r["u_dE_cm1"] or 0.0) ** 2 for r in rows if np.isfinite(r["u_dE_cm1"])))
    tot_a = math.sqrt(sum((r["u_alpha_rel_pct"] or 0.0) ** 2 for r in rows if np.isfinite(r["u_alpha_rel_pct"])))
    if bf.get("chi2_red") and bf["chi2_red"] > 1:                    # scatter beyond the error bars (Birge ratio): the statistical part is scaled
        scale = math.sqrt(bf["chi2_red"])
        extra = u_stat * (scale - 1)
        rows.append({"component": "excess scatter around the Boltzmann curve (Birge ratio)", "u_dE_cm1": u_stat * math.sqrt(scale ** 2 - 1), "u_alpha_rel_pct": a_stat * math.sqrt(scale ** 2 - 1),
                     "note": f"reduced χ² = {bf['chi2_red']:.1f}: errors scaled by √χ² = {scale:.2f}"})
        tot = math.sqrt(tot ** 2 + (u_stat * math.sqrt(scale ** 2 - 1)) ** 2)
        tot_a = math.sqrt(tot_a ** 2 + (a_stat * math.sqrt(scale ** 2 - 1)) ** 2)
        _ = extra
    # temperature uncertainties: calibration (accuracy of the curve, same processing assumed) and resolution (one reading)
    r_all = monte_carlo(T, D, s_stat * math.sqrt(max(bf.get("chi2_red") or 1.0, 1.0)), sT_rand, sT_sys, 2500, seed + 1, offset)
    T_cal = np.array([_robust_sd(r_all["T_pred"][:, i]) for i in range(len(T))])
    Sr = 100 * abs(bf["dE_eV"]) / (KB_EV * T ** 2)
    T_res = (s_stat / D) / (Sr / 100)
    return {"rows": rows, "u_dE_cm1": tot, "u_alpha_rel_pct": tot_a, "T_cal": T_cal, "T_res": T_res, "u_Sr_rel": tot / max(abs(bf["dE_cm1"]), 1e-12)}


def between_runs(T, D, run_ids, offset=False):
    """Energy gap of every independent run (a file / series with at least 3 temperatures). Returns {runs: [{run, n, dE_cm1}], sd_cm1, n_runs}."""
    from mfs import pl as PL
    T, D = np.asarray(T, float), np.asarray(D, float)
    run_ids = np.asarray(run_ids)
    out = []
    for r in dict.fromkeys(run_ids.tolist()):
        m = run_ids == r
        if len(set(np.round(T[m], 1))) < 3:
            continue
        b = PL.fit_boltzmann(T[m], D[m], offset=offset)
        if b.get("status") == "ok":
            out.append({"run": str(r), "n": int(m.sum()), "dE_cm1": float(b["dE_cm1"]), "alpha": float(b["alpha"])})
    sd = float(np.std([o["dE_cm1"] for o in out], ddof=1)) if len(out) >= 2 else None
    return {"runs": out, "sd_cm1": sd, "n_runs": len(out)}


# ------------------------------------------------------------------------------------------------------------------------------ experiment
def hysteresis(T, D, order, dE_eV, sD=None):
    """Heating against cooling: the measurement order gives a heating branch (temperature rising) and a cooling branch (falling). Where both cover the same
    temperatures, the ratio of the cooling branch is compared with the heating branch (interpolated in ln D against 1/T) and turned into an apparent
    temperature difference dT = d(ln D) kB T^2 / dE. Returns {status, max_dT, rows} or {status: 'one direction only'}."""
    T, D = np.asarray(T, float), np.asarray(D, float)
    o = np.argsort(np.asarray(order, float))
    T, D = T[o], D[o]
    sD = None if sD is None else np.asarray(sD, float)[o]
    if len(T) < 4:
        return {"status": "too few points"}
    step = np.sign(np.diff(T))
    dirn = np.concatenate([[step[0]], step])
    up = dirn > 0
    dn = dirn < 0
    if up.sum() < 2 or dn.sum() < 2:
        return {"status": "one direction only"}
    Th, Dh = T[up], D[up]
    oh = np.argsort(Th)
    Th, Dh = Th[oh], Dh[oh]
    rows = []
    for t, d, k in zip(T[dn], D[dn], np.nonzero(dn)[0]):
        if not (Th.min() <= t <= Th.max()):
            continue
        lnh = np.interp(1.0 / t, (1.0 / Th)[::-1], np.log(Dh)[::-1])
        dT = (math.log(d) - lnh) * KB_EV * t ** 2 / abs(dE_eV) if dE_eV else float("nan")
        e = (sD[k] / d) * KB_EV * t ** 2 / abs(dE_eV) if (sD is not None and dE_eV) else float("nan")
        rows.append({"T (K)": float(t), "apparent ΔT cooling − heating (K)": float(dT), "± (K)": float(e)})
    if not rows:
        return {"status": "the branches do not overlap"}
    mx = max(abs(r["apparent ΔT cooling − heating (K)"]) for r in rows)
    return {"status": "ok", "max_dT": mx, "rows": rows}


def laser_heating(P, D, alpha, dE_eV, C=0.0, sD=None):
    """Temperature read from the ratio at several excitation powers (same set temperature): linear in power -> the heating per unit power and the temperature at
    zero power. Returns {status, slope_K_per_unit, slope_err, T0, rows}."""
    P, D = np.asarray(P, float), np.asarray(D, float)
    ok = np.isfinite(P) & np.isfinite(D) & (D > C)
    P, D = P[ok], D[ok]
    if len(P) < 3:
        return {"status": "needs at least 3 powers"}
    Tm = predict_T(D, alpha, dE_eV, C)
    if not np.all(np.isfinite(Tm)):
        return {"status": "a ratio lies outside the calibration (no temperature)"}
    A = np.vstack([P, np.ones_like(P)]).T
    coef, *_ = np.linalg.lstsq(A, Tm, rcond=None)
    res = Tm - A @ coef
    s2 = float(np.sum(res ** 2) / max(len(P) - 2, 1))
    cov = np.linalg.inv(A.T @ A) * max(s2, 1e-12)
    return {"status": "ok", "slope": float(coef[0]), "slope_err": float(math.sqrt(cov[0, 0])), "T0": float(coef[1]), "T0_err": float(math.sqrt(cov[1, 1])),
            "rows": [{"power": float(p), "T from the ratio (K)": float(t)} for p, t in zip(P, Tm)], "heating_at_max": float(coef[0] * P.max())}


# ------------------------------------------------------------------------------------------------------------------------------ known-answer self-check
def simulate_from_fit(fit, T, dE_eV, alpha, counts=True, pedestal=None, seed=0, inverse=False):
    """A series like the measured one with a KNOWN Boltzmann ratio: the fitted line shapes and backgrounds of every spectrum are kept, the two groups are
    rescaled so that I1 / I2 = alpha exp(-dE / kB T) (their sum unchanged) and noise like the measured one is added (Poisson-like for counts, with the
    removed pedestal; otherwise the point-to-point noise of the spectrum)."""
    from mfs import plratio as PR
    rng = np.random.default_rng(seed)
    g1, g2 = fit.get("groups") or ([0], [1])
    out = []
    for i, (row, t) in enumerate(zip(fit["rows"], T)):
        la = np.asarray(row.get("line_areas", row["areas"]), float)
        comps = row["components"]
        s1, s2 = float(np.sum(la[g1])), float(np.sum(la[g2]))
        r = alpha * math.exp(-dE_eV / (KB_EV * t))
        if inverse:                                                # the ratio of the fit is I1/I2 but the thermometer uses I2/I1
            r = 1.0 / r
        tot = s1 + s2
        k1 = (tot * r / (1 + r)) / s1 if s1 > 0 else 1.0
        k2 = (tot / (1 + r)) / s2 if s2 > 0 else 1.0
        y = np.array(row["background"], float).copy()
        for k, c in enumerate(comps):
            y += np.asarray(c, float) * (k1 if k in g1 else k2 if k in g2 else 1.0)
        if counts:
            ped = float(pedestal[i]) if pedestal is not None else 0.0
            var = np.maximum(y + ped, 1.0)
            noisy = y + rng.standard_normal(len(y)) * np.sqrt(var)
        else:
            var = np.full(len(y), PR._noise(row["y"]) ** 2)
            noisy = y + rng.standard_normal(len(y)) * np.sqrt(var)
        out.append({"T": float(t), "x": np.asarray(row["x"], float), "y": noisy, "var": var, "truth": r})
    return out


def self_check(fit, T, dE_eV, alpha, counts=True, pedestal=None, n=6, seed=0, progress=None, inverse=False):
    """Known-answer test of the analysis on series simulated from the user's own fit (see simulate_from_fit): the same band model is fitted, the ratio and its error
    taken, the Boltzmann law fitted. Returns {truth_cm1, runs: [{dE_cm1, err_cm1, z}], bias_pct, z_mean, z_sd, coverage}: unbiased analysis -> bias near 0;
    honest error bars -> z spread near 1 and about 95 % of the runs within 2 error bars."""
    from mfs import pl as PL
    from mfs import plratio as PR
    allb = list(fit["bands"]) + list(fit.get("extra", []))
    runs = []
    for k in range(n):
        sim = simulate_from_fit(fit, T, dE_eV, alpha, counts, pedestal, seed + k, inverse)
        f2 = PR.fit_global(sim, allb, fit["window"], fit.get("bg_order", 1), fit.get("drift", "none"), counts, 60, fit.get("groups"))
        D = np.array([r["ratio"] for r in f2["rows"]])
        sD = np.array([r["ratio_err"] for r in f2["rows"]])
        if inverse:
            sD, D = sD / D ** 2, 1.0 / D
        b = PL.fit_boltzmann(T, D, sigma=sD)
        if b.get("status") == "ok":
            z = (b["dE_cm1"] - dE_eV * EV_TO_CM) / b["dE_err_cm1"] if b["dE_err_cm1"] > 0 else float("nan")
            runs.append({"dE_cm1": b["dE_cm1"], "err_cm1": b["dE_err_cm1"], "z": z})
        if progress:
            progress((k + 1) / n)
    truth = dE_eV * EV_TO_CM
    if not runs:
        return {"truth_cm1": truth, "runs": [], "status": "no run could be fitted"}
    de = np.array([r["dE_cm1"] for r in runs])
    z = np.array([r["z"] for r in runs])
    return {"status": "ok", "truth_cm1": truth, "runs": runs, "bias_pct": 100 * (float(np.mean(de)) - truth) / truth if truth else float("nan"),
            "bias_sigma": float(np.mean(z)) * math.sqrt(len(z)) if len(z) else float("nan"), "z_mean": float(np.mean(z)), "z_sd": float(np.std(z, ddof=1)) if len(z) > 1 else float("nan"),
            "coverage": float(np.mean(np.abs(z) <= 2))}


# ------------------------------------------------------------------------------------------------------------------------------ record
def _jsonable(v):
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, np.ndarray):
        return _jsonable(v.tolist())
    if isinstance(v, (np.floating, float)):
        return None if not np.isfinite(v) else float(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (str, int, bool)) or v is None:
        return v
    return str(v)


def analysis_record(version, files, settings, results, checks=None):
    """JSON text with everything needed to reproduce a result: program version, date, input files with their MD5 fingerprints, every setting, the results and the
    quality checks."""
    rec = {"program": f"Reticula {version}", "module": "Photoluminescence", "created": _dt.datetime.now().isoformat(timespec="seconds"),
           "files": files, "settings": settings, "results": results, "quality_checks": checks or []}
    return json.dumps(_jsonable(rec), indent=2, ensure_ascii=False)
