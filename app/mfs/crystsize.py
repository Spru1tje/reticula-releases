"""v2.15: crystallite size and microstrain from a powder pattern - automatic and with honest uncertainties.

Per peak a pseudo-Voigt fit gives the position and FWHM; the instrumental broadening is removed (Gaussian / Lorentzian / Voigt rule); then
  * Scherrer          D = K lambda / (beta cos theta)                         (one value per peak, mean and spread over the peaks)
  * Williamson-Hall   beta cos(theta) = K lambda / D + 4 eps sin(theta)       (straight line: intercept -> size, slope -> microstrain)
  * Halder-Wagner     (beta / tan theta)^2 = (K lambda / D) (beta / (tan theta sin theta)) + 16 eps^2   (weights the peaks differently)
beta in radians. The errors come from the weighted straight-line fits (weights = the fit error of every FWHM). Everything states when it is
not trustworthy (too few peaks, negative strain, size above ~100 nm = the limit of a laboratory diffractometer)."""
from __future__ import annotations

import math

import numpy as np
from scipy.signal import find_peaks, savgol_filter

from mfs import xrd as X

SIZE_LIMIT_NM = 100.0
TYPICAL_LAB_FWHM = 0.05          # deg 2theta of a well-aligned laboratory diffractometer with Cu radiation (typical, not measured)


def auto_peaks(x, y, lo=None, hi=None, n_max=12, min_prom=0.04):
    """Positions (deg) of the clearest peaks in [lo, hi], strongest first (peak height relative to the range in the window)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    lo = x.min() if lo is None else lo
    hi = x.max() if hi is None else hi
    m = (x >= lo) & (x <= hi)
    xs, ys = x[m], y[m]
    if len(xs) < 15:
        return []
    step = float(np.median(np.diff(xs)))
    win = max(5, int(round(0.2 / step)) | 1)
    ysm = savgol_filter(ys, min(win, len(ys) - (1 - len(ys) % 2)), 2) if len(ys) > win else ys
    rng = float(ysm.max() - ysm.min())
    if rng <= 0:
        return []
    idx, pr = find_peaks(ysm, prominence=min_prom * rng, distance=max(3, int(0.3 / step)))
    order = np.argsort(pr["prominences"])[::-1][:n_max]
    return sorted(float(xs[idx[i]]) for i in order)


def _wls(x, y, sy):
    """Weighted straight line y = a + b x with the errors of a and b and R^2."""
    x, y, sy = (np.asarray(v, float) for v in (x, y, sy))
    w = 1.0 / np.maximum(sy, 1e-12) ** 2
    S, Sx, Sy = w.sum(), (w * x).sum(), (w * y).sum()
    Sxx, Sxy = (w * x * x).sum(), (w * x * y).sum()
    d = S * Sxx - Sx ** 2
    if d <= 0 or len(x) < 3:
        return None
    b = (S * Sxy - Sx * Sy) / d
    a = (Sy - b * Sx) / S
    res = y - a - b * x
    chi2 = float((w * res ** 2).sum())
    dof = len(x) - 2
    scale = max(1.0, chi2 / dof) if dof > 0 else 1.0           # inflate the errors when the line does not fit within the fit errors
    ea, eb = math.sqrt(Sxx / d * scale), math.sqrt(S / d * scale)
    ss_tot = float((w * (y - (w * y).sum() / S) ** 2).sum())
    return {"a": a, "b": b, "ea": ea, "eb": eb, "r2": 1.0 - chi2 / ss_tot if ss_tot > 0 else float("nan")}


def analyse(x, y, wavelength=X.CU_KA1, inst=None, K=0.9, shape="voigt", peaks=None, window=None, lo=None, hi=None, max_peaks=12):
    """Complete crystallite-size analysis of ONE pattern. `inst` = callable 2theta -> instrumental FWHM (deg); default = typical laboratory value.
    Returns a dict: peaks (list of dicts), scherrer, wh, hw (each dict or None), notes (list of str)."""
    inst = inst or (lambda tt: TYPICAL_LAB_FWHM)
    notes: list[str] = []
    guess = list(peaks) if peaks else auto_peaks(x, y, lo, hi, max_peaks)
    rows = []
    for g in guess:
        w = window or 0.8
        try:
            r = X.fit_peak(x, y, g, w)
        except Exception as e:
            rows.append({"guess": g, "ok": False, "note": str(e)[:80]})
            continue
        bi = float(inst(r["center"]))
        bs = X.instrument_correction(r["fwhm"], bi, shape)
        ok = np.isfinite(bs) and bs > 0
        row = {"guess": g, "ok": bool(ok), "two_theta": r["center"], "fwhm_obs": r["fwhm"], "fwhm_err": r["fwhm_err"], "fwhm_inst": bi,
               "beta": bs, "eta": r["eta"], "height": r["height"]}
        if ok:
            th = math.radians(r["center"] / 2)
            row["D_scherrer_nm"] = K * wavelength / (math.radians(bs) * math.cos(th)) / 10
            # error of beta from the error of the observed FWHM (the instrumental value is taken as exact)
            db = r["fwhm_err"] * (r["fwhm"] / bs if shape == "gaussian" else 1.0)
            row["beta_err"] = float(db)
            row["D_err_nm"] = row["D_scherrer_nm"] * db / bs
            row["note"] = ""
        else:
            row["note"] = "observed width ≤ instrumental width: not resolvable"
        rows.append(row)
    good = [r for r in rows if r["ok"]]
    out = {"peaks": rows, "scherrer": None, "wh": None, "hw": None, "notes": notes, "K": K, "wavelength": wavelength, "shape": shape}
    if not good:
        notes.append("No usable peak: the peaks are as narrow as the instrument (size > ~100 nm) or none could be fitted.")
        return out
    D = np.array([r["D_scherrer_nm"] for r in good])
    out["scherrer"] = {"D_mean": float(D.mean()), "D_sd": float(D.std(ddof=1)) if len(D) > 1 else float("nan"), "n": len(D),
                       "D_weighted": float(np.sum(D / np.maximum([r["D_err_nm"] for r in good], 1e-9) ** 2) / np.sum(1 / np.maximum([r["D_err_nm"] for r in good], 1e-9) ** 2))}
    if out["scherrer"]["D_mean"] > SIZE_LIMIT_NM:
        notes.append(f"Mean size {out['scherrer']['D_mean']:.0f} nm is above ~{SIZE_LIMIT_NM:.0f} nm: the peaks are close to the instrumental width, treat the size as a lower limit.")
    if len(good) >= 3:
        th = np.radians([r["two_theta"] for r in good]) / 2
        beta = np.radians([r["beta"] for r in good])
        sb = np.radians([r["beta_err"] for r in good])
        # Williamson-Hall: y = beta cos(theta), x = 4 sin(theta)
        f = _wls(4 * np.sin(th), beta * np.cos(th), sb * np.cos(th))
        if f:
            Dwh = K * wavelength / f["a"] / 10 if f["a"] > 0 else float("nan")
            out["wh"] = {"x": 4 * np.sin(th), "y": beta * np.cos(th), "sy": sb * np.cos(th), "intercept": f["a"], "slope": f["b"], "r2": f["r2"],
                         "D_nm": Dwh, "D_err_nm": Dwh * f["ea"] / f["a"] if f["a"] > 0 else float("nan"), "strain": f["b"], "strain_err": f["eb"]}
            if not f["a"] > 0:
                notes.append("Williamson–Hall: the intercept is ≤ 0 – the width is dominated by strain or the instrument; no size from this fit.")
            if f["b"] < 0:
                notes.append("Williamson–Hall: negative slope (negative strain is not physical) – the peaks are broadened by size only; use the Scherrer value.")
        # Halder-Wagner: y = (beta / tan theta)^2, x = beta / (tan theta sin theta)
        t, s = np.tan(th), np.sin(th)
        yy = (beta / t) ** 2
        xx = beta / (t * s)
        syy = 2 * beta / t ** 2 * sb
        g2 = _wls(xx, yy, syy)
        if g2 and g2["b"] > 0:
            Dhw = K * wavelength / g2["b"] / 10
            eps2 = g2["a"] / 16.0
            eps = math.sqrt(eps2) if eps2 > 0 else 0.0
            out["hw"] = {"x": xx, "y": yy, "intercept": g2["a"], "slope": g2["b"], "r2": g2["r2"], "D_nm": Dhw, "D_err_nm": Dhw * g2["eb"] / g2["b"],
                         "strain": eps, "strain_err": (g2["ea"] / 16.0) / (2 * eps) if eps > 0 else float("nan")}
    else:
        notes.append(f"Only {len(good)} usable peak(s): Williamson–Hall and Halder–Wagner need ≥ 3 peaks; the Scherrer value is given.")
    if out["wh"] and out["scherrer"]:
        a, b = out["scherrer"]["D_mean"], out["wh"]["D_nm"]
        if np.isfinite(b) and b > 0 and (max(a, b) / min(a, b) > 2):
            notes.append("Scherrer and Williamson–Hall sizes differ by more than a factor 2: strain / anisotropic broadening is likely – check the fit lines.")
    return out


def summary_line(res: dict) -> str:
    """One readable sentence with the main result."""
    s, w = res.get("scherrer"), res.get("wh")
    if not s:
        return "No size could be determined."
    t = f"Scherrer: D = {s['D_mean']:.1f}" + (f" ± {s['D_sd']:.1f}" if np.isfinite(s["D_sd"]) else "") + f" nm ({s['n']} peaks)"
    if w and np.isfinite(w["D_nm"]):
        t += f"; Williamson–Hall: D = {w['D_nm']:.1f} ± {w['D_err_nm']:.1f} nm, ε = {w['strain'] * 100:.3f} ± {w['strain_err'] * 100:.3f} %"
    return t + "."


def draw_wh(ax, results: dict, colors=None):
    """Williamson–Hall plot of several samples: points with error bars and the fitted line. results = {name: analyse(...)}."""
    colors = colors or {}
    for k, (nm, r) in enumerate(results.items()):
        w = r.get("wh")
        if not w:
            continue
        c = colors.get(nm) or f"C{k}"
        ax.errorbar(w["x"], w["y"], yerr=w["sy"], fmt="o", ms=4, color=c, capsize=2, lw=0.8, label=nm)
        xx = np.linspace(0, float(np.max(w["x"])) * 1.05, 20)
        ax.plot(xx, w["intercept"] + w["slope"] * xx, "-", color=c, lw=1.0)
    ax.set_xlabel(r"4 sin $\theta$")
    ax.set_ylabel(r"$\beta$ cos $\theta$ (rad)")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
