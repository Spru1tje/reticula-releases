"""PXRD analysis: peak fitting (pseudo-Voigt), Scherrer size with instrument-broadening correction (Gaussian,
Lorentzian or Voigt/Langford), crystallinity index against a reference, and peak-shift tracking (with d-spacing and
the lattice-parameter change for cubic phases)."""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import curve_fit

CU_KA1 = 1.5406


def pv(x, a, c, w, eta, b0, b1):
    g = np.exp(-4 * math.log(2) * ((x - c) / w) ** 2)
    l = 1 / (1 + ((x - c) / (w / 2)) ** 2)
    return a * (eta * l + (1 - eta) * g) + b0 + b1 * (x - c)


def fit_peak(x, y, guess, window=0.8):
    """Pseudo-Voigt + linear background in guess ± window/… Returns centre, FWHM, eta, area, height, errors."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = (x > guess - window) & (x < guess + window)
    xs, ys = x[m], y[m]
    if len(xs) < 8:
        raise ValueError("too few points around the peak")
    i = int(np.argmax(ys))
    p0 = [ys.max() - ys.min(), xs[i], 0.15, 0.5, ys.min(), 0.0]
    lb = [0, guess - window, 1e-3, 0, -np.inf, -np.inf]; ub = [np.inf, guess + window, 2 * window, 1, np.inf, np.inf]
    p, cov = curve_fit(pv, xs, ys, p0=p0, bounds=(lb, ub), maxfev=20000)
    e = np.sqrt(np.diag(cov))
    a, c, w, eta = p[:4]
    if abs(c - guess) > 0.9 * window or a < 3 * np.std(ys - pv(xs, *p)):
        raise ValueError(f"no peak near {guess:g}° (fit ran to the window edge / below noise)")
    area = a * w * (eta * math.pi / 2 + (1 - eta) * math.sqrt(math.pi / (4 * math.log(2))))
    return {"center": c, "center_err": e[1], "fwhm": w, "fwhm_err": e[2], "eta": eta, "height": a, "area": area, "params": p,
            "x": xs, "fit": pv(xs, *p)}


def instrument_correction(beta_obs, beta_inst, shape="gaussian"):
    """Sample broadening from observed and instrumental FWHM (°): gaussian β² = βo² − βi², lorentzian β = βo − βi,
    'voigt' (Langford/de Keijser approximation for pseudo-Voigt profiles): β = sqrt((βo − βi)·sqrt(βo² − βi²))."""
    bo, bi = float(beta_obs), float(beta_inst)
    if bo <= bi:
        return float("nan")
    if shape == "lorentzian":
        return bo - bi
    if shape == "voigt":
        return math.sqrt((bo - bi) * math.sqrt(bo ** 2 - bi ** 2))
    return math.sqrt(bo ** 2 - bi ** 2)


def scherrer_size(center_deg, fwhm_deg, inst_fwhm_deg=0.05, K=0.9, wavelength=CU_KA1, shape="gaussian"):
    b = instrument_correction(fwhm_deg, inst_fwhm_deg, shape)
    if not np.isfinite(b):
        return {"size_nm": float("nan"), "beta_sample": b, "note": "observed width ≤ instrumental width – size not resolvable"}
    D = K * wavelength / (math.radians(b) * math.cos(math.radians(center_deg / 2))) / 10
    return {"size_nm": D, "beta_sample": b, "note": ""}


def instrument_fwhm_from_standard(x, y, peaks, window=0.6):
    """Caglioti-free shortcut: fit the given peaks of a line-profile standard (LaB6 SRM 660, Si) and return a linear
    FWHM(2θ) interpolator."""
    pos, fw = [], []
    for g in peaks:
        try:
            r = fit_peak(x, y, g, window); pos.append(r["center"]); fw.append(r["fwhm"])
        except Exception:
            continue
    pos, fw = np.asarray(pos), np.asarray(fw)
    if len(pos) == 0:
        raise ValueError("no standard peak could be fitted")
    if len(pos) == 1:
        return lambda tt: float(fw[0]), pos, fw
    c = np.polyfit(pos, fw, 1)
    return (lambda tt: float(np.polyval(c, tt))), pos, fw


def crystallinity_index(x, y, x_ref, y_ref, windows, background="linear"):
    """CI (%) = Σ integrated (background-subtracted) intensity of the chosen reflections in the sample / same in the
    reference × 100. Both patterns should be measured and normalised the same way (same mass/holder or normalised to an
    internal standard); windows = [(lo, hi), …] in °2θ."""
    def area(xx, yy, lo, hi):
        m = (xx >= lo) & (xx <= hi)
        xs, ys = xx[m], yy[m]
        if len(xs) < 3:
            return 0.0
        base = np.interp(xs, [xs[0], xs[-1]], [ys[:3].mean(), ys[-3:].mean()]) if background == "linear" else 0 * xs
        return float(np.trapezoid(np.clip(ys - base, 0, None), xs))
    A_s = sum(area(np.asarray(x), np.asarray(y), lo, hi) for lo, hi in windows)
    A_r = sum(area(np.asarray(x_ref), np.asarray(y_ref), lo, hi) for lo, hi in windows)
    return {"CI_pct": 100 * A_s / A_r if A_r else float("nan"), "A_sample": A_s, "A_ref": A_r,
            "per_window": [(lo, hi, area(np.asarray(x), np.asarray(y), lo, hi), area(np.asarray(x_ref), np.asarray(y_ref), lo, hi)) for lo, hi in windows]}


def d_spacing(two_theta, wavelength=CU_KA1):
    return wavelength / (2 * np.sin(np.radians(np.asarray(two_theta, float) / 2)))


def track_shifts(traces, peaks, window=0.6, ref_index=0, hkl=None, wavelength=CU_KA1):
    """Fit every peak in every trace; shift Δ2θ vs the reference trace, d-spacing and (cubic, if hkl given) a = d·√(h²+k²+l²)."""
    rows = []
    ref = {}
    for i, (nm, x, y) in enumerate(traces):
        for j, g in enumerate(peaks):
            try:
                r = fit_peak(x, y, g, window)
                c, ce = r["center"], r["center_err"]
            except Exception:
                c, ce, r = float("nan"), float("nan"), {"fwhm": float("nan")}
            d = float(d_spacing(c, wavelength)) if np.isfinite(c) else float("nan")
            row = {"sample": nm, "peak (°2θ, guess)": g, "2θ (°)": c, "± (°)": ce, "FWHM (°)": r["fwhm"], "d (Å)": d}
            if hkl and j < len(hkl) and hkl[j]:
                h = hkl[j]; row["hkl"] = "".join(str(v) for v in h); row["a cubic (Å)"] = d * math.sqrt(sum(v * v for v in h))
            rows.append(row)
            if i == ref_index:
                ref[j] = c
    for r_ in rows:
        j = peaks.index(r_["peak (°2θ, guess)"])
        r_["Δ2θ vs ref (°)"] = r_["2θ (°)"] - ref.get(j, float("nan"))
    return rows
