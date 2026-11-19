"""v2.15.1: photoluminescence extras - spectral-response correction, temperature uncertainty from REPEAT measurements, power dependence
(slope n of log I vs log P) and CIE 1931 colour coordinates (x, y, CCT) of emission spectra."""
from __future__ import annotations

import io
import math
import re

import numpy as np


# ------------------------------------------------------------------------------------------------- spectral correction
def read_correction(src) -> tuple[np.ndarray, np.ndarray]:
    """Spectral-response correction file: two columns (wavelength nm, correction factor), any separator, '#' comments. Returns (wl, factor)."""
    raw = src.read() if hasattr(src, "read") else src
    txt = raw.decode("utf-8", "replace") if isinstance(raw, (bytes, bytearray)) else str(raw)
    rows = []
    for ln in txt.splitlines():
        ln = ln.strip()
        if not ln or ln[0] in "#%;":
            continue
        ln = ln.replace(";", " ").replace("	", " ")
        parts = [q.strip() for q in ln.split(",")] if "," in ln and ln.count(",") == 1 and " " not in ln else ln.split()
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except (ValueError, IndexError):
            continue
    if len(rows) < 3:
        raise ValueError("correction file: need at least 3 lines 'wavelength factor'")
    a = np.array(sorted(rows))
    if np.any(a[:, 1] <= 0):
        raise ValueError("correction factors must be positive")
    return a[:, 0], a[:, 1]


def apply_correction(wl, y, cal_wl, cal_f):
    """Multiplies the spectrum by the correction factor at each wavelength (linear interpolation; constant outside the calibrated range)."""
    f = np.interp(np.asarray(wl, float), cal_wl, cal_f)
    return np.asarray(y, float) * f


# ------------------------------------------------------------------------------------------------- repeats
def repeat_uncertainty(T, D, sr_pct, tol=0.6):
    """Spectra measured several times at the same temperature (within `tol` K): per temperature the mean ratio, its standard deviation and the
    temperature uncertainty sigma_T = (std(D) / mean(D)) / Sr that the scatter between the repeats implies. `sr_pct` = callable T -> Sr (% per K).
    Returns a list of dicts for the temperatures with >= 2 repeats (empty when there are none)."""
    T = np.asarray(T, float)
    D = np.asarray(D, float)
    order = np.argsort(T)
    groups, cur = [], [order[0]] if len(T) else []
    for i in order[1:]:
        if abs(T[i] - T[cur[0]]) <= tol:
            cur.append(i)
        else:
            groups.append(cur)
            cur = [i]
    if cur:
        groups.append(cur)
    out = []
    for g in groups:
        if len(g) < 2:
            continue
        d = D[g]
        sd = float(np.std(d, ddof=1))
        m = float(np.mean(d))
        tt = float(np.mean(T[g]))
        sr = float(sr_pct(tt)) / 100.0
        out.append({"T": tt, "n": len(g), "mean": m, "std": sd, "rel_pct": 100 * sd / m if m else float("nan"),
                    "sigma_T": (sd / m) / sr if m and sr else float("nan")})
    return out


# ------------------------------------------------------------------------------------------------- power dependence
_UNITS = {"nw": 1e-9, "uw": 1e-6, "µw": 1e-6, "μw": 1e-6, "mw": 1e-3, "w": 1.0, "kw": 1e3}


def parse_power(text):
    """Excitation power / power density from a name such as '250mW', 'P=1.5 W', '12_mW_cm2', '0.8W/cm2' -> number in the unit written
    (converted to W if the unit is a plain power; kept as is for densities). None when nothing is found."""
    m = re.search(r"(\d+(?:[.,]\d+)?)\s*[_ ]?\s*(nW|uW|µW|μW|mW|W|kW)(?![A-Za-z])", str(text), re.I)
    if not m:
        return None
    return float(m.group(1).replace(",", ".")) * _UNITS[m.group(2).lower()]


def fit_power(P, I):
    """log10 I = n log10 P + c (n = number of photons for upconversion; ~1 for a linear process). Returns dict(n, n_err, c, R2, f)."""
    P = np.asarray(P, float)
    I = np.asarray(I, float)
    ok = np.isfinite(P) & np.isfinite(I) & (P > 0) & (I > 0)
    P, I = P[ok], I[ok]
    if len(P) < 3:
        return {"status": "needs at least 3 powers"}
    x, y = np.log10(P), np.log10(I)
    A = np.vstack([x, np.ones_like(x)]).T
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    fit = A @ coef
    dof = max(len(x) - 2, 1)
    s2 = float(np.sum((y - fit) ** 2) / dof)
    cov = np.linalg.inv(A.T @ A) * s2
    ss = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - float(np.sum((y - fit) ** 2)) / ss if ss > 0 else float("nan")
    out = {"status": "ok", "n": float(coef[0]), "n_err": float(math.sqrt(cov[0, 0])), "c": float(coef[1]), "R2": r2,
           "f": lambda p: 10 ** (coef[1] + coef[0] * np.log10(np.asarray(p, float))), "n_points": int(len(x))}
    # saturation: slope of the upper half of the powers clearly lower than that of the lower half
    if len(x) >= 6:
        o = np.argsort(x)
        h = len(x) // 2
        n_lo = np.polyfit(x[o][:h + 1], y[o][:h + 1], 1)[0]
        n_hi = np.polyfit(x[o][h:], y[o][h:], 1)[0]
        out["slope_low"], out["slope_high"] = float(n_lo), float(n_hi)
        out["saturating"] = bool(n_hi < 0.7 * n_lo and n_lo > 0)
    return out


# ------------------------------------------------------------------------------------------------- CIE 1931
def _g(x, mu, s1, s2):
    s = np.where(x < mu, s1, s2)
    return np.exp(-0.5 * ((x - mu) / s) ** 2)


def cmf(wl):
    """CIE 1931 2° colour-matching functions (x̄, ȳ, z̄): multi-lobe Gaussian fit of Wyman, Sloan & Shirley (2013), accurate to a few 10⁻³."""
    wl = np.asarray(wl, float)
    xb = 1.056 * _g(wl, 599.8, 37.9, 31.0) + 0.362 * _g(wl, 442.0, 16.0, 26.7) - 0.065 * _g(wl, 501.1, 20.4, 26.2)
    yb = 0.821 * _g(wl, 568.8, 46.9, 40.5) + 0.286 * _g(wl, 530.9, 16.3, 31.1)
    zb = 1.217 * _g(wl, 437.0, 11.8, 36.0) + 0.681 * _g(wl, 459.0, 26.0, 13.8)
    return xb, yb, zb


def cie_xy(wl, I, photon_counts=True, lo=380.0, hi=780.0):
    """Chromaticity (x, y), CCT (K, McCamy) and the uncorrected tristimulus values of an emission spectrum. `photon_counts`: the intensity is
    photons per wavelength interval → converted to energy (× 1/λ) first. Only the 380–780 nm part of the spectrum is used."""
    wl = np.asarray(wl, float)
    I = np.asarray(I, float)
    m = (wl >= lo) & (wl <= hi)
    w, v = wl[m], I[m]
    if len(w) < 5 or not np.any(v > 0):
        return {"x": float("nan"), "y": float("nan"), "CCT": float("nan")}
    e = v / w if photon_counts else v
    xb, yb, zb = cmf(w)
    X, Y, Z = (float(np.trapezoid(e * c, w)) for c in (xb, yb, zb))
    s = X + Y + Z
    x, y = X / s, Y / s
    n = (x - 0.3320) / (0.1858 - y) if abs(0.1858 - y) > 1e-9 else float("nan")
    cct = 449.0 * n ** 3 + 3525.0 * n ** 2 + 6823.3 * n + 5520.33
    return {"x": float(x), "y": float(y), "CCT": float(cct) if 1500 < cct < 25000 else float("nan"), "X": X, "Y": Y, "Z": Z}


def spectral_locus(step=5.0):
    wl = np.arange(420.0, 681.0, step)
    xb, yb, zb = cmf(wl)
    s = xb + yb + zb
    return wl, xb / s, yb / s


def draw_cie(ax, points, labels=None, colors=None, locus=True):
    """CIE 1931 chromaticity diagram (spectral locus + purple line) with the given points [(x, y), ...]."""
    if locus:
        wl, lx, ly = spectral_locus(2.0)
        ax.plot(np.append(lx, lx[0]), np.append(ly, ly[0]), "-", color="#444444", lw=0.8)
        for w in (460, 480, 500, 520, 540, 560, 580, 600, 620):
            i = int(np.argmin(abs(wl - w)))
            ax.plot(lx[i], ly[i], "|", color="#444444", ms=4)
            ax.annotate(str(w), (lx[i], ly[i]), xytext=(-8 if w < 520 else 4, 2), textcoords="offset points", fontsize=5, color="#444444")
    for k, (x, y) in enumerate(points):
        ax.plot(x, y, "o", ms=4, color=(colors[k] if colors else "#D55E00"), mec="k", mew=0.3, label=(labels[k] if labels else None))
    ax.set_xlim(0, 0.8)
    ax.set_ylim(0, 0.9)
    ax.set_xlabel("CIE x")
    ax.set_ylabel("CIE y")
    ax.set_aspect("equal")


# ------------------------------------------------------------------------------------------------- v2.18.4: correction from a calibration lamp, coverage
def correction_from_lamp(meas_wl, meas_I, cert_wl, cert_I, smooth_nm=5.0):
    """Spectral-response correction from a calibration lamp: factor = certified spectrum / measured spectrum (both on the measured grid), lightly smoothed,
    normalised to 1 at its median. The lamp must be measured with the same detector, grating, slits and filters as the samples. Returns (wl, factor)."""
    mw = np.asarray(meas_wl, float)
    mI = np.asarray(meas_I, float)
    o = np.argsort(np.asarray(cert_wl, float))
    cw, cI = np.asarray(cert_wl, float)[o], np.asarray(cert_I, float)[o]
    m = (mw >= cw.min()) & (mw <= cw.max()) & (mI > 0.02 * np.nanmax(mI))
    if m.sum() < 5:
        raise ValueError("the lamp measurement and the certified spectrum do not overlap enough (or the measured lamp signal is too weak)")
    w, mi = mw[m], mI[m]
    f = np.interp(w, cw, cI) / mi
    step = float(np.median(np.diff(w))) if len(w) > 2 else 1.0
    k = max(1, int(round(smooth_nm / max(step, 1e-9))))
    if k > 1 and len(f) > 2 * k:
        f = np.convolve(np.pad(f, k, mode="edge"), np.ones(2 * k + 1) / (2 * k + 1), mode="same")[k:-k]
    f = f / np.median(f)
    return w, f


def correction_coverage(data_wl_min, data_wl_max, cal_wl):
    """Fraction of the data range that the correction covers (outside it the factor is held constant - not a real correction)."""
    lo, hi = max(data_wl_min, float(np.min(cal_wl))), min(data_wl_max, float(np.max(cal_wl)))
    span = max(data_wl_max - data_wl_min, 1e-9)
    return max(0.0, hi - lo) / span


# ------------------------------------------------------------------------------------------------- v2.18.4: dominant wavelength and colour purity
WHITE = {"E (equal energy, 0.3333, 0.3333)": (1 / 3, 1 / 3), "D65 (0.3127, 0.3290)": (0.3127, 0.3290), "C (0.3101, 0.3162)": (0.3101, 0.3162)}


def dominant_wavelength(x, y, white=(1 / 3, 1 / 3)):
    """Dominant wavelength (nm) and colour purity (%) of a chromaticity (x, y) relative to a white point: the ray from the white point through the colour meets the
    spectral locus at the dominant wavelength; purity = distance colour-white / distance locus-white. When the ray meets the purple line the complementary
    wavelength is returned with complementary=True (purple colours have no dominant wavelength)."""
    if not (np.isfinite(x) and np.isfinite(y)):
        return {"dominant_nm": float("nan"), "purity_pct": float("nan"), "complementary": False}
    wl = np.arange(380.0, 701.0, 0.5)
    xb, yb, zb = cmf(wl)
    s = xb + yb + zb
    lx, ly = xb / s, yb / s
    wx, wy = white
    dx, dy = x - wx, y - wy
    if abs(dx) + abs(dy) < 1e-9:
        return {"dominant_nm": float("nan"), "purity_pct": 0.0, "complementary": False}

    def hit(dirx, diry):
        """the crossing with the locus at or beyond the colour itself (t >= 1): the blue end of the locus folds back, so an earlier crossing is not the exit"""
        best, far = None, None
        for i in range(len(wl) - 1):
            ax_, ay_, bx_, by_ = lx[i], ly[i], lx[i + 1], ly[i + 1]
            ex, ey = bx_ - ax_, by_ - ay_
            den = dirx * ey - diry * ex
            if abs(den) < 1e-15:
                continue
            t = ((ax_ - wx) * ey - (ay_ - wy) * ex) / den
            u = ((ax_ - wx) * diry - (ay_ - wy) * dirx) / den
            if t > 0 and 0 <= u <= 1:
                cand = (t, wl[i] + u * (wl[i + 1] - wl[i]))
                if t >= 0.999 and (best is None or t < best[0]):
                    best = cand
                if far is None or t > far[0]:
                    far = cand
        return best or far
    h = hit(dx, dy)
    if h is not None:
        t, lam = h
        return {"dominant_nm": float(lam), "purity_pct": float(100 * min(1.0, 1.0 / t)), "complementary": False}
    h2 = hit(-dx, -dy)                          # purple: complementary wavelength; purity relative to the purple line
    px, py = lx[0], ly[0]
    qx, qy = lx[-1], ly[-1]
    ex, ey = qx - px, qy - py
    den = dx * ey - dy * ex
    t = ((px - wx) * ey - (py - wy) * ex) / den if abs(den) > 1e-15 else float("nan")
    return {"dominant_nm": float(h2[1]) if h2 else float("nan"), "purity_pct": float(100 * min(1.0, 1.0 / t)) if t == t and t > 0 else float("nan"), "complementary": True}
