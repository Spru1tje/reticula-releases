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
