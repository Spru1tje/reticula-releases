"""Data processing and analysis: baselines, smoothing, normalisation, PXRD simulation,
BET (Rouquerol), t-plot, TGA/DTG, NMR integration, kinetics and isotherm fits."""
from __future__ import annotations
import math
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.signal import savgol_filter, find_peaks
from scipy.optimize import curve_fit

# --------------------------------------------------------------------------- generic processing
def als_baseline(y, lam=1e5, p=0.01, niter=10):
    """Asymmetric least squares baseline (Eilers & Boelens 2005)."""
    y = np.asarray(y, float)
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    DD = lam * D.dot(D.transpose())
    w = np.ones(L)
    z = y
    for _ in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        z = spsolve(sparse.csc_matrix(W + DD), w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def rolling_ball_min(y, window=101):
    s = pd.Series(y)
    base = s.rolling(window, center=True, min_periods=1).min().rolling(window, center=True, min_periods=1).mean()
    return base.values


def poly_baseline(x, y, order=2, anchors=None):
    if anchors is None:
        # iterative: fit, keep points below fit
        mask = np.ones_like(y, bool)
        for _ in range(20):
            c = np.polyfit(x[mask], y[mask], order)
            b = np.polyval(c, x)
            mask = y <= b + 1e-12
        return b
    xa = np.array(anchors); ya = np.interp(xa, x, y)
    return np.polyval(np.polyfit(xa, ya, min(order, len(xa) - 1)), x)


def smooth(y, method="Savitzky-Golay", window=11, order=3):
    y = np.asarray(y, float)
    if method == "none" or window < 3:
        return y
    window = int(window) | 1
    if method == "Savitzky-Golay":
        return savgol_filter(y, min(window, len(y) - (1 - len(y) % 2)), min(order, window - 1))
    if method == "moving average":
        return pd.Series(y).rolling(window, center=True, min_periods=1).mean().values
    return y


def normalise(x, y, mode="max", region=None):
    y = np.asarray(y, float)
    if mode == "none":
        return y
    if mode == "max":
        return y / np.nanmax(y)
    if mode == "min-max":
        return (y - np.nanmin(y)) / (np.nanmax(y) - np.nanmin(y))
    if mode == "area":
        return y / np.trapezoid(np.abs(y), x) * (1 if x[-1] > x[0] else -1)
    if mode == "peak in region" and region:
        m = (x >= min(region)) & (x <= max(region))
        return y / np.nanmax(y[m])
    return y


def process(x, y, baseline="none", smooth_method="none", window=11, norm="max", lam=1e5, p=0.01, region=None):
    y = np.asarray(y, float); x = np.asarray(x, float)
    b = np.zeros_like(y)
    if baseline == "ALS":
        b = als_baseline(y, lam, p)
    elif baseline == "rolling minimum":
        b = rolling_ball_min(y, max(11, len(y) // 30))
    elif baseline == "polynomial":
        b = poly_baseline(x, y, 3)
    elif baseline == "linear (ends)":
        b = np.interp(x, [x[0], x[-1]], [np.mean(y[:5]), np.mean(y[-5:])])
    y2 = smooth(y - b, smooth_method, window)
    return normalise(x, y2, norm, region)


def transmittance_to_absorbance(y):
    y = np.asarray(y, float)
    if np.nanmax(y) > 2:   # percent
        y = y / 100
    return -np.log10(np.clip(y, 1e-6, None))


# --------------------------------------------------------------------------- PXRD
def simulate_pxrd(cif_bytes_or_path, two_theta=(3, 50), wavelength="CuKa1", fwhm=0.12, step=0.01,
                  substitute: dict | None = None):
    """Simulate a powder pattern from a CIF with pymatgen (kinematic, Lorentz-polarisation, Debye-Waller = 0),
    broadened with pseudo-Voigt profiles. `substitute` e.g. {'Zr': 'Hf'} for the Hf analogue."""
    import warnings, tempfile, os
    warnings.filterwarnings("ignore")
    from pymatgen.core import Structure
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    if isinstance(cif_bytes_or_path, (bytes, bytearray)):
        with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as t:
            t.write(cif_bytes_or_path); p = t.name
    else:
        p = cif_bytes_or_path
    s = Structure.from_file(p)
    if substitute:
        s.replace_species(substitute)
    pat = XRDCalculator(wavelength=wavelength).get_pattern(s, two_theta_range=two_theta)
    x = np.arange(two_theta[0], two_theta[1], step)
    y = np.zeros_like(x)
    eta = 0.5
    sig = fwhm / (2 * math.sqrt(2 * math.log(2)))
    for pos, inten in zip(pat.x, pat.y):
        g = np.exp(-0.5 * ((x - pos) / sig) ** 2)
        l = 1 / (1 + ((x - pos) / (fwhm / 2)) ** 2)
        y += inten * (eta * l + (1 - eta) * g)
    sticks = pd.DataFrame({"two_theta": pat.x, "I": pat.y, "hkl": [", ".join("".join(str(i) for i in h["hkl"]) for h in hk) for hk in pat.hkls],
                           "d": pat.d_hkls})
    return pd.DataFrame({"x": x, "y": y / y.max()}), sticks, s.composition.reduced_formula


def scherrer(x, y, peak_guess, K=0.9, wavelength=1.5406, instr_fwhm=0.05):
    m = (x > peak_guess - 1) & (x < peak_guess + 1)
    xs, ys = x[m], y[m]
    def pv(xx, a, c, w, b):
        return a * (0.5 / (1 + ((xx - c) / (w / 2)) ** 2) + 0.5 * np.exp(-4 * math.log(2) * ((xx - c) / w) ** 2)) + b
    p, _ = curve_fit(pv, xs, ys, p0=[ys.max() - ys.min(), xs[np.argmax(ys)], 0.2, ys.min()])
    beta = math.radians(math.sqrt(max(p[2] ** 2 - instr_fwhm ** 2, 1e-8)))
    D = K * wavelength / (beta * math.cos(math.radians(p[1] / 2))) / 10
    return {"center": p[1], "fwhm": p[2], "size_nm": D}


# --------------------------------------------------------------------------- N2 sorption
CM3STP_TO_MOL = 1 / 22414.0
N2_CROSS = 0.162e-18  # m2
LIQ_N2 = 0.0015468     # cm3 liquid per cm3 STP


def bet_rouquerol(p, q, min_points=4, pmax=0.35, allow_three=True):
    """BET with Rouquerol criteria; if fewer than `min_points` points lie in the criterion-1 range (typical for
    steep Type-I MOF isotherms measured without low-pressure doses), retry with 3 points and flag it
    (result['few_points'] = True) instead of returning nothing."""
    res = _bet_rouquerol(p, q, min_points=min_points, pmax=pmax)
    if res is None and allow_three and min_points > 3:
        res = _bet_rouquerol(p, q, min_points=3, pmax=pmax)
        if res is not None:
            res["few_points"] = True
    elif res is not None:
        res["few_points"] = False
    return res


def bet_range_spread(res, rel_tol=0.0):
    """How much does the BET area depend on the chosen fitting range? Uses every window the Rouquerol search evaluated.
    Returns {'n_valid', 'min', 'median', 'max', 'spread_pct'} over the windows that satisfy ALL Rouquerol criteria
    (or over all windows with criterion 1 and 3 when none does, flagged 'strict': False); None without windows.
    spread_pct = (max - min) / median × 100 - report it next to the area (BETSI-style sensitivity)."""
    if not res or not res.get("windows"):
        return None
    w = res["windows"]; strict = [x for x in w if x["ok"]]
    use = strict or w
    a = np.array([x["area"] for x in use], float)
    med = float(np.median(a))
    return {"n_valid": len(strict), "n_windows": len(w), "strict": bool(strict), "min": float(a.min()), "median": med, "max": float(a.max()),
            "spread_pct": float(100 * (a.max() - a.min()) / med) if med > 0 else float("nan")}


def isotherm_quality(p, q, branch=None, tol=0.01):
    """Physical-consistency checks on an N2 77 K isotherm (adsorption branch). Returns a list of plain-text warnings.
    - the uptake must not decrease with increasing p/p0 (a drop of more than `tol` x max uptake is flagged);
    - the first point above p/p0 0.01 means the micropore-filling step is not resolved (BET then has very few points)."""
    p = np.asarray(p, float); q = np.asarray(q, float)
    o = np.argsort(p); p, q = p[o], q[o]
    out = []
    if len(p) < 3:
        return out
    run = np.maximum.accumulate(q)
    drop = run - q
    k = int(np.argmax(drop))
    if drop[k] > tol * q.max():
        i0 = int(np.argmax(q[:k + 1]))
        out.append(f"Adsorption uptake DECREASES from {q[i0]:.1f} cm³ g⁻¹ at p/p₀ {p[i0]:.2f} to {q[k]:.1f} cm³ g⁻¹ at p/p₀ {p[k]:.2f} "
                   f"(−{drop[k]:.1f} cm³ g⁻¹, {100 * drop[k] / q[i0]:.0f} %). Physically the amount adsorbed cannot fall with rising pressure: "
                   "this points to a measurement artefact (free-space / He-calibration error, Dewar level or bath-temperature drift, "
                   "a small leak, or too little sample in the tube). The data are plotted as measured – check the instrument's validation report.")
    if p.min() > 0.01:
        out.append(f"First adsorption point at p/p₀ = {p.min():.3f}: the micropore-filling step (p/p₀ < 0.01) is not resolved, so only "
                   "a few points fall in the Rouquerol BET range. Add low-pressure doses (e.g. p/p₀ 10⁻⁶–10⁻²) for a robust BET area.")
    return out


def _bet_rouquerol(p, q, min_points=4, pmax=0.35):
    """BET with the four Rouquerol consistency criteria (Rouquerol, Llewellyn & Rouquerol 2007; ISO 9277:2022).
    Returns dict with area, C, n_m, range, R2, criteria flags and the full linear-plot data.
    v2.13: every fitting window [i..j] is evaluated with running sums per start point (vectorised over j) instead of one
    least-squares call per window - same results (tests compare with the former implementation), ~15x faster."""
    p = np.asarray(p, float); q = np.asarray(q, float)
    o = np.argsort(p); p, q = p[o], q[o]
    m = (p > 0) & (p < pmax)
    p, q = p[m], q[m]
    rq = q * (1 - p)
    # criterion 1: n(1-p) increases -> restrict to p <= p at max of n(1-p)
    imax = int(np.argmax(rq)) if len(rq) else -1
    cand = []
    if imax >= 0:
        with np.errstate(divide="ignore", invalid="ignore"):
            yb_all = p / (q * (1 - p))
        rising = np.diff(rq) > 0                              # rising[k]: n(1-p) increases from point k to k+1
        bad = np.concatenate([[0], np.cumsum(~rising)])       # bad[k] = number of non-increasing steps before point k
    for i in range(0, imax + 1):
        if i + min_points - 1 > imax:
            break
        x = p[i:imax + 1]; y = yb_all[i:imax + 1]
        dx = x - x[0]; dy = y - y[0]                          # shifted for numerical stability
        n = np.arange(1, len(x) + 1, dtype=float)
        Sx, Sy = np.cumsum(dx), np.cumsum(dy)
        Sxx, Sxy, Syy = np.cumsum(dx * dx), np.cumsum(dx * dy), np.cumsum(dy * dy)
        with np.errstate(divide="ignore", invalid="ignore"):
            den = n * Sxx - Sx ** 2
            sl = (n * Sxy - Sx * Sy) / den
            ic = y[0] + (Sy - sl * Sx) / n - sl * x[0]
            ss = n * Syy - Sy ** 2
            R2 = np.where(ss > 0, (n * Sxy - Sx * Sy) ** 2 / (den * ss), 0.0)
        qmin, qmax = np.minimum.accumulate(q[i:imax + 1]), np.maximum.accumulate(q[i:imax + 1])
        for k in range(min_points - 1, len(x)):
            s_, b_ = float(sl[k]), float(ic[k])
            if not (b_ > 0 and s_ > 0):
                continue  # criterion 2: C > 0
            j = i + k
            C = s_ / b_ + 1
            nm = 1 / (s_ + b_)
            c3 = bool(qmin[k] <= nm <= qmax[k])
            p_nm = float(np.interp(nm, q[i:j + 1], p[i:j + 1])) if c3 else float("nan")
            c3 = not np.isnan(p_nm)
            # criterion 4: 1/(sqrt(C)+1) ~ p_nm within 20 % (often 10 %)
            p_th = 1 / (math.sqrt(C) + 1)
            c4 = bool(c3 and abs(p_th - p_nm) / p_nm < 0.2)
            mono = bool(bad[j] - bad[i] == 0)
            cand.append(dict(i=i, j=j, slope=s_, intercept=b_, C=C, n_m=nm, R2=float(R2[k]), c1=mono, c2=True, c3=c3, c4=c4,
                             p_nm=p_nm, p_th=p_th, n=j - i + 1))
    if not cand:
        return None
    valid = [c for c in cand if c["c1"] and c["c3"] and c["c4"]] or [c for c in cand if c["c1"] and c["c3"]] or cand
    # prefer highest R2, then more points
    best = max(valid, key=lambda c: (round(c["R2"], 4), c["n"]))
    area = best["n_m"] * CM3STP_TO_MOL * 6.02214076e23 * N2_CROSS
    for c in cand:                                       # v2.13: area of every window, for the range-sensitivity view
        c["area"] = c["n_m"] * CM3STP_TO_MOL * 6.02214076e23 * N2_CROSS
    best["windows"] = [dict(i=c["i"], j=c["j"], p_lo=float(p[c["i"]]), p_hi=float(p[c["j"]]), n=c["n"], area=float(c["area"]), C=float(c["C"]),
                            R2=float(c["R2"]), ok=bool(c["c1"] and c["c3"] and c["c4"])) for c in cand]
    best.update(area=area, p_lo=p[best["i"]], p_hi=p[best["j"]],
                all_criteria=bool(best["c1"] and best["c3"] and best["c4"]),
                plot_p=p[: imax + 1], plot_y=p[: imax + 1] / (q[: imax + 1] * (1 - p[: imax + 1])),
                rq_p=p, rq=rq)
    return best


def bet_fixed_range(p, q, p_lo=0.05, p_hi=0.30):
    """Classical BET on a FIXED p/p0 window (what instrument software does by default, typically 0.05-0.30).
    Shown only for comparison with the reported value - the app's own result is bet_rouquerol()."""
    p = np.asarray(p, float); q = np.asarray(q, float)
    m = (p >= p_lo * 0.999) & (p <= p_hi * 1.001)
    if m.sum() < 3:
        return None
    pp, qq = p[m], q[m]
    y = pp / (qq * (1 - pp))
    s, b = np.polyfit(pp, y, 1)
    nm = 1 / (s + b)
    C = s / b + 1 if b != 0 else np.inf
    R2 = 1 - ((y - (s * pp + b)) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    area = nm * CM3STP_TO_MOL * 6.02214076e23 * N2_CROSS
    return {"area": area, "C": C, "n_m": nm, "R2": R2, "n": int(m.sum()), "p_lo": float(pp.min()), "p_hi": float(pp.max()),
            "intercept_negative": b <= 0}


def t_plot(p, q, t_range=(3.5, 5.0)):
    """Harkins-Jura t-plot: t (Å) = [13.99 / (0.034 - log10 p)]^0.5. Returns micropore volume (cm3/g) and external area."""
    p = np.asarray(p, float); q = np.asarray(q, float)
    m = (p > 1e-4) & (p < 0.8)
    t = np.sqrt(13.99 / (0.034 - np.log10(p[m])))
    qq = q[m]
    s = (t >= t_range[0]) & (t <= t_range[1])
    if s.sum() < 3:
        return None
    slope, icpt = np.polyfit(t[s], qq[s], 1)
    return {"V_micro": icpt * LIQ_N2, "S_ext": slope * 15.47, "t": t, "q": qq, "slope": slope, "intercept": icpt}


def total_pore_volume(p, q, p0=0.95):
    p = np.asarray(p, float); q = np.asarray(q, float)
    o = np.argsort(p)
    return float(np.interp(p0, p[o], q[o]) * LIQ_N2)


# --------------------------------------------------------------------------- TGA
def tga_analysis(T, w, smooth_window=21):
    """Normalise, compute DTG (-dw/dT, % °C-1), find steps and the dehydroxylated plateau."""
    T = np.asarray(T, float); w = np.asarray(w, float)
    o = np.argsort(T); T, w = T[o], w[o]
    # resample on uniform grid
    Tg = np.linspace(T.min(), T.max(), max(400, min(4000, len(T))))
    wg = np.interp(Tg, T, w)
    ws = smooth(wg, "Savitzky-Golay", smooth_window, 2)
    dtg = -np.gradient(ws, Tg)
    pk, props = find_peaks(dtg, prominence=max(dtg.max() * 0.05, 1e-4))
    W = lambda t: float(np.interp(t, Tg, ws))
    T_plat, info = flat_plateau(Tg, dtg, pk)
    return {"T": Tg, "w": ws, "dtg": dtg, "peaks": Tg[pk], "W200": W(200) if Tg.max() > 200 else None,
            "T_plateau": T_plat, "W_plateau": W(T_plat) if T_plat else None, "W_end": float(ws[-20:].mean()),
            "W_initial": float(ws[0]), "plateau_info": info}


def flat_plateau(Tg, dtg, peak_idx, half_width=25.0, slope_tol=None):
    """v2.17: the FLATTEST stretch of the TGA curve between the solvent / dehydroxylation steps and the first linker-combustion step.

    Window: from the end of the last step below ~300 °C (or 200 °C) up to the start of the first DTG peak above 300 °C (or 480 °C).
    Within it, every temperature is scored by the mean |DTG| over +/- half_width °C (the local slope of the mass curve): the plateau is the
    point with the smallest local slope, and it must be flat in absolute terms (the local slope below slope_tol, default 0.05 %/°C).
    Returns (T_plateau or None, info dict with the window, the local slope and a flatness verdict)."""
    Tg = np.asarray(Tg, float); dtg = np.asarray(dtg, float)
    peaks = np.asarray(peak_idx, int)
    # the linker combustion = the largest DTG peak above 380 °C (or nothing: the window then ends at 480 °C)
    lo, hi = 150.0, 480.0
    if len(peaks):
        T_pk = Tg[peaks]
        big = [i for i, t in zip(peaks, T_pk) if t >= 380]
        if big:
            L = max(big, key=lambda i: dtg[i])
            hi = min(hi, float(Tg[L]) - 2 * half_width)
        before = [Tg[i] for i in peaks if Tg[i] < hi]
        if before:
            lo = max(lo, max(before) + 1.5 * half_width)     # the last step below the window (dehydroxylation / solvent) is over
    if hi - lo < 2 * half_width:
        return None, {"window": (float(lo), float(hi)), "slope_pct_per_C": None, "flat": False,
                      "reason": "the steps leave no flat stretch (window narrower than 2 x the averaging width)"}
    # the plateau must lie INSIDE the window: the edges of a window are the steps themselves, not a plateau
    m = (Tg >= lo + half_width) & (Tg <= hi - half_width)
    if m.sum() < 5:
        return None, {"window": (float(lo), float(hi)), "slope_pct_per_C": None, "flat": False, "reason": "no plateau inside the window"}
    dx = float(np.median(np.diff(Tg)))
    k = max(1, int(round(half_width / dx)))
    absd = np.abs(dtg)
    local = np.convolve(absd, np.ones(2 * k + 1) / (2 * k + 1), mode="same")
    idx = np.nonzero(m)[0]
    j = idx[int(np.argmin(local[idx]))]
    tol = 0.05 if slope_tol is None else slope_tol
    flat = bool(local[j] <= tol)
    return float(Tg[j]), {"window": (float(lo), float(hi)), "slope_pct_per_C": float(local[j]), "flat": flat,
                          "half_width": half_width, "reason": "" if flat else "the mass still changes here (no clear plateau in this window)"}


def renormalise_to_end(w, W_end):
    return 100 * np.asarray(w) / W_end


# --------------------------------------------------------------------------- NMR
def integrate(x, y, lo, hi, baseline=True):
    m = (x >= min(lo, hi)) & (x <= max(lo, hi))
    xs, ys = x[m], y[m]
    if len(xs) < 3:
        return 0.0
    o = np.argsort(xs); xs, ys = xs[o], ys[o]
    if baseline:
        ys = ys - np.interp(xs, [xs[0], xs[-1]], [ys[0], ys[-1]])
    return float(np.trapezoid(ys, xs))


def nmr_concentration(I, nH, I_std, nH_std=9, c_std=1.0):
    """qNMR: [X] = (I_X / nH_X) / (I_std / nH_std) * [std]"""
    return (I / nH) / (I_std / nH_std) * c_std


# --------------------------------------------------------------------------- kinetics
def first_order_fit(t, c, sigma=None, offset=False):
    """C(t) = C0 exp(-k t) (+ C_inf). Returns k, k_err (1 SD from covariance), t_half +- err, C0, R2, model."""
    t = np.asarray(t, float); c = np.asarray(c, float)
    if offset:
        f = lambda tt, c0, k, cinf: (c0 - cinf) * np.exp(-k * tt) + cinf
        p0 = [c.max(), 1 / max(t.max() / 3, 1e-9), c.min() * 0.5]
    else:
        f = lambda tt, c0, k: c0 * np.exp(-k * tt)
        p0 = [c.max(), 1 / max(t.max() / 3, 1e-9)]
    p, cov = curve_fit(f, t, c, p0=p0, sigma=sigma, absolute_sigma=sigma is not None, maxfev=20000)
    perr = np.sqrt(np.diag(cov))
    k, kerr = p[1], perr[1]
    res = c - f(t, *p)
    R2 = 1 - (res ** 2).sum() / ((c - c.mean()) ** 2).sum()
    th = math.log(2) / k
    return {"k": k, "k_err": kerr, "t_half": th, "t_half_err": th * kerr / k, "C0": p[0], "params": p, "R2": R2,
            "model": lambda tt: f(np.asarray(tt, float), *p), "cinf": p[2] if offset else None}


def langmuir_fit(ce, qe):
    """q = q_max K C / (1 + K C)"""
    ce = np.asarray(ce, float); qe = np.asarray(qe, float)
    f = lambda C, qm, K: qm * K * C / (1 + K * C)
    p, cov = curve_fit(f, ce, qe, p0=[qe.max() * 1.2, 1 / max(np.median(ce), 1e-9)], maxfev=20000)
    e = np.sqrt(np.diag(cov))
    res = qe - f(ce, *p)
    R2 = 1 - (res ** 2).sum() / ((qe - qe.mean()) ** 2).sum()
    return {"q_max": p[0], "q_max_err": e[0], "K": p[1], "K_err": e[1], "R2": R2, "model": lambda C: f(np.asarray(C, float), *p)}


def freundlich_fit(ce, qe):
    f = lambda C, KF, n: KF * C ** (1 / n)
    p, cov = curve_fit(f, ce, qe, p0=[qe.max() / max(ce.max(), 1e-9) ** 0.5, 2], maxfev=20000)
    e = np.sqrt(np.diag(cov))
    return {"K_F": p[0], "K_F_err": e[0], "n": p[1], "n_err": e[1], "model": lambda C: f(np.asarray(C, float), *p)}


def eyring(T_K, k_s):
    """ln(k/T) = ln(kB/h) + ΔS‡/R - ΔH‡/(RT)"""
    R, kB, h = 8.314462618, 1.380649e-23, 6.62607015e-34
    X = 1 / np.asarray(T_K, float); Y = np.log(np.asarray(k_s, float) / np.asarray(T_K, float))
    (s, b), cov = np.polyfit(X, Y, 1, cov=True) if len(X) > 2 else (np.polyfit(X, Y, 1), np.zeros((2, 2)))
    dH = -s * R; dS = (b - math.log(kB / h)) * R
    return {"dH_kJ": dH / 1000, "dS_J": dS, "dH_err": math.sqrt(cov[0, 0]) * R / 1000, "dS_err": math.sqrt(cov[1, 1]) * R}


# --------------------------------------------------------------------------- titration
def equivalence_points(V, pH, n=3, smooth_window=15):
    """Equivalence volumes = maxima of the smoothed first derivative dpH/dV (Klet 2016 used Lorentzian fits of the same curve)."""
    V = np.asarray(V, float); pH = np.asarray(pH, float)
    o = np.argsort(V); V, pH = V[o], pH[o]
    d1 = np.gradient(smooth(pH, "Savitzky-Golay", smooth_window, 2), V)
    pk, pr = find_peaks(d1, prominence=np.nanmax(d1) * 0.05, distance=max(3, len(V) // 40))
    if len(pk) == 0:
        return []
    top = pk[np.argsort(pr["prominences"])[::-1][:n]]
    return sorted(V[top].tolist())


def pka_half_volume(V, pH, V_eq_prev, V_eq):
    """pKa = pH at the half-equivalence volume between two equivalence points."""
    return float(np.interp((V_eq_prev + V_eq) / 2, V, pH))
