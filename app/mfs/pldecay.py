r"""v2.18.4: luminescence decays (lifetimes) - time-correlated single-photon counting (TCSPC) and multichannel-scaling (MCS) histograms, or any decay curve.

Models (intensity I(t), background B):
  tail fit            I(t) = sum_i a_i exp(-(t - t0) / tau_i) + B          for t after the peak (the instrument response is short compared with tau)
  rise and decay      I(t) = a [exp(-(t - t0)/tau_d) - exp(-(t - t0)/tau_r)] + B   (upconversion, energy transfer: the emitting level is fed with a delay)
  reconvolution       I(t) = IRF(t - s) (*) sum_i a_i exp(-t / tau_i) + B      with a measured instrument response function (IRF) and a shift s
Weights: photon counts -> 1 / max(counts, 1) (Neyman chi-squared; slightly biased below about 10 counts per channel); otherwise the noise of the tail.
The number of exponentials is chosen by the information criterion; a higher number must clearly improve the fit AND give separated lifetimes (ratio > 1.5)
with every amplitude above 2 % - otherwise the simpler model is kept (General practice: report the simplest model that leaves a random residual).
Averages: amplitude-weighted tau_amp = sum a_i tau_i / sum a_i (rates, energy transfer), intensity-weighted tau_int = sum a_i tau_i^2 / sum a_i tau_i (what most
instruments call the average lifetime). Lifetime thermometry: tau(T) = tau0 / (1 + C exp(-Ea / kB T)) (Mott-Seitz), Sr = |d tau / d T| / tau."""
from __future__ import annotations

import math
import re

import numpy as np
from scipy.optimize import least_squares, nnls

UNITS = {"ps": 1e-3, "ns": 1.0, "us": 1e3, "µs": 1e3, "μs": 1e3, "ms": 1e6, "s": 1e9}


def _unit_from(text):
    t = str(text)
    for u in ("ps", "ns", "µs", "μs", "us", "ms"):
        if re.search(r"(?<![A-Za-z])" + re.escape(u) + r"(?![A-Za-z])", t):
            return "µs" if u in ("us", "μs") else u
    if re.search(r"\(\s*s\s*\)|\[\s*s\s*\]|(?<![A-Za-z])sec(?:onds?)?(?![A-Za-z])", t):
        return "s"
    return None


def read_decay(src, name=""):
    """Decay curves from text / CSV / Excel: first column time, then one column per decay (header = names). The time unit is read from the header
    ('Time (ns)', 'us', 'ms'); without one it is 'ns'. Returns [{name, T_text, t, y, unit}]."""
    from mfs import pl as PL
    raw = src if isinstance(src, (bytes, bytearray)) else src.read()
    rows = PL._rows_from_bytes(bytes(raw))
    data = [i for i, r in enumerate(rows) if len(r) >= 2 and PL._is_num(r[0]) and PL._is_num(r[1])]
    if len(data) < 10:
        raise ValueError("No decay found: expected a time column followed by one or more count columns.")
    lens = [len(rows[i]) for i in data]
    ncol = max(set(lens), key=lens.count)
    use = [i for i in data if len(rows[i]) == ncol and all(PL._is_num(c) for c in rows[i][:ncol])]
    arr = np.array([[PL._f(c) for c in rows[i][:ncol]] for i in use], float)
    head = " ".join(" ".join(r) for r in rows[:use[0]])
    header = rows[use[0] - 1] if use[0] > 0 else None
    unit = _unit_from(head) or "ns"
    stem = re.sub(r"\.[A-Za-z0-9]{2,5}$", "", re.split(r"[\\/]", name)[-1]) if name else "decay"
    out = []
    o = np.argsort(arr[:, 0])
    for k in range(1, ncol):
        lab = header[k] if header is not None and len(header) >= ncol and header[k] and not PL._is_num(header[k]) else (stem if ncol == 2 else f"{stem} #{k}")
        out.append({"name": lab, "T_text": lab if lab != stem else stem, "t": arr[o, 0], "y": arr[o, k], "unit": unit})
    return out


def _weights(y, counts):
    if counts:
        return 1.0 / np.maximum(y, 1.0)
    tail = y[int(0.8 * len(y)):]
    d = np.diff(tail) if len(tail) > 4 else np.diff(y)
    s = float(np.median(np.abs(d - np.median(d))) / 0.6745 / math.sqrt(2)) if len(d) > 3 else float(np.std(y))
    return np.full(len(y), 1.0 / max(s, 1e-12) ** 2)


def _basis(t, taus, model, t0, irf=None, shift=0.0, tr=None):
    cols = []
    if model == "rise":
        tt = np.clip(t - t0, 0, None)
        cols.append(np.where(t >= t0, np.exp(-tt / taus[0]) - np.exp(-tt / max(tr, 1e-12)), 0.0))
        return np.stack(cols, 1)
    if model == "reconv":
        dt = float(np.median(np.diff(t)))
        tg = np.arange(len(t)) * dt
        ir = np.interp(t - shift, t, irf, left=0.0, right=0.0)
        ir = ir / max(ir.sum(), 1e-300)
        for tau in taus:
            e = np.exp(-tg / tau)
            c = np.convolve(ir, e)[: len(t)]
            cols.append(c)
        return np.stack(cols, 1)
    tt = t - t0
    for tau in taus:
        cols.append(np.exp(-tt / tau))
    return np.stack(cols, 1)


def _solve_lin(Bm, y, w, bg=True):
    M = np.hstack([Bm, np.ones((len(y), 1))]) if bg else Bm
    sw = np.sqrt(w)
    c, _ = nnls(M * sw[:, None], y * sw, maxiter=5000)
    return c, M @ c


def fit_decay(t, y, n_exp=1, model="tail", irf=None, t_start=None, t_end=None, counts=True, bg=True):
    """One decay fit. model: 'tail' | 'rise' | 'reconv' (needs irf on the same time grid). Returns a result dict (taus, errors, amplitudes, averages, chi2_red, dw,
    fit curve, residuals on the fitted range) or {'status': reason}."""
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    ipk = int(np.argmax(y))
    if model == "tail":
        start = t_start if t_start is not None else t[min(ipk + max(2, len(t) // 200), len(t) - 5)]
        t0 = float(t[ipk])                                   # amplitudes refer to the peak (the end of the excitation), not to the start of the fitted range
    elif model == "rise":
        on = np.nonzero(y > 0.05 * y.max())[0]
        t0 = float(t[on[0]]) if len(on) else float(t[0])
        start = t_start if t_start is not None else t0
    else:
        if irf is None:
            return {"status": "reconvolution needs an instrument response (IRF)"}
        on = np.nonzero(np.asarray(irf, float) > 0.01 * float(np.max(irf)))[0]
        start = t_start if t_start is not None else float(t[on[0]] if len(on) else t[0])
        t0 = float(t[0])
    end = t_end if t_end is not None else t[-1]
    m = (t >= start) & (t <= end)
    if m.sum() < 8 + 2 * n_exp:
        return {"status": "too few points in the fitted range"}
    ts, ys = t[m], y[m]
    w = _weights(ys, counts)
    span = float(ts[-1] - ts[0])
    dt = float(np.median(np.diff(ts)))
    # start lifetimes: log-spaced inside the window
    if model == "rise":
        p0 = [math.log(span / 4), math.log(max(span / 40, dt)), t0 - 2 * dt]
        lb, ub = [math.log(dt / 2)] * 2 + [float(t[0])], [math.log(span * 20)] * 2 + [t0 + dt]
        start = float(t[0]) if t_start is None else start
        m = (t >= start) & (t <= end)
        ts, ys = t[m], y[m]
        w = _weights(ys, counts)
    else:
        p0 = list(np.log(np.geomspace(span / 20, span / 2, n_exp))) if n_exp > 1 else [math.log(span / 5)]
        lb, ub = [math.log(dt / 4)] * n_exp, [math.log(span * 20)] * n_exp
    if model == "reconv":
        p0.append(0.0)
        lb.append(-10 * dt)
        ub.append(10 * dt)
    irf_full = None if irf is None else np.asarray(irf, float)

    def basis(taus, sh, tr):
        if model == "reconv":
            return _basis(t, taus, "reconv", 0.0, irf_full, sh)[m]
        if model == "rise":
            return _basis(ts, taus, "rise", sh, None, 0.0, tr)
        return _basis(ts, taus, "tail", t0)

    def unpack(p):
        if model == "rise":
            return [math.exp(p[0])], p[2], math.exp(p[1])
        if model == "reconv":
            return [math.exp(v) for v in p[:-1]], p[-1], None
        return [math.exp(v) for v in p], 0.0, None

    def resid(p):
        taus, sh, tr = unpack(p)
        Bm = basis(taus, sh, tr)
        c, f = _solve_lin(Bm, ys, w, bg)
        return (ys - f) * np.sqrt(w)
    try:
        sol = least_squares(resid, p0, bounds=(lb, ub), x_scale=1.0, max_nfev=400)
    except Exception as e:
        return {"status": f"fit failed: {e}"}
    taus, sh, tr = unpack(sol.x)
    Bm = basis(taus, sh, tr)
    c, f = _solve_lin(Bm, ys, w, bg)
    amps = c[: Bm.shape[1]]
    B = float(c[-1]) if bg else 0.0
    r = (ys - f)
    chi2 = float(np.sum(w * r * r) / max(len(ys) - len(sol.x) - len(c), 1))
    dw = float(np.sum(np.diff(r) ** 2) / max(np.sum(r ** 2), 1e-300))
    # errors of the nonlinear parameters (lifetimes) from the Jacobian of the reduced problem, scaled by chi2
    try:
        J = sol.jac
        cov = np.linalg.inv(J.T @ J) * max(chi2, 1.0)
        se = np.sqrt(np.clip(np.diag(cov), 0, None))
    except Exception:
        se = np.full(len(sol.x), np.nan)
    if model == "rise":
        tau_list = [taus[0]]
        tau_err = [taus[0] * se[0]]
        rise = tr
        rise_err = tr * se[1]
    else:
        n_t = len(taus)
        tau_list = list(taus)
        tau_err = [tau * s for tau, s in zip(taus, se[:n_t])]
        rise = rise_err = None
    o = np.argsort(tau_list)
    tau_list = [tau_list[i] for i in o]
    tau_err = [tau_err[i] for i in o]
    amps = np.asarray(amps, float)[o] if model != "rise" else np.asarray(amps, float)
    a_sum = float(np.sum(amps))
    frac_a = amps / a_sum if a_sum > 0 else amps * np.nan
    ai_ti = amps * np.array(tau_list)
    frac_i = ai_ti / max(float(np.sum(ai_ti)), 1e-300)
    tau_amp = float(np.sum(ai_ti) / a_sum) if a_sum > 0 else float("nan")
    tau_int = float(np.sum(amps * np.array(tau_list) ** 2) / max(np.sum(ai_ti), 1e-300))
    flag = "good"
    if chi2 > 1.3 and counts or dw < 1.5:
        flag = "check"
    if chi2 > 3 and counts or dw < 1.0:
        flag = "poor"
    return {"status": "ok", "model": model, "n_exp": len(tau_list), "taus": tau_list, "tau_err": tau_err, "amps": amps.tolist(), "frac_amp": frac_a.tolist(),
            "frac_int": frac_i.tolist(), "tau_amp": tau_amp, "tau_int": tau_int, "rise": rise, "rise_err": rise_err, "B": B, "shift": sh if model == "reconv" else None, "onset": sh if model == "rise" else None,
            "chi2_red": chi2, "dw": dw, "flag": flag, "t": ts, "y": ys, "fit": f, "wres": r * np.sqrt(w), "t0": t0,
            "bic": float(len(ys) * math.log(max(np.sum(w * r * r), 1e-300) / len(ys)) + (len(sol.x) + len(c)) * math.log(len(ys)))}


def fit_auto(t, y, model="tail", irf=None, t_start=None, t_end=None, counts=True, max_exp=3):
    """Fits 1 .. max_exp exponentials and keeps the simplest one that is not clearly beaten (BIC lower by more than 10 AND separated lifetimes AND every amplitude
    fraction above 2 %). Returns (best, all fits)."""
    if model == "rise":
        r = fit_decay(t, y, 1, "rise", irf, t_start, t_end, counts)
        return r, [r]
    fits = []
    best = None
    for n in range(1, max_exp + 1):
        r = fit_decay(t, y, n, model, irf, t_start, t_end, counts)
        fits.append(r)
        if r.get("status") != "ok":
            continue
        if best is None:
            best = r
            continue
        taus = sorted(r["taus"])
        sep = all(b / a > 1.5 for a, b in zip(taus, taus[1:]))
        big = all(f > 0.02 for f in r["frac_amp"])
        if r["bic"] < best["bic"] - 10 and sep and big:
            best = r
        else:
            break
    return best or fits[0], fits


def draw_decay(ax, rax, r, unit="ns", log=True, colour="#0072B2"):
    """Decay with the fit (log scale) and the weighted residuals below (on rax, or offset on the same axis when rax is None)."""
    ax.plot(r["t"], r["y"], ".", ms=1.6, color="#555555", label="measured")
    ax.plot(r["t"], r["fit"], "-", lw=1.0, color=colour, label="fit")
    if log:
        ax.set_yscale("log")
        ax.set_ylim(max(np.min(r["y"][r["y"] > 0]) * 0.7 if np.any(r["y"] > 0) else 0.5, 0.5), np.max(r["y"]) * 1.5)
    ax.set_ylabel("Counts" if log else "Intensity")
    ax.legend(frameon=False, fontsize=6)
    if rax is not None:
        rax.plot(r["t"], r["wres"], "-", lw=0.5, color="#555555")
        rax.axhline(0, lw=0.4, color="k")
        rax.set_ylabel("w. res.")
        rax.set_xlabel(f"Time ({unit})")
    else:
        ax.set_xlabel(f"Time ({unit})")


def fit_lifetime_T(T, tau):
    """tau(T) = tau0 / (1 + C exp(-Ea / kB T)) - the Mott-Seitz law has the form of the thermal-quenching law, so its fit is reused. Adds the relative sensitivity
    Sr(T) = 100 % |d tau / d T| / tau."""
    from mfs import pl as PL
    f = PL.fit_quench(T, tau)
    if f.get("status") != "ok":
        return f
    T = np.asarray(T, float)
    h = 0.01
    dt = (f["f"](T + h) - f["f"](T - h)) / (2 * h)
    f["Sr"] = 100 * np.abs(dt) / f["f"](T)
    return f
