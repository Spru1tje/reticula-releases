"""Kinetics toolbox: first-order fits with 95 % confidence intervals, initial rates, Arrhenius and Eyring analysis
from several temperatures, TOF per open site, and pD (pH)–rate profiles (sigmoidal or bell-shaped, giving pKa)."""
from __future__ import annotations
import math
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

R_GAS, KB, H_PLANCK = 8.314462618, 1.380649e-23, 6.62607015e-34


def _tcrit(dof, level=0.95):
    return float(stats.t.ppf(0.5 + level / 2, max(dof, 1)))


def first_order_ci(t, c, offset=False, level=0.95):
    """C(t) = C0 exp(-k t) (+ C_inf), with the parameter SE from the covariance and the CI = t(dof) × SE."""
    t = np.asarray(t, float); c = np.asarray(c, float)
    if offset:
        f = lambda tt, c0, k, ci: (c0 - ci) * np.exp(-k * tt) + ci
        p0 = [c.max(), 1 / max(t.max() / 3, 1e-9), c.min() * 0.5]
    else:
        f = lambda tt, c0, k: c0 * np.exp(-k * tt)
        p0 = [c.max(), 1 / max(t.max() / 3, 1e-9)]
    p, cov = curve_fit(f, t, c, p0=p0, maxfev=20000)
    se = np.sqrt(np.diag(cov)); dof = len(t) - len(p); tc = _tcrit(dof, level)
    res = c - f(t, *p); R2 = 1 - (res ** 2).sum() / ((c - c.mean()) ** 2).sum()
    k, kse = p[1], se[1]
    th = math.log(2) / k
    return {"k_err": kse, "t_half_err": th * kse / k, "cinf": p[2] if offset else None,
            "k": k, "k_se": kse, "k_ci": tc * kse, "k_lo": k - tc * kse, "k_hi": k + tc * kse, "C0": p[0], "R2": R2, "n": len(t),
            "dof": dof, "level": level, "t_half": th, "t_half_ci": th * tc * kse / k, "params": p, "cov": cov,
            "model": lambda tt: f(np.asarray(tt, float), *p), "band": _band(f, p, cov, tc)}


def _band(f, p, cov, tc):
    """Confidence band of the fitted curve (delta method)."""
    def band(tt):
        tt = np.asarray(tt, float); y0 = f(tt, *p); J = []
        for i in range(len(p)):
            dp = np.zeros_like(p); h = 1e-6 * max(abs(p[i]), 1e-12); dp[i] = h
            J.append((f(tt, *(p + dp)) - y0) / h)
        J = np.array(J).T
        s = np.sqrt(np.einsum("ij,jk,ik->i", J, cov, J))
        return y0 - tc * s, y0 + tc * s
    return band


def initial_rate(t, c, frac=0.1, n_min=3, level=0.95):
    """Initial rate from a straight line through the points up to `frac` conversion (at least n_min points):
    r0 = -dC/dt at t -> 0, with the CI of the slope."""
    t = np.asarray(t, float); c = np.asarray(c, float); o = np.argsort(t); t, c = t[o], c[o]
    c0 = c[0]; conv = 1 - c / c0 if c0 else np.zeros_like(c)
    m = conv <= frac
    m[:n_min] = True
    tt, cc = t[m], c[m]
    r = stats.linregress(tt, cc)
    tc = _tcrit(len(tt) - 2, level)
    return {"r0": -r.slope, "r0_se": r.stderr, "r0_ci": tc * r.stderr, "intercept": r.intercept, "n": int(m.sum()),
            "t_max": float(tt.max()), "R2": r.rvalue ** 2, "frac": frac}


def arrhenius(T_C, k, k_err=None, level=0.95):
    """ln k = ln A - Ea/(R T). Weighted by 1/σ(ln k)² if k_err is given. Ea (kJ/mol) and ln A with CI."""
    T = np.asarray(T_C, float) + 273.15; k = np.asarray(k, float)
    X, Y = 1 / T, np.log(k)
    w = None if k_err is None else (np.asarray(k, float) / np.asarray(k_err, float)) ** 2
    (s, b), cov = _wlin(X, Y, w)
    tc = _tcrit(len(X) - 2, level)
    Ea = -s * R_GAS / 1000
    return {"Ea_kJ": Ea, "Ea_ci": tc * math.sqrt(cov[0, 0]) * R_GAS / 1000, "lnA": b, "lnA_ci": tc * math.sqrt(cov[1, 1]),
            "A": math.exp(b), "X": X, "Y": Y, "slope": s, "intercept": b, "R2": _r2(X, Y, s, b), "n": len(X)}


def eyring(T_C, k_s, k_err=None, level=0.95):
    """ln(k/T) = ln(kB/h) + ΔS‡/R − ΔH‡/(R T); k in s⁻¹. ΔH‡ (kJ/mol), ΔS‡ (J/mol/K), ΔG‡(298 K) with CIs."""
    T = np.asarray(T_C, float) + 273.15; k = np.asarray(k_s, float)
    X, Y = 1 / T, np.log(k / T)
    w = None if k_err is None else (k / np.asarray(k_err, float)) ** 2
    (s, b), cov = _wlin(X, Y, w)
    tc = _tcrit(len(X) - 2, level)
    dH = -s * R_GAS; dS = (b - math.log(KB / H_PLANCK)) * R_GAS
    return {"dH_kJ": dH / 1000, "dH_ci": tc * math.sqrt(cov[0, 0]) * R_GAS / 1000, "dS_J": dS, "dS_ci": tc * math.sqrt(cov[1, 1]) * R_GAS,
            "dG298_kJ": (dH - 298.15 * dS) / 1000, "X": X, "Y": Y, "slope": s, "intercept": b, "R2": _r2(X, Y, s, b), "n": len(X)}


def _wlin(X, Y, w=None):
    if len(X) < 2:
        raise ValueError("need at least 2 temperatures")
    if len(X) == 2:
        p = np.polyfit(X, Y, 1); return p, np.zeros((2, 2))
    if w is None:
        return np.polyfit(X, Y, 1, cov=True)
    return np.polyfit(X, Y, 1, w=np.sqrt(w), cov="unscaled") if np.all(np.isfinite(w)) else np.polyfit(X, Y, 1, cov=True)


def _r2(X, Y, s, b):
    r = Y - (s * X + b); return float(1 - (r ** 2).sum() / max(((Y - Y.mean()) ** 2).sum(), 1e-300))


def tof(k_obs_per_s, substrate_mM, volume_mL, catalyst_mg, open_sites_per_node, M_node_g_mol, k_err=None):
    """Turnover frequency per open site (s⁻¹): initial rate (k_obs·C0, mol/s) / open sites (mol).
    open sites (mol) = catalyst mass / M(node formula) × open sites per node."""
    n_sub = substrate_mM * 1e-3 * volume_mL * 1e-3                 # mol substrate
    rate = k_obs_per_s * n_sub                                     # mol/s at t -> 0
    n_os = catalyst_mg * 1e-3 / M_node_g_mol * open_sites_per_node   # mol open sites
    t_ = rate / n_os
    return {"TOF_s": t_, "TOF_h": t_ * 3600, "TOF_err_s": (t_ * k_err / k_obs_per_s) if k_err else None,
            "n_open_sites_umol": n_os * 1e6, "rate_umol_s": rate * 1e6}


def pd_profile_sigmoid(pD, kmax, pKa, kmin=0.0):
    """Single ionisation: active form = deprotonated. k = kmin + (kmax-kmin)/(1+10^(pKa-pD))."""
    return kmin + (kmax - kmin) / (1 + 10 ** (pKa - np.asarray(pD, float)))


def pd_profile_bell(pD, kmax, pKa1, pKa2):
    """Two ionisations (active = singly deprotonated): k = kmax / (1 + 10^(pKa1-pD) + 10^(pD-pKa2))."""
    pD = np.asarray(pD, float)
    return kmax / (1 + 10 ** (pKa1 - pD) + 10 ** (pD - pKa2))


def fit_pd_profile(pD, k, model="sigmoidal", k_err=None, level=0.95, pD_to_pH=False):
    """Fit k(pD). pD_to_pH: report pKa(H2O) ≈ pKa(D2O) − 0.4 too (approximate solvent isotope correction; the glass
    electrode reading correction pD = pH* + 0.4 must be applied by the user to the x values)."""
    pD = np.asarray(pD, float); k = np.asarray(k, float)
    sig = None if k_err is None else np.asarray(k_err, float)
    if model.startswith("sig"):
        f = pd_profile_sigmoid; p0 = [k.max(), float(np.median(pD)), max(k.min(), 0)]
        names = ["k_max", "pKa", "k_min"]
    else:
        f = pd_profile_bell; i = int(np.argmax(k)); p0 = [k.max() * 1.5, pD[i] - 1.5, pD[i] + 1.5]
        names = ["k_max", "pKa1", "pKa2"]
    p, cov = curve_fit(f, pD, k, p0=p0, sigma=sig, absolute_sigma=sig is not None, maxfev=50000)
    se = np.sqrt(np.diag(cov)); tc = _tcrit(len(pD) - len(p), level)
    res = k - f(pD, *p); R2 = 1 - (res ** 2).sum() / ((k - k.mean()) ** 2).sum()
    out = {"model": model, "params": dict(zip(names, p)), "ci": dict(zip(names, tc * se)), "R2": R2, "n": len(pD),
           "curve": lambda x: f(np.asarray(x, float), *p)}
    return out
