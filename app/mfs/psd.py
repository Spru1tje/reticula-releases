"""Pore-size distributions computed from an N2 (77 K) isotherm.

Methods (say which one is used on every figure / table):
  * "DFT kernel (pyGAPS, carbon slit)": non-negative least-squares fit of the openly licensed pyGAPS kernel
    DFT-N2-77K-carbon-slit (MIT) with 2nd-difference (Tikhonov) regularisation. Only for slit pores / carbons;
    for MOFs it is an approximation (other wall potential, other pore shape).
  * "HK slit" / "SF cylinder": Horvath-Kawazoe (1983) slit and Saito-Foley (1991) cylinder micropore models with the
    Cheng-Yang correction; equations and parameters ported from pyGAPS (MIT, psd_micro.py / models_hk.py).
  * "BJH": Barrett-Joyner-Halenda (1951), cylindrical pores, Kelvin equation (hemispherical meniscus on the
    desorption or adsorption branch) + Harkins-Jura thickness. Valid for mesopores (2-50 nm).
No openly licensed NLDFT/QSDFT kernel for oxide / cylindrical pores exists; see references/kernels/README.txt."""
from __future__ import annotations
import os
import numpy as np
from scipy import constants, optimize

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KERNEL_DIR = os.path.join(HERE, "references", "kernels")
KERNELS = {"DFT kernel (pyGAPS, N2 77 K, carbon slit)": "DFT-N2-77K-carbon-slit.csv"}
N2 = {"molecular_diameter": 0.3, "polarizability": 1.7403e-3, "magnetic_susceptibility": 3.6e-8, "surface_density": 6.71e18,
      "liquid_density": 0.806, "molar_mass": 28.0134}
HK_MATERIALS = {   # pyGAPS models_hk.py (MIT)
    "Carbon (HK)": {"molecular_diameter": 0.34, "polarizability": 1.02e-3, "magnetic_susceptibility": 1.35e-7, "surface_density": 3.845e19},
    "AlSi oxide ion (SF)": {"molecular_diameter": 0.276, "polarizability": 2.5e-3, "magnetic_susceptibility": 1.3e-8, "surface_density": 1.315e19},
}
VSTP = 22414.0      # cm3 STP per mol
_CACHE = {}


def load_kernel(name):
    fn = KERNELS.get(name, name)
    path = fn if os.path.isabs(fn) else os.path.join(KERNEL_DIR, fn)
    if path in _CACHE:
        return _CACHE[path]
    a = np.genfromtxt(path, delimiter=",")
    w = a[0, 1:]; p = a[1:, 0]; K = a[1:, 1:]
    _CACHE[path] = (p, w, K)
    return _CACHE[path]


def dft_fit(p, q, kernel=list(KERNELS)[0], lam=None, p_range=(1e-6, 0.995)):
    """p: relative pressure, q: cm3 STP g-1 (adsorption branch). Returns dict(w nm, dVdw cm3 g-1 nm-1, V_cum, fit (cm3 STP/g at p),
    p, q, residual_rel, lam, method)."""
    pk, w, K = load_kernel(kernel)
    p = np.asarray(p, float); q = np.asarray(q, float)
    m = (p >= max(p_range[0], pk.min())) & (p <= min(p_range[1], pk.max())) & np.isfinite(q)
    p, q = p[m], q[m]
    o = np.argsort(p); p, q = p[o], q[o]
    n = q / VSTP * 1e3                                          # mmol g-1
    lp, lk = np.log(p), np.log(pk)
    A = np.column_stack([np.interp(lp, lk, K[:, j]) for j in range(len(w))])   # mmol/g per cm3/g
    # 2nd-difference regularisation, lambda scaled to the data
    D = np.diff(np.eye(len(w)), 2, axis=0)
    if lam is None:
        lam = 0.02 * np.linalg.norm(A, 2)
    from scipy.optimize import nnls
    x, _ = nnls(np.vstack([A, lam * D]), np.concatenate([n, np.zeros(D.shape[0])]), maxiter=50 * len(w))
    fit = A @ x * VSTP / 1e3
    dw = np.gradient(w)
    return {"w": w, "dVdw": x / dw, "V_cum": np.cumsum(x), "V": x, "fit": fit, "p": p, "q": q, "lam": lam,
            "residual_rel": float(np.linalg.norm(fit - q) / np.linalg.norm(q)), "method": kernel}


# ------------------------------------------------------------------ Horvath-Kawazoe / Saito-Foley (ported from pyGAPS, MIT)
def _kirkwood_muller(ads, mat):
    p_a, m_a = ads["polarizability"] * 1e-27, ads["magnetic_susceptibility"] * 1e-27
    p_m, m_m = mat["polarizability"] * 1e-27, mat["magnetic_susceptibility"] * 1e-27
    mc2 = constants.electron_mass * constants.speed_of_light ** 2
    return 1.5 * mc2 * p_a * m_a, 6 * mc2 * p_a * p_m / (p_a / m_a + p_m / m_m)


def hk_psd(p, q, geometry="slit", material="Carbon (HK)", T=77.355, cheng_yang=True, p_max=0.1):
    """Micropore PSD. geometry 'slit' (Horvath-Kawazoe) or 'cylinder' (Saito-Foley). q in cm3 STP g-1."""
    ads, mat = N2, HK_MATERIALS[material]
    p = np.asarray(p, float); q = np.asarray(q, float)
    m = (p > 0) & (p <= p_max) & np.isfinite(q)
    p, q = p[m], q[m]; o = np.argsort(p); p, q = p[o], q[o]
    n = q / VSTP * 1e3                                          # mmol/g
    d_a, d_m = ads["molecular_diameter"], mat["molecular_diameter"]
    a_a, a_m = _kirkwood_muller(ads, mat)
    d = (d_a + d_m) / 2
    NRT = constants.Avogadro / constants.gas_constant / T
    nsum = ads["surface_density"] * a_a + mat["surface_density"] * a_m
    if geometry == "slit":
        sig = 0.8583742 * d; s4, s10 = sig ** 4 / 3, sig ** 10 / 9
        cc = NRT * nsum / (sig * 1e-9) ** 4; ct = s10 / d ** 9 - s4 / d ** 3

        def pot(L):
            return cc / (L - 2 * d) * (s4 / (L - d) ** 3 - s10 / (L - d) ** 9 + ct)
        bound, geo = 2 * d, 1
    else:
        cc = 0.75 * np.pi * NRT * nsum / (d * 1e-9) ** 4
        ak, bk = [1.0], [1.0]
        for k in range(1, 2000):
            ak.append(((-4.5 - k) / k) ** 2 * ak[-1]); bk.append(((-1.5 - k) / k) ** 2 * bk[-1])

        def pot(L):
            r = d / L; r4 = r ** 4; r10 = 0.65625 * r ** 10
            s = r10 - r4
            for k in range(1, int(L * 25)):
                s += (1 / (k + 1) * (1 - r) ** (2 * k)) * (ak[k] * r10 - bk[k] * r4)
            return cc * s
        bound, geo = d, 2
    cov = n / (n.max() * 1.01)
    L_out = []
    for pp, cv in zip(p, cov):
        corr = (1 + 1 / cv * np.log(1 - cv)) if cheng_yang else 0.0
        f = lambda L: (np.exp(pot(L) - corr) - pp) ** 2
        r = optimize.minimize_scalar(f, method="bounded", bounds=(bound * 1.0001, 50))
        L_out.append(r.x)
        if r.x > 10 / geo:
            break
    L_out = np.asarray(L_out)
    w = (L_out - d_m) if geometry == "slit" else (2 * L_out - d_m)
    k = len(w)
    V = n[:k] * ads["molar_mass"] / ads["liquid_density"] / 1000      # cm3 liquid / g
    o = np.argsort(w); w, V = w[o], np.maximum.accumulate(V[o])
    # resample the cumulative curve on a regular grid before differentiating: the raw w(p) points crowd together near
    # the top of the micropore range, and point-to-point differences there give spikes
    grid = np.arange(w[0], w[-1], 0.01)
    if len(grid) < 5:
        grid = np.linspace(w[0], w[-1], 5)
    Vg = np.interp(grid, w, V)
    win = 5; Vs = np.convolve(np.pad(Vg, win // 2, mode="edge"), np.ones(win) / win, mode="valid")
    dV = np.gradient(Vs, grid)
    return {"w": grid, "dVdw": np.clip(dV, 0, None), "V_cum": Vg, "method": ("Horvath–Kawazoe slit" if geometry == "slit" else "Saito–Foley cylinder")
            + f" ({material}, {'Cheng–Yang corr.' if cheng_yang else 'original'})"}


# ------------------------------------------------------------------ BJH
def harkins_jura(p):
    return 0.1 * np.sqrt(13.99 / (0.034 - np.log10(p)))           # nm


def kelvin_radius(p, T=77.355, gamma=8.88e-3, Vm=34.68e-6):
    """Kelvin core radius (nm), hemispherical meniscus, wetting (cos θ = 1)."""
    return -2 * gamma * Vm / (constants.gas_constant * T * np.log(p)) * 1e9


def bjh_psd(p, q, branch="desorption", p_lo=0.12, p_hi=0.99):
    """BJH (cylinder) with the Harkins–Jura thickness. p, q of one branch; returns dV/dw (cm3 g-1 nm-1) vs pore width (nm)."""
    p = np.asarray(p, float); q = np.asarray(q, float)
    m = (p >= p_lo) & (p <= p_hi) & np.isfinite(q)
    p, q = p[m], q[m]
    o = np.argsort(p)[::-1]; p, q = p[o], q[o]                     # from high to low pressure
    Vl = q / VSTP * N2["molar_mass"] / N2["liquid_density"]         # cm3 liquid / g
    rk = kelvin_radius(p); t = harkins_jura(p); rp = rk + t
    Vp, wp, sumA = [], [], 0.0            # sumA = Σ 2V/r of the pores emptied so far (cm3 g-1 nm-1; × Δt (nm) -> cm3 g-1)
    for i in range(len(p) - 1):
        dV = Vl[i] - Vl[i + 1]
        dt = t[i] - t[i + 1]
        r_avg = 0.5 * (rp[i] + rp[i + 1]); rk_avg = 0.5 * (rk[i] + rk[i + 1])
        Rn = (r_avg / (rk_avg + dt)) ** 2
        v = max(0.0, Rn * (dV - dt * sumA))
        Vp.append(v); wp.append(2 * r_avg)
        sumA += 2 * v / r_avg
    wp, Vp = np.asarray(wp), np.asarray(Vp)
    o = np.argsort(wp); wp, Vp = wp[o], Vp[o]
    edges = np.concatenate([[wp[0] * 0.95], 0.5 * (wp[1:] + wp[:-1]), [wp[-1] * 1.05]])
    dw = np.diff(edges)
    return {"w": wp, "dVdw": Vp / dw, "V_cum": np.cumsum(Vp), "method": f"BJH ({branch}, cylinder, Harkins–Jura)"}
