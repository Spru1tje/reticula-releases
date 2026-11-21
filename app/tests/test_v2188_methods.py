"""v2.18.8: new non-isothermal kinetics (Kissinger, Ozawa-Flynn-Wall) and Sips isotherm, each tested on data with a known answer."""
import numpy as np

from mfs import kinetics_nonisothermal as KN
from mfs import isotherm_models as IM

R = 8.314462618


def _tp_for_kissinger(betas, ea, a_pre):
    """Solve ln(beta/T^2) = ln(A R/Ea) - Ea/(R T) for T by bisection (exact Kissinger peak temperatures)."""
    out = []
    for b in betas:
        lo, hi = 200.0, 2000.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            f = np.log(b / mid ** 2) - (np.log(a_pre * R / ea) - ea / (R * mid))
            lo, hi = (mid, hi) if f > 0 else (lo, mid)
        out.append(0.5 * (lo + hi))
    return np.array(out)


def test_kissinger_recovers_known_activation_energy():
    ea = 120_000.0
    betas = np.array([2.0, 5.0, 10.0, 20.0]) / 60.0          # K/s
    tp = _tp_for_kissinger(betas, ea, a_pre=1e12)
    res = KN.kissinger(betas, tp)
    assert res["status"] == "ok"
    assert abs(res["Ea_kJ"] - 120.0) < 0.05
    assert res["R2"] > 0.999999


def test_ozawa_flynn_wall_recovers_known_activation_energy():
    ea = 90_000.0
    betas = np.array([2.0, 5.0, 10.0, 20.0])                  # K/min, any unit works
    # exact OFW relation: log10(beta) = c - 0.4567 Ea / (R T)  ->  choose c, then T
    c = 3.0
    t_alpha = 0.4567 * ea / (R * (c - np.log10(betas)))
    res = KN.ozawa_flynn_wall(betas, t_alpha)
    assert res["status"] == "ok"
    assert abs(res["Ea_kJ"] - 90.0) < 0.05


def test_kissinger_needs_three_rates():
    assert KN.kissinger([5.0, 10.0], [400.0, 410.0])["status"].startswith("needs")


def test_sips_recovers_known_parameters():
    qm, k, n = 4.0, 0.8, 0.7
    p = np.linspace(0.05, 10.0, 25)
    q = IM.sips(p, qm, k, n)
    res = IM.fit_sips(p, q)
    assert res["status"] == "ok"
    assert abs(res["qm"] - qm) < 1e-4
    assert abs(res["K"] - k) < 1e-4
    assert abs(res["n"] - n) < 1e-4
    assert res["R2"] > 0.999999


def test_sips_reduces_to_langmuir_when_n_is_one():
    p = np.linspace(0.1, 5.0, 20)
    q = IM.sips(p, 3.0, 1.2, 1.0)
    langmuir = 3.0 * (1.2 * p) / (1 + 1.2 * p)
    assert np.allclose(q, langmuir)


def test_sips_too_few_points():
    assert IM.fit_sips([1.0, 2.0], [1.0, 1.5])["status"].startswith("needs")
