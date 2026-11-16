"""v2.13: the vectorised Rouquerol window search must give the same result as the original one-lstsq-per-window code."""
import math, os, sys
import numpy as np
import pytest
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from mfs import analysis as A, readers as R  # noqa: E402
CM3STP_TO_MOL, N2_CROSS = A.CM3STP_TO_MOL, A.N2_CROSS


# ---- the v2.12 implementation, kept here as the reference ----
def _bet_reference(p, q, min_points=4, pmax=0.35):
    """BET with the four Rouquerol consistency criteria (Rouquerol, Llewellyn & Rouquerol 2007; ISO 9277:2022).
    Returns dict with area, C, n_m, range, R2, criteria flags and the full linear-plot data."""
    p = np.asarray(p, float); q = np.asarray(q, float)
    o = np.argsort(p); p, q = p[o], q[o]
    m = (p > 0) & (p < pmax)
    p, q = p[m], q[m]
    rq = q * (1 - p)
    # criterion 1: n(1-p) increases -> restrict to p <= p at max of n(1-p)
    imax = int(np.argmax(rq))
    best = None
    cand = []
    for i in range(0, imax + 1):
        for j in range(i + min_points - 1, imax + 1):
            pp, qq = p[i:j + 1], q[i:j + 1]
            yb = pp / (qq * (1 - pp))
            A = np.vstack([pp, np.ones_like(pp)]).T
            (s, b), res, *_ = np.linalg.lstsq(A, yb, rcond=None)
            if b <= 0 or s <= 0:
                continue  # criterion 2: C > 0
            C = s / b + 1
            nm = 1 / (s + b)
            pred = A @ [s, b]
            ss = ((yb - yb.mean()) ** 2).sum()
            R2 = 1 - ((yb - pred) ** 2).sum() / ss if ss > 0 else 0
            # criterion 3: p at n_m within range
            p_nm = np.interp(nm, qq, pp) if qq.min() <= nm <= qq.max() else np.nan
            c3 = not np.isnan(p_nm)
            # criterion 4: 1/(sqrt(C)+1) ~ p_nm within 20 % (often 10 %)
            p_th = 1 / (math.sqrt(C) + 1)
            c4 = c3 and abs(p_th - p_nm) / p_nm < 0.2
            mono = bool(np.all(np.diff(qq * (1 - pp)) > 0))
            cand.append(dict(i=i, j=j, slope=s, intercept=b, C=C, n_m=nm, R2=R2, c1=mono, c2=True, c3=c3, c4=c4,
                             p_nm=p_nm, p_th=p_th, n=j - i + 1))
    if not cand:
        return None
    valid = [c for c in cand if c["c1"] and c["c3"] and c["c4"]] or [c for c in cand if c["c1"] and c["c3"]] or cand
    # prefer highest R2, then more points
    best = max(valid, key=lambda c: (round(c["R2"], 4), c["n"]))
    area = best["n_m"] * CM3STP_TO_MOL * 6.02214076e23 * N2_CROSS
    best.update(area=area, p_lo=p[best["i"]], p_hi=p[best["j"]],
                all_criteria=bool(best["c1"] and best["c3"] and best["c4"]),
                plot_p=p[: imax + 1], plot_y=p[: imax + 1] / (q[: imax + 1] * (1 - p[: imax + 1])),
                rq_p=p, rq=rq)
    best['_ncand'] = len(cand)
    return best


def _same(a, b):
    assert (a is None) == (b is None)
    if a is None:
        return
    assert (a["i"], a["j"], a["n"]) == (b["i"], b["j"], b["n"])
    for k in ("area", "C", "n_m", "slope", "intercept", "R2"):
        assert a[k] == pytest.approx(b[k], rel=1e-7, abs=1e-9), k
    for k in ("c1", "c2", "c3", "c4", "all_criteria"):
        assert bool(a[k]) == bool(b[k]), k
    assert len(a["windows"]) == b["_ncand"]


def _examples():
    out = []
    d = os.path.join(HERE, "example_data", "Sorption")
    for f in sorted(os.listdir(d)):
        if f.endswith((".xlsx", ".csv")):
            iso = R.read_sorption(os.path.join(d, f), f)["iso"]; a = iso[iso.branch == "ads"]
            out.append((f, a.p_rel.values, a.q.values))
    return out


@pytest.mark.parametrize("name,p,q", _examples(), ids=[e[0] for e in _examples()])
def test_example_isotherms_identical(name, p, q):
    for mp in (4, 3):
        _same(A._bet_rouquerol(p, q, min_points=mp), _bet_reference(p, q, min_points=mp))


def test_random_noisy_isotherms_identical():
    rng = np.random.default_rng(42)
    for t in range(40):
        n = int(rng.integers(8, 120))
        p = np.sort(np.concatenate([np.logspace(-6, -2, n // 2), rng.uniform(0.01, 0.4, n - n // 2)]))
        nm_, C = rng.uniform(20, 400), 10 ** rng.uniform(0.7, 3.8)
        q = nm_ * C * p / ((1 - p) * (1 + (C - 1) * p)) * (1 + rng.uniform(0, 0.02) * rng.normal(size=p.size))
        q = np.abs(q) + 1e-3
        _same(A._bet_rouquerol(p, q), _bet_reference(p, q))


def test_degenerate_inputs():
    assert A._bet_rouquerol(np.array([0.1, 0.2]), np.array([1.0, 2.0])) is None
    assert A._bet_rouquerol(np.array([]), np.array([])) is None
    assert A.bet_rouquerol(np.array([0.05, 0.1, 0.2]), np.array([10.0, 14.0, 17.0]))["few_points"] is True


def test_speed_is_reasonable():
    import time
    p = np.concatenate([np.logspace(-6, -2, 60), np.linspace(0.011, 0.34, 90)])
    q = 100 * 3000 * p / ((1 - p) * (1 + 2999 * p))
    t = time.perf_counter(); A.bet_rouquerol(p, q); assert time.perf_counter() - t < 0.5     # was 0.44 s before v2.13 (CI machines vary)
