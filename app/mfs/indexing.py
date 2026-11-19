"""Automatic indexing of a powder pattern with an unknown phase: finds the unit cell from the peak positions.

v2.18.4 (rewritten for a higher and more honest M20; benchmark on 31 real structures in tools/benchmark_index.py):
  * every candidate cell is REFINED by robust least squares on the indexed lines, together with a zero shift (2θ offset) - the figure of merit is computed with the refined cell;
  * lattice centring is tested (P, I, F, A, B, C, R): a centred lattice has fewer allowed lines, so the true conventional cell of NaCl, CeO2, anatase … gets its due M20;
  * the M20 counts the indexed lines (lines that do not fit - an impurity - are listed as 'unindexed' instead of destroying the figure of merit);
  * systems: cubic, tetragonal, hexagonal / trigonal (hexagonal axes), orthorhombic, and (optionally, slower) monoclinic.
Quality: de Wolff M20 (> 20 = reliable, 10-20 = probable, 4-10 = possible, < 4 = unreliable) and the number of unindexed lines.
Physical limits: every cell edge >= 2.4 Å and a volume >= 30 Å³. The cell is refined afterwards with Le Bail (Phase.from_cell): a good cell gives a good Le Bail fit."""
from __future__ import annotations

import itertools
import math

import numpy as np
from scipy.optimize import least_squares

_CUBIC_N = [n for n in range(1, 60) if any(a * a + b * b + c * c == n for a in range(8) for b in range(8) for c in range(8))]
_H = np.arange(0, 8)
_GH, _GK, _GL = (v.ravel() for v in np.meshgrid(np.arange(-7, 8), np.arange(-7, 8), _H, indexing="ij"))
MIN_EDGE, MIN_VOLUME, MAX_VOLUME = 2.4, 30.0, 15000.0
CENTRINGS = {
    "P": lambda h, k, l: np.ones_like(h, bool),
    "I": lambda h, k, l: (h + k + l) % 2 == 0,
    "F": lambda h, k, l: ((h % 2 == k % 2) & (k % 2 == l % 2)),
    "A": lambda h, k, l: (k + l) % 2 == 0,
    "B": lambda h, k, l: (h + l) % 2 == 0,
    "C": lambda h, k, l: (h + k) % 2 == 0,
    "R": lambda h, k, l: (-h + k + l) % 3 == 0,
}
_ALLOWED = {"cubic": "PIF", "tetragonal": "PI", "hexagonal": "PR", "orthorhombic": "PICFAB", "monoclinic": "PC"}


def _qs(two_theta, wavelength):
    return (2 * np.sin(np.radians(np.asarray(two_theta, float) / 2)) / wavelength) ** 2


def _tt_of_q(q, wavelength):
    return 2 * np.degrees(np.arcsin(np.clip(wavelength * np.sqrt(np.maximum(q, 0)) / 2, 0, 1)))


def _volume(system, p):
    if system == "cubic":
        return p[0] ** 3
    if system == "tetragonal":
        return p[0] ** 2 * p[1]
    if system == "hexagonal":
        return p[0] ** 2 * p[1] * math.sqrt(3) / 2
    if system == "monoclinic":
        return p[0] * p[1] * p[2] * math.sin(math.radians(p[3]))
    return p[0] * p[1] * p[2]


def _q_of(system, p, H, K, L):
    if system == "cubic":
        return (H ** 2 + K ** 2 + L ** 2) / p[0] ** 2
    if system == "tetragonal":
        return (H ** 2 + K ** 2) / p[0] ** 2 + L ** 2 / p[1] ** 2
    if system == "hexagonal":
        return 4.0 / 3.0 * (H ** 2 + H * K + K ** 2) / p[0] ** 2 + L ** 2 / p[1] ** 2
    if system == "monoclinic":                      # unique axis b: p = (a, b, c, beta in degrees)
        s = math.sin(math.radians(p[3])); c_ = math.cos(math.radians(p[3]))
        return (H ** 2 / p[0] ** 2 + L ** 2 / p[2] ** 2 - 2 * H * L * c_ / (p[0] * p[2])) / s ** 2 + K ** 2 / p[1] ** 2
    return H ** 2 / p[0] ** 2 + K ** 2 / p[1] ** 2 + L ** 2 / p[2] ** 2


def _hkl_grid(system, p, qmax):
    hm = [int(math.ceil(v * math.sqrt(qmax))) + 1 for v in ((p[0], p[0], p[-1]) if system in ("cubic", "tetragonal", "hexagonal") and len(p) < 3
                                                             else (p[0], p[1] if len(p) > 1 else p[0], p[2] if len(p) > 2 else p[-1]))]
    if system == "cubic":
        hm = [int(math.ceil(p[0] * math.sqrt(qmax))) + 1] * 3
    elif system in ("tetragonal", "hexagonal"):
        hm = [int(math.ceil(p[0] * math.sqrt(qmax))) + 1] * 2 + [int(math.ceil(p[1] * math.sqrt(qmax))) + 1]
    elif system == "monoclinic":
        hm = [int(math.ceil(p[0] * math.sqrt(qmax) / math.sin(math.radians(p[3])))) + 1, int(math.ceil(p[1] * math.sqrt(qmax))) + 1,
              int(math.ceil(p[2] * math.sqrt(qmax) / math.sin(math.radians(p[3])))) + 1]
    H, K, L = np.meshgrid(np.arange(-hm[0], hm[0] + 1), np.arange(-hm[1], hm[1] + 1), np.arange(0, hm[2] + 1), indexing="ij")
    return H.ravel(), K.ravel(), L.ravel()


def _lines(system, p, qmax, centring="P"):
    """All unique q = 1/d² of the cell up to qmax for the lattice centring (no other absence rules)."""
    H, K, L = _hkl_grid(system, p, qmax)
    q = _q_of(system, p, H, K, L)
    ok = (q > 1e-9) & (q <= qmax * 1.02) & CENTRINGS[centring](H, K, L)
    return np.unique(np.round(q[ok], 7))


def _calc_q(system, p, qmax):
    return _lines(system, p, qmax, "P")


def _nearest(qc, qo):
    j = np.searchsorted(qc, qo)
    left, right = np.clip(j - 1, 0, len(qc) - 1), np.clip(j, 0, len(qc) - 1)
    dl, dr = np.abs(qc[left] - qo), np.abs(qc[right] - qo)
    return np.where(dl <= dr, left, right), np.minimum(dl, dr)


def _m_n(q_obs, q_calc):
    """de Wolff M(n) for the observed peaks q_obs (sorted) against the calculated lines (every observed line counts)."""
    qo = np.sort(q_obs)
    _, dq = _nearest(q_calc, qo)
    n_calc = int(np.sum(q_calc <= qo[-1] * 1.001))
    return qo[-1] / (2 * max(float(np.mean(dq)), 1e-9) * max(n_calc, 1)), dq


def _m_robust(qo, qc, tol_q):
    """M(n) over the INDEXED lines (|ΔQ| <= tol_q); lines that do not fit are counted as unindexed. Returns (M, n_unindexed, mean|ΔQ|)."""
    _, dq = _nearest(qc, qo)
    ok = dq <= tol_q
    n_ok = int(ok.sum())
    if n_ok < 4:
        return 0.0, int(len(qo) - n_ok), float("nan")
    qmax = float(qo[ok].max())
    n_calc = max(int(np.sum(qc <= qmax * 1.001)), 1)
    m = qmax / (2 * max(float(np.mean(dq[ok])), 1e-9) * n_calc)
    return m * (n_ok / len(qo)) ** 2, int(len(qo) - n_ok), float(np.mean(dq[ok]))


def _tol_q(q, wavelength, tt_tol):
    s = np.sqrt(np.clip(q, 1e-12, None)) * wavelength / 2
    return 2 * s * np.sqrt(np.clip(1 - s * s, 0, 1)) / wavelength ** 2 * math.radians(tt_tol) + 1e-7


def _bounds(system, p):
    if system == "cubic":
        return [(MIN_EDGE, 90)]
    if system in ("tetragonal", "hexagonal"):
        return [(MIN_EDGE, 90)] * 2
    if system == "monoclinic":
        return [(MIN_EDGE, 90)] * 3 + [(90.0, 135.0)]
    return [(MIN_EDGE, 90)] * 3


def _grid_for(system, p, qmax, centring, margin=1.12):
    """Fixed (h, k, l) list (centring applied) that covers the cell p up to qmax, also when the lengths grow by `margin` during the refinement."""
    pg = np.array(p, float)
    n_len = 1 if system == "cubic" else 2 if system in ("tetragonal", "hexagonal") else 3
    pg[:n_len] *= margin
    H, K, L = _hkl_grid(system, pg, qmax * margin * margin)
    sel = CENTRINGS[centring](H, K, L)
    return H[sel], K[sel], L[sel]


def _basis(system, H, K, L):
    """Design matrix of Q = Σ x_j b_j(hkl) for the reciprocal parameters x of the system."""
    if system == "cubic":
        return (H ** 2 + K ** 2 + L ** 2)[:, None]
    if system == "tetragonal":
        return np.stack([H ** 2 + K ** 2, L ** 2], 1)
    if system == "hexagonal":
        return np.stack([4.0 / 3.0 * (H ** 2 + H * K + K ** 2), L ** 2], 1)
    if system == "monoclinic":
        return np.stack([H ** 2, K ** 2, L ** 2, H * L], 1)
    return np.stack([H ** 2, K ** 2, L ** 2], 1)


def _x_to_params(system, x):
    """Reciprocal parameters -> cell edges (and beta); None when the metric is not positive."""
    x = np.asarray(x, float)
    if np.any(x[: 3 if system in ("orthorhombic", "monoclinic") else len(x)] <= 0):
        return None
    if system == "cubic":
        return np.array([1 / math.sqrt(x[0])])
    if system == "tetragonal":
        return np.array([1 / math.sqrt(x[0]), 1 / math.sqrt(x[1])])
    if system == "hexagonal":
        return np.array([1 / math.sqrt(x[0]), 1 / math.sqrt(x[1])])
    if system == "monoclinic":
        A, B, C, D = x
        if abs(D) >= 2 * math.sqrt(A * C):
            return None
        cb = -D / (2 * math.sqrt(A * C))
        s2 = 1 - cb * cb
        return np.array([1 / math.sqrt(A * s2), 1 / math.sqrt(B), 1 / math.sqrt(C * s2), math.degrees(math.acos(cb))])
    return np.array([1 / math.sqrt(x[0]), 1 / math.sqrt(x[1]), 1 / math.sqrt(x[2])])


def _params_to_x(system, p):
    p = np.asarray(p, float)
    if system == "cubic":
        return np.array([1 / p[0] ** 2])
    if system == "tetragonal":
        return np.array([1 / p[0] ** 2, 1 / p[1] ** 2])
    if system == "hexagonal":
        return np.array([1 / p[0] ** 2, 1 / p[1] ** 2])
    if system == "monoclinic":
        cb, sb2 = math.cos(math.radians(p[3])), math.sin(math.radians(p[3])) ** 2
        return np.array([1 / (p[0] ** 2 * sb2), 1 / p[1] ** 2, 1 / (p[2] ** 2 * sb2), -2 * cb / (p[0] * p[2] * sb2)])
    return np.array([1 / p[0] ** 2, 1 / p[1] ** 2, 1 / p[2] ** 2])


def _grid_for(system, p, qmax, centring, margin=1.12):
    """Fixed (h, k, l) list (centring applied) that covers the cell p up to qmax, also when the lengths grow by `margin` during the refinement."""
    pg = np.array(p, float)
    n_len = 1 if system == "cubic" else 2 if system in ("tetragonal", "hexagonal") else 3
    pg[:n_len] *= margin
    H, K, L = _hkl_grid(system, pg, qmax * margin * margin)
    sel = CENTRINGS[centring](H, K, L)
    return H[sel], K[sel], L[sel]


def refine_cell(system, p0, tt_obs, wavelength, centring="P", tt_tol=0.08, zero0=0.0, fit_zero=True, rounds=4):
    """Refinement of the cell (and a zero shift) on the lines that are indexed, by iterated linear least squares: Q = Σ x_j b_j(hkl) is linear in the reciprocal parameters
    and a zero shift z enters as q_obs(2θ − z) ≈ q_obs(2θ) − z·dq/d2θ. The assignment (nearest calculated line) is repeated, with Huber weights against impurity lines.
    Returns (params, zero, indexed mask, rms 2θ residual (°), calculated lines q)."""
    tt_obs = np.asarray(tt_obs, float)
    p = np.array(p0, float)
    z = zero0
    qmax = float(_qs(tt_obs.max() + 0.2, wavelength))
    H, K, L = _grid_for(system, p, qmax, centring)
    Bm = _basis(system, H, K, L)
    qo0 = _qs(tt_obs, wavelength)
    g = 2.0 / wavelength ** 2 * np.sin(np.radians(tt_obs)) * math.pi / 180.0          # dq / d(2θ) per degree
    x = _params_to_x(system, p)
    for rnd in range(rounds):
        qo = qo0 - z * g
        qall = Bm @ x
        ok = qall > 1e-9
        idx_all = np.nonzero(ok)[0]
        order = idx_all[np.argsort(qall[ok])]
        qs = qall[order]
        j = np.searchsorted(qs, qo)
        left, right = np.clip(j - 1, 0, len(qs) - 1), np.clip(j, 0, len(qs) - 1)
        pick = np.where(np.abs(qs[left] - qo) <= np.abs(qs[right] - qo), left, right)
        dq = np.abs(qs[pick] - qo)
        tolq = _tol_q(qo, wavelength, tt_tol)
        mask = dq <= tolq
        if mask.sum() < (len(x) + (1 if fit_zero else 0) + 1):
            break
        rows = order[pick[mask]]
        A = Bm[rows]
        if fit_zero:
            A = np.hstack([A, g[mask][:, None]])
        yv = qo0[mask]
        w = 1.0 / np.maximum(g[mask], 1e-9)                                              # equal weight in 2θ
        wr = np.ones(len(yv))
        sol = None
        for _ in range(3):                                                               # Huber re-weighting
            sw = (w * wr)[:, None]
            try:
                sol, *_ = np.linalg.lstsq(A * sw, yv * sw[:, 0], rcond=None)
            except np.linalg.LinAlgError:
                sol = None
                break
            res_tt = (A @ sol - yv) * w
            s = max(1.4826 * float(np.median(np.abs(res_tt))), 0.01)
            wr = np.where(np.abs(res_tt) <= 1.5 * s, 1.0, 1.5 * s / np.maximum(np.abs(res_tt), 1e-12))
        if sol is None:
            break
        xn = sol[: len(x)]
        pn = _x_to_params(system, xn)
        if pn is None:
            break
        x = xn
        p = pn
        if fit_zero:
            z = float(np.clip(sol[len(x)], -0.3, 0.3))
    qall = Bm @ x
    qc = np.unique(np.round(qall[(qall > 1e-9) & (qall <= qmax * 1.05)], 7))
    qo = qo0 - z * g
    idx, dq = _nearest(qc, qo)
    mask = dq <= _tol_q(qo, wavelength, tt_tol)
    rms = float(np.sqrt(np.mean((_tt_of_q(qc[idx[mask]], wavelength) - _tt_of_q(qo[mask], wavelength)) ** 2))) if mask.any() else float("nan")
    return p, z, mask, rms, qc

def _quality(m):
    return "reliable" if m >= 20 else "probable" if m >= 10 else "possible" if m >= 4 else "unreliable"


_TBL = None


def _table():
    """Small hkl table (h -7..7, k -7..7, l 0..7) with the design matrices per system, for the vectorised pre-test of the seeds."""
    global _TBL
    if _TBL is None:
        ok = ~((_GH == 0) & (_GK == 0) & (_GL == 0))
        H, K, L = _GH[ok].astype(float), _GK[ok].astype(float), _GL[ok].astype(float)
        _TBL = {s: _basis(s, H, K, L) for s in ("cubic", "tetragonal", "hexagonal", "orthorhombic", "monoclinic")}
    return _TBL


def _prefilter(system, P, qo8, tol8, min_hit):
    """Vectorised: which seeds P (N x n) explain at least `min_hit` of the lowest lines qo8 (each within tol8)? Returns the boolean mask."""
    if len(P) == 0:
        return np.zeros(0, bool)
    Bt = _table()[system]
    ok_all = np.zeros(len(P), bool)
    for s0 in range(0, len(P), 4000):
        X = np.array([_params_to_x(system, p) for p in P[s0:s0 + 4000]])
        Q = Bt @ X.T                                            # (T, n)
        hits = np.zeros(X.shape[0], int)
        for qi, ti in zip(qo8, tol8):
            hits += (np.abs(Q - qi).min(axis=0) <= ti)
        ok_all[s0:s0 + 4000] = hits >= min_hit
    return ok_all


def _seeds(system, q_obs, n_use, stage):
    """Candidate parameter sets (cell edges) for one crystal system. stage 1 = small indices (fast), stage 2 = larger indices (thorough)."""
    out = []
    if system == "cubic":
        for N in _CUBIC_N[: (16 if stage == 1 else 26)]:
            for k in range(min(5, n_use)):
                out.append([math.sqrt(N / q_obs[k])])
    elif system in ("tetragonal", "hexagonal"):
        f = 1.0 if system == "tetragonal" else 4.0 / 3.0
        nm = 5 if stage == 1 else 7
        ms = sorted({a * a + b * b for a in range(0, nm) for b in range(0, nm)}) if system == "tetragonal" else \
            sorted({a * a + a * b + b * b for a in range(0, nm) for b in range(0, nm)})
        nl = 7 if stage == 1 else 9
        for k1, k2 in itertools.combinations(range(min(7, n_use)), 2):
            for m1, l1, m2, l2 in itertools.product(ms[:7 if stage == 1 else 10], range(0, nl), ms[:7 if stage == 1 else 10], range(0, nl)):
                if (m1, l1) == (m2, l2) or (m1 == 0 and l1 == 0) or (m2 == 0 and l2 == 0):
                    continue
                det = f * m1 * l2 * l2 - f * m2 * l1 * l1
                if abs(det) < 1e-9:
                    continue
                Ac = (q_obs[k1] * l2 * l2 - q_obs[k2] * l1 * l1) / det
                Cc = (f * m1 * q_obs[k2] - f * m2 * q_obs[k1]) / det
                if Ac > 0 and Cc > 0:
                    out.append([1 / math.sqrt(Ac), 1 / math.sqrt(Cc)])
    elif system == "orthorhombic":
        rng_ = range(0, 2) if stage == 1 else range(0, 3)
        trip = [(h, k, l) for h in rng_ for k in rng_ for l in rng_ if (h, k, l) != (0, 0, 0)]
        A_all = np.array([[[a_[0] ** 2, a_[1] ** 2, a_[2] ** 2], [b_[0] ** 2, b_[1] ** 2, b_[2] ** 2], [c_[0] ** 2, c_[1] ** 2, c_[2] ** 2]]
                          for a_, b_, c_ in itertools.product(trip, trip, trip)], float)
        A_all = A_all[np.abs(np.linalg.det(A_all)) > 1e-6]
        for ks in itertools.combinations(range(min(6 if stage == 1 else 7, n_use)), 3):
            sol = np.linalg.solve(A_all, np.broadcast_to(q_obs[list(ks)], (len(A_all), 3))[..., None])[..., 0]
            good = np.all(sol > 0, axis=1)
            for s_ in sol[good]:
                out.append(list(1 / np.sqrt(s_)))
    elif system == "monoclinic":
        out = [list(s) for s in _monoclinic_seeds(q_obs, n_use, 0, 0, None)]
    return out


def _unique_rows(P, step=0.04):
    seen, keep = set(), []
    for i, p in enumerate(P):
        key = tuple(round(v / step) for v in p[: (1 if len(p) == 1 else 3 if len(p) >= 3 else len(p))]) + ((round(p[3]) ,) if len(p) == 4 else ())
        if key in seen:
            continue
        seen.add(key)
        keep.append(i)
    return keep


def index(two_theta, wavelength=1.540593, tt_tol=0.08, systems=("cubic", "tetragonal", "hexagonal", "orthorhombic"), max_solutions=6, n_fit=20, refine=True, centring=True,
          max_refine=80, progress=None, thorough=None):
    """Returns the best cells, best first: [{system, cell (a, b, c, alpha, beta, gamma), centring, zero_shift, M20, n_unindexed, n_used, volume, quality, note, rms_2theta}].
    two_theta: peak positions (°, Kα1). progress: optional callable(fraction, text). The search runs in stages: stage 1 (small Miller indices, fast) and, when it finds no
    reliable cell, stage 2 (larger indices; monoclinic when it is in `systems` or `thorough`)."""
    prog = progress or (lambda f, t: None)
    tt = np.sort(np.asarray(two_theta, float))
    if len(tt) < 6:
        return []
    q_obs = _qs(tt, wavelength)
    n_use = min(len(tt), n_fit)
    qo = q_obs[:n_use]
    ttu = tt[:n_use]
    n8 = min(8, n_use)
    tt_seed = tt_tol * 1.6                                      # the seed search is loose (zero shift and noise); the refinement tightens it
    tol8 = _tol_q(qo[:n8], wavelength, tt_seed)
    min_hit = max(5, n8 - 3)
    quota = {"cubic": 100, "tetragonal": 160, "hexagonal": 160, "orthorhombic": 250, "monoclinic": 250}
    sols_all = []

    def run_stage(stage, systs, f0, f1):
        sols = []
        for si, system in enumerate(systs):
            prog(f0 + (f1 - f0) * si / max(len(systs), 1), f"stage {stage}: {system} cells")
            P = _seeds(system, q_obs, n_use, stage)
            if not P:
                continue
            P = np.array([p for p in P if len(p) and all(np.isfinite(p)) and not any(v < MIN_EDGE for v in p[: 3 if system != "cubic" else 1])
                          and MIN_VOLUME <= _volume(system, p) <= MAX_VOLUME and not any(v > 90 for v in p[:3])], float)
            if len(P) == 0:
                continue
            keep = _unique_rows(P, 0.03)
            P = P[keep]
            ok = _prefilter(system, P, qo[:n8], tol8, min_hit)
            P = P[ok]
            if len(P) == 0:
                continue
            # score survivors by the number of explained lines, then refine the best ones
            sc = []
            tol_seed = _tol_q(qo, wavelength, tt_seed)
            cens = _ALLOWED[system] if centring else "P"
            for p in P:
                best_m, best_n = 0.0, 0
                for cen in cens:
                    qc = _lines(system, p, float(qo[-1]), cen)
                    if len(qc) == 0:
                        continue
                    _, dq = _nearest(qc, qo)
                    good = dq <= tol_seed
                    n_g = int(good.sum())
                    if n_g < 4:
                        continue
                    n_calc = max(int(np.sum(qc <= float(qo[good].max()) * 1.001)), 1)
                    m_seed = float(qo[good].max()) / (2 * max(float(np.mean(dq[good])), 1e-9) * n_calc) * (n_g / len(qo)) ** 2   # line density matters: a dense cell explains anything
                    if m_seed > best_m:
                        best_m, best_n = m_seed, n_g
                sc.append((best_n, best_m))
            order = sorted(range(len(P)), key=lambda i: -sc[i][1])
            kept, seen = [], set()
            for i in order:
                if sc[i][0] < max(6, int(0.7 * n_use)):
                    continue
                key = tuple(round(v / 0.05) for v in P[i][: 3 if system != "cubic" else 1])
                if key in seen:
                    continue
                seen.add(key)
                kept.append(P[i])
                if len(kept) >= quota[system]:
                    break
            for p in kept:
                best = None
                for cen in (_ALLOWED[system] if centring else "P"):
                    try:
                        pr, z, mask, rms, qc = refine_cell(system, p, ttu, wavelength, cen, tt_tol, 0.0, True) if refine else (np.array(p), 0.0, np.ones(len(ttu), bool), 0.0, _lines(system, p, float(qo[-1]), cen))
                    except Exception:
                        continue
                    if any(v < MIN_EDGE for v in pr[: 3 if system != "cubic" else 1]) or not (MIN_VOLUME <= _volume(system, pr) <= MAX_VOLUME):
                        continue
                    qo_z = _qs(ttu - z, wavelength)
                    m, n_un, dqm = _m_robust(qo_z, qc, _tol_q(qo_z, wavelength, min(tt_tol, max(0.03, 2.5 * (rms if rms == rms else 0.03)))))
                    # the centring that leaves FEWER lines unindexed wins (a centring that forbids an observed line is wrong, however good its M20 over the other lines looks);
                    # among equals the more constrained one (higher M20) wins
                    if best is None or n_un < best[4] or (n_un == best[4] and m > best[0] * 1.001):
                        best = (m, cen, pr, z, n_un, rms)
                if best is None:
                    continue
                m, cen, pr, z, n_un, rms = best
                if m < 1.5 or n_un > max(3, int(0.3 * n_use)):
                    continue
                sols.append({"system": system, "params": tuple(float(v) for v in pr), "M20": float(m), "n_unindexed": int(n_un), "centring": cen, "zero": float(z), "rms": float(rms)})
        return sols

    base_sys = [s for s in ("cubic", "tetragonal", "hexagonal", "orthorhombic") if s in systems]
    sols_all += run_stage(1, base_sys, 0.02, 0.45)

    def reliable(sols):
        return any(s["M20"] >= 20 and s["n_unindexed"] <= max(2, 0.25 * n_use) for s in sols)

    want_mono = "monoclinic" in systems or bool(thorough)
    if not reliable(sols_all):
        prog(0.5, "no reliable cell yet: thorough search (larger indices)")
        sols_all += run_stage(2, base_sys, 0.5, 0.75)
    if (want_mono or (thorough is None and not reliable(sols_all) and len(tt) >= 10)) and "monoclinic" not in systems and thorough is None:
        pass                                                      # monoclinic only when asked for (it is slower): systems=(..., "monoclinic") or thorough=True
    if want_mono and not reliable(sols_all) and len(tt) >= 10:
        prog(0.78, "monoclinic cells")
        sols_all += run_stage(2, ["monoclinic"], 0.78, 0.97)
    sols = sols_all
    # symmetrise: an orthorhombic cell with two equal edges is tetragonal, with three equal edges cubic; a tetragonal c/a = sqrt(2) is the F/I-centred cubic lattice
    for s in sols:
        if s["system"] == "orthorhombic":
            e = sorted(s["params"])
            if abs(e[0] / e[1] - 1) < 0.004 and abs(e[1] / e[2] - 1) < 0.004:
                s["system"], s["params"], s["note"] = "cubic", (float(np.mean(e)),), "orthorhombic cell with three equal edges: cubic"
            else:
                for a_i, b_i, c_i in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
                    if abs(s["params"][a_i] / s["params"][b_i] - 1) < 0.004:
                        s["system"], s["params"], s["note"] = "tetragonal", (float(np.mean([s["params"][a_i], s["params"][b_i]])), float(s["params"][c_i])), "orthorhombic cell with two equal edges: tetragonal"
                        break
    for s in sols:
        if s["system"] == "hexagonal" and s["centring"] == "R" and abs(s["params"][1] / s["params"][0] - math.sqrt(6)) < 0.012:
            s["system"], s["params"], s["centring"], s["note"] = "cubic", (s["params"][1] / math.sqrt(3),), "F", "cubic F (the rhombohedral setting of the face-centred cubic lattice)"
    for s in sols:
        if s["system"] == "tetragonal":
            r = s["params"][1] / s["params"][0]
            if abs(r - math.sqrt(2)) < 0.004 and s["centring"] == "I":
                s["system"], s["params"], s["centring"], s["note"] = "cubic", (s["params"][1],), "F", "cubic F (the body-centred tetragonal cell with c/a = √2 is the face-centred cubic lattice)"
    sols.sort(key=lambda s: -s["M20"])
    out = []
    for s in sols:
        p = s["params"]
        if s["system"] == "cubic":
            cell = (p[0], p[0], p[0], 90.0, 90.0, 90.0)
        elif s["system"] == "tetragonal":
            cell = (p[0], p[0], p[1], 90.0, 90.0, 90.0)
        elif s["system"] == "hexagonal":
            cell = (p[0], p[0], p[1], 90.0, 90.0, 120.0)
        elif s["system"] == "monoclinic":
            cell = (p[0], p[1], p[2], 90.0, p[3], 90.0)
        else:
            cell = (p[0], p[1], p[2], 90.0, 90.0, 90.0)
        vol = _volume(s["system"], p) if s["system"] not in ("cubic",) or len(p) == 1 else p[0] ** 3
        # the same lattice described by a super- or sub-cell of a solution already kept: keep the one that indexes MORE lines (then the higher M20)
        dup = False
        for oi, o in enumerate(out):
            r = vol / o["volume"]
            near_int = (abs(r - round(r)) < 0.03 * max(r, 1) and round(r) >= 1) or (abs(1 / r - round(1 / r)) < 0.03 and round(1 / r) >= 1)
            if near_int and abs(math.log(r)) < math.log(8.5) and (abs(r - 1) < 0.02 or o["M20"] >= s["M20"] * 0.6):
                if s["n_unindexed"] + 1 < o["n_unindexed"]:        # the new one explains clearly more lines: it replaces the kept one
                    out[oi] = None
                    continue
                dup = True
                break
        out = [o for o in out if o is not None]
        if dup:
            continue
        out.append({"system": s["system"], "cell": tuple(round(v, 5) for v in cell), "centring": s["centring"], "zero_shift": round(s["zero"], 4), "M20": round(s["M20"], 1),
                    "n_unindexed": s["n_unindexed"], "n_used": n_use, "volume": round(vol, 2), "quality": _quality(s["M20"]), "note": s.get("note", ""),
                    "rms_2theta": round(s["rms"], 4)})
    rank = {"cubic": 0, "hexagonal": 1, "tetragonal": 1, "orthorhombic": 2, "monoclinic": 3}
    if out:                                               # among solutions of comparable quality the HIGHEST symmetry first (a cubic cell and its sub-cell descriptions index the same lines)
        top = max(o["M20"] for o in out if o["n_unindexed"] <= max(3, 0.2 * o["n_used"])) if any(o["n_unindexed"] <= max(3, 0.2 * o["n_used"]) for o in out) else max(o["M20"] for o in out)
        best_un = min(o["n_unindexed"] for o in out)
        out.sort(key=lambda s: (s["n_unindexed"] > max(3, 0.2 * s["n_used"]), not (s["M20"] >= 0.3 * top and s["n_unindexed"] <= best_un + 1),
                                rank[s["system"]] if (s["M20"] >= 0.3 * top and s["n_unindexed"] <= best_un + 1) else 0, -s["M20"], s["volume"]))
    prog(1.0, "done")
    return out[:max_solutions]


def _monoclinic_seeds(q_obs, n_use, wavelength, tt_tol, prog=None):
    """Candidate monoclinic cells (unique axis b) from four of the six lowest lines with small Miller indices: Q = A h² + B k² + C l² + D hl solved exactly (batched),
    positive-definite solutions only. The refinement then decides. The number of index combinations is limited."""
    trip = [(h, k, l) for h in range(-1, 3) for k in range(0, 3) for l in range(0, 3) if (h, k, l) != (0, 0, 0) and not (k == 0 and l == 0 and h < 0)]
    T = np.array(trip, float)
    M = np.stack([T[:, 0] ** 2, T[:, 1] ** 2, T[:, 2] ** 2, T[:, 0] * T[:, 2]], 1)
    n = len(M)
    seeds = []
    base = min(6, n_use)
    I1 = np.array(np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij")).reshape(3, -1)
    for ci, ks in enumerate(itertools.combinations(range(base), 4)):
        q4 = q_obs[list(ks)]
        for i0 in range(n):
            A = np.stack([np.broadcast_to(M[i0], (I1.shape[1], 4)), M[I1[0]], M[I1[1]], M[I1[2]]], 1)
            det = np.linalg.det(A)
            ok = np.abs(det) > 1e-4
            if not ok.any():
                continue
            sol = np.linalg.solve(A[ok], np.broadcast_to(q4, (int(ok.sum()), 4))[..., None])[..., 0]
            A_, B_, C_, D_ = sol.T
            good = (A_ > 0) & (B_ > 0) & (C_ > 0) & (np.abs(D_) < 2 * np.sqrt(np.maximum(A_ * C_, 1e-12)) * 0.98)
            for a, b, c, d in sol[good]:
                p = _x_to_params("monoclinic", [a, b, c, d])
                if p is not None and 90.0 <= p[3] <= 130.0:
                    seeds.append(tuple(p))
            if len(seeds) > 6000:
                return seeds
    return seeds