"""v2.16: automatic indexing of a powder pattern with an unknown phase: finds the unit cell from the peak positions.

Systems: cubic, tetragonal, hexagonal, orthorhombic (Miller indices are assigned to the lowest-angle peaks, the cell follows from them, every other peak is tested).
Monoclinic and triclinic are NOT searched (they need far better data and a much larger search - use a dedicated program such as DICVOL / TREOR / McMaille).
Quality: the de Wolff figure of merit M20 (> 10 = a reliable solution, 4-10 = possible, < 4 = unreliable) and the number of unindexed lines.
Physical limits: every cell edge >= 2.4 Å and a volume >= 30 Å³ (smaller 'cells' index anything and give absurd M20). The solution is refined afterwards with
Le Bail (Phase.from_cell): a good cell gives a good Le Bail fit."""
from __future__ import annotations

import itertools
import math

import numpy as np

_CUBIC_N = [n for n in range(1, 40) if any(a * a + b * b + c * c == n for a in range(7) for b in range(7) for c in range(7))]
_H = np.arange(0, 8)
_GH, _GK, _GL = (v.ravel() for v in np.meshgrid(np.arange(-7, 8), np.arange(-7, 8), _H, indexing="ij"))
MIN_EDGE, MIN_VOLUME, MAX_VOLUME = 2.4, 30.0, 15000.0


def _qs(two_theta, wavelength):
    return (2 * np.sin(np.radians(np.asarray(two_theta, float) / 2)) / wavelength) ** 2


def _volume(system, p):
    if system == "cubic":
        return p[0] ** 3
    if system == "tetragonal":
        return p[0] ** 2 * p[1]
    if system == "hexagonal":
        return p[0] ** 2 * p[1] * math.sqrt(3) / 2
    return p[0] * p[1] * p[2]


def _q_of(system, p, H, K, L):
    if system == "cubic":
        return (H ** 2 + K ** 2 + L ** 2) / p[0] ** 2
    if system == "tetragonal":
        return (H ** 2 + K ** 2) / p[0] ** 2 + L ** 2 / p[1] ** 2
    if system == "hexagonal":
        return 4.0 / 3.0 * (H ** 2 + H * K + K ** 2) / p[0] ** 2 + L ** 2 / p[1] ** 2
    return H ** 2 / p[0] ** 2 + K ** 2 / p[1] ** 2 + L ** 2 / p[2] ** 2


def _calc_q(system, p, qmax):
    """All unique q = 1/d² of the cell up to qmax (no extinction rules)."""
    hm = [int(math.ceil(v * math.sqrt(qmax))) + 1 for v in (p if len(p) == 3 else (p[0], p[0], p[-1]))]
    H, K, L = np.meshgrid(np.arange(-hm[0], hm[0] + 1), np.arange(-hm[1], hm[1] + 1), np.arange(0, hm[2] + 1), indexing="ij")
    q = _q_of(system, p, H.ravel(), K.ravel(), L.ravel())
    q = q[(q > 1e-9) & (q <= qmax * 1.02)]
    return np.unique(np.round(q, 7))


def _m_n(q_obs, q_calc):
    """de Wolff M(n) for the observed peaks q_obs (sorted) against the calculated lines."""
    qo = np.sort(q_obs)
    j = np.searchsorted(q_calc, qo)
    left, right = np.clip(j - 1, 0, len(q_calc) - 1), np.clip(j, 0, len(q_calc) - 1)
    dq = np.minimum(np.abs(q_calc[left] - qo), np.abs(q_calc[right] - qo))
    n_calc = int(np.sum(q_calc <= qo[-1] * 1.001))
    return qo[-1] / (2 * max(float(np.mean(dq)), 1e-9) * max(n_calc, 1)), dq


def index(two_theta, wavelength=1.540593, tt_tol=0.06, systems=("cubic", "tetragonal", "hexagonal", "orthorhombic"), max_solutions=6, n_fit=20):
    """Returns the best cells, best first: [{system, cell (a, b, c, alpha, beta, gamma), M20, n_unindexed, n_used, volume, quality}]."""
    tt = np.sort(np.asarray(two_theta, float))
    if len(tt) < 6:
        return []
    q_obs = _qs(tt, wavelength)
    n_use = min(len(tt), n_fit)
    qo = q_obs[:n_use]
    qmax = float(qo[-1])
    n8 = min(8, n_use)

    def qtol(q):
        s = np.sqrt(np.clip(q, 1e-12, None)) * wavelength / 2
        return 2 * s * np.sqrt(np.clip(1 - s * s, 0, 1)) / wavelength ** 2 * math.radians(tt_tol) + 1e-7

    tol8 = qtol(qo[:n8]) * 1.5
    sols = []

    def consider(system, p):
        if any(v < MIN_EDGE for v in p) or not (MIN_VOLUME <= _volume(system, p) <= MAX_VOLUME) or any(v > 90 for v in p):
            return
        # fast pre-test: the lowest peaks must be explained by the small hkl table
        qt = np.unique(np.round(_q_of(system, p, _GH, _GK, _GL), 6))
        qt = qt[qt > 1e-9]
        j = np.searchsorted(qt, qo[:n8])
        left, right = np.clip(j - 1, 0, len(qt) - 1), np.clip(j, 0, len(qt) - 1)
        d8 = np.minimum(np.abs(qt[left] - qo[:n8]), np.abs(qt[right] - qo[:n8]))
        if np.sum(d8 <= tol8) < max(5, n8 - 2):
            return
        qc = _calc_q(system, p, qmax)
        if len(qc) == 0:
            return
        j = np.searchsorted(qc, qo)
        left, right = np.clip(j - 1, 0, len(qc) - 1), np.clip(j, 0, len(qc) - 1)
        dq = np.minimum(np.abs(qc[left] - qo), np.abs(qc[right] - qo))
        n_ok = int(np.sum(dq <= qtol(qo) * 1.5))
        if n_ok < max(6, int(0.75 * n_use)):
            return
        m, _ = _m_n(qo, qc)
        sols.append({"system": system, "params": tuple(float(v) for v in p), "M20": float(m), "n_unindexed": int(n_use - n_ok)})

    if "cubic" in systems:
        for N in _CUBIC_N[:12]:
            for k in range(min(4, n_use)):
                consider("cubic", [math.sqrt(N / q_obs[k])])
    for system in [s for s in ("tetragonal", "hexagonal") if s in systems]:
        f = 1.0 if system == "tetragonal" else 4.0 / 3.0
        ms = sorted({a * a + b * b for a in range(0, 5) for b in range(0, 5)} - {0}) if system == "tetragonal" else \
            sorted({a * a + a * b + b * b for a in range(0, 5) for b in range(0, 5)} - {0})
        for k1, k2 in itertools.combinations(range(min(6, n_use)), 2):
            for m1, l1, m2, l2 in itertools.product(ms[:7], range(0, 5), ms[:7], range(0, 5)):
                if (m1, l1) == (m2, l2):
                    continue
                det = f * m1 * l2 * l2 - f * m2 * l1 * l1
                if abs(det) < 1e-9:
                    continue
                Ac = (q_obs[k1] * l2 * l2 - q_obs[k2] * l1 * l1) / det
                Cc = (f * m1 * q_obs[k2] - f * m2 * q_obs[k1]) / det
                if Ac > 0 and Cc > 0:
                    consider(system, [1 / math.sqrt(Ac), 1 / math.sqrt(Cc)])
    reliable = any(s["M20"] >= 10 for s in sols)
    if "orthorhombic" in systems and len(tt) >= 8 and not reliable:
        trip = [(h, k, l) for h in range(0, 3) for k in range(0, 3) for l in range(0, 3) if (h, k, l) != (0, 0, 0)]
        for ks in itertools.combinations(range(min(6, n_use)), 3):
            qs3 = q_obs[list(ks)]
            for a_, b_, c_ in itertools.product(trip, trip, trip):
                A = np.array([[a_[0] ** 2, a_[1] ** 2, a_[2] ** 2], [b_[0] ** 2, b_[1] ** 2, b_[2] ** 2], [c_[0] ** 2, c_[1] ** 2, c_[2] ** 2]], float)
                if abs(np.linalg.det(A)) < 1e-6:
                    continue
                sol = np.linalg.solve(A, qs3)
                if np.all(sol > 0):
                    consider("orthorhombic", list(1 / np.sqrt(sol)))
    # a tetragonal cell with c/a = sqrt(2) (or 1/sqrt(2)) is the F- / I-centred CUBIC lattice described with a smaller cell: report the cubic one
    for s in sols:
        if s["system"] == "tetragonal":
            r = s["params"][1] / s["params"][0]
            if abs(r - math.sqrt(2)) < 0.004:
                s["system"], s["params"], s["note"] = "cubic", (s["params"][1],), "cubic (F-centred setting of a smaller tetragonal cell)"
            elif abs(r - 1 / math.sqrt(2)) < 0.004:
                s["system"], s["params"], s["note"] = "cubic", (s["params"][0],), "cubic (F-centred setting of a smaller tetragonal cell)"
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
        else:
            cell = (p[0], p[1], p[2], 90.0, 90.0, 90.0)
        vol = _volume(s["system"], p)
        # a super- or sub-cell of a solution already kept (volume ratio close to an integer) is the same lattice: skip it
        dup = False
        for o in out:
            r = vol / o["volume"]
            if abs(r - round(r)) < 0.03 * max(r, 1) and round(r) >= 1 or abs(1 / r - round(1 / r)) < 0.03 and round(1 / r) >= 1:
                if abs(math.log(r)) < math.log(8.5):
                    if (abs(r - 1) < 0.02) or o["M20"] >= s["M20"] * 0.6:
                        dup = True
                        break
        if dup:
            continue
        out.append({"system": s["system"], "cell": tuple(round(v, 5) for v in cell), "M20": round(s["M20"], 1), "n_unindexed": s["n_unindexed"], "n_used": n_use,
                    "volume": round(vol, 2), "quality": "reliable" if s["M20"] >= 10 else "possible" if s["M20"] >= 4 else "unreliable", "note": s.get("note", "")})
    rank = {"cubic": 0, "hexagonal": 1, "tetragonal": 1, "orthorhombic": 2}
    out.sort(key=lambda s: (-(round(s["M20"] / 5)), s["n_unindexed"], rank[s["system"]], s["volume"]))
    return out[:max_solutions]
