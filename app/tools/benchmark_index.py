"""Developer tool: benchmark of the powder indexing on real crystal structures (COD, tests/data/cod_bench). For every structure a realistic peak list is made
(Cu Kα1, peak-position noise, a zero shift, 12 % of the weak lines missing, one impurity line) and the indexing is asked for the cell.
Usage: python tools/benchmark_index.py [--sigma 0.015] [--seed 1] [--only name,name]"""
import argparse
import glob
import os
import sys
import time
import warnings

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
warnings.filterwarnings("ignore")
from mfs import indexing as IX  # noqa: E402

LAM = 1.540593


def true_cell(cif):
    from pymatgen.core import Structure
    s = Structure.from_file(cif)
    L = s.lattice
    return s, (L.a, L.b, L.c, L.alpha, L.beta, L.gamma)


def make_peaks(structure, rng, sigma=0.015, zero=None, miss=0.12, n_max=25, impurity=True):
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    p = XRDCalculator(wavelength=LAM).get_pattern(structure, two_theta_range=(5, 90))
    tt, I = np.array(p.x), np.array(p.y)
    order = np.argsort(tt)
    tt, I = tt[order], I[order]
    # merge lines closer than 0.12 deg (they are not resolved in a laboratory pattern)
    keep = [0]
    for k in range(1, len(tt)):
        if tt[k] - tt[keep[-1]] < 0.12:
            if I[k] > I[keep[-1]]:
                keep[-1] = k
        else:
            keep.append(k)
    tt, I = tt[keep], I[keep]
    sel = (I > 0.8) & (rng.random(len(tt)) > np.where(I < 5, miss * 2, miss * 0.3))
    tt = tt[sel]
    z = rng.uniform(-0.03, 0.03) if zero is None else zero
    tt = tt + z + rng.normal(0, sigma, len(tt))
    tt = np.sort(tt)[:n_max]
    if impurity and len(tt) > 8:
        tt = np.sort(np.append(tt, tt[4] + rng.uniform(0.6, 1.4)))[: n_max + 1]
    return tt, z


def volume(cell):
    a, b, c, al, be, ga = cell
    ca, cb, cg = (np.cos(np.radians(v)) for v in (al, be, ga))
    return a * b * c * np.sqrt(max(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg, 0))


def same_lattice(found, truth, tol=0.012):
    """volume within tol and, when the angles are right-angled, the sorted edges within tol. Returns 'exact', 'cell-multiple' or ''."""
    vf, vt = volume(found), volume(truth)
    if abs(vf / vt - 1) < tol * 2:
        if all(abs(x - 90) < 0.5 for x in found[3:]) and all(abs(x - 90) < 0.5 for x in truth[3:]):
            ef, et = sorted(found[:3]), sorted(truth[:3])
            return "exact" if all(abs(a / b - 1) < tol for a, b in zip(ef, et)) else ""
        return "exact"
    r = vf / vt
    if (abs(r - round(r)) < 0.02 * round(r) and 1 < round(r) <= 6) or (abs(1 / r - round(1 / r)) < 0.02 and 1 < round(1 / r) <= 6):
        return "cell-multiple"
    return ""


def run(index_fn, sigma=0.015, seed=1, only=None, verbose=True, **kw):
    rng = np.random.default_rng(seed)
    rows = []
    for cif in sorted(glob.glob(os.path.join(HERE, "tests", "data", "cod_bench", "*.cif"))):
        name = os.path.splitext(os.path.basename(cif))[0]
        if only and name not in only:
            continue
        try:
            s, cell = true_cell(cif)
        except Exception:
            continue
        tt, z = make_peaks(s, rng, sigma)
        if len(tt) < 8:
            rows.append((name, "few peaks", None, None, 0.0)); continue
        t0 = time.time()
        try:
            sols = index_fn(tt, LAM, **kw)
        except Exception as e:
            sols = []
            print(name, "ERROR", repr(e)[:100])
        dt = time.time() - t0
        found, rank, m20, kind = None, None, None, ""
        for r, sol in enumerate(sols[:6]):
            k = same_lattice(sol["cell"], cell)
            if k == "exact":
                found, rank, m20, kind = sol, r + 1, sol["M20"], k; break
            if k and not found:
                found, rank, m20, kind = sol, r + 1, sol["M20"], k
        best = sols[0]["M20"] if sols else None
        rows.append((name, kind or "MISSED", rank, m20, dt, best, len(tt), [round(v, 3) for v in cell[:3]]))
        if verbose:
            print(f"{name:18s} {kind or 'MISSED':14s} rank={rank} M20(true)={m20} best M20={best}  n={len(tt)}  {dt:5.1f}s  cell={[round(v, 2) for v in cell[:3]]}")
    ok = sum(1 for r in rows if r[1] == "exact")
    mult = sum(1 for r in rows if r[1] == "cell-multiple")
    top1 = sum(1 for r in rows if r[1] in ("exact", "cell-multiple") and r[2] == 1)
    m = [r[3] for r in rows if r[3] is not None]
    print(f"\nexact {ok}/{len(rows)}, cell multiple {mult}, first-ranked {top1}; mean M20 of the true cell {np.mean(m) if m else float('nan'):.1f}, "
          f"median {np.median(m) if m else float('nan'):.1f}; time {sum(r[4] for r in rows):.0f} s")
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigma", type=float, default=0.015)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    run(IX.index, a.sigma, a.seed, set(a.only.split(",")) if a.only else None)
