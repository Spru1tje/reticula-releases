r"""v2.18.4: charge flipping (Oszlanyi & Suto 2004) for powder data, with space-group symmetry, normalised amplitudes and iterative re-partitioning of
overlapping reflections (Baerlocher, McCusker & Palatinus 2007).

The method in words. A crystal structure is a periodic electron density rho(x). From a powder pattern only the AMPLITUDES |F(h)| of the Fourier coefficients
are known (from a Le Bail fit); the phases are lost. Charge flipping recovers them by a cycle between the two spaces on a grid:
  1. start with random phases and the observed amplitudes;  2. Fourier transform to the density;  3. FLIP the sign of every density point below a small
  threshold delta (a fraction of the standard deviation of the density) - real atoms are positive peaks on an empty background, so this pushes the density
  towards 'atoms only';  4. transform back; keep the new phases but the observed amplitudes; repeat. The cycle converges to an atomic density when the
  amplitudes are good enough (no atom positions are assumed).
Powder-specific parts: reflections with equal d-spacing cannot be told apart, so the amplitude of a group is first shared equally and, from the cycle's own
density, re-partitioned between the members; the amplitudes are normalised (E-values) so that the density is sharp; the space group (when known) is imposed
on the structure factors in every cycle, which fixes the origin and removes most of the freedom.

All numpy / scipy; one trial takes seconds, several trials with different random phases run in parallel threads."""
from __future__ import annotations

import math
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy import fft as sfft
from scipy import ndimage


def _smooth_n(n):
    """smallest 5-smooth integer >= n (fast FFT sizes)"""
    while True:
        m = n
        for p in (2, 3, 5):
            while m % p == 0:
                m //= p
        if m == 1:
            return n
        n += 1


def _grid_n(nmin, mult):
    n = max(int(nmin), mult)
    while n % mult != 0 or _smooth_n(n) != n:
        n += 1
    return n


def reciprocal_metric(cell):
    from mfs.rietveld import _cell_metric
    return np.linalg.inv(_cell_metric(cell))


def reflections_in_sphere(cell, dmin):
    """All hkl (full sphere, without 000) with d >= dmin; returns hkl (M,3) int, d (M,)"""
    a, b, c = cell[:3]
    Gs = reciprocal_metric(cell)
    hm = [int(math.ceil(v / dmin)) + 1 for v in (a, b, c)]
    h, k, l = np.meshgrid(*[np.arange(-m, m + 1) for m in hm], indexing="ij")
    hkl = np.stack([h.ravel(), k.ravel(), l.ravel()], 1)
    inv_d2 = np.einsum("ki,ij,kj->k", hkl, Gs, hkl)
    ok = (inv_d2 <= 1.0 / dmin ** 2) & (inv_d2 > 1e-12)
    return hkl[ok], 1.0 / np.sqrt(inv_d2[ok])


def normalise_amplitudes(s2, F2, m, window=None):
    """Normalised structure factors E^2 = F^2 / <F^2>(s): the running mean over the d-spacing groups (weighted by multiplicity m) is divided out.
    s2 = sin^2(theta)/lambda^2 of every group, F2 = mean |F|^2 per reflection of the group."""
    o = np.argsort(s2)
    n = len(o)
    w = window or max(7, n // 8)
    F2o, mo = np.asarray(F2, float)[o], np.asarray(m, float)[o]
    mean = np.empty(n)
    for i in range(n):
        lo, hi = max(0, i - w // 2), min(n, i + w // 2 + 1)
        mean[i] = np.sum(F2o[lo:hi] * mo[lo:hi]) / max(np.sum(mo[lo:hi]), 1e-12)
    E2 = np.empty(n)
    E2[o] = F2o / np.maximum(mean, 1e-12 * max(F2o.max(), 1e-12))
    return E2


class Problem:
    """Everything the charge-flipping cycle needs: grid, reflections, amplitudes, symmetry."""

    def __init__(self, cell, hkl, group, G, dmin, sym=None, forbidden=None, oversample=2.4, f000_free=True, max_grid=128):
        """hkl (M,3) int: all reflections of the sphere (both signs). group (M,) int: index of the d-spacing group of each reflection (-1 = not observed).
        G (K,): the total |F|^2 (or E^2 times the number of members) of every group. sym = (R (k,3,3), t (k,3)) point operations or None;
        forbidden (M,) bool: systematically absent reflections (their amplitude is set to zero)."""
        self.cell = tuple(cell)
        self.hkl = np.asarray(hkl, int)
        self.group = np.asarray(group, int)
        self.dmin = float(dmin)
        M = len(self.hkl)
        self.M = M
        hm = np.abs(self.hkl).max(axis=0)
        mult = [1, 1, 1]
        if sym is not None:                                          # the grid must carry the translations of the operations (1/2, 1/3, 1/4, 1/6, 1/8, 1/12, ...)
            for ax in range(3):
                den = 1
                for tt in sym[1][:, ax]:
                    fr = _fraction_denominator(tt)
                    den = den * fr // math.gcd(den, fr)
                mult[ax] = den
        nmin = [int(math.ceil(oversample * max(h, 1))) * 2 // 2 + 1 for h in hm]
        nmin = [max(2 * int(h) + 1, int(oversample * h)) for h in hm]
        n = [_grid_n(nmin[i], mult[i]) for i in range(3)]
        if sym is not None:                                          # axes that are mixed by the rotations need the same n
            R_all = sym[0]
            changed = True
            while changed:
                changed = False
                for ax1 in range(3):
                    for ax2 in range(ax1 + 1, 3):
                        mixed = np.any(np.abs(R_all[:, ax1, ax2]) > 0.5) or np.any(np.abs(R_all[:, ax2, ax1]) > 0.5)
                        if mixed and n[ax1] != n[ax2]:
                            v = max(n[ax1], n[ax2])
                            lcm = mult[ax1] * mult[ax2] // math.gcd(mult[ax1], mult[ax2])
                            v = _grid_n(v, lcm)
                            n[ax1] = n[ax2] = v
                            changed = True
        self.n = tuple(int(v) for v in n)
        if np.prod(self.n) > max_grid ** 3:
            raise ValueError(f"grid {self.n} too large")
        self.idx = tuple((self.hkl[:, i] % self.n[i]) for i in range(3))
        self.flat = np.ravel_multi_index(self.idx, self.n)
        self.forbidden = np.zeros(M, bool) if forbidden is None else np.asarray(forbidden, bool)
        obs = (self.group >= 0) & ~self.forbidden                 # the group total is shared between the ALLOWED members only
        self.obs = obs
        G = np.asarray(G, float)
        self.G = G
        self.members = np.bincount(self.group[obs], minlength=len(G)).astype(float)
        self.A0 = np.zeros(M)
        self.A0[obs] = np.sqrt(np.maximum(G[self.group[obs]], 0.0) / np.maximum(self.members[self.group[obs]], 1.0))
        self.sym = sym
        self._build_symmetry()

    def _build_symmetry(self):
        self.orb = None
        if self.sym is None or len(self.sym[0]) <= 1:
            return
        R, t = self.sym
        key = lambda h: ((h[:, 0] + 512) * 1024 + (h[:, 1] + 512)) * 1024 + (h[:, 2] + 512)
        k0 = key(self.hkl)
        order = np.argsort(k0)
        ks = k0[order]
        pos = np.empty((len(R), self.M), int)
        ph = np.empty((len(R), self.M), complex)
        valid = np.zeros((len(R), self.M), bool)
        for g, (r, tt) in enumerate(zip(R, t)):
            ri = np.rint(np.linalg.inv(r)).astype(int)
            hp = self.hkl @ ri                                       # h' = R^-T h (row form)
            kk = key(hp)
            j = np.searchsorted(ks, kk)
            j = np.clip(j, 0, self.M - 1)
            ok = ks[j] == kk
            pos[g] = np.where(ok, order[j], np.arange(self.M))
            valid[g] = ok
            ph[g] = np.exp(2j * np.pi * (hp @ tt))
        self.orb = (pos, ph, valid)

    def symmetrise(self, f):
        if self.orb is None:
            return f
        pos, ph, valid = self.orb
        acc = np.zeros_like(f)
        cnt = np.zeros(len(f))
        for g in range(len(pos)):
            acc += np.where(valid[g], ph[g] * f[pos[g]], 0.0)
            cnt += valid[g]
        return acc / np.maximum(cnt, 1.0)


def _fraction_denominator(x, maxden=24):
    x = float(x) % 1.0
    for d in (1, 2, 3, 4, 6, 8, 12, 24):
        if abs(x * d - round(x * d)) < 1e-3:
            return d
    return 24


def run_trial(P: Problem, seed=0, n_iter=300, delta=1.1, repartition=True, f000=0.0, progress=None, stop_event=None):
    """One charge-flipping run; returns a dict {F (M,) complex of the final structure factors, R (history), R_final, rho (map from the final F), seed}."""
    rng = np.random.default_rng(seed)
    n = P.n
    M = P.M
    f = P.A0 * np.exp(2j * np.pi * rng.random(M))
    f = P.symmetrise(f)
    ph = np.exp(1j * np.angle(f + 1e-30))
    f = P.A0 * ph
    # Friedel: F(-h) = conj F(h): build the partner map once
    key = lambda h: ((h[:, 0] + 512) * 1024 + (h[:, 1] + 512)) * 1024 + (h[:, 2] + 512)
    k0 = key(P.hkl)
    o = np.argsort(k0)
    j = np.clip(np.searchsorted(k0[o], key(-P.hkl)), 0, M - 1)
    partner = o[j]
    sel = np.nonzero(key(P.hkl) <= key(-P.hkl))[0]               # one of each pair
    A = P.A0.copy()
    grid = np.zeros(n, complex)
    Rh = []
    unobs = ~P.obs & ~P.forbidden
    has_groups = repartition and P.members.max() > 1
    for it in range(n_iter):
        if stop_event is not None and stop_event.is_set():
            break
        grid[:] = 0
        f[sel] = f[sel]
        f[partner[sel]] = np.conj(f[sel])
        grid.reshape(-1)[P.flat] = f
        grid.reshape(-1)[0] = f000
        rho = sfft.ifftn(grid, workers=1).real
        sg = float(rho.std())
        flipped = np.where(rho < delta * sg, -rho, rho)
        Fc = sfft.fftn(flipped, workers=1)
        fc = Fc.reshape(-1)[P.flat]
        fc = P.symmetrise(fc)
        f000 = float(Fc.reshape(-1)[0].real) if it > 5 else f000
        amp = np.abs(fc)
        # R factor of the observed amplitudes (best scale)
        o_ = P.obs & ~P.forbidden
        s = float(np.sum(A[o_] * amp[o_]) / max(np.sum(amp[o_] ** 2), 1e-30))
        R = float(np.sum(np.abs(A[o_] - s * amp[o_])) / max(np.sum(A[o_]), 1e-30))
        Rh.append(R)
        if progress is not None and it % 20 == 0:
            progress(it / n_iter, R)
        # re-partition of overlapping reflections (members of one d-spacing group share the group total)
        if has_groups and it >= 25 and it % 4 == 0:
            a2 = (s * amp) ** 2
            grp = P.group
            tot_c = np.bincount(grp[P.obs], weights=a2[P.obs], minlength=len(P.G))
            share = np.where(tot_c[grp] > 0, a2 / np.maximum(tot_c[grp], 1e-30), 1.0 / np.maximum(P.members[grp], 1.0))
            share = 0.5 * share + 0.5 / np.maximum(P.members[grp], 1.0)
            A_new = np.sqrt(np.maximum(P.G[grp] * share, 0.0))
            A = np.where(P.obs & ~P.forbidden, A_new, 0.0)
            # members that are related by symmetry (or +-h) must share one amplitude: average the squared amplitude over the symmetry orbit
            if P.orb is not None:
                pos, _, valid = P.orb
                a2o = A ** 2
                acc = np.zeros_like(a2o)
                cnt = np.zeros_like(a2o)
                for g in range(len(pos)):
                    acc += np.where(valid[g], a2o[pos[g]], 0.0)
                    cnt += valid[g]
                A = np.sqrt(acc / np.maximum(cnt, 1.0))
            else:
                A = np.sqrt(0.5 * (A ** 2 + A[partner] ** 2))
        # new structure factors: observed amplitudes with the new phases, unobserved ones free
        new = np.where(np.abs(fc) > 1e-30, fc / np.maximum(np.abs(fc), 1e-30), 1.0)
        f = np.where(P.obs & ~P.forbidden, A * new, 0.0)
        f = np.where(unobs, fc, f)
        f[P.forbidden] = 0.0
    grid[:] = 0
    f[partner[sel]] = np.conj(f[sel])
    grid.reshape(-1)[P.flat] = f
    grid.reshape(-1)[0] = f000
    rho = sfft.ifftn(grid, workers=1).real
    tail = Rh[-10:] if len(Rh) >= 10 else Rh
    return {"F": f, "R": Rh, "R_final": float(np.mean(tail)) if tail else 1.0, "rho": rho, "seed": seed, "delta": delta, "A": A}


def run_trials(P: Problem, n_trials=12, n_iter=300, threads=None, seed0=0, deltas=(1.0, 1.3, 0.8), progress=None, stop_event=None):
    """Several charge-flipping runs with different random phases and thresholds, in parallel threads. Returns the results sorted by R (best first)."""
    import os
    threads = threads or max(1, min(n_trials, ((os.cpu_count() or 2) // 2) or 1))     # FFTs release the GIL; hyper-threads add little, so physical cores only
    done = [0]

    def job(k):
        r = run_trial(P, seed=seed0 + k, n_iter=n_iter, delta=deltas[k % len(deltas)], stop_event=stop_event)
        done[0] += 1
        if progress:
            progress(done[0] / n_trials, f"trial {done[0]} of {n_trials}: R = {r['R_final'] * 100:.1f} %")
        return r
    with ThreadPoolExecutor(max_workers=threads) as ex:
        res = list(ex.map(job, range(n_trials)))
    res.sort(key=lambda r: r["R_final"])
    return res


def peaks(rho, cell, n_max=200, min_rel=0.08, min_dist=0.55):
    """Maxima of the periodic density map. Returns [(x, y, z, height)] (fractional), strongest first; peaks closer than min_dist (A) are merged."""
    n = rho.shape
    mx = ndimage.maximum_filter(rho, size=3, mode="wrap")
    cand = np.argwhere((rho == mx) & (rho > min_rel * rho.max()))
    if len(cand) == 0:
        return []
    out = []
    for i, j, k in cand:
        # parabolic refinement per axis
        fr = []
        for ax, (p, npt) in enumerate(zip((i, j, k), n)):
            idx = [(i, j, k)] * 3
            m1 = list((i, j, k)); p1 = list((i, j, k))
            m1[ax] = (p - 1) % npt
            p1[ax] = (p + 1) % npt
            y0, ym, yp = rho[i, j, k], rho[tuple(m1)], rho[tuple(p1)]
            den = ym - 2 * y0 + yp
            d = 0.5 * (ym - yp) / den if abs(den) > 1e-12 else 0.0
            fr.append(((p + max(-0.5, min(0.5, d))) / npt) % 1.0)
        out.append((fr[0], fr[1], fr[2], float(rho[i, j, k])))
    out.sort(key=lambda p: -p[3])
    from mfs.rietveld import _cell_metric
    Gm = _cell_metric(cell)
    keep = []
    for p in out:
        ok = True
        for q in keep:
            d = np.array(p[:3]) - np.array(q[:3])
            d -= np.round(d)
            if float(d @ Gm @ d) ** 0.5 < min_dist:
                ok = False
                break
        if ok:
            keep.append(p)
        if len(keep) >= n_max:
            break
    return keep
