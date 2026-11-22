"""v2.15: robust band ratios for luminescence thermometry - a GLOBAL fit of the whole temperature series instead of integrating windows.

Why: with overlapping, wavy, dirty or doubled bands an integration window cannot separate the bands (the overlap is cut at an arbitrary
valley, and the result changes with the window). Here every spectrum of the series is described by the SAME bands:
  * every band is a pseudo-Voigt with unit area; its position and width may drift LINEARLY with temperature (shared by all spectra);
  * the areas and the (polynomial) background of every spectrum are linear parameters, solved exactly in every step (variable projection);
  * extra 'nuisance' bands (a shoulder, a third peak) can be added so that they do not contaminate the two bands of the ratio;
  * the simplest model (shared shape) and the drifting model are compared with the BIC and the better one is used.
The ratio is area1 / area2 with a standard error that includes the covariance of the two areas; a robustness check refits with other windows
and the other model and adds the spread as a systematic error. Quality flags per spectrum (R², reduced χ², residual autocorrelation)
tell the user where not to trust a point. Nothing here depends on the window the user would have drawn."""
from __future__ import annotations

import contextlib
import functools
import math
import os
import threading

import numpy as np
from scipy.optimize import least_squares, nnls
from scipy.signal import find_peaks

LN2 = math.log(2.0)


_BLAS_LOCK = threading.Lock()
_BLAS_DEPTH = [0]


def single_blas(fn):
    """v2.19.5: the fits work on many SMALL matrices (QR of 400 x 10, SVD of the Jacobian); a multi-threaded BLAS made each call up to 5x slower
    (thread hand-over costs more than the work). The outermost fit call limits BLAS to one thread (threadpoolctl) and restores it afterwards;
    without threadpoolctl nothing changes."""
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        with _BLAS_LOCK:
            outer = _BLAS_DEPTH[0] == 0
            _BLAS_DEPTH[0] += 1
        try:
            ctx = contextlib.nullcontext()
            if outer and os.environ.get("MFS_PL_BLAS_THREADS") != "keep":
                try:
                    from threadpoolctl import threadpool_limits
                    ctx = threadpool_limits(limits=1, user_api="blas")
                except Exception:
                    ctx = contextlib.nullcontext()
            with ctx:
                return fn(*args, **kwargs)
        finally:
            with _BLAS_LOCK:
                _BLAS_DEPTH[0] -= 1
    return wrapper


_POOL = None
_POOL_LOCK = threading.Lock()
_IN_WORKER = threading.local()


def _parallel_ok():
    return os.environ.get("MFS_NO_PARALLEL") != "1" and (os.cpu_count() or 1) > 1 and not getattr(_IN_WORKER, "on", False)


def _pool():
    global _POOL
    with _POOL_LOCK:
        if _POOL is None:
            from concurrent.futures import ThreadPoolExecutor
            _POOL = ThreadPoolExecutor(max_workers=max(2, min(8, os.cpu_count() or 2)), thread_name_prefix="pl-fit")
        return _POOL


def _call(fn, args, kwargs):
    _IN_WORKER.on = True
    try:
        return True, fn(*args, **kwargs)
    except Exception as e:                                               # the caller decides (robust_ratio skips a failed alternative)
        return False, e
    finally:
        _IN_WORKER.on = False


def run_independent(calls):
    """v2.19.5: independent fits [(fn, args, kwargs), ...] at the same time on several CPU cores (threads: NumPy / LAPACK release the GIL).
    Every fit is independent, so the results are identical to running them one after the other. Returns [(ok, result_or_exception)].
    MFS_NO_PARALLEL=1 (or one core, or a call from inside a worker) runs them in series."""
    if len(calls) < 2 or not _parallel_ok():
        out = []
        for fn, a, k in calls:
            try:
                out.append((True, fn(*a, **k)))
            except Exception as e:
                out.append((False, e))
        return out
    pool = _pool()
    futs = [pool.submit(_call, fn, a, k) for fn, a, k in calls]
    return [f.result() for f in futs]


def _warm():
    """MFS_PL_WARM=0 switches the v2.19.5 warm starts off (comparisons)."""
    return os.environ.get("MFS_PL_WARM") != "0" and not _fd_mode()


def _start_for(prev_params, prev_bands, new_bands):
    """Parameter values of an earlier fit for the bands that are also in a new band list (matched by their start guess), renamed to the new
    band indices - a warm start for the next fit. Bands that are new keep their own start values."""
    if not prev_params or prev_bands is None:
        return None
    key = lambda b: (round(float(b[0]), 6), round(float(b[1]), 6))     # noqa: E731
    old = {key(b): k for k, b in enumerate(prev_bands)}
    out = {}
    for j, b in enumerate(new_bands):
        k = old.get(key(b))
        if k is None:
            continue
        for pfx in "cswge":
            v = prev_params.get(f"{pfx}{k}")
            if v is not None:
                out[f"{pfx}{j}"] = float(v)
    return out or None


def _fd_mode():
    """MFS_PL_FD_JAC=1: the pre-2.19.5 finite-difference Jacobian and ratio gradient (for comparisons and tests)."""
    return os.environ.get("MFS_PL_FD_JAC") == "1"


def _noise(y):
    """Noise level of a spectrum from the point-to-point differences (robust: a smooth band hardly contributes)."""
    d = np.diff(np.asarray(y, float))
    if len(d) < 3:
        return 1.0
    return float(np.median(np.abs(d - np.median(d))) / 0.6745 / math.sqrt(2.0)) or float(np.std(d) / math.sqrt(2.0)) or 1.0


def pv_unit(x, c, w, eta):
    """Unit-area pseudo-Voigt (per nm) with FWHM w and Lorentzian share eta."""
    z = (x - c) / w
    g = (2.0 * math.sqrt(LN2 / math.pi) / w) * np.exp(-4.0 * LN2 * z * z)
    l = (2.0 / (math.pi * w)) / (1.0 + 4.0 * z * z)
    return eta * l + (1.0 - eta) * g


def pv_matrix(x, c, w, eta):
    """Same as pv_unit for several bands at once: returns an array (len(x), number of bands)."""
    x = np.asarray(x, float)[:, None]
    c = np.asarray(c, float)[None, :]
    w = np.asarray(w, float)[None, :]
    eta = np.asarray(eta, float)[None, :]
    z = (x - c) / w
    g = (2.0 * math.sqrt(LN2 / math.pi) / w) * np.exp(-4.0 * LN2 * z * z)
    l = (2.0 / (math.pi * w)) / (1.0 + 4.0 * z * z)
    return eta * l + (1.0 - eta) * g


# ------------------------------------------------------------------------------------------------------------- start values
def auto_bands(x, y, lo=None, hi=None, n_extra=0):
    """The two strongest separate bands (+ optionally further ones) of a (mean) spectrum: [(centre, fwhm_guess), ...] sorted by centre.
    The first two are the ratio pair (I1 = shorter wavelength). Returns None when fewer than two bands are found."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    lo = x.min() if lo is None else lo
    hi = x.max() if hi is None else hi
    m = (x >= lo) & (x <= hi)
    xs, ys = x[m], y[m]
    if len(xs) < 15:
        return None
    from scipy.signal import savgol_filter
    step = float(np.median(np.diff(xs)))
    win = max(5, int(round(2.0 / step)) | 1)
    ysm = savgol_filter(ys, min(win, len(ys) - (1 - len(ys) % 2)), 2) if len(ys) > win else ys
    rng = float(ysm.max() - ysm.min())
    if rng <= 0:
        return None
    # candidates: maxima of the spectrum AND of the negative second derivative (that reveals shoulders of overlapping bands)
    w2 = max(7, int(round(0.06 * len(ys))) | 1)
    d2 = savgol_filter(ys, min(w2, len(ys) - (1 - len(ys) % 2)), 3, deriv=2)
    cand = list(find_peaks(ysm, prominence=0.03 * rng)[0])
    pk2, pr2 = find_peaks(-d2, prominence=0.04 * float(np.max(-d2) - np.min(-d2)))
    cand += [int(i) for i in pk2]
    cand = sorted(set(int(i) for i in cand))
    merged = []
    for i in cand:                                                # merge candidates closer than ~1.5 nm
        if merged and xs[i] - xs[merged[-1]] < max(1.5, 3 * step):
            if ysm[i] > ysm[merged[-1]]:
                merged[-1] = i
        else:
            merged.append(i)
    base = float(np.percentile(ysm, 5))
    merged = [i for i in merged if ysm[i] - base > 0.06 * rng]
    if len(merged) < 2:
        return None
    strongest = sorted(merged, key=lambda i: ysm[i], reverse=True)[:2 + max(0, n_extra)]
    chosen = sorted(strongest)
    out = []
    for p_ in chosen:
        h = ysm[p_]
        half = h - (h - base) / 2.0
        i0 = p_
        while i0 > 0 and ysm[i0] > half:
            i0 -= 1
        i1 = p_
        while i1 < len(ysm) - 1 and ysm[i1] > half:
            i1 += 1
        w = max(float(xs[i1] - xs[i0]), 3 * step)
        out.append((float(xs[p_]), w))
    # a shoulder has no half-height crossing of its own: never guess a width larger than the distance to the nearest other band x 1.2
    for k, (c, w) in enumerate(out):
        gap = min([abs(c - o) for j, (o, _) in enumerate(out) if j != k] or [w])
        out[k] = (c, float(min(w, 1.2 * gap)))
    return out


def default_window(bands, x):
    """Fit window: the bands plus 3 FWHM on each side (at least as wide as the data allow)."""
    lo = min(c - 3.0 * w for c, w in bands)
    hi = max(c + 3.0 * w for c, w in bands)
    return float(max(lo, np.min(x))), float(min(hi, np.max(x)))


# ------------------------------------------------------------------------------------------------------------------ the fit
class _Model:
    def __init__(self, spectra, bands, window, bg_order, drift_pos, drift_width, counts):
        self.nb = len(bands)
        self.bg_order = int(bg_order)
        self.drift_pos, self.drift_width = bool(drift_pos), bool(drift_width)
        self.T = np.array([s["T"] for s in spectra], float)
        self.T0 = float(np.mean(self.T))
        self.Ts = max(float(np.ptp(self.T)), 1.0)                      # temperatures enter as (T - T0) / range
        lo, hi = window
        self.data = []
        self.var = []
        for s in spectra:
            x = np.asarray(s["x"], float)
            y = np.asarray(s["y"], float)
            m = (x >= lo) & (x <= hi) & np.isfinite(y)
            self.data.append((x[m], y[m]))
            # v2.18.4: after a background subtraction the Poisson variance is that of the RAW counts (signal + background): spectra may carry it as "var"
            self.var.append(np.asarray(s["var"], float)[m] if s.get("var") is not None else y[m])
        self.counts = counts
        self.w = [1.0 / np.maximum(v, 1.0) if counts else np.ones_like(y) / max(_noise(y), 1e-12) ** 2 for (x, y), v in zip(self.data, self.var)]
        self.guess = list(bands)
        # parameter layout per band: c0, [s], w0, [g], eta
        self.names, p0, lb, ub = [], [], [], []
        allc = sorted(c for c, _ in bands)
        for k, (c, w) in enumerate(bands):
            gap = min([abs(c - o) for o in allc if o != c] or [w * 4])
            dc = max(min(gap * 0.45, 3.0 * w), 0.5 * w)
            self.names += [f"c{k}"]; p0 += [c]; lb += [c - dc]; ub += [c + dc]
            if drift_pos:
                self.names += [f"s{k}"]; p0 += [0.0]; lb += [-0.5 * w]; ub += [0.5 * w]          # nm over the whole temperature range
            self.names += [f"w{k}"]; p0 += [w]; lb += [0.25 * w]; ub += [4.0 * w]
            if drift_width:
                self.names += [f"g{k}"]; p0 += [0.0]; lb += [-0.6 * w]; ub += [0.6 * w]
            self.names += [f"e{k}"]; p0 += [0.5]; lb += [0.0]; ub += [1.0]
        self.p0, self.lb, self.ub = np.array(p0), np.array(lb), np.array(ub)

        # v2.19.5: everything the evaluation needs is computed once (it was rebuilt for every spectrum in every function call)
        self.n_par = len(self.p0)
        self.tf = (self.T - self.T0) / self.Ts
        self.sw = [np.sqrt(w) for w in self.w]
        self.yw = [sw * y for sw, (_, y) in zip(self.sw, self.data)]
        self.Bs = []
        for x, _ in self.data:
            xn = (x - x.mean()) / max(np.ptp(x), 1e-9) * 2 if len(x) else x
            self.Bs.append(np.stack([np.polynomial.chebyshev.chebval(xn, [0] * j + [1]) for j in range(self.bg_order + 1)], 1))
        ic, is_, iw, ig, ie = [], [], [], [], []
        i = 0
        for _k in range(self.nb):                                       # same layout as self.names: c, [s], w, [g], eta per band
            ic.append(i); i += 1
            is_.append(i if self.drift_pos else -1); i += 1 if self.drift_pos else 0
            iw.append(i); i += 1
            ig.append(i if self.drift_width else -1); i += 1 if self.drift_width else 0
            ie.append(i); i += 1
        self._ic, self._is, self._iw, self._ig, self._ie = (np.array(v, int) for v in (ic, is_, iw, ig, ie))
        self._pband = np.zeros(self.n_par, int)                         # band of every parameter
        for idx in (self._ic, self._is, self._iw, self._ig, self._ie):
            for k, p in enumerate(idx):
                if p >= 0:
                    self._pband[p] = k
        self._cache = {}                                                 # th -> per-spectrum evaluation (least_squares asks f and J at the same th)
        # v2.19.5: one common wavelength grid (the usual case: one instrument) -> every spectrum is handled in ONE batched array operation
        x0 = self.data[0][0] if self.data else np.zeros(0)
        self.shared = bool(self.data) and len(x0) > 0 and all(len(x) == len(x0) and np.array_equal(x, x0) for x, _ in self.data)
        if self.shared:
            self.X = x0
            self.SW = np.stack(self.sw)
            self.YW = np.stack(self.yw)
            self.Y = np.stack([y for _, y in self.data])
            self.B0 = self.Bs[0]

    # ---------------------------------------------------------------------------------------------------- band shapes
    def _bp(self, th, t):
        c = th[self._ic] + (th[self._is] * t if self.drift_pos else 0.0)
        wr = th[self._iw] + (th[self._ig] * t if self.drift_width else 0.0)
        return c, np.maximum(wr, 1e-3), th[self._ie], wr > 1e-3

    def band_params(self, th, ti):
        """(centres, widths, etas) of all bands at the temperature index ti."""
        c, w, e, _ = self._bp(np.asarray(th, float), self.tf[ti])
        return [float(v) for v in c], [float(v) for v in w], [float(v) for v in e]

    def design(self, th, ti):
        c, w, e, _ = self._bp(np.asarray(th, float), self.tf[ti])
        return pv_matrix(self.data[ti][0], c, w, e), self.Bs[ti]

    # ---------------------------------------------------------------------------------------------------- linear step
    def _linear(self, A, ti, try_free=True):
        """Areas and background of one spectrum for fixed shapes: the non-negative least-squares solution.
        v2.19.5: solved by QR; NNLS only when an area would be negative (when the free solution has all areas >= 0 it IS the NNLS solution).
        Returns (amp, bg, Q, R, active) with Q, R the QR factors of the weighted design of the active bands + background."""
        B, sw, yw, nb = self.Bs[ti], self.sw[ti], self.yw[ti], self.nb
        M = np.hstack([A, B]) * sw[:, None]
        try:
            if not try_free:
                raise ValueError("free solution already known to fail")
            Q, R = np.linalg.qr(M)
            d = np.abs(np.diag(R))
            if d.size and d.min() > 1e-11 * d.max():
                coef = np.linalg.solve(R, Q.T @ yw)
                if np.all(np.isfinite(coef)) and np.all(coef[:nb] >= 0.0):
                    return coef[:nb], coef[nb:], Q, R, np.arange(nb)
        except (np.linalg.LinAlgError, ValueError):
            pass
        y = self.data[ti][1]
        scale = max(float(np.abs(y).max()), 1e-12)
        c, _ = nnls(np.hstack([A, B, -B]) * sw[:, None], y / scale * sw, maxiter=5000)
        amp = c[:nb] * scale
        bg = (c[nb:nb + B.shape[1]] - c[nb + B.shape[1]:]) * scale
        active = np.flatnonzero(amp > 0)
        Q, R = np.linalg.qr(np.hstack([A[:, active], B]) * sw[:, None])
        return amp, bg, Q, R, active

    def _bp_all(self, th):
        """Band parameters of all spectra at once: arrays (n_spectra, n_bands)."""
        t = self.tf[:, None]
        c = th[self._ic][None, :] + (th[self._is][None, :] * t if self.drift_pos else 0.0)
        wr = th[self._iw][None, :] + (th[self._ig][None, :] * t if self.drift_width else 0.0) + np.zeros_like(t)
        c = c + np.zeros_like(t)
        return c, np.maximum(wr, 1e-3), np.broadcast_to(th[self._ie][None, :], c.shape), wr > 1e-3

    def _eval_batch(self, th):
        """All spectra in one batched QR (common grid). Spectra whose free solution has a negative area, or a singular design, are redone one by one
        with NNLS. Returns the per-spectrum list (as _eval) and the batched arrays for the vectorised derivatives (or None)."""
        C, W, E, _ = self._bp_all(th)
        z = (self.X[None, :, None] - C[:, None, :]) / W[:, None, :]
        G = (2.0 * math.sqrt(LN2 / math.pi) / W)[:, None, :] * np.exp(-4.0 * LN2 * z * z)
        L = (2.0 / (math.pi * W))[:, None, :] / (1.0 + 4.0 * z * z)
        P = E[:, None, :] * L + (1.0 - E[:, None, :]) * G                       # (S, n, nb)
        S = P.shape[0]
        Bb = np.broadcast_to(self.B0[None, :, :], (S,) + self.B0.shape)
        M = np.concatenate([P, Bb], axis=2) * self.SW[:, :, None]
        nb = self.nb
        Q, R = np.linalg.qr(M)
        d = np.abs(np.diagonal(R, axis1=1, axis2=2))
        good = (d.min(axis=1) > 1e-11 * d.max(axis=1))
        coef = np.full((S, M.shape[2]), np.nan)
        if good.any():
            try:
                coef[good] = np.linalg.solve(R[good], np.einsum("snk,sn->sk", Q[good], self.YW[good])[:, :, None])[:, :, 0]
            except np.linalg.LinAlgError:
                good[:] = False
        good &= np.all(np.isfinite(coef), axis=1) & np.all(coef[:, :nb] >= 0.0, axis=1)
        ev = []
        all_active = np.arange(nb)
        cf = np.where(good[:, None], coef, 0.0)
        Rres = self.SW * (self.Y - np.einsum("snk,sk->sn", P, cf[:, :nb]) - cf[:, nb:] @ self.B0.T)
        for ti in range(S):
            A = P[ti]
            if good[ti]:
                ev.append((coef[ti, :nb], coef[ti, nb:], A, Q[ti], R[ti], all_active, Rres[ti]))
            else:
                amp, bg, Qi, Ri, active = self._linear(A, ti, try_free=False)
                r = self.sw[ti] * (self.data[ti][1] - A @ amp - self.B0 @ bg)
                ev.append((amp, bg, A, Qi, Ri, active, r))
        # the batched pieces are kept also when some spectra needed NNLS: the Jacobian then still builds every derivative in one go
        batch = {"z": z, "G": G, "L": L, "P": P, "W": W, "E": E, "Q": Q, "R": R, "coef": coef, "good": good,
                 "r": Rres if bool(good.all()) else np.stack([e[6] for e in ev])}
        return ev, batch

    def _eval(self, th, with_batch=False):
        th = np.asarray(th, float)
        key = th.tobytes()
        hit = self._cache.get(key)
        if hit is None:
            if self.shared:
                ev, batch = self._eval_batch(th)
            else:
                ev, batch = [], None
                for ti in range(len(self.data)):
                    A, B = self.design(th, ti)
                    amp, bg, Q, R, active = self._linear(A, ti)
                    r = self.sw[ti] * (self.data[ti][1] - A @ amp - B @ bg)
                    ev.append((amp, bg, A, Q, R, active, r))
            if len(self._cache) > 6:
                self._cache.clear()
            hit = self._cache[key] = (ev, batch)
        return hit if with_batch else hit[0]

    def _batch_derivs(self, th, batch, amp=None):
        """Vectorised version of _derivs for all spectra (common grid, every band active): D, DA (S, n, n_par)."""
        z, G, L, P, W, E = batch["z"], batch["G"], batch["L"], batch["P"], batch["W"], batch["E"]
        _, _, _, free = self._bp_all(th)
        den = 1.0 + 4.0 * z * z
        Ee = E[:, None, :]
        dPz = Ee * (L * (-8.0 * z / den)) + (1.0 - Ee) * (G * (-8.0 * LN2 * z))
        Wn = W[:, None, :]
        dPc = -dPz / Wn
        dPw = -(P + z * dPz) / Wn * free[:, None, :]
        dPe = L - G
        S, n, _ = P.shape
        sw = self.SW[:, :, None]
        D = np.empty((S, n, self.n_par))
        D[:, :, self._ic] = dPc * sw
        D[:, :, self._iw] = dPw * sw
        D[:, :, self._ie] = dPe * sw
        t = self.tf[:, None, None]
        if self.drift_pos:
            D[:, :, self._is] = dPc * t * sw
        if self.drift_width:
            D[:, :, self._ig] = dPw * t * sw
        amp = batch["coef"][:, :self.nb] if amp is None else amp
        DA = D * amp[:, self._pband][:, None, :]
        return D, DA

    def solve(self, th, ti):
        amp, bg, A, _, _, _, _ = self._eval(th)[ti]
        return amp, bg, A, self.Bs[ti]

    def residual(self, th):
        return np.concatenate([e[6] for e in self._eval(th)])

    # ---------------------------------------------------------------------------------------------------- derivatives
    def _derivs(self, th, ti, ev):
        """Pieces of the exact derivative of the projected residual of spectrum ti (variable projection, Golub & Pereyra 1973):
        D[:, p] = weighted d(column of band k)/d(theta_p), DA = D * a_k, and the position of band k among the active columns."""
        t = self.tf[ti]
        x, sw = self.data[ti][0], self.sw[ti]
        amp, _, _, Q, R, active, r = ev
        c, w, e, free = self._bp(np.asarray(th, float), t)
        z = (x[:, None] - c[None, :]) / w[None, :]
        z2 = z * z
        G = (2.0 * math.sqrt(LN2 / math.pi) / w)[None, :] * np.exp(-4.0 * LN2 * z2)
        den = 1.0 + 4.0 * z2
        L = (2.0 / (math.pi * w))[None, :] / den
        P = e[None, :] * L + (1.0 - e[None, :]) * G
        dPz = e[None, :] * (L * (-8.0 * z / den)) + (1.0 - e[None, :]) * (G * (-8.0 * LN2 * z))
        dPc = -dPz / w[None, :]
        dPw = -(P + z * dPz) / w[None, :] * free[None, :]
        dPe = L - G
        D = np.zeros((len(x), self.n_par))
        jpos = np.full(self.n_par, -1)
        ak = np.zeros(self.n_par)
        for j, k in enumerate(active):
            cols = [(self._ic[k], dPc[:, k]), (self._iw[k], dPw[:, k]), (self._ie[k], dPe[:, k])]
            if self.drift_pos:
                cols.append((self._is[k], t * dPc[:, k]))
            if self.drift_width:
                cols.append((self._ig[k], t * dPw[:, k]))
            for p, v in cols:
                D[:, p] = v * sw
                jpos[p] = j
                ak[p] = amp[k]
        return D, D * ak[None, :], jpos, Q, R, r

    @staticmethod
    def _rinv(R):
        try:
            Ri = np.linalg.inv(R)
            if np.all(np.isfinite(Ri)):
                return Ri
        except (np.linalg.LinAlgError, ValueError):
            pass
        return np.linalg.pinv(R)

    def jacobian(self, th):
        """Exact Jacobian of the projected residual: J_p = -(P_perp dM/dtheta_p a + (M+)^T (dM/dtheta_p)^T r). One evaluation instead of n_par + 1."""
        th = np.asarray(th, float)
        ev, batch = self._eval(th, with_batch=True)
        if batch is not None and bool(batch["good"].all()):              # all spectra at once
            D, DA = self._batch_derivs(th, batch)
            Q, R, r = batch["Q"], batch["R"], batch["r"]
            Qt = np.swapaxes(Q, 1, 2)
            J = DA - Q @ (Qt @ DA)
            U = Q @ np.swapaxes(self._rinv_batch(R), 1, 2)
            J += U[:, :, self._pband] * np.einsum("snp,sn->sp", D, r)[:, None, :]
            return -J.reshape(-1, self.n_par)
        if batch is not None:                                            # mixed: derivatives in one go, projection per group
            good = batch["good"]
            AMP = np.stack([e[0] for e in ev])
            D, DA = self._batch_derivs(th, batch, AMP)
            Jb = np.empty_like(D)
            g = np.flatnonzero(good)
            if len(g):
                Q, R, r = batch["Q"][g], batch["R"][g], batch["r"][g]
                Jg = DA[g] - Q @ (np.swapaxes(Q, 1, 2) @ DA[g])
                U = Q @ np.swapaxes(self._rinv_batch(R), 1, 2)
                Jg += U[:, :, self._pband] * np.einsum("snp,sn->sp", D[g], r)[:, None, :]
                Jb[g] = -Jg
            for ti in np.flatnonzero(~good):
                _, _, _, Qi, Ri, active, ri = ev[ti]
                pos = {int(k): j for j, k in enumerate(active)}
                jpos = np.array([pos.get(int(k), -1) for k in self._pband])
                m = jpos >= 0
                Di = D[ti] * m[None, :]
                Ji = DA[ti] - Qi @ (Qi.T @ DA[ti])
                if m.any() and Ri.size:
                    Ji[:, m] += (Qi @ self._rinv(Ri).T)[:, jpos[m]] * (Di[:, m].T @ ri)[None, :]
                Jb[ti] = -Ji
            blocks = list(Jb)
        else:
            blocks = []
        for ti in (range(len(self.data)) if batch is None else ()):
            D, DA, jpos, Q, R, r = self._derivs(th, ti, ev[ti])
            J = DA - Q @ (Q.T @ DA)
            m = jpos >= 0
            if m.any() and R.size:
                U = Q @ self._rinv(R).T
                J[:, m] += U[:, jpos[m]] * (D[:, m].T @ r)[None, :]
            blocks.append(-J)
        J = np.vstack(blocks)
        # A line that NNLS holds at zero area in a spectrum has an exactly zero shape derivative there, so the optimiser could never move it
        # back into place (tried: difference columns only for lines that are off in EVERY spectrum - the fits then need more steps, slower). For the parameters of such lines the column is the forward difference of the pre-2.19.5 code (same step, same bound rule):
        # it sees that a small shift switches the line on again. One batched residual per parameter (cheaper than spectrum by spectrum when the
        # line is off in many spectra); usually 0-2 lines, so the saving stays.
        inactive = sorted({int(k) for e in ev for k in range(self.nb) if k not in set(int(a) for a in e[5])})
        if inactive:
            r0 = np.concatenate([e[6] for e in ev])
            for p in np.flatnonzero(np.isin(self._pband, inactive)):
                h = 1e-4 * (1.0 if th[p] >= 0 else -1.0) * max(1.0, abs(th[p]))
                if th[p] + h > self.ub[p] or th[p] + h < self.lb[p]:
                    h = -h
                tp = th.copy()
                tp[p] += h
                J[:, p] = (self.residual(tp) - r0) / h
        return J

    @staticmethod
    def _rinv_batch(R):
        try:
            Ri = np.linalg.inv(R)
            if np.all(np.isfinite(Ri)):
                return Ri
        except (np.linalg.LinAlgError, ValueError):
            pass
        return np.linalg.pinv(R)

    def ratio_grads(self, th, g1, g2):
        """d(ratio)/d(theta) of every spectrum: array (n_spectra, n_par). Vectorised on a common grid (the spectra that needed NNLS one by one).
        For a line held at zero area the column is the central difference of the pre-2.19.5 code (see jacobian); one pair of series
        evaluations per parameter serves every spectrum in which that line is off."""
        th = np.asarray(th, float)
        ev, batch = self._eval(th, with_batch=True)
        S = len(self.data)
        out = np.full((S, self.n_par), np.nan)
        good = batch["good"] if batch is not None else np.zeros(S, bool)
        g = np.flatnonzero(good)
        if len(g):
            AMP = np.stack([e[0] for e in ev])
            D, DA = self._batch_derivs(th, batch, AMP)
            Q, R, r = batch["Q"][g], batch["R"][g], batch["r"][g]
            Ri = self._rinv_batch(R)
            dcoef = -Ri @ (np.swapaxes(Q, 1, 2) @ DA[g])
            dcoef += (Ri @ np.swapaxes(Ri, 1, 2))[:, :, self._pband] * np.einsum("snp,sn->sp", D[g], r)[:, None, :]
            da = dcoef[:, :self.nb, :]
            amp = AMP[g]
            s1, s2 = amp[:, g1].sum(1), amp[:, g2].sum(1)
            with np.errstate(divide="ignore", invalid="ignore"):
                og = (da[:, g1, :].sum(1) * s2[:, None] - s1[:, None] * da[:, g2, :].sum(1)) / (s2 ** 2)[:, None]
            og[s2 <= 0] = np.nan
            out[g] = og
        for ti in np.flatnonzero(~good):
            out[ti] = self.ratio_grad(th, ti, g1, g2)
        off = [set(range(self.nb)) - set(int(a) for a in e[5]) for e in ev]
        need = sorted({k for o in off for k in o})
        if need:
            for p in np.flatnonzero(np.isin(self._pband, need)):
                k = int(self._pband[p])
                h = 1e-5 * max(1.0, abs(th[p]))
                tp, tm = th.copy(), th.copy()
                tp[p] += h
                tm[p] -= h
                evp, evm = self._eval(tp), self._eval(tm)
                for ti in range(S):
                    if k not in off[ti]:
                        continue
                    ap, am = evp[ti][0], evm[ti][0]
                    rp = float(np.sum(ap[g1])) / float(np.sum(ap[g2])) if np.sum(ap[g2]) > 0 else np.nan
                    rm = float(np.sum(am[g1])) / float(np.sum(am[g2])) if np.sum(am[g2]) > 0 else np.nan
                    out[ti, p] = (rp - rm) / (2 * h)
            self._eval(th)                                               # keep the solution in the cache
        return out

    def ratio_grad(self, th, ti, g1, g2):
        """d(ratio)/d(theta) of spectrum ti (ratio = summed areas of group g1 / g2), exact: da/dtheta_p = -M+ dM/dtheta_p a + (M^T M)^-1 (dM/dtheta_p)^T r."""
        ev = self._eval(th)[ti]
        amp, active = ev[0], ev[5]
        D, DA, jpos, Q, R, r = self._derivs(th, ti, ev)
        Ri = self._rinv(R)
        dcoef = -Ri @ (Q.T @ DA)
        m = jpos >= 0
        if m.any():
            dcoef[:, m] += (Ri @ Ri.T)[:, jpos[m]] * (D[:, m].T @ r)[None, :]
        da = np.zeros((self.nb, self.n_par))
        da[active] = dcoef[:len(active)]
        s1, s2 = float(np.sum(amp[g1])), float(np.sum(amp[g2]))
        if s2 <= 0:
            return np.full(self.n_par, np.nan)
        return (da[g1].sum(0) * s2 - s1 * da[g2].sum(0)) / s2 ** 2


def _group_ratio(amp, cov, g1, g2):
    """Ratio of the summed areas of two groups of bands and its standard error with the full covariance of the areas."""
    s1, s2 = float(np.sum(amp[g1])), float(np.sum(amp[g2]))
    v1 = float(np.sum(cov[np.ix_(g1, g1)]))
    v2 = float(np.sum(cov[np.ix_(g2, g2)]))
    c12 = float(np.sum(cov[np.ix_(g1, g2)]))
    r = s1 / s2 if s2 > 0 else float("nan")
    if s1 > 0 and s2 > 0:
        rel2 = v1 / s1 ** 2 + v2 / s2 ** 2 - 2 * c12 / (s1 * s2)
        return r, abs(r) * math.sqrt(max(rel2, 0.0)), s1, s2, math.sqrt(max(v1, 0.0)), math.sqrt(max(v2, 0.0))
    return r, float("nan"), s1, s2, float("nan"), float("nan")


@single_blas
def fit_global(spectra, bands, window=None, bg_order=1, drift="auto", counts=False, max_nfev=80, groups=None, start=None):
    """Global fit of the whole series. spectra = [{'T': K, 'x': nm, 'y': counts}], bands = [(centre, fwhm_guess), ...] with the ratio pair FIRST
    (I1 = bands[0], I2 = bands[1]). drift = 'none' | 'drift' | 'auto' (BIC decides). Returns the result dict (see keys below).
    v2.18.4: groups = ([band indices of I1], [band indices of I2]) - the ratio of SUMMED areas, for bands that are split into several lines (Stark
    components of a lanthanide manifold): every line has its own shape, the ratio uses the sum of each group. Default ([0], [1]).
    v2.19.5: start = {parameter name: value} of an earlier fit (warm start; see _start_for)."""
    window = window or default_window(bands, np.concatenate([np.asarray(s["x"], float) for s in spectra]))
    if drift == "auto":
        if _warm():
            # v2.19.5: the drifting model CONTAINS the shared one; started from the shared solution (drifts 0) it can only do better - the BIC
            # comparison then compares the two models, not two convergences
            a = fit_global(spectra, bands, window, bg_order, "none", counts, max_nfev, groups, start=start)
            b = fit_global(spectra, bands, window, bg_order, "drift", counts, max_nfev, groups, start=a["params"])
        else:
            (oka, a), (okb, b) = run_independent([(fit_global, (spectra, bands, window, bg_order, "none", counts, max_nfev, groups), {}),
                                                  (fit_global, (spectra, bands, window, bg_order, "drift", counts, max_nfev, groups), {})])
            if not oka:
                raise a
            if not okb:
                raise b
        best = a if a["bic"] <= b["bic"] + 2.0 else b                    # prefer the simpler model unless the drift model is clearly better
        best["model_compare"] = {"shared shape": a["bic"], "drifting shape": b["bic"]}
        best["other_model_fit"] = b if best is a else a                 # v2.19.5: robust_ratio reuses it instead of fitting it again
        return best
    mdl = _Model(spectra, bands, window, bg_order, drift == "drift", drift == "drift", counts)
    th0 = mdl.p0.copy()

    fd = _fd_mode()
    xs = np.maximum(np.abs(th0) * 0.1, 0.05)
    if start and not fd:
        eps = 1e-6 * (mdl.ub - mdl.lb)
        for i_, nm_ in enumerate(mdl.names):
            if nm_ in start and np.isfinite(start[nm_]):
                th0[i_] = float(np.clip(start[nm_], mdl.lb[i_] + eps[i_], mdl.ub[i_] - eps[i_]))
    if fd:                                                               # pre-2.19.5 path (comparisons, tests)
        sol = least_squares(mdl.residual, th0, jac="2-point", bounds=(mdl.lb, mdl.ub), x_scale=xs, diff_step=1e-4, max_nfev=max_nfev,
                            ftol=1e-9, xtol=1e-10, loss="soft_l1", f_scale=3.0)       # robust: one dirty spectrum must not bend the shared shape
        th = sol.x
    else:
        # v2.19.5: (1) exact analytic Jacobian of the projected residual (was: finite differences = n_par + 1 evaluations of every spectrum per step);
        # (2) the fit moves the OFFSETS from the start values: the trust region then starts at one scale unit instead of the size of the parameters
        # themselves (centres ~500 nm), so the first step cannot throw a line far away (it used to widen one line over its neighbour and switch that
        # neighbour off - a poor local minimum the old code sometimes ended in, e.g. a shared-shape fit with a reduced chi-square of 58 instead of 1.0).
        sol = least_squares(lambda d: mdl.residual(th0 + d), np.zeros_like(th0), jac=lambda d: mdl.jacobian(th0 + d), bounds=(mdl.lb - th0, mdl.ub - th0),
                            x_scale=xs, max_nfev=max_nfev, ftol=1e-9, xtol=1e-10, loss="soft_l1", f_scale=3.0)
        th = th0 + sol.x
    n_pts = sum(len(d[0]) for d in mdl.data)
    n_lin = len(mdl.data) * (mdl.nb + mdl.bg_order + 1)
    n_par = len(th) + n_lin
    rss = float(np.sum(sol.fun ** 2))
    bic = n_pts * math.log(max(rss, 1e-300) / n_pts) + n_par * math.log(n_pts)
    # v2.18.5: covariance of the NON-LINEAR shape parameters (centres, widths, eta, drifts) from the Jacobian of the projected residual.
    # Their uncertainty propagates into every ratio through its sensitivity to them; without this term the ratio errors were 10-20 % too small.
    try:
        J_ = sol.jac
        cov_th = np.linalg.pinv(J_.T @ J_) * (rss / max(n_pts - n_par, 1))
    except Exception:
        cov_th = None
    g1_, g2_ = (list(groups[0]), list(groups[1])) if groups else ([0], [1])
    grads = mdl.ratio_grads(th, g1_, g2_) if (cov_th is not None and not fd) else None     # v2.19.5: exact, all spectra at once

    def ratio_at(th_, ti_):
        a_, _, _, _ = mdl.solve(th_, ti_)
        s1_, s2_ = float(np.sum(a_[g1_])), float(np.sum(a_[g2_]))
        return s1_ / s2_ if s2_ > 0 else float("nan")
    rows = []
    for ti, (x, y) in enumerate(mdl.data):
        amp, bg, A, B = mdl.solve(th, ti)
        fit = A @ amp + B @ bg
        res = y - fit
        w = mdl.w[ti]
        dof = max(len(y) - mdl.nb - mdl.bg_order - 1 - 1, 1)
        s2 = float(np.sum(w * res ** 2) / dof)
        Mfull = np.hstack([A, B]) * np.sqrt(w)[:, None]
        try:
            cov = np.linalg.inv(Mfull.T @ Mfull) * s2 * 1.0
        except np.linalg.LinAlgError:
            cov = np.linalg.pinv(Mfull.T @ Mfull) * s2
        # amplitudes were solved on a scaled problem; cov above is already in data units because Mfull uses the real columns
        var = np.diag(cov)[:mdl.nb]
        g1, g2 = (list(groups[0]), list(groups[1])) if groups else ([0], [1])
        r, r_err, a0, a1, e0, e1 = _group_ratio(amp, cov[:mdl.nb, :mdl.nb], g1, g2)
        if cov_th is not None and np.isfinite(r_err):                 # add the shape-parameter term: (dR/dtheta)^T Cov(theta) (dR/dtheta)
            if fd:
                grad = np.zeros(len(th))
                for k in range(len(th)):
                    h_ = 1e-5 * max(1.0, abs(th[k]))
                    tp_, tm_ = th.copy(), th.copy()
                    tp_[k] += h_; tm_[k] -= h_
                    grad[k] = (ratio_at(tp_, ti) - ratio_at(tm_, ti)) / (2 * h_)
            else:                                                     # v2.19.5: exact derivative (was 2 x n_par extra solves per spectrum)
                grad = grads[ti]
            var_shape = float(grad @ cov_th @ grad)
            r_err = math.sqrt(r_err ** 2 + max(var_shape, 0.0))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - float(np.sum(res ** 2)) / ss_tot if ss_tot > 0 else float("nan")
        rr = res - res.mean()
        dw = float(np.sum(np.diff(res) ** 2) / max(np.sum(res ** 2), 1e-300))         # Durbin-Watson: 2 = uncorrelated residuals
        cc, ww, ee = mdl.band_params(th, ti)
        flag = "good"
        if r2 < 0.98 or dw < 0.8 or (counts and s2 > 3.0) or not np.isfinite(r_err):
            flag = "check"
        if r2 < 0.9 or dw < 0.4 or (counts and s2 > 8.0) or a0 <= 0 or a1 <= 0:
            flag = "poor"
        areas_out = amp.copy()
        if groups:                                                  # the two ratio 'areas' are the group sums; the single lines stay in 'line_areas'
            areas_out = np.concatenate([[a0, a1], amp])
        rows.append({"T": float(mdl.T[ti]), "areas": areas_out, "line_areas": amp.copy(), "area_err": np.sqrt(np.maximum(var, 0)), "group_err": (e0, e1),
                     "ratio": float(r), "ratio_err": float(r_err),
                     "r2": r2, "chi2_red": s2, "dw": dw, "flag": flag, "x": x, "y": y, "fit": fit, "background": B @ bg,
                     "components": [A[:, k] * amp[k] for k in range(mdl.nb)], "centres": cc, "widths": ww, "etas": ee})
    esd = {}
    try:
        J = sol.jac
        cov_th = np.linalg.inv(J.T @ J) * (rss / max(n_pts - n_par, 1))
        esd = {n: float(math.sqrt(max(cov_th[i, i], 0))) for i, n in enumerate(mdl.names)}
    except Exception:
        pass
    return {"rows": rows, "params": dict(zip(mdl.names, map(float, th))), "esd": esd, "model": "drifting shape" if mdl.drift_pos else "shared shape",
            "bic": float(bic), "window": tuple(map(float, window)), "n_bands": mdl.nb, "bg_order": mdl.bg_order, "rss": rss, "T0": mdl.T0, "Ts": mdl.Ts,
            "groups": (list(groups[0]), list(groups[1])) if groups else None, "drift": "drift" if mdl.drift_pos else "none",
            "nfev": int(sol.nfev), "status": int(sol.status), "fit_bands": [tuple(map(float, b)) for b in bands], "counts": bool(counts)}


def window_integration(spectra, w1, w2, baseline="linear"):
    """The classic method for comparison: area in two fixed windows above the straight line through each window's ends."""
    I1, I2 = [], []
    for s in spectra:
        x = np.asarray(s["x"], float)
        y = np.asarray(s["y"], float)
        out = []
        for lo, hi in (w1, w2):
            m = (x >= min(lo, hi)) & (x <= max(lo, hi))
            if m.sum() < 3:
                out.append(0.0)
                continue
            xs, v = x[m], y[m]
            if baseline == "linear":
                v = v - np.interp(xs, [xs[0], xs[-1]], [v[0], v[-1]])
            out.append(float(np.trapezoid(v, xs)))
        I1.append(out[0]); I2.append(out[1])
    return np.array(I1), np.array(I2)


@single_blas
def robust_ratio(spectra, bands, window=None, bg_order=1, counts=False, rel_windows=(0.8, 1.25), main=None, bg_alternatives=None):
    """Main fit + refits with narrower / wider windows and the other shape model. Returns the main result with 'ratio_sys' (spread between the
    alternatives per spectrum) and 'ratio_total' = √(stat² + sys²).
    v2.18.4: bg_alternatives = [(label, spectra)] - the same series after OTHER reasonable background subtractions (see plbackground.alternatives); each is refitted
    with the same bands, window and shape model, and the spread of the ratio over them is the background part of the systematic error ('ratio_sys_bg';
    'ratio_sys_fit' is the part from windows and shape model)."""
    main = main or fit_global(spectra, bands, window, bg_order, "auto", counts)
    bands = list(main.get("bands", bands)) + list(main.get("extra", []))
    groups = main.get("groups")
    lo, hi = main["window"]
    mid, half = (lo + hi) / 2, (hi - lo) / 2
    alts = []
    other = "drift" if main["model"] == "shared shape" else "none"
    allx = np.concatenate([np.asarray(s["x"], float) for s in spectra])
    win_list = [(max(mid - half * f, allx.min()), min(mid + half * f, allx.max())) for f in rel_windows]
    win_labels = [f"fit range x{f}" for f in rel_windows]
    if main.get("ranges"):                                   # manifold: the two ranges must stay whole; only the margin outside them changes (7 or 25 nm instead of 15)
        r1, r2 = main["ranges"]
        a, b = min(r1[0], r2[0]), max(r1[1], r2[1])
        win_list = [(max(a - pad, allx.min()), min(b + pad, allx.max())) for pad in (7.0, 25.0)]
        win_labels = ["margin 7 nm around the ranges", "margin 25 nm around the ranges"]
    used_labels = []
    same = "none" if main["model"] == "shared shape" else "drift"
    win_calls, win_labs = [], []
    for w, lab in zip(win_list, win_labels):
        try:
            bb = [bd for bd in bands if w[0] <= bd[0] <= w[1]] if groups else bands
            if groups:
                g1 = [k for k, (c, _) in enumerate(bb) if min(main["ranges"][0]) <= c <= max(main["ranges"][0])] if main.get("ranges") else groups[0]
                g2 = [k for k, (c, _) in enumerate(bb) if min(main["ranges"][1]) <= c <= max(main["ranges"][1])] if main.get("ranges") else groups[1]
                gg = (g1, g2)
            else:
                gg = None
        except Exception:
            continue
        win_calls.append((fit_global, (spectra, bb, w, bg_order, same, counts, 80, gg), {"start": _start_for(main.get("params"), main.get("fit_bands", bands), bb) if _warm() else None}))
        win_labs.append(lab)
    # v2.19.5: the other shape model was already fitted by the 'auto' choice; it is reused when it is EXACTLY the fit that would be made here
    # (same lines, window, groups, background order, weights, and converged in fewer than 80 evaluations = the same path as a refit capped at 80)
    of = main.pop("other_model_fit", None)
    gl = lambda g: None if g is None else [list(map(int, x)) for x in g]      # noqa: E731
    reuse_other = bool(of is not None and of.get("drift") == other and of.get("fit_bands") == [tuple(map(float, b)) for b in bands]
                       and tuple(of.get("window", ())) == tuple(map(float, main["window"])) and gl(of.get("groups")) == gl(groups)
                       and of.get("bg_order") == int(bg_order) and of.get("counts") == bool(counts) and of.get("status", 0) > 0 and of.get("nfev", 10 ** 9) < 80)
    other_calls = [] if reuse_other else [(fit_global, (spectra, bands, main["window"], bg_order, other, counts, 80, groups),
                                           {"start": _start_for(main.get("params"), main.get("fit_bands", bands), bands) if _warm() else None})]
    bg_list = list(bg_alternatives or [])
    st0 = _start_for(main.get("params"), main.get("fit_bands", bands), bands) if _warm() else None
    bg_calls = [(fit_global, (sp_alt, bands, main["window"], bg_order, same, counts, 40, groups), {"start": st0}) for _, sp_alt in bg_list]
    results = run_independent(win_calls + other_calls + bg_calls)          # v2.19.5: all alternatives at the same time
    for (ok, r), lab in zip(results[:len(win_calls)], win_labs):
        if ok:
            alts.append(r)
            used_labels.append(lab)
    if reuse_other:
        alts.append(of)
    elif other_calls and results[len(win_calls)][0]:
        alts.append(results[len(win_calls)][1])
    R = np.array([r["ratio"] for r in main["rows"]])
    stat = np.array([r["ratio_err"] for r in main["rows"]])
    if alts:
        A = np.array([[r["ratio"] for r in a["rows"]] for a in alts])
        sys_ = np.sqrt(np.nanmean((A - R[None, :]) ** 2, axis=0))
    else:
        sys_ = np.zeros_like(R)
    sys_bg = np.zeros_like(R)
    bg_rows, bg_labels = [], []
    for (lab, _), (ok, a) in zip(bg_list, results[len(win_calls) + len(other_calls):]):
        if ok and len(a["rows"]) == len(R):
            bg_rows.append([r["ratio"] for r in a["rows"]])
            bg_labels.append(lab)
    if bg_rows:
        Ab = np.array(bg_rows)
        sys_bg = np.sqrt(np.nanmean((Ab - R[None, :]) ** 2, axis=0))
    main["ratio_sys_fit"] = sys_
    main["ratio_sys_bg"] = sys_bg
    main["bg_alt_labels"] = bg_labels
    main["bg_alt_ratios"] = np.array(bg_rows) if bg_rows else None
    main["fit_alt_ratios"] = A if alts else None                      # v2.18.4: kept so that the Boltzmann fit can be repeated on every alternative
    main["fit_alt_labels"] = used_labels + (["other peak model"] if len(alts) > len(used_labels) else [])
    main["ratio_sys"] = np.sqrt(sys_ ** 2 + sys_bg ** 2)
    main["ratio_total"] = np.sqrt(np.nan_to_num(stat) ** 2 + main["ratio_sys"] ** 2)
    main["n_alternatives"] = len(alts) + len(bg_rows)
    return main


def suggest_extra_band(res, min_sigma=4.0, exclude=()):
    """Looks for a clear structure left in the residual (summed over all spectra) and proposes an extra band: (centre, fwhm) or None."""
    xs = np.unique(np.round(np.concatenate([r["x"] for r in res["rows"]]), 3))
    acc = np.zeros_like(xs)
    cnt = np.zeros_like(xs)
    for r in res["rows"]:
        d = r["y"] - r["fit"]
        s = max(float(np.std(d)), 1e-12)
        acc += np.interp(xs, r["x"], d / s)
        cnt += 1
    z = acc / max(cnt.max(), 1) * math.sqrt(len(res["rows"]))
    step = float(np.median(np.diff(xs))) if len(xs) > 1 else 1.0
    zp = np.concatenate([[z.min() - 1.0], z, [z.min() - 1.0]])              # pad: a structure at the edge of the window is a peak too
    pk, pr = find_peaks(zp, height=min_sigma, distance=max(3, int(1.0 / max(step, 1e-6))))
    pk = pk - 1
    pk = pk[(pk >= 0) & (pk < len(z))]
    if not len(pk):
        return None
    pk = [q for q in pk if not any(abs(xs[q] - c) < 0.7 * w for c, w in exclude)]       # not inside the cores of the bands of the ratio
    if not len(pk):
        return None
    i = pk[int(np.argmax(z[pk]))]
    half = z[i] / 2
    a = i
    while a > 0 and z[a] > half:
        a -= 1
    b = i
    while b < len(z) - 1 and z[b] > half:
        b += 1
    return float(xs[i]), float(max(xs[b] - xs[a], 3 * step))


@single_blas
def auto_fit(spectra, counts=False, lo=None, hi=None, bg_order=1, max_extra=3, bands=None):
    """Complete automatic analysis: find the two bands (also hidden shoulders), fit the series, then add up to `max_extra` nuisance bands
    while a clear structure remains in the residual AND the BIC improves. Returns the fit result with 'bands' (centre, fwhm) and 'extra' (list)."""
    x = np.concatenate([np.asarray(s["x"], float) for s in spectra])
    grid = np.asarray(spectra[0]["x"], float)
    mean = np.mean([np.interp(grid, np.asarray(s["x"], float), np.asarray(s["y"], float)) for s in spectra], axis=0)
    bands = list(bands) if bands else auto_bands(grid, mean, lo, hi)
    if not bands or len(bands) < 2:
        return None
    bands = bands[:2]
    win = default_window(bands, x)
    if lo is not None and hi is not None:                          # the user's / full spectral range is the fit window (tails and neighbours matter)
        win = (float(max(lo, x.min())), float(min(hi, x.max())))
    best = fit_global(spectra, bands, win, bg_order, "auto", counts)
    extra = []
    for _ in range(max_extra):
        sug = suggest_extra_band(best, exclude=bands + extra)
        if sug is None or any(abs(sug[0] - c) < 0.6 * w for c, w in bands + extra):
            break
        trial = fit_global(spectra, bands + extra + [sug], win, bg_order, "auto", counts,
                           start=_start_for(best.get("params"), bands + extra, bands + extra + [sug]) if _warm() else None)
        if trial["bic"] < best["bic"] - 10.0:
            best, extra = trial, extra + [sug]
        else:
            break
    best["bands"], best["extra"] = bands, extra
    return best


# ------------------------------------------------------------------------------------------------------------- v2.18.4: crowded spectra (manifolds)
def auto_lines(x, y, lo=None, hi=None, max_lines=16, min_rel=0.02):
    """ALL lines of a (mean) spectrum between lo and hi: maxima and shoulders (minima of the second derivative), merged when closer than ~1.5 nm,
    weaker than `min_rel` of the strongest dropped. Returns [(centre, fwhm_guess)] sorted by centre (at most max_lines, strongest kept)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    lo = x.min() if lo is None else lo
    hi = x.max() if hi is None else hi
    m = (x >= lo) & (x <= hi)
    xs, ys = x[m], y[m]
    if len(xs) < 15:
        return []
    from scipy.signal import savgol_filter
    step = float(np.median(np.diff(xs)))
    win = max(5, int(round(2.0 / step)) | 1)
    ysm = savgol_filter(ys, min(win, len(ys) - (1 - len(ys) % 2)), 2) if len(ys) > win else ys
    base = float(np.percentile(ysm, 5))
    rng = float(ysm.max() - base)
    if rng <= 0:
        return []
    w2 = max(7, int(round(2.5 / step)) | 1)
    d2 = savgol_filter(ys, min(w2, len(ys) - (1 - len(ys) % 2)), 3, deriv=2)
    cand = list(find_peaks(ysm, prominence=min_rel * rng)[0]) + list(find_peaks(-d2, prominence=0.02 * float(np.ptp(d2)))[0])
    cand = sorted(set(int(i) for i in cand if ysm[int(i)] - base > min_rel * rng))
    merged = []
    for i in cand:
        if merged and xs[i] - xs[merged[-1]] < max(1.5, 2.5 * step):
            if ysm[i] > ysm[merged[-1]]:
                merged[-1] = i
        else:
            merged.append(i)
    merged = sorted(sorted(merged, key=lambda i: -ysm[i])[:max_lines])
    out = []
    for k, p_ in enumerate(merged):
        gaps = [abs(xs[p_] - xs[q]) for q in merged if q != p_]
        h = ysm[p_]
        half = h - (h - base) / 2.0
        i0, i1 = p_, p_
        while i0 > 0 and ysm[i0] > half:
            i0 -= 1
        while i1 < len(ysm) - 1 and ysm[i1] > half:
            i1 += 1
        wd = max(float(xs[i1] - xs[i0]), 2.5 * step)
        if gaps:
            wd = min(wd, 1.1 * min(gaps))
        out.append((float(xs[p_]), float(max(wd, 2.5 * step))))
    return out


@single_blas
def fit_manifolds(spectra, range1, range2, counts=False, lo=None, hi=None, bg_order=1, max_extra=4, lines=None):
    """Ratio of two MANIFOLDS: every line found between lo and hi gets its own pseudo-Voigt (shape shared by the series, optionally drifting); the lines
    whose centre lies in range1 are summed into I1, those in range2 into I2, the others are nuisance lines. Lines left in the residual are added while the
    BIC improves. Returns the fit result with 'bands' (all lines), 'groups', 'extra' = [] and 'mode' = 'manifold', or None."""
    grid = np.asarray(spectra[0]["x"], float)
    mean = np.mean([np.interp(grid, np.asarray(s["x"], float), np.asarray(s["y"], float)) for s in spectra], axis=0)
    xall = np.concatenate([np.asarray(s["x"], float) for s in spectra])
    lo = float(max(lo if lo is not None else xall.min(), xall.min()))
    hi = float(min(hi if hi is not None else xall.max(), xall.max()))
    lines = list(lines) if lines else auto_lines(grid, mean, lo, hi)

    def grp(ls):
        g1 = [k for k, (c, _) in enumerate(ls) if min(range1) <= c <= max(range1)]
        g2 = [k for k, (c, _) in enumerate(ls) if min(range2) <= c <= max(range2)]
        return g1, g2
    g1, g2 = grp(lines)
    if not g1 or not g2:
        return None
    best = fit_global(spectra, lines, (lo, hi), bg_order, "auto", counts, 120, (g1, g2))
    for _ in range(max_extra):
        sug = suggest_extra_band(best, min_sigma=4.0, exclude=[(c, 0.45 * w) for c, w in lines])       # lines can be closer than a width: only the cores are excluded
        if sug is None or any(abs(sug[0] - c) < 0.3 * w for c, w in lines):
            break
        trial_lines = sorted(lines + [sug])
        t1, t2 = grp(trial_lines)
        trial = fit_global(spectra, trial_lines, (lo, hi), bg_order, "auto", counts, 120, (t1, t2),
                           start=_start_for(best.get("params"), lines, trial_lines) if _warm() else None)
        if trial["bic"] < best["bic"] - 10.0:
            best, lines, g1, g2 = trial, trial_lines, t1, t2
        else:
            break
    best["bands"], best["extra"], best["groups"], best["mode"] = lines, [], (g1, g2), "manifold"
    best["ranges"] = (tuple(sorted(range1)), tuple(sorted(range2)))
    return best


def flag_share(fit):
    """Fraction of the spectra of a fit result that are not 'good'."""
    rows = fit.get("rows") or []
    return float(np.mean([r["flag"] != "good" for r in rows])) if rows else 1.0
