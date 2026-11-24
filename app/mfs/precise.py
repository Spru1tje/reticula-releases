"""v2.21: precise sphere mode for small particles (sub-20 nm on the usual TEM scale).

The circle fit of mfs.spheres uses the visible edge points (ray casting). On a small particle there are only a few
pixels on the edge, so its radius rests on a handful of noisy samples. This module fits a physical model to ALL pixels
of a window around each circle: a disc with an anti-aliased, blurred edge,

    model(x, y) = b + A * 0.5 * erfc((d - r) / (sqrt(2) * s)),   d = distance from the centre (cx, cy),

with a background b, a contrast amplitude A, the radius r (the 50 % point of the edge, independent of any threshold)
and the edge width s (camera blur and halo). The fit is a bounded non-linear least-squares problem per particle.

The radius returned is the 50 % point, so no inward correction is needed; that is what the blur parameter does.
Measured on synthetic particles with known radius: see tests/test_v221_precise.py.
"""
from __future__ import annotations

import functools
import math
import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy import ndimage as ndi
from scipy import special
from scipy.optimize import least_squares


@functools.lru_cache(maxsize=256)
def _kernel(s_q):
    """1-D Gaussian kernel of width s_q (px), radius and normalisation as scipy.ndimage.gaussian_filter (truncate 4)."""
    radius = int(4.0 * s_q + 0.5)
    x = np.arange(-radius, radius + 1, dtype=float)
    k = np.exp(-0.5 * x * x / (s_q * s_q))
    return k / k.sum()


def _blur(img, s):
    """Separable Gaussian blur with a cached kernel (the width is quantised to 0.02 px: the fit cannot tell the difference)."""
    k = _kernel(round(float(s) / 0.02) * 0.02)
    out = ndi.correlate1d(img, k, axis=0, mode="nearest")
    return ndi.correlate1d(out, k, axis=1, mode="nearest")


def _thickness(dist2, r):
    """Projected thickness of a sphere of radius r at squared distance dist2 from its centre (0 outside)."""
    return 2.0 * np.sqrt(np.clip(r * r - dist2, 0.0, None))


def _model(p, yy, xx, y0, x0, shape):
    """Sphere image: b + A * (1 - exp(-k t)), t = projected thickness, convolved with a Gaussian of width s (pixels).
    Evaluated on the window (shape) with the blur done on a padded grid, so the edge of the image does not leak in."""
    x_c, y_c, r, k, A, b, s = p
    t = _thickness((yy - y_c) ** 2 + (xx - x_c) ** 2, r)
    img = A * (1.0 - np.exp(-k * t))
    if s > 0.05:
        img = _blur(img, s)
    return b + img


def fit_disc(c, cx, cy, r0, s0=1.0, pad=None):
    """Fit one absorbing sphere to the contrast image c around (cx, cy). The model is the blurred projected sphere
    (dome-shaped absorption), so the 50 % radius is physical and does not depend on a threshold.
    Returns dict(cx, cy, r, k, A, b, s, rms, ok) or None."""
    H, W = c.shape
    r0 = float(max(r0, 1.0))
    pad = int(math.ceil(3.0 * s0 + 4)) if pad is None else int(pad)
    R = int(math.ceil(r0 + pad))
    y0, y1 = max(0, int(round(cy)) - R), min(H, int(round(cy)) + R + 1)
    x0, x1 = max(0, int(round(cx)) - R), min(W, int(round(cx)) + R + 1)
    if y1 - y0 < 6 or x1 - x0 < 6:
        return None
    yy, xx = np.mgrid[y0:y1, x0:x1].astype(float)
    z = c[y0:y1, x0:x1].astype(float)
    if z.size < 36:
        return None
    inside = np.hypot(yy - cy, xx - cx) <= 0.6 * r0
    outside = np.hypot(yy - cy, xx - cx) >= r0 + pad - 2
    b0 = float(np.median(z[outside])) if outside.any() else 0.0
    A0 = max(float(np.median(z[inside]) - b0) if inside.any() else 1.0, 1e-6)

    def resid(p):
        return (_model(p, yy, xx, y0, x0, z.shape) - z).ravel()

    # start: absorption k with the core at about A0; edge width s0
    p0 = np.array([cx, cy, r0, 0.5 / r0, A0, b0, s0])
    lo = np.array([cx - 0.5 * r0, cy - 0.5 * r0, 0.5 * r0, 1e-4, 0.0, -np.inf, 0.2])
    hi = np.array([cx + 0.5 * r0, cy + 0.5 * r0, 1.6 * r0, 10.0, np.inf, np.inf, 3.0])
    p0 = np.clip(p0, lo + 1e-9, hi - 1e-9)
    sol = least_squares(resid, p0, bounds=(lo, hi), x_scale=np.array([1, 1, 1, 0.1, A0, A0, 0.3]), max_nfev=300)
    x_c, y_c, r, k, A, b, s = sol.x
    rms = float(np.sqrt(np.mean(sol.fun ** 2)))
    ok = bool(sol.success and A > 0.05 * A0 and np.isfinite(rms))
    return {"cx": float(x_c), "cy": float(y_c), "r": float(r), "k": float(k), "A": float(A), "b": float(b),
            "s": float(s), "rms": rms, "ok": ok}


def refine_circles(c, circles, s0=1.2):
    """Replace the circle fits of spheres.resolve() by the model fit where it converged; others keep their values.
    Returns the number refined. Each circle dict is updated in place (cx, cy, r, rms, precise=True)."""
    n_ok = 0
    for o in circles:
        f = fit_disc(c, o["cx"], o["cy"], o.get("r_edge", o["r"]), s0=s0)
        if f is None or not f["ok"]:
            o["precise"] = False
            continue
        # keep the fit only when it stays close to the ray fit (a large jump means the model found another object)
        r_old = o.get("r_edge", o["r"])
        if math.hypot(f["cx"] - o["cx"], f["cy"] - o["cy"]) > 0.5 * r_old or not (0.5 * r_old < f["r"] < 1.6 * r_old):
            o["precise"] = False
            continue
        o.update(cx=f["cx"], cy=f["cy"], r_edge=f["r"], r=f["r"], rms_model=f["rms"], edge_s=f["s"], precise=True)
        n_ok += 1
    return n_ok


def fit_disc_cell(c, cx, cy, r0, nearest, owner, s0=1.0, tol=1e-8, max_nfev=300):
    """As fit_disc, but only the pixels of this particle's own cell (cell_mask on the whole image) enter the fit, so
    neighbours do not bias the radius. Returns the same dict or None."""
    H, W = c.shape
    r0 = float(max(r0, 1.0))
    pad = int(math.ceil(3.0 * s0 + 4))
    R = int(math.ceil(r0 + pad))
    y0, y1 = max(0, int(round(cy)) - R), min(H, int(round(cy)) + R + 1)
    x0, x1 = max(0, int(round(cx)) - R), min(W, int(round(cx)) + R + 1)
    if y1 - y0 < 6 or x1 - x0 < 6:
        return None
    yy, xx = np.mgrid[y0:y1, x0:x1].astype(float)
    z = c[y0:y1, x0:x1].astype(float)
    cell = nearest[y0:y1, x0:x1] == owner
    inside = (np.hypot(yy - cy, xx - cx) <= 0.6 * r0) & cell
    outside = (np.hypot(yy - cy, xx - cx) >= r0 + pad - 2) & cell
    if inside.sum() < 4 or outside.sum() < 4:
        return None
    b0 = float(np.median(z[outside]))
    A0 = max(float(np.median(z[inside]) - b0), 1e-6)
    # the model is evaluated on the full window (so the blur sees the true neighbourhood); only the cell pixels count
    def resid(p):
        m = _model(p, yy, xx, y0, x0, z.shape)
        return ((m - z) * cell).ravel()
    p0 = np.array([cx, cy, r0, 0.5 / r0, A0, b0, s0])
    lo = np.array([cx - 0.5 * r0, cy - 0.5 * r0, 0.5 * r0, 1e-4, 0.0, -np.inf, 0.2])
    hi = np.array([cx + 0.5 * r0, cy + 0.5 * r0, 1.6 * r0, 10.0, np.inf, np.inf, 3.0])
    p0 = np.clip(p0, lo + 1e-9, hi - 1e-9)
    sol = least_squares(resid, p0, bounds=(lo, hi), x_scale=np.array([1, 1, 1, 0.1, A0, A0, 0.3]), max_nfev=max_nfev, ftol=tol, xtol=tol, gtol=tol)
    x_c, y_c, r, k, A, b, s = sol.x
    npx = int(cell.sum())
    rms = float(np.sqrt(np.mean(sol.fun[cell.ravel()] ** 2))) if npx else float("inf")
    ok = bool(sol.success and A > 0.05 * A0 and np.isfinite(rms))
    return {"cx": float(x_c), "cy": float(y_c), "r": float(r), "k": float(k), "A": float(A), "b": float(b),
            "s": float(s), "rms": rms, "ok": ok, "rel_rms": rms / max(A, 1e-9)}


def refine_labels(c, lab, info, s0=1.0):
    """Precise sizes for the particles of an outline label image (the candidates of the mask segmentation).
    Each particle is fitted with the blurred absorbing-sphere model on the pixels of its own Voronoi cell. Returns
    (labels, info, summary). Particles whose fit fails keep kind 'lowconf' (excluded by default); touching particles
    are fitted separately; particles touching the image edge are 'cut' and not measured."""
    H, W = lab.shape
    ids = [i for i, v in info.items() if isinstance(v, dict) and isinstance(i, int) and v.get("kind") == "particle"]
    sls = ndi.find_objects(lab)
    seeds = np.zeros((H, W), int)
    cent = {}
    for i in ids:
        sl = sls[i - 1] if i - 1 < len(sls) else None
        if sl is None:
            continue
        sub = lab[sl] == i
        yy_s, xx_s = np.nonzero(sub)
        cy = sl[0].start + yy_s.mean(); cx = sl[1].start + xx_s.mean()
        cent[i] = (cx, cy, math.sqrt(sub.sum() / math.pi))
        seeds[int(round(cy)), int(round(cx))] = i
    _, (iy, ix) = ndi.distance_transform_edt(seeds == 0, return_indices=True)
    nearest = seeds[iy, ix]                                   # Voronoi cell of every pixel
    out = np.zeros_like(lab)
    new_info = {}
    n_ok = n_cut = n_low = 0
    jobs = []
    for i, (cx, cy, r0) in cent.items():
        # the particle is cut when its outline reaches the image border
        if cx - r0 < 1.5 or cy - r0 < 1.5 or cx + r0 > W - 2.5 or cy + r0 > H - 2.5:
            n_cut += 1
            continue
        jobs.append((i, cx, cy, r0))

    def _run(job):
        i, cx, cy, r0 = job
        return fit_disc_cell(c, cx, cy, r0, nearest, i, s0=s0)

    # the particles are independent: fit them in parallel (NumPy/SciPy release the GIL for the array work)
    workers = max(1, min(8, os.cpu_count() or 1))
    if workers > 1 and len(jobs) > 8:
        with ThreadPoolExecutor(max_workers=workers) as ex:
            fits = list(ex.map(_run, jobs))
    else:
        fits = [_run(j) for j in jobs]
    for (i, cx, cy, r0), f in zip(jobs, fits):
        good = bool(f and f["ok"] and 0.5 * r0 < f["r"] < 1.6 * r0 and math.hypot(f["cx"] - cx, f["cy"] - cy) < 0.5 * r0)
        if not good:
            n_low += 1
            new_info[i] = {"kind": "lowconf", "solidity": info.get(i, {}).get("solidity", 1.0), "precise": False}
            out[lab == i] = i
            continue
        n_ok += 1
        # the disc is drawn only in its window (a full-image grid per particle was most of the time before)
        rr = int(math.ceil(f["r"])) + 2
        y0, y1 = max(0, int(f["cy"]) - rr), min(H, int(f["cy"]) + rr + 1)
        x0, x1 = max(0, int(f["cx"]) - rr), min(W, int(f["cx"]) + rr + 1)
        yy, xx = np.mgrid[y0:y1, x0:x1]
        disc = np.hypot(yy - f["cy"], xx - f["cx"]) <= f["r"]
        win = out[y0:y1, x0:x1]
        win[disc & (nearest[y0:y1, x0:x1] == i) & (win == 0)] = i
        new_info[i] = {"kind": "particle", "solidity": 1.0, "precise": True, "d_px": 2 * f["r"], "rms": f["rms"],
                       "edge_s": f["s"], "circle": (f["cx"], f["cy"], f["r"]), "conf": 1.0}
    summary = {"mode": "precise", "active": True, "n_precise": n_ok, "n_lowconf": n_low, "n_edge_excluded": n_cut,
               "n_candidates": len(cent)}
    info_out = {k: v for k, v in info.items() if not isinstance(k, int)}
    info_out.update(new_info)
    return out, info_out, summary
