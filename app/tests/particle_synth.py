"""Synthetic TEM-like images with particles of KNOWN size and shape (for accuracy tests of the particle sizer).
Dark particles on a bright, unevenly lit support; camera blur + noise; the true area of every particle is computed on a 4x
super-sampled grid, so the true equivalent diameter is exact."""
import math

import numpy as np
from scipy import ndimage as ndi

SHAPES = ("disc", "ellipse", "square", "hexagon", "rod", "blob")


def _mask(shape, cx, cy, size, rng, ss=4):
    """Boolean mask on the ss-times finer grid of one particle. size = equivalent diameter (px)."""
    H, W = shape
    y, x = np.mgrid[0:H * ss, 0:W * ss] / ss
    X, Y = x - cx, y - cy
    ang = rng.uniform(0, math.pi)
    Xr, Yr = X * math.cos(ang) + Y * math.sin(ang), -X * math.sin(ang) + Y * math.cos(ang)
    R = size / 2
    return Xr, Yr, R


def make(shape_name, n=24, size=56, size_cv=0.06, hw=(900, 900), blur=2.0, noise=3.0, seed=0, ss=4):
    """Returns (image uint8, [dict(cx, cy, d_true (equivalent circle diameter, px), shape)])."""
    rng = np.random.default_rng(seed)
    H, W = hw
    y, x = np.mgrid[0:H * ss, 0:W * ss] / ss
    mask = np.zeros((H * ss, W * ss), bool)
    truth, centres = [], []
    tries = 0
    while len(truth) < n and tries < 5000:
        tries += 1
        d = size * rng.normal(1.0, size_cv)
        cx, cy = rng.uniform(1.0 * d, W - 1.0 * d), rng.uniform(1.0 * d, H - 1.0 * d)
        if any(math.hypot(cx - a, cy - b) < 1.45 * (d + e) / 2 + 12 for a, b, e in centres):
            continue
        ang = rng.uniform(0, math.pi)
        X, Y = x - cx, y - cy
        Xr, Yr = X * math.cos(ang) + Y * math.sin(ang), -X * math.sin(ang) + Y * math.cos(ang)
        R = d / 2
        if shape_name == "disc":
            m = Xr ** 2 + Yr ** 2 <= R ** 2
        elif shape_name == "ellipse":                       # aspect 1.6 at equal area
            a, b = R * math.sqrt(1.6), R / math.sqrt(1.6)
            m = (Xr / a) ** 2 + (Yr / b) ** 2 <= 1
        elif shape_name == "square":                        # equal area: side = d*sqrt(pi)/2
            h = d * math.sqrt(math.pi) / 4
            m = (np.abs(Xr) <= h) & (np.abs(Yr) <= h)
        elif shape_name == "hexagon":
            a = R * math.sqrt(2 * math.pi / (3 * math.sqrt(3)))        # circumradius for equal area
            m = (np.abs(Yr) <= a * math.sqrt(3) / 2) & (np.abs(Xr) * math.sqrt(3) / 2 + np.abs(Yr) / 2 <= a * math.sqrt(3) / 2)
        elif shape_name == "rod":                           # stadium 3:1 at equal area
            Lr = 3.0
            w = R * 2 * math.sqrt(1 / (Lr - 1 + math.pi / 4)) / 2 * 1.0
            w = math.sqrt(math.pi * R * R / ((Lr - 1) + math.pi / 4)) / 1.0 / 2 * 1.0
            half = (Lr - 1) * w
            m = ((np.abs(Xr) <= half) & (np.abs(Yr) <= w)) | ((np.abs(Xr) - half) ** 2 + Yr ** 2 <= w * w)
        elif shape_name == "blob":                          # random Fourier outline, 10-15 % radial irregularity
            th = np.linspace(0, 2 * np.pi, 361)
            r = 1 + sum(rng.normal(0, 0.06 / k) * np.cos(k * th + rng.uniform(0, 6.28)) for k in (2, 3, 4, 5))
            area = 0.5 * np.trapezoid(r ** 2, th)
            r = r / math.sqrt(area / math.pi) * R            # equal-area scaling
            T = np.arctan2(Yr, Xr) % (2 * np.pi)
            m = np.hypot(Xr, Yr) <= np.interp(T, th, r)
        else:
            raise ValueError(shape_name)
        if (m & mask).any():
            continue
        mask |= m
        centres.append((cx, cy, d))
        truth.append({"cx": cx, "cy": cy, "d_true": 2 * math.sqrt(m.sum() / ss / ss / math.pi), "shape": shape_name})
    fine = mask.astype(float)
    part = fine.reshape(H, ss, W, ss).mean((1, 3))                  # coverage per pixel
    bg = 205 + 18 * (np.mgrid[0:H, 0:W][1] / W - 0.5) + 10 * (np.mgrid[0:H, 0:W][0] / H - 0.5)   # uneven illumination
    core = ndi.gaussian_filter(rng.normal(0, 1, (H, W)), 6) * 12                                      # particle texture
    img = bg * (1 - part) + (75 + core) * part
    img = ndi.gaussian_filter(img, blur) + rng.normal(0, noise, img.shape)
    return np.clip(img, 0, 255).astype(np.uint8), truth


def evaluate(img, truth, lab, f=1, info=None):
    """Match detections to truth by centroid; returns list of (shape, d_true, d_est) for the clean (non-aggregate) detections."""
    import math as _m
    out = []
    ids = [i for i in np.unique(lab) if i > 0]
    cm = ndi.center_of_mass(lab > 0, lab, ids) if ids else []
    areas = ndi.sum(lab > 0, lab, ids) if ids else []
    for t in truth:
        best, bd = None, 1e9
        for i, (cy, cx), A in zip(ids, cm, areas):
            dd = _m.hypot(cx * f - t["cx"], cy * f - t["cy"])
            if dd < bd:
                bd, best = dd, (i, A)
        if best is not None and bd < 0.3 * t["d_true"]:
            i, A = best
            if isinstance(info, dict) and isinstance(info.get(i), dict) and info[i].get("kind") in ("aggregate",):
                continue
            ci = info.get(i) if isinstance(info, dict) else None
            d = ci["d_px"] * f if isinstance(ci, dict) and ci.get("d_px") else 2 * _m.sqrt(A * f * f / _m.pi)
            out.append((t["shape"], t["d_true"], d))
    return out
