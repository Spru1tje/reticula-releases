"""v2.12 sphere model for the particle sizer: circle fitting on the real particle edge, splitting of touching
doubles / triples, overlapping (occluded) spheres and a per-particle confidence.

Pipeline (per connected cluster of the foreground mask of the adaptive segmentation):
  1. candidate centres from a gradient-direction Hough transform of the particle contrast (edge pixels vote along
     their gradient for centres at radius rmin..rmax) plus the distance-map maxima;
  2. each candidate is refined by casting rays from its centre; on every ray the sub-pixel position where the contrast
     falls to `level` × the particle's core contrast is the edge. Rays that meet a neighbour before reaching the
     background (contact / overlap) give no edge point; a robust least-squares circle is fitted to the free edge points
     (so an occluded sphere is fitted on its visible arcs only); the edge is then moved inward by `delta_px` (blur /
     halo of the camera, calibrated on the real TEM series, see NOTES v2.12);
  3. greedy selection, best confidence first: a circle is kept when it explains enough new foreground area and its
     centre is not a duplicate. Overlap lenses that are darker than both spheres (double thickness) confirm an
     occlusion; a heavily overlapping second circle without a darker lens is rejected (one particle, not two);
  4. confidence = visible-arc term × fit term × coverage term (0–1). Below `conf_min` the particle is kept but
     classified 'lowconf' (excluded from the histogram by default, orange outline, click to include).
Only numpy / scipy (no scikit-image), works on the working grid of segment_adaptive."""
from __future__ import annotations
import math
import numpy as np
from scipy import ndimage as ndi

# calibrated defaults (cross-validated on the 15 usable SiO2-NH2 TEM images 12k–120k, see v2.12 NOTES)
LEVEL = 0.35         # edge = where the contrast falls to LEVEL × core contrast (CV optimum 0.25–0.4)
DELTA_PX = 2.5       # inward edge correction, px of the working grid (camera blur + halo; CV on 12k–120k)
CONF_MIN = 0.5
N_RAYS = 96


def roundness(lab, info, L, min_n=3):
    """Is this a population of round particles? Uses the clean isolated singles: mean |r(θ)/r̄ − 1| of their outline.
    Returns {'n', 'roundness' (lower = rounder), 'round' (bool), 'aspect'}."""
    from mfs import shapeprior as SP
    A1 = L * L
    areas = np.bincount(lab.ravel())
    for max_contact, sol_min in ((0.15, 0.9), (0.35, 0.95)):    # isolated singles first; few of them → clean touching ones
        rs, asp = _round_stats(lab, info, A1, areas, max_contact, sol_min, SP)
        if len(rs) >= min_n:
            break
    if len(rs) < min_n:
        return {"n": len(rs), "roundness": None, "round": False, "aspect": None,
                "reason": f"only {len(rs)} clean isolated particle(s) – shape not decided"}
    r = float(np.median(rs)); a = float(np.median(asp))
    return {"n": len(rs), "roundness": r, "aspect": a, "round": bool(r < 0.05 and a < 1.15)}


def _round_stats(lab, info, A1, areas, max_contact, sol_min, SP):
    rs, asp = [], []
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        if sl is None or not isinstance(info.get(i), dict) or info[i].get("kind") != "particle":
            continue
        if not (0.4 * A1 <= areas[i] <= 1.6 * A1) or info[i].get("solidity", 0) < sol_min or not SP._isolated(lab, i, sl, max_contact):
            continue
        sub = lab[sl] == i
        p, _ = SP.radial_profile(sub)
        rs.append(float(np.mean(np.abs(p / p.mean() - 1))))
        yy, xx = np.nonzero(sub); ang = np.linspace(0, np.pi, 18, endpoint=False)
        P = np.stack([xx, yy], 1) @ np.stack([np.cos(ang), np.sin(ang)]); w = P.max(0) - P.min(0) + 1
        asp.append(w.max() / w.min())
    return rs, asp


def fit_circle(x, y, w=None, iters=4):
    """Robust algebraic (Kåsa) circle fit with Tukey-like reweighting. Returns cx, cy, r, rms."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    w = np.ones_like(x) if w is None else np.asarray(w, float)
    cx = cy = r = np.nan
    for _ in range(iters):
        A = np.stack([x, y, np.ones_like(x)], 1) * w[:, None]
        b = (x * x + y * y) * w
        try:
            s, *_ = np.linalg.lstsq(A, b, rcond=None)
        except np.linalg.LinAlgError:
            break
        cx, cy = s[0] / 2, s[1] / 2; r = math.sqrt(max(s[2] + cx * cx + cy * cy, 1e-9))
        res = np.hypot(x - cx, y - cy) - r
        sc = 1.4826 * np.median(np.abs(res)) + 0.3
        w = np.clip(1 - (res / (3 * sc)) ** 2, 0, 1) ** 2
        if w.sum() < 5:
            break
    res = np.hypot(x - cx, y - cy) - r
    good = w > 0.05
    rms = float(np.sqrt(np.mean(res[good] ** 2))) if good.any() else float("inf")
    return float(cx), float(cy), float(r), rms


def hough(cs, m, rmin, rmax, edge_q=0.6):
    """Gradient-direction Hough accumulator of circle centres (votes summed over rmin..rmax)."""
    gy, gx = np.gradient(cs)
    mag = np.hypot(gx, gy)
    ring = ndi.binary_dilation(m, iterations=3) & ~ndi.binary_erosion(m, iterations=max(2, int(0.35 * rmin)))
    thr = np.quantile(mag[ring], edge_q) if ring.any() else 0
    yy, xx = np.nonzero(ring & (mag > thr))
    if not len(yy):
        return np.zeros(cs.shape)
    ux, uy = gx[yy, xx] / (mag[yy, xx] + 1e-12), gy[yy, xx] / (mag[yy, xx] + 1e-12)   # towards the particle (contrast up)
    wv = mag[yy, xx]
    H, W = cs.shape
    acc = np.zeros(H * W)
    for r in np.arange(rmin, rmax + 1, 1.0):
        cy = np.rint(yy + uy * r).astype(int); cx = np.rint(xx + ux * r).astype(int)
        ok = (cy >= 0) & (cy < H) & (cx >= 0) & (cx < W)
        acc += np.bincount(cy[ok] * W + cx[ok], weights=wv[ok], minlength=H * W)
    return ndi.gaussian_filter(acc.reshape(H, W), max(1.0, 0.06 * rmin))


def ray_edges(c, cx, cy, r0, core, level=LEVEL, n=N_RAYS, rmin_f=0.45, rmax_f=1.55, back_f=0.35):
    """Sub-pixel edge points on n rays from (cx, cy): first outward crossing of level·core after the core.
    Returns (ex, ey, free_mask[n], kinds[n]) with kind 0 = free edge, 1 = contact / overlap (no background reached),
    2 = outside the image."""
    H, W = c.shape
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    step = 0.33
    rr = np.arange(rmin_f * r0, rmax_f * r0, step)
    X = cx + np.outer(np.cos(th), rr); Y = cy + np.outer(np.sin(th), rr)
    out = (X < 0) | (X > W - 1) | (Y < 0) | (Y > H - 1)
    P = ndi.map_coordinates(c, [np.clip(Y, 0, H - 1), np.clip(X, 0, W - 1)], order=1, mode="nearest")
    lv = level * core
    ex = np.full(n, np.nan); ey = np.full(n, np.nan); kinds = np.ones(n, int)
    nb = int(back_f * r0 / step) + 1
    for k in range(n):
        p = P[k]
        below = p < lv
        j = np.argmax(below) if below.any() else -1
        if j <= 0:
            if out[k].any():
                kinds[k] = 2
            continue
        if out[k, : j + 1].any():
            kinds[k] = 2; continue
        seg = p[j: j + nb]
        if len(seg) < nb or out[k, j: j + nb].any():
            kinds[k] = 2 if out[k, j: j + nb].any() else 1
            if kinds[k] == 2:
                continue
        if (seg > lv).any() or seg.min() > 0.5 * lv:      # rises again (neighbour) or never reaches the background
            kinds[k] = 1; continue
        t = (p[j - 1] - lv) / ((p[j - 1] - p[j]) or 1e-9)
        rj = rr[j - 1] + t * step
        ex[k] = cx + math.cos(th[k]) * rj; ey[k] = cy + math.sin(th[k]) * rj; kinds[k] = 0
    return ex, ey, kinds == 0, kinds


def _core(c, cx, cy, r):
    H, W = c.shape
    y0, y1 = int(max(0, cy - r)), int(min(H, cy + r + 1)); x0, x1 = int(max(0, cx - r)), int(min(W, cx + r + 1))
    if y1 <= y0 or x1 <= x0:
        return 0.0
    yy, xx = np.mgrid[y0:y1, x0:x1]
    d = np.hypot(yy - cy, xx - cx) <= 0.4 * r
    v = c[y0:y1, x0:x1][d]
    return float(np.percentile(v, 60)) if v.size else 0.0


def refine(c, cx, cy, r0, level=LEVEL, n=N_RAYS, iters=3):
    """Ray-cast + robust circle fit around a candidate. Returns dict or None."""
    best = None
    for _ in range(iters):
        core = _core(c, cx, cy, r0)
        if core <= 0:
            return best
        ex, ey, free, kinds = ray_edges(c, cx, cy, r0, core, level, n)
        nf = int(free.sum())
        if nf < max(8, 0.15 * n):
            return best if best else {"cx": cx, "cy": cy, "r": r0, "rms": float("inf"), "arc": nf / n, "core": core,
                                      "n_free": nf, "contact": float((kinds == 1).mean()), "fitted": False}
        fx, fy, fr, rms = fit_circle(ex[free], ey[free])
        if not np.isfinite(fr) or fr < 0.4 * r0 or fr > 1.8 * r0 or math.hypot(fx - cx, fy - cy) > 0.6 * r0:
            return best
        best = {"cx": fx, "cy": fy, "r": fr, "rms": rms, "arc": nf / n, "core": core, "n_free": nf,
                "contact": float((kinds == 1).mean()), "fitted": True}
        if math.hypot(fx - cx, fy - cy) < 0.3 and abs(fr - r0) < 0.3:
            break
        cx, cy, r0 = fx, fy, fr
    return best


def confidence(cand, r_typ, cov):
    """0–1: visible-arc term × fit term × coverage term × size plausibility."""
    arc = min(1.0, cand["arc"] / 0.55)
    fit = math.exp(-((cand["rms"] / max(cand["r"], 1.0)) / 0.045) ** 2) if cand.get("fitted") else 0.25
    cv = float(np.clip((cov - 0.72) / 0.2, 0, 1))
    sz = math.exp(-(math.log(max(cand["r"], 1e-3) / r_typ) / 0.25) ** 2) if r_typ else 1.0
    return float(arc * fit * cv * (0.4 + 0.6 * sz))


def _disc(shape, cx, cy, r, sl=None):
    H, W = shape
    y0, y1 = int(max(0, math.floor(cy - r - 1))), int(min(H, math.ceil(cy + r + 2)))
    x0, x1 = int(max(0, math.floor(cx - r - 1))), int(min(W, math.ceil(cx + r + 2)))
    yy, xx = np.mgrid[y0:y1, x0:x1]
    return (slice(y0, y1), slice(x0, x1)), (np.hypot(yy - cy, xx - cx) <= r)


def resolve(c, m, L, level=LEVEL, delta_px=DELTA_PX, conf_min=CONF_MIN, n_rays=N_RAYS, min_new=0.3, kappa=0.32,
            lens_dark=1.06, max_clusters=2000):
    """Fit circles to every cluster of the mask. Returns (circles: list of dicts, unexplained mask)."""
    r_est = L / 2.0
    rmin, rmax = 0.6 * r_est, 1.4 * r_est
    cs = ndi.gaussian_filter(c, max(0.8, 0.03 * L))
    acc = hough(cs, m, rmin, rmax)
    mf = ndi.binary_fill_holes(m)
    lab, n = ndi.label(mf)
    sls = ndi.find_objects(lab)
    dist = ndi.distance_transform_edt(mf)
    cands = []
    # Hough peaks + distance-map maxima (the latter catches particles with weak edges)
    for src, arr, sz in (("hough", acc, int(0.5 * r_est) | 1), ("dist", dist, int(0.6 * r_est) | 1)):
        pk = (arr == ndi.maximum_filter(arr, size=max(3, sz))) & mf & (dist > 0.25 * r_est) & (arr > 0)
        yy, xx = np.nonzero(pk)
        if src == "hough" and len(yy):
            v = arr[yy, xx]; keep = v > 0.12 * np.percentile(v, 95)
            yy, xx = yy[keep], xx[keep]
        for y_, x_ in zip(yy, xx):
            r0 = float(np.clip(dist[y_, x_], 0.7 * r_est, 1.3 * r_est)) if src == "dist" else r_est
            cands.append((src, float(x_), float(y_), r0))
    # grid seeds over the clusters (every particle gets a start point even without a Hough / distance peak; the ray fit
    # converges to the nearest sphere and duplicates are merged below)
    gs = max(3, int(0.55 * r_est))
    yy, xx = np.mgrid[gs // 2: m.shape[0]: gs, gs // 2: m.shape[1]: gs]
    ok = mf[yy, xx] & (dist[yy, xx] > 0.3 * r_est)
    for y_, x_ in zip(yy[ok], xx[ok]):
        cands.append(("grid", float(x_), float(y_), r_est))
    # refine all candidates
    fitted = []
    seen = []
    for src, x_, y_, r0 in cands:
        f_ = refine(c, x_, y_, r0, level, n_rays)
        if not f_ or not np.isfinite(f_["r"]):
            continue
        if any(math.hypot(f_["cx"] - a, f_["cy"] - b) < 1.5 and abs(f_["r"] - rr) < 1.5 for a, b, rr in seen):
            continue
        seen.append((f_["cx"], f_["cy"], f_["r"]))
        f_["src"] = src
        fitted.append(f_)
    if not fitted:
        return [], mf
    rs = np.array([f_["r"] for f_ in fitted if f_.get("fitted") and f_["arc"] > 0.5])
    r_typ = float(np.median(rs)) if len(rs) >= 3 else r_est
    for f_ in fitted:
        sl, d = _disc(m.shape, f_["cx"], f_["cy"], f_["r"])
        cov = float(mf[sl][d].mean()) if d.any() else 0.0
        f_["cov"] = cov
        f_["conf0"] = confidence(f_, r_typ, cov)
    fitted.sort(key=lambda f_: -(f_["conf0"] + 0.05 * f_["arc"]))
    taken = np.zeros(m.shape, bool)
    out = []
    for f_ in fitted:
        if f_["r"] < 0.5 * r_typ or f_["r"] > 1.6 * r_typ or f_["cov"] < 0.6:
            continue
        cx, cy, r = f_["cx"], f_["cy"], f_["r"]
        if any(math.hypot(cx - o["cx"], cy - o["cy"]) < kappa * (r + o["r"]) for o in out):
            continue
        sl, d = _disc(m.shape, cx, cy, r)
        new = float((d & mf[sl] & ~taken[sl]).sum()) / max(1.0, d.sum())
        if new < min_new:
            continue
        # heavy overlap with an accepted circle: needs a darker lens (double thickness) to count as occlusion
        occl = []
        ok = True
        for o in out:
            dd = math.hypot(cx - o["cx"], cy - o["cy"])
            if dd < 0.92 * (r + o["r"]):
                sl2, d2 = _disc(m.shape, o["cx"], o["cy"], o["r"])
                full = np.zeros(m.shape, bool); full[sl] = d
                full2 = np.zeros(m.shape, bool); full2[sl2] = d2
                lens = full & full2
                if lens.sum() >= 6:
                    cl = float(np.median(cs[lens])); c1 = float(np.median(cs[full & ~full2])) if (full & ~full2).any() else cl
                    c2 = float(np.median(cs[full2 & ~full])) if (full2 & ~full).any() else cl
                    ratio = cl / max(1e-6, max(c1, c2))
                    if ratio >= lens_dark:
                        occl.append((o["id"], ratio))
                    elif dd < 0.55 * (r + o["r"]):
                        ok = False; break
        if not ok:
            continue
        f_ = dict(f_, id=len(out) + 1, occluded_by=occl)
        out.append(f_)
        taken[sl] |= d
    for o in out:
        conf = o["conf0"]
        if o["occluded_by"]:
            conf = min(1.0, conf * 1.1)
        o["conf"] = float(conf)
        o["r_edge"] = o["r"]
        o["r"] = max(1.0, o["r"] - delta_px)
        o["kind"] = "particle" if conf >= conf_min else "lowconf"
    rest = mf & ~ndi.binary_dilation(taken, iterations=2)
    return out, rest


def to_labels(circles, shape, rest=None, min_rest=0.35):
    """Label image: each pixel of a disc goes to the circle with the smallest d/r (overlap lenses split between the
    two spheres). Big unexplained remainders are returned as extra 'aggregate' labels."""
    lab = np.zeros(shape, int); best = np.full(shape, np.inf)
    for o in circles:
        sl, d = _disc(shape, o["cx"], o["cy"], o["r"])
        yy, xx = np.mgrid[sl]
        q = np.hypot(yy - o["cy"], xx - o["cx"]) / o["r"]
        upd = d & (q < best[sl])
        lab[sl][upd] = o["id"]; best[sl][upd] = q[upd]
    extra = []
    if rest is not None and circles:
        a1 = math.pi * float(np.median([o["r"] for o in circles])) ** 2
        rl, n = ndi.label(rest & (lab == 0))
        sz = np.bincount(rl.ravel())
        nxt = max(o["id"] for o in circles) + 1
        for j in range(1, n + 1):
            if sz[j] >= min_rest * a1:
                lab[rl == j] = nxt; extra.append(nxt); nxt += 1
    return lab, extra
