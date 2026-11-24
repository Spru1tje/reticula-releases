"""v2.21: elongated (rod / worm-like) particles on TEM micrographs: length and width per particle.

Method
  1. Contrast against the local film level: dark particles are darker than a Gaussian-smoothed neighbourhood
     (sigma `local_sigma` px, of the order of the rod width), normalised by that level.
  2. Hysteresis threshold on the contrast, at Otsu level t: seeds above 1.1 t, growth down to 0.7 t; holes filled.
  3. Objects smaller than `min_area_px` are dropped. Each object is measured from its pixel moments: the length is the
     full major-axis length and the width the full minor-axis length of the ellipse with the same second moments.
  4. Objects in the rim of the carbon-film holes (bright regions) are flagged 'film edge' and not counted.
  5. The scale comes from the scale bar (the longest bright horizontal bar in the bottom 12 %), or from `nm_per_px`.

Limits (see the methods text): no ground truth is available here, so accuracy is not validated; chains of touching rods
are measured as one object (flagged by length/width ratio); film texture can give small false objects.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import ndimage as ndi


def _solidity(xs, ys):
    """Area of the object divided by the area of its convex hull (1 = convex, low = irregular / clumped)."""
    from scipy.spatial import ConvexHull, QhullError
    pts = np.unique(np.column_stack([xs, ys]).astype(float), axis=0)
    if len(pts) < 3:
        return 1.0
    try:
        hull = ConvexHull(pts).volume
    except QhullError:
        return 1.0
    return float(min(1.0, len(xs) / max(hull, 1.0)))


def bar_nm_per_px(im, bar_nm):
    """nm per pixel from the longest bright horizontal bar in the bottom 12 % of the image (bar_nm: its label, e.g. 100)."""
    H, W = im.shape
    band = np.asarray(im)[int(H * 0.88):] > 235
    lab, _ = ndi.label(band)
    best = 0
    for sl in ndi.find_objects(lab):
        if sl is None:
            continue
        w = sl[1].stop - sl[1].start
        h = sl[0].stop - sl[0].start
        if w > 150 and h < 15:
            best = max(best, w)
    if best == 0:
        return None, 0
    return float(bar_nm) / best, best


def _otsu(x):
    from mfs.particles import otsu
    return float(otsu(x))


def segment_rods(im, nm_per_px, local_sigma=10.0, min_area_px=12, hole_level=230, edge_zone_px=10):
    """Returns (labels, objects). objects: list of dicts with id, area_px, length_nm, width_nm, aspect, centre, film_edge."""
    from mfs.particles import hysteresis
    a = np.asarray(im, float)
    loc = ndi.gaussian_filter(a, local_sigma)
    con = (loc - ndi.gaussian_filter(a, 1.0)) / (loc + 1e-6)
    t = _otsu(con)
    mask = ndi.binary_fill_holes(hysteresis(con, 0.7 * t, 1.1 * t))
    mask = ndi.binary_opening(mask, iterations=1)
    lab, n = ndi.label(mask)
    areas = np.bincount(lab.ravel())
    keep = np.zeros(n + 1, bool); keep[1:] = areas[1:] >= min_area_px
    lab = np.where(keep[lab], lab, 0)
    lab, n = ndi.label(lab > 0)
    # carbon-film holes (bright) and a band around them: the dark rim of the film is not a particle
    holes = ndi.gaussian_filter(a, 2.0) > hole_level
    holes = ndi.binary_opening(holes, iterations=3)
    rim = ndi.binary_dilation(holes, iterations=edge_zone_px) & ~holes
    objs = []
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        if sl is None:
            continue
        ys, xs = np.nonzero(lab[sl] == i)
        ys = ys + sl[0].start; xs = xs + sl[1].start
        A = ys.size
        cov = np.cov(np.vstack([xs, ys]))
        ev = np.sort(np.linalg.eigvalsh(cov))[::-1] if A > 2 else np.array([0.0, 0.0])
        major = 4.0 * math.sqrt(max(ev[0], 0.0))          # full length of the ellipse with the same second moments
        minor = 4.0 * math.sqrt(max(ev[1], 0.0))
        in_rim = float(rim[ys, xs].mean())
        sol = _solidity(xs, ys)
        objs.append({"id": i, "area_px": int(A), "length_nm": major * nm_per_px, "width_nm": max(minor, 1e-9) * nm_per_px,
                     "aspect": major / max(minor, 1e-9), "solidity": sol, "cy": float(ys.mean()), "cx": float(xs.mean()),
                     "film_edge": bool(in_rim > 0.5)})
    return lab, objs


def flag_agglomerates(objs, nm_per_px, min_big_nm2=150.0, rel_big=6.0):
    """v2.21: mark agglomerates (clumps of particles) so they are not counted as rods.
    An object is an agglomerate when it is large AND compact (aspect < 2.5 or solidity < 0.85), or when it is very
    irregular (solidity < 0.6). 'Large' = more than max(150 nm^2, 6 x the median area of the image), so the rule adapts
    to the particle size of each micrograph. Single rods are long and thin (aspect >= 2.5) and stay in."""
    if not objs:
        return objs
    areas = np.array([o["area_px"] * nm_per_px ** 2 for o in objs])
    big_nm2 = max(min_big_nm2, rel_big * float(np.median(areas)))
    for o, a in zip(objs, areas):
        large_compact = a > big_nm2 and (o["aspect"] < 2.5 or o["solidity"] < 0.85)
        o["agglomerate"] = bool(large_compact or o["solidity"] < 0.6)
    return objs


def summary(objs):
    """Median and inter-quartile range of length and width (only objects that are not at the film edge)."""
    good = [o for o in objs if not o["film_edge"] and not o.get("agglomerate", False)]
    if not good:
        return {"n": 0}
    L = np.array([o["length_nm"] for o in good]); W = np.array([o["width_nm"] for o in good])
    return {"n": len(good), "n_film_edge": sum(o["film_edge"] for o in objs),
            "n_agglomerate": sum(o.get("agglomerate", False) for o in objs),
            "length_nm_median": float(np.median(L)), "length_nm_q25": float(np.percentile(L, 25)), "length_nm_q75": float(np.percentile(L, 75)),
            "width_nm_median": float(np.median(W)), "width_nm_q25": float(np.percentile(W, 25)), "width_nm_q75": float(np.percentile(W, 75)),
            "aspect_median": float(np.median([o["aspect"] for o in good]))}


def overlay(im, lab, objs, path, crop=None):
    """Save a picture of the micrograph with the accepted objects (red) and the film-edge objects (blue) for checking."""
    from PIL import Image
    edge = (lab > 0) & (lab != ndi.grey_erosion(lab, size=3))
    rgb = np.stack([np.asarray(im)] * 3, -1).astype(np.uint8)
    film = {o["id"] for o in objs if o["film_edge"]}
    rgb[edge & ~np.isin(lab, list(film))] = [255, 0, 0]
    rgb[edge & np.isin(lab, list(film))] = [0, 120, 255]
    agg = {o["id"] for o in objs if o.get("agglomerate")}
    rgb[edge & np.isin(lab, list(agg))] = [0, 200, 0]          # agglomerates (not counted) in green
    if crop is not None:
        y0, y1, x0, x1 = crop
        rgb = rgb[y0:y1, x0:x1]
    Image.fromarray(rgb).save(path)
    return path
