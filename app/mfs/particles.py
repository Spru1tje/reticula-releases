"""Particle sizing on SEM/TEM images (no scikit-image needed): Otsu threshold, morphological clean-up, connected
components, optional splitting of touching particles (distance-transform markers + nearest-marker assignment inside
each blob – a fast watershed-like split), per-particle size; log-normal fit of the size distribution."""
from __future__ import annotations
import math
import numpy as np
from scipy import ndimage as ndi, stats


def otsu(img):
    h, e = np.histogram(img[np.isfinite(img)], 256)
    c = (e[:-1] + e[1:]) / 2; w = h.cumsum(); w2 = w[-1] - w
    m = (h * c).cumsum(); mt = m[-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        var = (mt * w / w[-1] - m) ** 2 / (w * w2)
    var[~np.isfinite(var)] = 0
    return float(c[int(np.argmax(var))])


def _downsample(img, max_px):
    f = max(1, int(math.ceil(max(img.shape) / max_px)))
    if f == 1:
        return img.astype(float), 1
    h, w = (img.shape[0] // f) * f, (img.shape[1] // f) * f
    return img[:h, :w].astype(float).reshape(h // f, f, w // f, f).mean(axis=(1, 3)), f


def segment(img, dark_particles=True, threshold=None, smooth=1.0, open_px=1, split=True, min_dist_px=None,
            clear_border=True, max_px=1024):
    """Returns labels (on the working grid), the down-sampling factor and the threshold used."""
    a, f = _downsample(np.asarray(img), max_px)
    if a.ndim == 3:
        a = a.mean(axis=2)
    if smooth:
        a = ndi.gaussian_filter(a, smooth)
    t = otsu(a) if threshold is None else threshold
    m = a < t if dark_particles else a > t
    if open_px:
        m = ndi.binary_opening(m, iterations=int(open_px))
    m = ndi.binary_fill_holes(m)
    lab, n = ndi.label(m)
    if split and n:
        d = ndi.distance_transform_edt(m)
        md = min_dist_px or max(2, int(np.percentile(d[m], 90) * 0.8)) if m.any() else 3
        peaks = (d == ndi.maximum_filter(d, size=2 * int(md) + 1)) & (d > 0.5 * md)
        mk, nm = ndi.label(ndi.binary_dilation(peaks, iterations=1))
        if nm:
            _, idx = ndi.distance_transform_edt(mk == 0, return_indices=True)
            near = mk[idx[0], idx[1]]
            # keep each pixel inside its own blob: markers from another blob are not allowed
            blob_of_marker = ndi.maximum(lab, mk, index=np.arange(1, nm + 1))
            ok = blob_of_marker[np.clip(near - 1, 0, nm - 1)] == lab
            out = np.where(m & ok & (near > 0), near, 0)
            # blobs without a marker keep their own label
            nomk = m & (out == 0)
            if nomk.any():
                l2, n2 = ndi.label(nomk); out = np.where(nomk, l2 + nm, out)
            lab = _merge_slivers(out, lab_blob=lab)
    if clear_border:
        edge = np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))
        lab = np.where(np.isin(lab, edge[edge > 0]), 0, lab)
    lab, n = ndi.label(lab > 0) if not split else (lab, int(lab.max()))
    return lab, f, t


def _merge_slivers(lab, lab_blob, neck=0.65):
    """Undo over-splitting: two pieces of one blob are merged again when the line they share is not a neck, i.e. longer
    than `neck` × the smaller piece's equivalent diameter (true touching particles meet at a short neck)."""
    lab = lab.copy()
    for bi, sl in enumerate(ndi.find_objects(lab_blob), start=1):
        if sl is None:
            continue
        sl = tuple(slice(max(0, s.start - 1), s.stop + 1) for s in sl)
        sub = lab[sl]; inb = lab_blob[sl] == bi
        changed = True
        while changed:
            changed = False
            ids = [i for i in np.unique(sub[inb]) if i > 0]
            if len(ids) < 2:
                break
            area = {i: int((sub == i).sum()) for i in ids}
            for i in sorted(ids, key=area.get):
                ring = ndi.binary_dilation(sub == i) & (sub != i) & (sub > 0) & inb
                if not ring.any():
                    continue
                nb, cnt = np.unique(sub[ring], return_counts=True)
                j, L = int(nb[np.argmax(cnt)]), int(cnt.max())
                dmin = 2 * math.sqrt(min(area[i], area[j]) / math.pi)
                if L > neck * dmin:
                    sub[sub == i] = j; changed = True; break
        lab[sl] = sub
    return lab


def measure(lab, pixel_nm, factor=1, min_area_nm2=0.0, max_area_nm2=None, info=None):
    """Per particle: area, equivalent circle diameter, max Feret, min Feret (caliper, 36 angles), aspect ratio.
    v2.12: particles with a fitted circle (info[i]['d_px'], sphere model) report the circle: ECD = Feret = fitted
    diameter (an occluded sphere is not shrunk to its visible part), plus 'confidence'."""
    px = pixel_nm * factor
    ids = [i for i in np.unique(lab) if i > 0]
    sl = ndi.find_objects(lab)
    rows = []
    ang = np.linspace(0, np.pi, 36, endpoint=False); U = np.stack([np.cos(ang), np.sin(ang)])
    for i in ids:
        s_ = sl[i - 1] if i - 1 < len(sl) else None
        if s_ is None:
            continue
        sub = lab[s_] == i
        A = sub.sum() * px * px
        if A < min_area_nm2 or (max_area_nm2 and A > max_area_nm2) or sub.sum() < 5:
            continue
        yy, xx = np.nonzero(sub)
        P = np.stack([xx, yy]).T.astype(float) @ U
        w = P.max(0) - P.min(0) + 1
        ci = info.get(int(i)) if isinstance(info, dict) else None
        if isinstance(ci, dict) and ci.get("d_px"):
            D_ = ci["d_px"] * px; cx_, cy_, _r = ci["circle"]
            rows.append({"id": int(i), "area (nm²)": math.pi * D_ * D_ / 4, "ECD (nm)": D_, "Feret max (nm)": D_,
                         "Feret min (nm)": D_, "aspect ratio": 1.0, "x (px)": cx_ * factor, "y (px)": cy_ * factor,
                         "confidence": float(ci.get("conf", 1.0))})
            continue
        rows.append({"id": int(i), "area (nm²)": A, "ECD (nm)": 2 * math.sqrt(A / math.pi), "Feret max (nm)": w.max() * px,
                     "Feret min (nm)": w.min() * px, "aspect ratio": w.max() / w.min(),
                     "x (px)": (xx.mean() + s_[1].start) * factor, "y (px)": (yy.mean() + s_[0].start) * factor})
    return rows


def lognormal_fit(d):
    d = np.asarray(d, float); d = d[np.isfinite(d) & (d > 0)]
    if len(d) < 3:
        raise ValueError("need at least 3 particles")
    s, loc, sc = stats.lognorm.fit(d, floc=0)
    mu = math.log(sc)
    return {"n": len(d), "mean": float(d.mean()), "sd": float(d.std(ddof=1)), "median": float(np.median(d)),
            "mu": mu, "sigma": s, "geo_mean": sc, "geo_sd": math.exp(s), "mode": math.exp(mu - s * s),
            "pdf": lambda x: stats.lognorm.pdf(x, s, 0, sc),
            "ks_p": float(stats.kstest(d, "lognorm", args=(s, 0, sc)).pvalue)}


def draw_hist(ax, d, fit=None, bins="auto", color="#0072B2", label=None, unit="nm", show_stats=True):
    d = np.asarray(d, float)
    n_, b_, _ = ax.hist(d, bins=bins, color=color, alpha=0.55, edgecolor="white", lw=0.4, label=label)
    if fit:
        xx = np.linspace(max(1e-9, b_[0] * 0.8), b_[-1] * 1.1, 300)
        ax.plot(xx, fit["pdf"](xx) * len(d) * (b_[1] - b_[0]), color=color, lw=1.0)
        if show_stats:
            ax.text(0.97, 0.95, f"{fit['mean']:.0f} ± {fit['sd']:.0f} {unit}\n(n = {fit['n']})\nlog-normal: median {fit['geo_mean']:.0f} {unit}, "
                    f"σ$_g$ {fit['geo_sd']:.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=6)
    ax.set_xlabel(f"Particle size ({unit})"); ax.set_ylabel("Count")


def overlay_outlines(ax, img, lab, factor=1, color="#FF3B30"):
    ax.imshow(img, cmap="gray"); edge = lab ^ ndi.binary_erosion(lab > 0) * lab
    e = (edge > 0)
    if factor > 1:
        e = np.kron(e, np.ones((factor, factor), bool))[: img.shape[0], : img.shape[1]]
    rgba = np.zeros(e.shape + (4,)); rgba[e] = list(int(color[i:i + 2], 16) / 255 for i in (1, 3, 5)) + [1]
    ax.imshow(rgba, interpolation="nearest"); ax.set_axis_off()


# ============================================================== v2.8: marker watershed, aggregates, support holes, edits
def _solidity(mask):
    yy, xx = np.nonzero(mask)
    if len(yy) < 8:
        return 1.0
    try:
        from scipy.spatial import ConvexHull
        pts = np.stack([xx, yy], 1).astype(float)
        pts = np.vstack([pts, pts + [1, 0], pts + [0, 1], pts + [1, 1]])        # pixel corners
        return float(min(1.0, len(yy) / ConvexHull(pts).volume))
    except Exception:
        return 1.0


def _peaks(d, mask, min_dist, h):
    """Markers = h-maxima of the distance map (maxima rising more than h above every saddle to a higher peak; grey-scale
    reconstruction by dilation), then markers closer than min_dist merged."""
    d = np.where(mask, d, 0.0)
    rec = np.clip(d - h, 0, None)
    for _ in range(int(d.max() / max(h, 0.5) * 4) + 50):          # reconstruction by dilation of (d − h) under d
        nxt = np.minimum(ndi.grey_dilation(rec, size=3), d)
        if np.array_equal(nxt, rec):
            break
        rec = nxt
    pk = (d - rec >= h - 1e-6) & mask & (d > 0)
    mk, n = ndi.label(pk, structure=np.ones((3, 3)))
    if n > 1 and min_dist > 1:
        cy, cx = np.array(ndi.center_of_mass(pk, mk, range(1, n + 1))).T
        from scipy.spatial import cKDTree
        par = list(range(n + 1))

        def f(a_):
            while par[a_] != a_:
                par[a_] = par[par[a_]]; a_ = par[a_]
            return a_
        for i, j in cKDTree(np.stack([cx, cy], 1)).query_pairs(min_dist):
            par[f(i + 1)] = f(j + 1)
        root = np.array([f(i) for i in range(n + 1)])
        mk = root[mk]; mk[~pk] = 0
        u = np.unique(mk[mk > 0]); remap = np.zeros(mk.max() + 1, int); remap[u] = np.arange(1, len(u) + 1)
        mk = remap[mk]; n = len(u)
    return mk, n


def estimate_size_px(mask):
    """Typical particle radius (px): median inscribed radius of the distance-map peaks (robust to clusters)."""
    d = ndi.distance_transform_edt(mask)
    if not mask.any():
        return 5.0
    r = 3.0
    for _ in range(3):                    # iterate: the peak search window follows the size estimate
        pk = (d == ndi.maximum_filter(d, size=2 * int(r) + 1)) & (d > max(2.0, 0.5 * r))
        v = d[pk]
        if not len(v):
            break
        r = float(np.median(v))
    return r


def _ws(cost, markers, mask):
    """Marker-controlled watershed (priority flood, 4-connected) restricted to `mask`; lower cost floods first."""
    import heapq
    H, W = mask.shape
    out = np.where(mask, markers, 0).astype(np.int64).ravel()
    c = np.asarray(cost, float).ravel(); mk = mask.ravel()
    heap = []
    seeds = np.nonzero(out > 0)[0]
    cnt = 0
    for p in seeds:                           # push the unlabelled neighbours of every seed pixel
        heap.append((c[p], cnt, p)); cnt += 1
    heapq.heapify(heap)
    done = out > 0
    while heap:
        _, _, p = heapq.heappop(heap)
        lp = out[p]; y, x = divmod(p, W)
        for q, ok in ((p - W, y > 0), (p + W, y < H - 1), (p - 1, x > 0), (p + 1, x < W - 1)):
            if ok and mk[q] and not done[q]:
                done[q] = True; out[q] = lp
                heapq.heappush(heap, (c[q], cnt, q)); cnt += 1
    return out.reshape(H, W).astype(int)


def segment_ws(img, dark_particles=True, threshold=None, smooth=1.0, open_px=1, min_dist_factor=0.7, h_factor=0.1,
               method="distance", gradient_weight=0.5, exclude_holes=True, clear_border=True, solidity_min=0.85,
               area_factor=2.0, max_px=1024):
    """Particles on a support (SEM/TEM): returns (labels, factor, threshold, info) with info = {id: {'kind':
    'particle'|'aggregate', 'solidity': s}, '_holes': n, '_r_est_px': r, '_min_dist_px': md}.
    1) threshold (Otsu; with dark particles the bright support holes are removed first and Otsu is taken on the rest),
    2) distance-transform watershed with h-maxima markers, min-distance = factor × the auto-estimated radius,
    method 'gradient' / 'combined' floods the Sobel gradient (facet edges) instead of / mixed with the distance map,
    3) regions with solidity < solidity_min or area > area_factor × median are re-split with half the min-distance;
    what stays non-convex is classed 'aggregate', 4) edge-touching objects removed."""
    a, f = _downsample(np.asarray(img), max_px)
    if a.ndim == 3:
        a = a.mean(axis=2)
    if smooth:
        a = ndi.gaussian_filter(a, smooth)
    holes = np.zeros(a.shape, bool); nh = 0
    if exclude_holes and dark_particles:
        th = otsu(a)
        hm = a > th + 0.5 * (a.max() - th) * 0.5
        hm = ndi.binary_opening(hm, iterations=2)
        lh, nh0 = ndi.label(hm)
        for i, sl in enumerate(ndi.find_objects(lh), start=1):     # round bright blobs = holes in the carbon film
            sub = lh[sl] == i
            if sub.sum() > 30 and _solidity(sub) > 0.8:
                holes[sl] |= sub; nh += 1
        holes = ndi.binary_dilation(holes, iterations=2)
    pool = a[~holes]
    t = (otsu(pool) if pool.size else otsu(a)) if threshold is None else threshold
    m = (a < t) if dark_particles else (a > t)
    m &= ~holes
    if open_px:
        m = ndi.binary_opening(m, iterations=int(open_px))
    m = ndi.binary_fill_holes(m)
    r = estimate_size_px(m); md = max(2.0, min_dist_factor * r)
    g = None
    if method in ("gradient", "combined"):
        g = np.hypot(ndi.sobel(a, 0), ndi.sobel(a, 1))
        g = g / (np.percentile(g[m], 99) + 1e-9) if m.any() else g

    # 2) whole-image watershed
    d = ndi.gaussian_filter(ndi.distance_transform_edt(m), max(0.5, 0.15 * r)) * m
    mk, n = _peaks(d, m, md, h_factor * r)
    dn = d / (d.max() + 1e-9)
    cost = -dn if g is None else (g if method == "gradient" else gradient_weight * g - (1 - gradient_weight) * dn)
    lab = _ws(cost, mk, m) if n else ndi.label(m)[0]
    rest = m & (lab == 0)
    if rest.any():
        l2, _ = ndi.label(rest); lab = np.where(rest, l2 + lab.max(), lab)
    lab = _merge_small(lab, 0.2 * math.pi * r * r)
    # 3) non-convex or oversized regions: local re-split with a smaller min-distance; else 'aggregate'
    areas = np.bincount(lab.ravel()).astype(float); areas[0] = 0
    med = float(np.median(areas[areas > 0])) if np.any(areas > 0) else 1.0
    info = {}
    nxt = int(lab.max()) + 1
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        if sl is None:
            continue
        sub = lab[sl] == i
        s_ = _solidity(sub)
        big = sub.sum() > area_factor * med
        if s_ >= solidity_min and sub.sum() <= 0.72 * area_factor * med:
            info[i] = {"kind": "particle", "solidity": s_}; continue
        dd = ndi.distance_transform_edt(sub)
        ok = False
        for fct, hf in ((0.6, 0.06), (0.45, 0.03)):
            mk2, n2 = _peaks(dd, sub, max(2.0, md * fct), hf * r)
            if n2 < 2:
                continue
            c2 = -dd / (dd.max() + 1e-9) if g is None else cost[sl]
            l2 = _ws(c2, mk2, sub)
            l2 = _merge_small(l2, 0.2 * math.pi * r * r)
            ids = [j for j in np.unique(l2) if j > 0]
            sols = [_solidity(l2 == j) for j in ids]
            areas2 = [(l2 == j).sum() for j in ids]
            if len(ids) > 1 and min(areas2) > 0.3 * med and (np.mean(sols) > s_ + 0.02 or big):
                blk = lab[sl]
                for j, sj in zip(ids, sols):
                    blk[l2 == j] = nxt
                    info[nxt] = {"kind": "particle" if sj >= solidity_min - 0.05 else "aggregate", "solidity": sj, "split": "re-split"}
                    nxt += 1
                ok = True; break
        if not ok:
            info[i] = {"kind": "aggregate" if (s_ < solidity_min or big) else "particle", "solidity": s_}
    if clear_border:
        edge = np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))
        lab = np.where(np.isin(lab, edge[edge > 0]), 0, lab)
        info = {k: v for k, v in info.items() if k not in set(edge.tolist())}
    present = set(np.unique(lab).tolist()) - {0}
    info = {k: v for k, v in info.items() if k in present}
    for k in present - set(info):
        info[k] = {"kind": "particle", "solidity": 1.0}
    info.update({"_holes": nh, "_r_est_px": r, "_min_dist_px": md, "_hole_mask": holes})
    return lab, f, t, info


def _merge_small(lab, min_area):
    """Pieces smaller than min_area (px) join the neighbour they share the longest boundary with (faceted particles can
    touch along a whole face, so the neck criterion of v2.6 is not used here)."""
    lab = lab.copy()
    for _ in range(3):
        areas = np.bincount(lab.ravel())
        small = [i for i in np.nonzero(areas)[0] if i > 0 and areas[i] < min_area]
        if not small:
            break
        for i in small:
            m = lab == i
            if not m.any():
                continue
            ring = ndi.binary_dilation(m) & ~m & (lab > 0)
            if ring.any():
                nb, cnt = np.unique(lab[ring], return_counts=True); lab[m] = nb[np.argmax(cnt)]
    return lab


def apply_edits(lab, edits):
    """Manual corrections on the working grid. edits: list of ('delete', (x, y)) | ('merge', (x1, y1), (x2, y2)) |
    ('split', (x1, y1), (x2, y2)) [a cut line through one particle] | ('add', (x, y), r) [a disc]."""
    lab = lab.copy()
    H, W = lab.shape
    at = lambda p: int(lab[min(H - 1, max(0, int(p[1]))), min(W - 1, max(0, int(p[0])))])
    for e in edits:
        kind = e[0]
        if kind == "delete":
            i = at(e[1])
            if i:
                lab[lab == i] = 0
        elif kind == "merge":
            i, j = at(e[1]), at(e[2])
            if i and j and i != j:
                lab[lab == j] = i
        elif kind == "split":
            (x1, y1), (x2, y2) = e[1], e[2]
            n = int(max(abs(x2 - x1), abs(y2 - y1)) * 2) + 2
            xs = np.linspace(x1, x2, n); ys = np.linspace(y1, y2, n)
            line = np.zeros_like(lab, bool)
            line[np.clip(ys.round().astype(int), 0, H - 1), np.clip(xs.round().astype(int), 0, W - 1)] = True
            line = ndi.binary_dilation(line)
            ids = [i for i in np.unique(lab[line]) if i > 0]
            for i in ids:
                reg = (lab == i) & ~line
                l2, n2 = ndi.label(reg)
                if n2 > 1:
                    nxt = int(lab.max()) + 1
                    lab[lab == i] = 0
                    for j in range(1, n2 + 1):
                        lab[l2 == j] = i if j == 1 else nxt + j - 2
                    # give the cut pixels back to the nearest piece
                    cut = line & (lab == 0) & (ndi.binary_dilation(reg))
                    if cut.any():
                        _, idx = ndi.distance_transform_edt(lab == 0, return_indices=True)
                        lab[cut] = lab[idx[0][cut], idx[1][cut]]
        elif kind == "add":
            (x, y), rr = e[1], e[2]
            yy, xx = np.ogrid[:H, :W]
            lab[((xx - x) ** 2 + (yy - y) ** 2 <= rr * rr) & (lab == 0)] = int(lab.max()) + 1
    return lab


TEAL, RED = (0.078, 0.722, 0.651), (0.898, 0.282, 0.302)
ORANGE = (0.976, 0.451, 0.086)            # v2.12: low-confidence sphere fit (excluded by default)


def overlay_labels(ax, img, lab, factor=1, info=None, numbers=False, fs=4, alpha=0.13, seed=3, excluded=None, lw_px=None):
    """v2.10 clean overlay: thin teal outline = included particle, dotted teal = split by the shape prior / sphere model,
    red dashed outline + light red fill = aggregate / excluded, orange dashed = low-confidence sphere fit (v2.12);
    numbers optional (off by default). img on the full grid."""
    ax.imshow(img, cmap="gray")
    n = int(lab.max())
    if n == 0:
        ax.set_axis_off(); return
    L = np.kron(lab, np.ones((factor, factor), int))[: img.shape[0], : img.shape[1]] if factor > 1 else lab
    info = info or {}
    ex = set(int(i) for i in (excluded or []))
    lc_ids = [k for k, v in info.items() if isinstance(v, dict) and v.get("kind") == "lowconf"]
    agg_ids = [k for k, v in info.items() if isinstance(v, dict) and v.get("kind") == "aggregate"] + [i for i in ex if i not in lc_ids]
    pri_ids = [k for k, v in info.items() if isinstance(v, dict) and v.get("split") in ("prior", "circle")]
    lcm = np.isin(L, lc_ids) & (L > 0)
    agg = np.isin(L, agg_ids) & (L > 0) & ~lcm; pri = np.isin(L, pri_ids) & ~agg & ~lcm
    w = lw_px or max(1, int(round(max(img.shape[:2]) / 1400)))
    edge = (L != ndi.grey_erosion(L, size=2 * w + 1)) & (L > 0)
    yy, xx = np.indices(L.shape)
    dash = ((xx + yy) // (6 * w)) % 2 == 0; dot = ((xx + yy) // (2 * w)) % 2 == 0
    rgba = np.zeros(L.shape + (4,))
    rgba[agg] = list(RED) + [alpha]
    rgba[lcm & ~edge] = list(ORANGE) + [alpha * 0.8]
    rgba[edge & lcm & dash] = list(ORANGE) + [0.95]
    rgba[edge & ~agg & ~pri & ~lcm] = list(TEAL) + [0.95]
    rgba[edge & pri & dot] = list(TEAL) + [0.95]
    rgba[edge & agg & dash] = list(RED) + [0.95]
    ax.imshow(rgba, interpolation="nearest")
    if numbers:
        for i, (cy, cx) in zip(range(1, n + 1), ndi.center_of_mass(L > 0, L, range(1, n + 1))):
            if np.isfinite(cx):
                kind = info.get(i, {}).get("kind", "particle") if isinstance(info.get(i), dict) else "particle"
                ax.text(cx, cy, str(i), color="white" if kind != "aggregate" else "#FFD1CF",
                        fontsize=fs, ha="center", va="center")
    ax.set_axis_off()


# ============================================================== v2.10: adaptive segmentation for uneven TEM/SEM contrast
def _nms(pts, score, dmin):
    """Greedy non-maximum suppression (strongest first, no chaining): indices of the kept points."""
    if not len(pts):
        return np.zeros(0, int)
    from scipy.spatial import cKDTree
    order = np.argsort(-np.asarray(score)); tree = cKDTree(pts)
    dead = np.zeros(len(pts), bool); keep = []
    for i in order:
        if dead[i]:
            continue
        keep.append(i)
        for j in tree.query_ball_point(pts[i], dmin):
            dead[j] = True
    return np.array(keep, int)


def flatten_background(x, fg_mask, sigma):
    """Smooth background from the pixels outside fg_mask (normalised Gaussian on a ¼ grid, so large clusters and
    uneven illumination do not bias it)."""
    k = max(1, int(sigma // 16))
    xs, ws = x[::k, ::k], (~fg_mask[::k, ::k]).astype(float)
    num = ndi.gaussian_filter(xs * ws, sigma / k); den = ndi.gaussian_filter(ws, sigma / k)
    bg = num / np.maximum(den, 1e-3)
    if den.max() < 1e-3:
        bg = np.full(xs.shape, float(np.median(x)))
    if k > 1:
        bg = ndi.zoom(bg, (x.shape[0] / bg.shape[0], x.shape[1] / bg.shape[1]), order=1)[: x.shape[0], : x.shape[1]]
    return bg


def contrast_map(a, dark_particles=True, smooth=1.0, iters=2):
    """Background-flattened particle contrast (0 = background, positive = particle), the noise σ of the background and
    the foreground mask by hysteresis: seeds above hi = max(Otsu, 8σ) grown down to lo = max(low_k·σ, 0.3·hi).
    Catches light-grey particles next to dark ones on an unevenly lit support."""
    x = ndi.gaussian_filter(np.asarray(a, float), smooth) if smooth else np.asarray(a, float)
    if not dark_particles:
        x = x.max() + x.min() - x
    fg = x < otsu(x)
    sig = max(x.shape) / 8.0
    for _ in range(iters):
        bgm = ndi.binary_dilation(fg, iterations=4)
        bg = flatten_background(x, bgm, sig)
        c = (bg - x) / np.maximum(np.abs(bg), 1.0)
        cb = c[~bgm] if (~bgm).sum() > 100 else c.ravel()
        med = float(np.median(cb)); sn = float(1.4826 * np.median(np.abs(cb - med))) or 1e-3
        fg = c > max(otsu(c), 8 * sn)                      # provisional mask for the next background estimate
    return c, sn


def hysteresis(c, lo, hi):
    lab, n = ndi.label(c > lo)
    keep = np.zeros(n + 1, bool); keep[np.unique(lab[c > hi])] = True; keep[0] = False
    return keep[lab]


def estimate_side_px(m):
    """Typical particle width (px) from the strongest distance-map maxima (one per particle via NMS)."""
    d = ndi.distance_transform_edt(m)
    pk = (d == ndi.maximum_filter(d, size=7)) & (d > 2.5)
    yy, xx = np.nonzero(pk)
    if not len(yy):
        return 10.0
    keep = _nms(np.stack([yy, xx], 1).astype(float), d[pk], 6)
    v = d[yy[keep], xx[keep]]
    v = v[v >= 0.3 * np.percentile(v, 95)]             # v2.12: ignore specks / film edges when the particles are large
    return float(2 * np.median(v))


def _pair_geometry(mask):
    """(solidity, aspect ratio from the principal axes) of a boolean mask."""
    yy, xx = np.nonzero(mask)
    if len(yy) < 6:
        return 0.0, 1.0
    P = np.stack([xx, yy], 1).astype(float)
    w, _ = np.linalg.eigh(np.cov((P - P.mean(0)).T))
    return _solidity(mask), float(math.sqrt(max(w[1], 1e-9) / max(w[0], 1e-9)))


def elongated_population(m, L, min_n=3, aspect_min=1.8):
    """v2.14: does the sample consist of ELONGATED particles (rods)? Looks at the isolated, convex single particles of the
    foreground mask (connected components of 0.4-10 L^2 with solidity >= 0.9): median aspect ratio of their principal axes.
    Squares / cubes / spheres give ~1, rods > 1.8. Returns {'elongated', 'aspect', 'n'}."""
    lm, n = ndi.label(m)
    asp = []
    A1 = L * L
    for i, sl in enumerate(ndi.find_objects(lm), start=1):
        if sl is None:
            continue
        sub_ = lm[sl] == i
        A = int(sub_.sum())
        if not (0.4 * A1 <= A <= 10 * A1):
            continue
        sol, a = _pair_geometry(sub_)
        if sol >= 0.9:
            asp.append(a)
    med = float(np.median(asp)) if len(asp) >= min_n else None
    return {"elongated": bool(med is not None and med >= aspect_min), "aspect": med, "n": len(asp)}


def merge_elongated(lab, c, L, solidity_min=0.93, aspect_min=1.3, saddle_min=0.7, max_area_f=14.0):
    """v2.14: an elongated particle (rod, platelet seen edge-on) is cut into pieces by the markers that are placed at the scale of the
    particle WIDTH. Two touching segments are merged when (1) their union is convex (solidity >= solidity_min) and clearly elongated
    (aspect >= aspect_min) - two touching round particles are never convex - and (2) the contrast along their interface shows no valley
    (>= saddle_min of the weaker core): a gap or overlap line between two particles has one. Repeated until nothing changes."""
    lab = lab.copy()
    A_max = max_area_f * L * L
    for _ in range(12):
        objs = ndi.find_objects(lab)
        ids = [i for i, sl in enumerate(objs, start=1) if sl is not None]
        if len(ids) < 2:
            break
        areas = np.bincount(lab.ravel(), minlength=lab.max() + 1)
        cores = {}
        best = None
        for i in ids:
            sl = objs[i - 1]
            pad = 2
            sl2 = (slice(max(0, sl[0].start - pad), sl[0].stop + pad), slice(max(0, sl[1].start - pad), sl[1].stop + pad))
            sub_i = lab[sl2] == i
            ring = ndi.binary_dilation(sub_i, iterations=1) & ~sub_i
            for j in np.unique(lab[sl2][ring]):
                if j <= i:
                    continue
                inter = ring & (lab[sl2] == j)
                n_int = int(inter.sum())
                if n_int < 3:
                    continue
                sub_j = lab[sl2] == j
                union = sub_i | sub_j
                if areas[i] + areas[j] > A_max:
                    continue
                if n_int < 0.3 * math.sqrt(min(areas[i], areas[j])):        # only a point contact: not one particle
                    continue
                sol, asp = _pair_geometry(union)
                if sol < solidity_min or asp < aspect_min:
                    continue
                cc = c[sl2]
                for k, sm in ((i, sub_i), (j, sub_j)):
                    if k not in cores:
                        dd = ndi.distance_transform_edt(sm)
                        inner = dd >= max(1.0, 0.4 * dd.max())
                        cores[k] = float(np.percentile(cc[inner], 60)) if inner.any() else 0.0
                saddle = float(np.percentile(cc[inter | (ndi.binary_dilation(inter) & union)], 50))
                if min(cores[i], cores[j]) <= 0 or saddle < saddle_min * min(cores[i], cores[j]):
                    continue
                if best is None or sol > best[0]:
                    best = (sol, i, j)
        if best is None:
            break
        _, i, j = best
        lab[lab == j] = i
    return lab


def refine_edges(lab, c, level=None, delta_px=None, grow_f=0.35, min_grow_px=4):
    """v2.14: put the outline of NON-round particles on the same calibrated edge as the sphere model does for round ones.
    The thresholded mask of the segmentation sits too far out (blur and halo): on the SiO2-NH2 TEM series the mask-based
    equivalent diameter was 8-18 % larger than the circle-fit one. Per particle: core contrast = 60th percentile of the contrast
    inside its inner part, edge = where the contrast falls to `level` x core (iso-contour, inside the particle's own
    Voronoi cell so neighbours do not merge), moved inward by `delta_px` (camera blur) - the constants of mfs.spheres.
    Returns a new label image (same ids); particles that would vanish keep their old outline."""
    from mfs import spheres as SPH
    level = SPH.LEVEL if level is None else float(level)
    delta_px = SPH.DELTA_PX if delta_px is None else float(delta_px)
    H, W = lab.shape
    if not (lab > 0).any():
        return lab
    _, (iy, ix) = ndi.distance_transform_edt(lab == 0, return_indices=True)
    near = lab[iy, ix]                                         # nearest particle of every pixel (Voronoi cells)
    out = np.zeros_like(lab)
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        if sl is None:
            continue
        r = math.sqrt(max(1.0, float((lab[sl] == i).sum())) / math.pi)
        pad = int(max(min_grow_px, grow_f * r)) + 2
        y0, y1 = max(0, sl[0].start - pad), min(H, sl[0].stop + pad); x0, x1 = max(0, sl[1].start - pad), min(W, sl[1].stop + pad)
        reg = lab[y0:y1, x0:x1] == i
        cc = c[y0:y1, x0:x1]
        dist = ndi.distance_transform_edt(reg)
        inner = dist >= max(1.0, 0.4 * dist.max())
        core = float(np.percentile(cc[inner], 60)) if inner.any() else 0.0
        if core <= 0:
            out[y0:y1, x0:x1][reg] = i; continue
        cell = near[y0:y1, x0:x1] == i
        allowed = cell & ndi.binary_dilation(reg, iterations=int(max(min_grow_px, grow_f * r)))
        iso = (cc >= level * core) & allowed
        l2, n2 = ndi.label(iso)
        if n2 == 0 or not (l2[reg] > 0).any():
            out[y0:y1, x0:x1][reg] = i; continue
        keep = np.bincount(l2[reg & (l2 > 0)]).argmax()      # the component that holds most of the old particle
        new = ndi.binary_fill_holes(l2 == keep)
        if delta_px > 0:
            d2 = ndi.distance_transform_edt(np.pad(new, 1))[1:-1, 1:-1]
            new = d2 > delta_px
        if new.sum() < 0.25 * reg.sum():                        # something went wrong: keep the original outline
            out[y0:y1, x0:x1][reg] = i; continue
        sub = out[y0:y1, x0:x1]
        sub[new & (sub == 0)] = i
    return out


def segment_adaptive(img, dark_particles=True, smooth=1.0, low_k=4.0, size_px=None, marker_sigma=0.25, nms_factor=0.6,
                     core_weight=1.0, dist_weight=0.5, merge_factor=0.3, solidity_min=0.8, area_factor=2.0,
                     clear_border=True, max_px=1536, threshold=None, shape_prior=None, prior_strength=0.7, prior_k=10.0,
                     sphere="auto", conf_min=None, sphere_level=None, sphere_delta=None, edge_refine=True, merge_pieces="auto",
                     annotations=False):
    """v2.10 default: background flattening + hysteresis threshold, cube-scale markers (LoG of the contrast at σ =
    marker_sigma·L plus the distance map of the locally-darker cores), greedy NMS at nms_factor·L, marker watershed on
    the contrast; regions > area_factor·L² or with solidity < solidity_min are 'aggregate' (excluded by default).
    Optional learned shape prior (shape_prior=True): see mfs.shapeprior.
    v2.12 sphere model (mfs.spheres): sphere="auto" fits circles on the real particle edge when the clean particles are
    round (spheres), True always, False never (v2.10 outlines). Circle particles carry info[i]['d_px'] / 'conf';
    'lowconf' particles (confidence < conf_min) start excluded. Returns (labels, factor, threshold, info)."""
    a, f = _downsample(np.asarray(img), max_px)
    if a.ndim == 3:
        a = a.mean(axis=2)
    c, sn = contrast_map(a, dark_particles, smooth)
    hi = max(otsu(c), 8 * sn); lo = max(low_k * sn, 0.3 * hi)
    if threshold is not None:                      # manual: a fraction of the contrast range (0–1) or 8-bit level
        hi = lo = float(threshold)
    m = ndi.binary_opening(hysteresis(c, lo, hi), iterations=1)
    L = float(size_px / f) if size_px else estimate_side_px(m)
    hl, _ = ndi.label(ndi.binary_fill_holes(m) & ~m)       # fill only small holes (not the gaps inside a ring of particles)
    ha = np.bincount(hl.ravel()); small = ha < 0.5 * L * L; small[0] = False
    m |= small[hl]
    d = ndi.distance_transform_edt(m)
    s = max(1.0, L * marker_sigma)
    B = -ndi.gaussian_laplace(c, s) * s * s
    cs = ndi.gaussian_filter(c, max(1.0, L / 10))
    core = ndi.binary_opening(m & (cs > ndi.gaussian_filter(cs, L * 0.75)), iterations=2)
    dc = np.minimum(ndi.distance_transform_edt(core), L / 2) / (L / 2)
    S = B / (np.percentile(B[m], 95) + 1e-9 if m.any() else 1) + core_weight * dc
    pk = (S == ndi.maximum_filter(S, size=int(L * 0.5) | 1)) & m & (d > 0.2 * L) & (S > 0.05)
    yy, xx = np.nonzero(pk)
    keep = _nms(np.stack([yy, xx], 1).astype(float), S[pk], nms_factor * L)
    mk = np.zeros(m.shape, int); mk[yy[keep], xx[keep]] = np.arange(1, len(keep) + 1)
    cost = -(cs / (np.percentile(cs[m], 95) + 1e-9 if m.any() else 1)) - dist_weight * np.minimum(d, L / 2) / (L / 2)
    lab = _ws(cost, mk, m) if len(keep) else ndi.label(m)[0]
    rest = m & (lab == 0)
    if rest.any():
        l2, _ = ndi.label(rest); lab = np.where(rest, l2 + lab.max(), lab)
    A1 = L * L
    lab = _merge_small(lab, merge_factor * A1)
    ep = elongated_population(m, L)
    do_merge = bool(merge_pieces is True or (merge_pieces == "auto" and ep["elongated"]))
    if do_merge:
        lab = merge_elongated(lab, c, L)                  # v2.14: rods that were cut into pieces become one particle again
    info = {}
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        if sl is None:
            continue
        sub = lab[sl] == i; A = int(sub.sum())
        if A < 0.15 * A1:                          # isolated specks (noise, dust)
            lab[sl][sub] = 0; continue
        so = _solidity(sub)
        # a big particle is an aggregate only when it is not convex: an elongated single particle (rod) is much larger than a square
        big_ok = do_merge and so >= 0.93
        info[i] = {"kind": "aggregate" if ((A > area_factor * A1 and not big_ok) or so < solidity_min) else "particle", "solidity": so}
    # v2.21 (opt-in, annotations=True): text and the scale bar of a micrograph are not particles
    if annotations:
        for i in find_annotations(a, lab):
            lab[lab == i] = 0
            if isinstance(info.get(i), dict):
                info[i]["kind"] = "annotation"
    info_extra = {"_r_est_px": L / math.sqrt(math.pi), "_side_px": L, "_min_dist_px": nms_factor * L, "_holes": 0,
                  "_noise": sn, "_lo": lo, "_hi": hi, "_contrast": c, "_elongated": ep, "_pieces_merged": do_merge}
    sph = None
    if sphere:
        from mfs import spheres as SPH
        rd = SPH.roundness(lab, info, L)
        precise_mode = sphere == "precise"                  # v2.21: model fit of every sphere (sub-20 nm); needs round particles
        sph = {"mode": "precise" if precise_mode else ("auto" if sphere == "auto" else ("on" if sphere else "off")), "roundness": rd,
               "active": bool(precise_mode or sphere is True or sphere == "on" or (sphere == "auto" and rd.get("round")))}
        if sph["active"] and precise_mode:
            # v2.21 precise mode: the outline candidates are fitted one by one with the blurred absorbing-sphere model
            # on their own Voronoi cell (mfs.precise.refine_labels). Edge particles are dropped there.
            from mfs import precise as PR
            lab, info, psum = PR.refine_labels(c, lab, info)
            sph.update(psum)
            info_extra["_sphere"] = sph
            clear_border = False
            shape_prior = False
        elif sph["active"]:
            lv = SPH.LEVEL if sphere_level is None else float(sphere_level)
            de = SPH.DELTA_PX if sphere_delta is None else float(sphere_delta)
            cm = SPH.CONF_MIN if conf_min is None else float(conf_min)
            circ, rest = SPH.resolve(c, m, L, level=lv, delta_px=de, conf_min=cm)
            H_, W_ = m.shape
            n_edge = 0
            if clear_border:
                # v2.21: a particle cut by the image edge is not measured. The circle fit alone can put a cut particle inside
                # the image (it sees only the visible arc), so the ray test is used too: any ray that leaves the image before
                # it finds the particle edge means the particle is cut.
                keep_c = []
                for o in circ:
                    inside = (o["cx"] - o["r_edge"] > 0.5 and o["cy"] - o["r_edge"] > 0.5
                              and o["cx"] + o["r_edge"] < W_ - 1.5 and o["cy"] + o["r_edge"] < H_ - 1.5)
                    cut = (not inside) or o.get("out_frac", 0.0) > 0.02
                    if cut:
                        n_edge += 1
                    else:
                        keep_c.append(o)
                circ = keep_c
            for j, o in enumerate(circ, start=1):
                o["id"] = j
            lab, extra = SPH.to_labels(circ, m.shape, rest)
            clus, _ = ndi.label(ndi.binary_fill_holes(m))
            per = np.bincount([clus[min(H_ - 1, max(0, int(round(o["cy"])))), min(W_ - 1, max(0, int(round(o["cx"]))))]
                               for o in circ] or [0])
            info = {}
            for o in circ:
                cl = clus[min(H_ - 1, max(0, int(round(o["cy"])))), min(W_ - 1, max(0, int(round(o["cx"]))))]
                info[o["id"]] = {"kind": o["kind"], "solidity": 1.0, "conf": o["conf"], "d_px": 2 * o["r"],
                                 "circle": (o["cx"], o["cy"], o["r"]), "arc": o["arc"], "rms": o["rms"],
                                 "occluded": bool(o["occluded_by"]), "split": "circle" if (cl and per[cl] > 1) else None}
            for e in extra:
                info[e] = {"kind": "aggregate", "solidity": 0.0}
            sph.update({"n": len(circ), "n_lowconf": sum(o["kind"] == "lowconf" for o in circ),
                        "n_split": sum(1 for v in info.values() if v.get("split") == "circle"),
                        "n_occluded": sum(o["occluded_by"] != [] for o in circ), "level": lv, "delta_px": de, "conf_min": cm,
                        "n_unexplained": len(extra), "n_edge_excluded": n_edge})
            info_extra["_sphere"] = sph
            shape_prior = False                    # the circle model replaces the learned-outline splitting
            clear_border = False                   # done above on the fitted circles
        else:
            info_extra["_sphere"] = sph
    if shape_prior:
        from mfs import shapeprior as SP
        lab, info, pr = SP.apply_prior(lab, info, c, L, strength=prior_strength, k=prior_k,
                                       solidity_min=solidity_min, area_factor=area_factor)
        info_extra["_prior"] = pr
    if clear_border:
        edge = np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))
        info_extra["_edge_excluded"] = int((edge > 0).sum())      # v2.21: reported in the module
        lab = np.where(np.isin(lab, edge[edge > 0]), 0, lab)
        info = {k: v for k, v in info.items() if k not in set(edge.tolist())}
    if edge_refine and not (sph and sph.get("active")):      # v2.14: calibrated edge for non-round particles (circles already have it)
        lab = refine_edges(lab, c)
        info_extra["_edge_refined"] = True
    present = set(np.unique(lab).tolist()) - {0}
    info = {k: v for k, v in info.items() if k in present}
    info.update(info_extra)
    return lab, f, hi, info


def find_annotations(a, lab):
    """Labels of micrograph annotations: the scale bar (a long, nearly black bar), the text in the band just above it,
    and dark panel letters in the outer 12 % of the image. Particles are grey (median well above half the background),
    so only objects far darker than the background are considered. Returns a set of label ids (empty without a bar)."""
    H, W = lab.shape
    bg = float(np.median(a))
    stats = []
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        if sl is None:
            continue
        m = lab[sl] == i
        med = float(np.median(a[sl][m]))
        h, w = m.shape
        stats.append({"i": i, "sl": sl, "med": med, "area": int(m.sum()), "elong": max(h, w) / max(1.0, min(h, w)),
                      "cy": (sl[0].start + sl[0].stop) / 2, "cx": (sl[1].start + sl[1].stop) / 2, "box": sl})
    bars = [s for s in stats if s["med"] < 0.25 * bg and s["elong"] >= 6 and s["area"] >= 60]
    out = set()
    if bars:
        bar = max(bars, key=lambda s: s["area"])
        out.add(bar["i"])
        by0, by1 = bar["sl"][0].start, bar["sl"][0].stop
        bx0, bx1 = bar["sl"][1].start, bar["sl"][1].stop
        bh = by1 - by0
        zone_y0, zone_y1 = by0 - 3 * bh - 25, by1 + 4
        for s in stats:
            if s["med"] < 0.7 * bg and zone_y0 <= s["cy"] <= by0 and bx0 - 15 <= s["cx"] <= bx1 + 15:
                out.add(s["i"])
    for s in stats:
        if s["med"] < 0.45 * bg and s["area"] >= 60 and (s["cx"] < 0.12 * W or s["cy"] < 0.12 * H):
            if bars:                                    # letters only when a scale bar confirms this is a micrograph panel
                out.add(s["i"])
    return out
