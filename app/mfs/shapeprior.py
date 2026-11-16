"""v2.10 optional learned shape prior for particle sizing.

From the confident, isolated single particles of the image the typical projected outline is learned (radial profile
r(θ) aligned on its dominant Fourier harmonic → corner count m, orientation; mean template), with the size
distribution and aspect ratio. Its weight grows with the number of examples and how well they agree:
    w = n / (n + k) · consistency · strength
(consistency = mean polar IoU of the aligned examples with the template). The prior then breaks up 'aggregate'
regions: rotated templates at the learned scale are matched (normalised cross-correlation on the contrast map, FFT)
inside the aggregate (and inside 'particles' larger than 1.6× the learned median area), the peaks seed a marker watershed and every piece is accepted only if its area lies inside the
learned size distribution and its outline matches the template (IoU ≥ 0.85 − 0.3·w). Accepted pieces are
particles flagged split='prior' (dotted outline in the picker); the rest stays one aggregate."""
from __future__ import annotations
import math
import numpy as np
from scipy import ndimage as ndi

NB = 72                                     # angular bins of the radial profile


def radial_profile(mask):
    """r(θ) (NB bins, px) of a binary mask about its centroid: the farthest mask pixel per angle bin (gaps filled)."""
    yy, xx = np.nonzero(mask)
    cy, cx = yy.mean(), xx.mean()
    dy, dx = yy - cy, xx - cx
    th = np.mod(np.arctan2(dy, dx), 2 * np.pi); r = np.hypot(dy, dx) + 0.5
    b = np.minimum((th / (2 * np.pi) * NB).astype(int), NB - 1)
    prof = np.zeros(NB); np.maximum.at(prof, b, r)
    if (prof == 0).any():
        good = prof > 0
        idx = np.arange(NB)
        prof = np.interp(idx, idx[good], prof[good], period=NB)
    return prof, (cy, cx)


def harmonics(prof):
    F = np.fft.rfft(prof / prof.mean())
    return F


def polar_iou(r1, r2):
    a, b = np.minimum(r1, r2), np.maximum(r1, r2)
    return float((a * a).sum() / max((b * b).sum(), 1e-9))


def template_mask(prof_unit, scale, angle_deg=0.0, pad=2):
    """Binary mask of the template outline r(θ)·scale rotated by angle_deg."""
    R = int(math.ceil(prof_unit.max() * scale)) + pad
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
    th = np.mod(np.arctan2(yy, xx) - math.radians(angle_deg), 2 * np.pi)
    b = th / (2 * np.pi) * NB
    i0 = np.floor(b).astype(int) % NB; i1 = (i0 + 1) % NB; t = b - np.floor(b)
    rr = (prof_unit[i0] * (1 - t) + prof_unit[i1] * t) * scale
    return np.hypot(yy, xx) <= rr


def contact_fraction(lab, i, sl):
    """Fraction of the particle's outer ring that touches another particle (0 = fully isolated)."""
    H, W = lab.shape
    s2 = tuple(slice(max(0, s.start - 2), min(n, s.stop + 2)) for s, n in zip(sl, (H, W)))
    sub = lab[s2]; me = sub == i
    ring = ndi.binary_dilation(me) & ~me
    return float((sub[ring] > 0).mean()) if ring.any() else 0.0


def _isolated(lab, i, sl, max_contact=0.3):
    """Confident single: touches no image edge and at most max_contact of its outline touches other particles."""
    H, W = lab.shape
    touches_edge = sl[0].start == 0 or sl[1].start == 0 or sl[0].stop >= H or sl[1].stop >= W
    return not touches_edge and contact_fraction(lab, i, sl) <= max_contact


def learn(lab, info, contrast, L, k=10.0, strength=1.0, max_examples=400, n_show=12):
    """Learn the shape prior from confident isolated single particles. Returns a dict (see module doc)."""
    A1 = L * L
    areas = np.bincount(lab.ravel())
    ex = []
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        if sl is None or not isinstance(info.get(i), dict) or info[i].get("kind") != "particle":
            continue
        A = areas[i]
        if not (0.5 * A1 <= A <= 1.8 * A1) or info[i].get("solidity", 0) < 0.9 or not _isolated(lab, i, sl):
            continue
        ex.append(i)
        if len(ex) >= max_examples:
            break
    out = {"n": len(ex), "k": float(k), "strength": float(strength), "examples": ex, "weight": 0.0, "consistency": 0.0}
    if len(ex) < 3:
        out["reason"] = f"only {len(ex)} isolated single particle(s) – prior not learned"
        return out
    sls = ndi.find_objects(lab)
    profs, F_all, ar, A = [], [], [], []
    for i in ex:
        sub = lab[sls[i - 1]] == i
        p, _ = radial_profile(sub); profs.append(p); F_all.append(harmonics(p)); A.append(sub.sum())
        yy, xx = np.nonzero(sub); ang = np.linspace(0, np.pi, 36, endpoint=False)
        P = np.stack([xx, yy], 1) @ np.stack([np.cos(ang), np.sin(ang)])
        w = P.max(0) - P.min(0) + 1; ar.append(w.max() / w.min())
    F_all = np.array(F_all)
    amp = np.abs(F_all[:, 2:9]).mean(0)                      # harmonics m = 2..8
    m = [3, 4, 6][int(np.argmax(amp[[1, 2, 4]]))]           # corner count of the projected outline
    phi = -np.angle(F_all[:, m]) / m                          # orientation of each example (rad)
    al = []
    for p, ph in zip(profs, phi):                             # align: rotate by −φ
        sh = ph / (2 * np.pi) * NB
        idx = (np.arange(NB) + sh) % NB
        i0 = np.floor(idx).astype(int) % NB; t = idx - np.floor(idx)
        q = p[i0] * (1 - t) + p[(i0 + 1) % NB] * t
        al.append(q / q.mean())
    al = np.array(al)
    tmpl = np.median(al, 0)
    tmpl = tmpl / math.sqrt((tmpl ** 2).mean())               # unit RMS radius → area = π·scale²
    ious = np.array([polar_iou(q / math.sqrt((q ** 2).mean()), tmpl) for q in al])
    cons = float(np.clip((ious.mean() - 0.6) / 0.35, 0, 1))  # IoU 0.6 → 0, 0.95 → 1
    z = np.exp(1j * m * phi); R_ = abs(z.mean())
    osd = math.degrees(math.sqrt(max(0.0, -2 * math.log(max(R_, 1e-9)))) / m)
    A = np.array(A, float); ecd = 2 * np.sqrt(A / np.pi)
    w = len(ex) / (len(ex) + k) * cons * strength
    rng = np.random.default_rng(0)
    show = sorted(rng.choice(ex, size=min(n_show, len(ex)), replace=False).tolist())
    out.update({"corners": m, "template": tmpl, "consistency": cons, "iou_mean": float(ious.mean()),
                "area_px": A, "area_p": np.percentile(A, [5, 50, 95]).tolist(), "ecd_px": ecd,
                "aspect_mean": float(np.mean(ar)), "aspect_sd": float(np.std(ar)),
                "orientation_sd_deg": osd, "orientations_deg": np.degrees(phi).tolist(),
                "harmonics": amp.tolist(), "weight": float(min(1.0, w)), "show": show,
                "crops": [_crop(contrast, sls[i - 1]) for i in show]})
    return out


def _crop(c, sl, pad=4):
    H, W = c.shape
    s2 = tuple(slice(max(0, s.start - pad), min(n, s.stop + pad)) for s, n in zip(sl, (H, W)))
    return np.asarray(c[s2], np.float32)


def template_scale(prior):
    """Scale (px) of the unit-RMS template giving the median learned area."""
    return math.sqrt(prior["area_p"][1] / math.pi)


def match_map(contrast, region, prior, angles):
    """Normalised cross-correlation of the region's contrast with the rotated templates; max over angles (+ argmax)."""
    from scipy.signal import fftconvolve
    sig = np.where(region, contrast, 0.0)
    sig = sig / (np.percentile(sig[region], 95) + 1e-9)
    best = np.full(sig.shape, -1.0); barg = np.zeros(sig.shape, int)
    sc = template_scale(prior)
    for j, ang in enumerate(angles):
        T = template_mask(prior["template"], sc, ang).astype(float)
        n = T.size; Tz = T - T.mean(); nt = math.sqrt((Tz ** 2).sum()) + 1e-9
        box = np.ones_like(T)
        s1 = fftconvolve(sig, box[::-1, ::-1], "same"); s2 = fftconvolve(sig * sig, box[::-1, ::-1], "same")
        var = np.maximum(s2 - s1 * s1 / n, 1e-9)
        num = fftconvolve(sig, Tz[::-1, ::-1], "same")
        ncc = num / (nt * np.sqrt(var))
        upd = ncc > best; best[upd] = ncc[upd]; barg[upd] = j
    return best, barg


def shape_score(piece, prior):
    """Polar IoU of a piece with the template at its best orientation (size-normalised)."""
    p, _ = radial_profile(piece)
    q = p / math.sqrt((p ** 2).mean())
    t = prior["template"]
    return max(polar_iou(np.roll(q, s), t) for s in range(0, NB, 1))


def apply_prior(lab, info, contrast, L, strength=0.7, k=10.0, solidity_min=0.8, area_factor=2.0,
                ncc_min=0.35, max_regions=400):
    """Break up 'aggregate' regions with the learned prior. Returns (labels, info, prior summary)."""
    from mfs.particles import _ws, _nms
    pr = learn(lab, info, contrast, L, k=k, strength=strength)
    pr.update({"regions_tried": 0, "regions_split": 0, "new_particles": 0, "applied": False})
    if pr["n"] < 3 or pr["weight"] < 0.05:
        return lab, info, pr
    pr["applied"] = True
    w = pr["weight"]; iou_min = 0.85 - 0.3 * w
    lo_a, hi_a = 0.7 * pr["area_p"][0], 1.3 * pr["area_p"][2]
    sym = 360.0 / pr["corners"]
    angles = np.arange(0, sym, 7.5 if sym <= 90 else 15.0)
    sc = template_scale(pr); D = 2 * sc
    lab = lab.copy(); nxt = int(lab.max()) + 1
    H, W = lab.shape
    sls = ndi.find_objects(lab)
    areas = np.bincount(lab.ravel())
    big = 1.6 * pr["area_p"][1]                               # aggregates + oversized 'particles' (two fused singles)
    aggs = [i for i, v in info.items() if isinstance(v, dict) and isinstance(i, int) and i < len(areas)
            and (v.get("kind") == "aggregate" or areas[i] > big)][:max_regions]
    for i in aggs:
        sl = sls[i - 1] if i - 1 < len(sls) else None
        if sl is None:
            continue
        pad = int(D)
        s2 = tuple(slice(max(0, s.start - pad), min(n, s.stop + pad)) for s, n in zip(sl, (H, W)))
        reg = lab[s2] == i
        if reg.sum() < 1.5 * pr["area_p"][1]:
            continue
        pr["regions_tried"] += 1
        sc_map, _ = match_map(contrast[s2], reg, pr, angles)
        pk = (sc_map == ndi.maximum_filter(sc_map, size=max(3, int(0.5 * D)) | 1)) & reg & (sc_map > ncc_min)
        yy, xx = np.nonzero(pk)
        keep = _nms(np.stack([yy, xx], 1).astype(float), sc_map[pk], 0.7 * D)
        if len(keep) < 2:
            continue
        mk = np.zeros(reg.shape, int); mk[yy[keep], xx[keep]] = np.arange(1, len(keep) + 1)
        dd = np.minimum(ndi.distance_transform_edt(reg), sc) / sc
        pieces = _ws(-(np.clip(sc_map, 0, 1) + 0.5 * dd), mk, reg)
        good = []
        for j in range(1, len(keep) + 1):
            pm = pieces == j
            a_ = pm.sum()
            if a_ < 10 or not (lo_a <= a_ <= hi_a):
                continue
            s_ = shape_score(pm, pr)
            if s_ >= iou_min:
                good.append((j, s_))
        if len(good) < 2 and not (len(good) == 1 and reg.sum() <= 2.2 * pr["area_p"][1]):
            continue
        blk = lab[s2]
        for j, s_ in good:
            blk[pieces == j] = nxt
            info[nxt] = {"kind": "particle", "solidity": 1.0, "split": "prior", "shape_score": float(s_)}
            nxt += 1
        pr["regions_split"] += 1; pr["new_particles"] += len(good)
        if not (blk == i).any():
            info.pop(i, None)
    return lab, info, pr
