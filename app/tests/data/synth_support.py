"""Synthetic SEM of dark faceted MOF-801 particles (300–600 nm) in clusters of 1–6 on holey carbon (bright round holes)."""
import numpy as np
from scipy import ndimage as ndi
from matplotlib.path import Path


def make(seed=1, H=900, W=1200, nm_px=6.0):
    rng = np.random.default_rng(seed)
    img = np.full((H, W), 135.0)
    yy, xx = np.mgrid[:H, :W]
    for _ in range(9):                                   # holes
        cx, cy, r = rng.uniform(0, W), rng.uniform(0, H), rng.uniform(50, 110)
        img[(xx - cx) ** 2 + (yy - cy) ** 2 < r * r] = 235
    truth = []; lab = np.zeros((H, W), int); k = 0
    pts = np.stack([xx.ravel(), yy.ravel()], 1)
    centres = []
    for c in range(16):                                  # clusters
        n = rng.choice([1, 1, 2, 3, 4, 5, 6])
        cx, cy = rng.uniform(80, W - 80), rng.uniform(80, H - 80)
        ang = rng.uniform(0, 2 * np.pi)
        for j in range(n):
            d_nm = np.exp(rng.normal(np.log(430), 0.18)); r = d_nm / 2 / nm_px
            if j:
                ang += rng.uniform(1.0, 2.2); px, py = centres[-1]
                cx2, cy2 = px + 1.75 * r * np.cos(ang), py + 1.75 * r * np.sin(ang)
            else:
                cx2, cy2 = cx, cy
            if any(np.hypot(cx2 - a, cy2 - b) < 1.7 * r for a, b in centres):
                continue
            if j and not (r < cx2 < W - r and r < cy2 < H - 80 - r):
                continue
            centres.append((cx2, cy2))
            rot = rng.uniform(0, np.pi / 3)
            poly = [(cx2 + r * 1.07 * np.cos(rot + q * np.pi / 3), cy2 + r * 1.07 * np.sin(rot + q * np.pi / 3)) for q in range(6)]
            m = Path(poly).contains_points(pts).reshape(H, W)
            k += 1; lab[m & (lab == 0)] = k
            img[m] = 55 + rng.uniform(-8, 8)
            # facet: one brighter internal triangle (top face)
            tri = [(cx2, cy2), poly[0], poly[1]]
            img[Path(tri).contains_points(pts).reshape(H, W) & m] += 18
    img = ndi.gaussian_filter(img, 1.2) + rng.normal(0, 7, img.shape)
    img[-70:] = 10; img[-45:-30, 40:240] = 250            # info strip with a scale bar
    ids = [i for i in range(1, k + 1) if (lab == i).any()]
    ecd = [2 * np.sqrt((lab == i).sum() / np.pi) * nm_px for i in ids]
    return np.clip(img, 0, 255).astype(np.uint8), lab, np.array(ecd), nm_px
