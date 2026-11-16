"""EDX map overlays (additive colour blend over an optional SE/BSE image) and line scans."""
from __future__ import annotations
import numpy as np
from scipy import ndimage as ndi


def _rgb(hexc):
    h = hexc.lstrip("#"); return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], float) / 255


def norm_map(m, pct=(1, 99.5), gamma=1.0):
    m = np.asarray(m, float)
    lo, hi = np.percentile(m, pct)
    return np.clip((m - lo) / ((hi - lo) or 1), 0, 1) ** gamma


def overlay(maps, colors, base=None, base_alpha=0.6, gains=None, pct=(1, 99.5), smooth=0.0):
    """maps: {element: 2D array}; colors: {element: '#rrggbb'}. Element colours are added (clipped) and laid over the
    grey base image, which is dimmed where the elements are strong so the colours stay saturated."""
    shp = next(iter(maps.values())).shape
    col = np.zeros(shp + (3,)); a = np.zeros(shp)
    for el, m in maps.items():
        m = np.asarray(m, float)
        if m.shape != shp:
            m = ndi.zoom(m, (shp[0] / m.shape[0], shp[1] / m.shape[1]), order=1)
        if smooth:
            m = ndi.gaussian_filter(m, smooth)
        v = np.clip(norm_map(m, pct) * (gains or {}).get(el, 1.0), 0, 1)
        col += v[..., None] * _rgb(colors.get(el, "#FFFFFF")); a = np.maximum(a, v)
    col = np.clip(col, 0, 1)
    if base is None:
        return col
    b = np.asarray(base, float)
    if b.ndim == 3:
        b = b.mean(axis=2)
    if b.shape != shp:
        b = ndi.zoom(b, (shp[0] / b.shape[0], shp[1] / b.shape[1]), order=1)
    b = norm_map(b, (0.5, 99.5))
    return np.clip((base_alpha * b * (1 - a))[..., None] + col, 0, 1)


def line_scan(m, p0, p1, width=3, n=None):
    """Mean intensity across a band of `width` px along p0→p1 ((x, y) in px). Returns distance (px) and profile."""
    (x0, y0), (x1, y1) = p0, p1
    L = float(np.hypot(x1 - x0, y1 - y0)); n = n or max(2, int(L))
    t = np.linspace(0, 1, n); xs, ys = x0 + t * (x1 - x0), y0 + t * (y1 - y0)
    nx, ny = -(y1 - y0) / (L or 1), (x1 - x0) / (L or 1)
    offs = np.linspace(-(width - 1) / 2, (width - 1) / 2, max(1, int(width)))
    prof = np.mean([ndi.map_coordinates(np.asarray(m, float), [ys + o * ny, xs + o * nx], order=1) for o in offs], axis=0)
    return t * L, prof
