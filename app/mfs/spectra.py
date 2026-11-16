"""Vibrational-spectra helpers (FTIR / Raman): normalisation to a reference band, peak picking for labels, difference
spectra on a common grid, zoomed insets."""
from __future__ import annotations
import numpy as np
from scipy.signal import find_peaks


def _sorted(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); o = np.argsort(x)
    return x[o], y[o]


def band_value(x, y, lo, hi, mode="height"):
    """Height (max above a linear chord between the window ends) or area of a band in [lo, hi]."""
    x, y = _sorted(x, y)
    m = (x >= min(lo, hi)) & (x <= max(lo, hi))
    if m.sum() < 3:
        raise ValueError(f"reference band {lo:g}–{hi:g} outside the spectrum")
    xs, ys = x[m], y[m]
    base = np.interp(xs, [xs[0], xs[-1]], [ys[0], ys[-1]])
    s = ys - base
    return float(np.trapezoid(s, xs)) if mode == "area" else float(s.max())


def normalise_to_band(x, y, lo, hi, mode="height", minimum=False):
    """Scale a spectrum so the reference band (e.g. ν(C=C) or the framework μ3-O band) has height/area 1.
    minimum=True for transmittance spectra (band = dip): uses 1 − T style depth."""
    y = np.asarray(y, float)
    yy = -y if minimum else y
    v = band_value(x, yy, lo, hi, mode)
    if not v or v <= 0:
        raise ValueError("reference band has zero intensity")
    if minimum:          # transmittance (% or fraction): top of the spectrum -> 1, band depth -> 1
        top = np.nanmax(y)
        return 1 + (y - top) / v
    return y / v


def pick_peaks(x, y, prominence=0.05, n_max=12, xlim=None, minima=False, min_sep=None):
    """Peak positions for labels: prominence is relative to the trace range inside xlim; the n_max most prominent,
    returned sorted by position. minima=True for transmittance dips."""
    x, y = _sorted(x, y)
    if xlim:
        m = (x >= min(xlim)) & (x <= max(xlim)); x, y = x[m], y[m]
    if len(x) < 5:
        return []
    s = -y if minima else y
    rng = np.nanmax(s) - np.nanmin(s) or 1.0
    dx = np.median(np.diff(x)) or 1.0
    dist = max(1, int((min_sep or 8 * dx) / dx))
    pk, pr = find_peaks(np.nan_to_num(s, nan=np.nanmin(s)), prominence=prominence * rng, distance=dist)
    o = np.argsort(pr["prominences"])[::-1][:n_max]
    return sorted(float(x[i]) for i in pk[o])


def difference(ref, sample, factor=1.0, grid=None):
    """sample − factor × ref on a common grid (overlap region, step = the coarser of the two)."""
    (xr, yr), (xs, ys) = _sorted(*ref), _sorted(*sample)
    lo, hi = max(xr.min(), xs.min()), min(xr.max(), xs.max())
    if hi <= lo:
        raise ValueError("spectra do not overlap")
    if grid is None:
        step = max(np.median(np.diff(xr)), np.median(np.diff(xs)))
        grid = np.arange(lo, hi, step)
    d = np.interp(grid, xs, ys) - factor * np.interp(grid, xr, yr)
    return grid, d


def auto_factor(ref, sample, lo, hi):
    """Least-squares subtraction factor inside a window where only the reference contributes."""
    g, _ = difference(ref, sample)
    (xr, yr), (xs, ys) = _sorted(*ref), _sorted(*sample)
    m = (g >= min(lo, hi)) & (g <= max(lo, hi))
    a, b = np.interp(g[m], xr, yr), np.interp(g[m], xs, ys)
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).sum() / max((a * a).sum(), 1e-300))


def label_peaks(ax, x, y, positions, offset=0.0, fmt="{:.0f}", fs=6, color="0.25", above=True, dy=0.03):
    """Write wavenumber labels next to the peak maxima (vertical text) on an already drawn (stacked) trace."""
    x, y = _sorted(x, y)
    out = []
    for p in positions:
        v = float(np.interp(p, x, y)) + offset
        t = ax.annotate(fmt.format(p), (p, v), xytext=(0, 3 if above else -3), textcoords="offset points", rotation=90,
                        ha="center", va="bottom" if above else "top", fontsize=fs, color=color)
        t.set_gid("mfs_peaklabel"); out.append(t)
    return out


def add_inset(ax, draw_fn, xlim, loc="upper left", size=0.38, ylim=None, box=True):
    """Zoomed inset (draw_fn(ax_inset) redraws the traces) + the indicated region on the parent axes."""
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    ia = inset_axes(ax, width=f"{100 * size:.0f}%", height=f"{100 * size:.0f}%", loc=loc, borderpad=0.8)
    draw_fn(ia)
    ia.set_xlim(*xlim)
    if ylim:
        ia.set_ylim(*ylim)
    else:
        ia.relim(); ia.autoscale_view(scalex=False)
    ia.tick_params(labelsize=5, length=2); ia.set_yticks([])
    for t in ia.texts[:]:
        t.remove()
    if ia.get_legend():
        ia.get_legend().remove()
    ia.set_xlabel(""); ia.set_ylabel(""); ia.set_title("")
    if box:
        ax.axvspan(min(xlim), max(xlim), color="0.85", alpha=0.35, lw=0, zorder=0).set_gid("mfs_inset_region")
    ia.set_gid("mfs_inset")
    return ia
