"""v2.15: vibrational band assignment (own, editable list) and simultaneous multi-band fit for FTIR / Raman.

The built-in ranges are the TYPICAL ranges from standard IR / Raman tables; they are a starting point, not proof - the assignment of a band
in a MOF depends on the framework, so the table can be edited and extended (the user's own rows take precedence and are kept in the project)."""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import curve_fit

# (low cm-1, high cm-1, assignment, technique 'ir' / 'raman' / 'both')
BUILTIN = [
    (3600, 3200, "ν(O–H) – hydroxyl / adsorbed water (H-bonded, broad)", "ir"),
    (3500, 3300, "ν(N–H) – amine / amide", "ir"),
    (3100, 3000, "ν(C–H) aromatic", "both"),
    (3000, 2850, "ν(C–H) aliphatic (CH₂ / CH₃)", "both"),
    (1740, 1680, "ν(C=O) carboxylic acid / ester (protonated, free linker or modulator)", "ir"),
    (1690, 1640, "ν(C=O) amide I (e.g. DMF)", "ir"),
    (1660, 1620, "δ(H₂O) – water bending", "ir"),
    (1620, 1540, "ν_as(COO⁻) – coordinated carboxylate", "both"),
    (1620, 1570, "ν(C=C) aromatic ring", "raman"),
    (1520, 1450, "ν(C=C) aromatic ring", "ir"),
    (1450, 1360, "ν_s(COO⁻) – coordinated carboxylate", "both"),
    (1300, 1150, "ν_as(S=O) – sulfonate / sulfone", "ir"),
    (1250, 1000, "ν(C–O), ν(C–N), in-plane δ(C–H)", "ir"),
    (1080, 1020, "ν_s(S=O) – sulfonate", "ir"),
    (1020, 990, "ring breathing (aromatic)", "raman"),
    (900, 700, "γ(C–H) out-of-plane aromatic", "ir"),
    (800, 600, "ν(M–O) / μ₃-OH – metal–oxygen (Zr–O, Hf–O cluster)", "both"),
    (600, 400, "ν(M–O) lattice / M–O–M", "both"),
]


def library(user_rows=None, technique="ir"):
    """The assignment table for a technique ('ir' or 'raman'): user rows first, then the built-in ones. Rows = (lo, hi, text, source)."""
    rows = []
    for r in user_rows or []:
        try:
            lo, hi = sorted((float(r["from (cm⁻¹)"]), float(r["to (cm⁻¹)"])))
            if str(r["assignment"]).strip():
                rows.append((hi, lo, str(r["assignment"]).strip(), "own"))
        except (KeyError, TypeError, ValueError):
            continue
    for hi, lo, txt, tech in BUILTIN:
        if tech in ("both", technique):
            rows.append((max(hi, lo), min(hi, lo), txt, "typical range"))
    return rows


def assign(positions, user_rows=None, technique="ir", max_candidates=3):
    """For each band position: the matching rows (the narrowest range around the position first; own rows before typical ranges)."""
    lib = library(user_rows, technique)
    out = []
    for p in positions:
        cand = [(0 if src == "own" else 1, hi - lo, txt, src, lo, hi) for hi, lo, txt, src in lib if lo <= p <= hi]
        cand.sort(key=lambda c: (c[0], c[1]))
        out.append({"position": float(p), "candidates": [(c[2], c[3], (c[4], c[5])) for c in cand[:max_candidates]],
                    "ambiguous": len(cand) > 1})
    return out


def _model(x, *p):
    """Sum of pseudo-Voigt bands (area, centre, fwhm, eta each) + linear background (b0, b1)."""
    n = (len(p) - 2) // 4
    y = p[-2] + p[-1] * (x - x.mean())
    for i in range(n):
        a, c, w, e = p[4 * i:4 * i + 4]
        g = (2 * math.sqrt(math.log(2) / math.pi) / w) * np.exp(-4 * math.log(2) * ((x - c) / w) ** 2)
        l = (2 / (math.pi * w)) / (1 + 4 * ((x - c) / w) ** 2)
        y = y + a * (e * l + (1 - e) * g)
    return y


def fit_bands(x, y, centers, lo=None, hi=None, max_fwhm=200.0, shape="pseudo-Voigt"):
    """Fits all `centers` at once (pseudo-Voigt, a shared linear background) in [lo, hi].
    Returns {'bands': [{center, fwhm, area, height, eta, center_err, fwhm_err, area_err}], 'x', 'y', 'fit', 'background', 'r2', 'rms'}.
    For absorption / Raman peaks that point up; use −y for transmittance dips."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    centers = sorted(float(c) for c in centers)
    lo = (min(centers) - 60.0) if lo is None else lo
    hi = (max(centers) + 60.0) if hi is None else hi
    m = (x >= lo) & (x <= hi)
    xs, ys = x[m], y[m]
    if len(xs) < 8 + 4 * len(centers):
        raise ValueError("too few points in the fit window for this number of bands")
    step = float(np.median(np.abs(np.diff(xs))))
    p0, lb, ub = [], [], []
    ylo = float(np.percentile(ys, 5))
    for c in centers:
        h = float(np.interp(c, xs, ys)) - ylo
        w0 = max(8.0, 4 * step)
        p0 += [max(h, 1e-9) * w0 * 1.2, c, w0, 0.5]
        lb += [0.0, c - 0.6 * w0 - 3 * step, 1.5 * step, 0.0]
        ub += [np.inf, c + 0.6 * w0 + 3 * step, max_fwhm, 1.0]
    p0 += [ylo, 0.0]
    lb += [-np.inf, -np.inf]
    ub += [np.inf, np.inf]
    p, cov = curve_fit(_model, xs, ys, p0=p0, bounds=(lb, ub), maxfev=40000)
    e = np.sqrt(np.clip(np.diag(cov), 0, None))
    fit = _model(xs, *p)
    bg = p[-2] + p[-1] * (xs - xs.mean())
    ss_res = float(np.sum((ys - fit) ** 2))
    ss_tot = float(np.sum((ys - ys.mean()) ** 2))
    bands = []
    for i in range(len(centers)):
        a, c, w, eta = p[4 * i:4 * i + 4]
        ea, ec, ew, _ = e[4 * i:4 * i + 4]
        height = a * (eta * 2 / (math.pi * w) + (1 - eta) * 2 * math.sqrt(math.log(2) / math.pi) / w)
        bands.append({"center": float(c), "fwhm": float(w), "area": float(a), "height": float(height), "eta": float(eta),
                      "center_err": float(ec), "fwhm_err": float(ew), "area_err": float(ea)})
    comps = [_model(xs, *(list(p[4 * i:4 * i + 4]) + [0.0, 0.0])) for i in range(len(centers))]
    return {"bands": bands, "x": xs, "y": ys, "fit": fit, "background": bg, "components": comps,
            "r2": 1 - ss_res / ss_tot if ss_tot > 0 else float("nan"), "rms": math.sqrt(ss_res / len(xs))}


def draw_fit(ax, res, colors=None, label=True, invert_x=True):
    """Data, sum of the bands, background and every band as a filled component."""
    xs = res["x"]
    ax.plot(xs, res["y"], "o", ms=2.0, mfc="none", mec="#222222", mew=0.4, label="data")
    ax.plot(xs, res["fit"], "-", color="#D62728", lw=1.0, label="fit")
    cols = colors or ["#0072B2", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#D55E00", "#7F7F7F", "#8C564B"]
    for i, comp in enumerate(res["components"]):
        ax.fill_between(xs, res["background"], res["background"] + comp, color=cols[i % len(cols)], alpha=0.35, lw=0)
        if label:
            b = res["bands"][i]
            ax.annotate(f"{b['center']:.0f}", (b["center"], float(np.interp(b["center"], xs, res["background"] + comp))), xytext=(0, 3),
                        textcoords="offset points", ha="center", fontsize=6)
    ax.plot(xs, res["background"], "--", color="#7F7F7F", lw=0.6)
    ax.set_xlim((xs.max(), xs.min()) if invert_x else (xs.min(), xs.max()))
