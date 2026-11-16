"""XPS peak fitting: Shirley (iterative) and Tougaard (universal 2-parameter cross-section) backgrounds,
GL (pseudo-Voigt) components, spin–orbit doublets with fixed splitting and area ratio, atomic-% with RSFs."""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import least_squares

# spin–orbit doublets: (lower-j label, higher-j label, splitting eV (BE of high-j minus... see note), area ratio minor/major)
# major = j = l + 1/2 (lower BE), minor = j = l − 1/2 at BE + Δ, area minor/major = l / (l + 1)
DOUBLETS = {
    "Zr 3d": {"major": "3d5/2", "minor": "3d3/2", "split": 2.43, "ratio": 2 / 3, "window": (176.0, 190.0)},
    "Hf 4f": {"major": "4f7/2", "minor": "4f5/2", "split": 1.71, "ratio": 3 / 4, "window": (13.0, 23.0)},
    "S 2p": {"major": "2p3/2", "minor": "2p1/2", "split": 1.18, "ratio": 1 / 2, "window": (160.0, 174.0)},
    "Hf 4d": {"major": "4d5/2", "minor": "4d3/2", "split": 10.9, "ratio": 2 / 3, "window": (205.0, 230.0)},
    "Ti 2p": {"major": "2p3/2", "minor": "2p1/2", "split": 5.7, "ratio": 1 / 2, "window": (452.0, 470.0)},
}
SINGLETS = {"C 1s": (280.0, 294.0), "O 1s": (526.0, 538.0), "N 1s": (396.0, 406.0), "F 1s": (682.0, 692.0)}
# Scofield photoionisation cross-sections, Al Kα 1486.7 eV, C 1s = 1 (whole doublet). Without transmission / IMFP
# correction – replace with the instrument library RSFs (CasaXPS, Avantage, Kratos) for quantitative work.
RSF_SCOFIELD = {"C 1s": 1.00, "N 1s": 1.80, "O 1s": 2.93, "F 1s": 4.43, "Na 1s": 8.52, "Si 2p": 0.82, "S 2p": 1.68,
                "Cl 2p": 2.29, "Ti 2p": 7.81, "Zr 3d": 7.04, "Hf 4f": 7.52, "Hf 4d": 11.53}


def shirley(x, y, tol=1e-6, maxit=60):
    """Iterative Shirley (Proctor–Sherwood). x may be in either order; endpoints = averages of 3 points."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    rev = x[0] > x[-1]
    if rev:
        x, y = x[::-1], y[::-1]
    lo, hi = y[:3].mean(), y[-3:].mean()      # low-BE end, high-BE end
    B = np.full_like(y, lo)
    for _ in range(maxit):
        s = np.clip(y - B, 0, None)
        cum = np.concatenate([[0], np.cumsum(0.5 * (s[1:] + s[:-1]) * np.diff(x))])
        tot = cum[-1] or 1e-300
        Bn = lo + (hi - lo) * cum / tot
        if np.max(np.abs(Bn - B)) < tol * max(abs(hi), 1):
            B = Bn; break
        B = Bn
    return B[::-1] if rev else B


def tougaard(x, y, B=2866.0, C=1643.0, fit_B=True):
    """Tougaard universal cross-section K(T) = B T / (C + T²)²; background at E (BE) = ∫_{E'<E} K(E − E') j(E') dE'.
    fit_B: B is scaled so the background meets the data at the high-BE end."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    rev = x[0] > x[-1]
    if rev:
        x, y = x[::-1], y[::-1]
    y0 = y - y[:3].mean()
    dx = np.gradient(x)
    bg = np.zeros_like(y)
    for i in range(len(x)):
        T = x[i] - x[:i]
        bg[i] = np.sum(T / (C + T ** 2) ** 2 * y0[:i] * dx[:i])
    if fit_B:
        end = (y0[-3:].mean()) / (bg[-3:].mean() or 1e-300)
        bg = bg * end
        Bv = B * 0 + end
    else:
        bg = B * bg; Bv = B
    out = bg + y[:3].mean()
    return (out[::-1] if rev else out), Bv


def gl(x, center, fwhm, eta, area):
    """Pseudo-Voigt normalised to unit area × area. eta = Lorentzian fraction (CasaXPS GL(30) ≈ eta 0.3)."""
    s = fwhm / (2 * math.sqrt(2 * math.log(2)))
    g = np.exp(-0.5 * ((x - center) / s) ** 2) / (s * math.sqrt(2 * math.pi))
    l = (fwhm / 2) / math.pi / ((x - center) ** 2 + (fwhm / 2) ** 2)
    return area * ((1 - eta) * g + eta * l)


def fit(x, y, components, background="shirley", eta=0.3, fwhm_bounds=(0.6, 3.5), fit_eta=False, tougaard_B=2866.0):
    """components: list of dicts {"name", "region" (key of DOUBLETS or SINGLETS or None), "center" (major BE guess),
    optional "fix_center"}; doublets get a minor peak at +split with area × ratio and the same FWHM.
    Returns per-component table, fitted curves, background, residual."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    if background == "shirley":
        bg = shirley(x, y)
    elif background == "tougaard":
        bg, _ = tougaard(x, y, tougaard_B)
    elif background == "linear":
        bg = np.interp(x, [x[0], x[-1]], [y[:3].mean(), y[-3:].mean()])
    else:
        bg = np.zeros_like(y)
    s = y - bg
    tot = abs(np.trapezoid(s, x)) or 1.0
    p0, lb, ub = [], [], []
    for c in components:
        p0 += [c["center"], 1.2, tot / len(components)]
        w = 0.05 if c.get("fix_center") else 1.5
        lb += [c["center"] - w, fwhm_bounds[0], 0]; ub += [c["center"] + w, fwhm_bounds[1], 10 * tot]
    if fit_eta:
        p0.append(eta); lb.append(0); ub.append(1)

    def model(p, xx=x):
        e = p[-1] if fit_eta else eta
        out, parts = np.zeros_like(xx), []
        for k, c in enumerate(components):
            cen, fw, ar = p[3 * k: 3 * k + 3]
            d = DOUBLETS.get(c.get("region") or "")
            pk = gl(xx, cen, fw, e, ar)
            if d:
                pk2 = gl(xx, cen + d["split"], fw, e, ar * d["ratio"])
                parts.append((pk, pk2)); out = out + pk + pk2
            else:
                parts.append((pk, None)); out = out + pk
        return out, parts

    r = least_squares(lambda p: model(p)[0] - s, p0, bounds=(lb, ub), max_nfev=20000)
    fitc, parts = model(r.x)
    e = r.x[-1] if fit_eta else eta
    rows = []
    for k, c in enumerate(components):
        cen, fw, ar = r.x[3 * k: 3 * k + 3]
        d = DOUBLETS.get(c.get("region") or "")
        total = ar * (1 + d["ratio"]) if d else ar
        rows.append({"component": c["name"], "region": c.get("region"), "BE major (eV)": cen,
                     "BE minor (eV)": cen + d["split"] if d else None, "FWHM (eV)": fw, "area (total)": total,
                     "split (eV)": d["split"] if d else None, "area ratio": d["ratio"] if d else None})
    tot_area = sum(rw["area (total)"] for rw in rows) or 1
    for rw in rows:
        rw["% of envelope"] = 100 * rw["area (total)"] / tot_area
    resid = s - fitc
    return {"rows": rows, "bg": bg, "fit": fitc + bg, "parts": parts, "resid": resid, "eta": e,
            "rms_rel": float(np.sqrt(np.mean(resid ** 2)) / (s.max() or 1)), "x": x, "y": y}


def atomic_percent(areas: dict, rsf: dict | None = None):
    """areas: {"Zr 3d": A, "O 1s": A, ...} (whole doublets). at.% = (A/RSF) / Σ(A/RSF) × 100."""
    rsf = rsf or RSF_SCOFIELD
    n = {k: v / rsf[k] for k, v in areas.items() if k in rsf and rsf[k]}
    t = sum(n.values()) or 1
    return {k: 100 * v / t for k, v in n.items()}


def guess_region(x):
    x = np.asarray(x, float); lo, hi = x.min(), x.max(); mid = 0.5 * (lo + hi)
    for k, d in list(DOUBLETS.items()) + [(k, {"window": w}) for k, w in SINGLETS.items()]:
        a, b = d["window"]
        if a <= mid <= b:
            return k
    return None


def draw_fit(ax, f, colors=None, show_resid=True, label_comp=True):
    x = f["x"]
    ax.plot(x, f["y"], "o", ms=1.8, mfc="none", mec="0.25", mew=0.5, label="data")
    ax.plot(x, f["bg"], "-", color="0.55", lw=0.7, label="background")
    cols = colors or ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
    for k, (rw, (p1, p2)) in enumerate(zip(f["rows"], f["parts"])):
        c = cols[k % len(cols)]
        ax.fill_between(x, f["bg"], f["bg"] + p1, color=c, alpha=0.35, lw=0)
        ax.plot(x, f["bg"] + p1, color=c, lw=0.7, label=rw["component"])
        if p2 is not None:
            ax.fill_between(x, f["bg"], f["bg"] + p2, color=c, alpha=0.18, lw=0)
            ax.plot(x, f["bg"] + p2, color=c, lw=0.7, ls="--")
    ax.plot(x, f["fit"], "-", color="k", lw=0.9, label="envelope")
    if show_resid:
        sc = (f["y"].max() - f["bg"].min()) * 0.08
        ax.plot(x, f["bg"].min() - sc + f["resid"], "-", color="0.6", lw=0.5)
    ax.set_xlim(max(x), min(x)); ax.set_xlabel("Binding energy (eV)"); ax.set_ylabel("Intensity (counts s$^{-1}$)")
