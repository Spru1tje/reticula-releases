"""v2.19.1: Vegard's law for solid solutions: composition from the lattice parameter, and the lattice parameter from the composition.

Vegard (1921): a(x) = x a_A + (1 - x) a_B, where x is the fraction of end member A and a_A, a_B are the lattice parameters of the
end members. Many solid solutions bend away from the straight line; the bowing parameter b gives
a(x) = x a_A + (1 - x) a_B - b x (1 - x). The bibliographic entry Denton and Ashcroft 1991, Phys. Rev. A 43, 3161
(doi:10.1103/PhysRevA.43.3161) was checked: authors, journal, volume and first page match. The bowed FORMULA in this
module is NOT confirmed against that paper (its abstract is not in the landing-page metadata). The PXRD panel uses b = 0 only.
Inverting the bowed form is a quadratic in x; the root between 0 and 1 is returned.
General practice: Vegard is an approximation. Check the bowing parameter against a reference series, and do not read a composition
from a lattice parameter that lies outside the range of the two end members."""
from __future__ import annotations

import math


def lattice_from_composition(x: float, a_A: float, a_B: float, bowing: float = 0.0) -> float:
    """Lattice parameter (same unit as a_A, a_B) of the solid solution with fraction x of end member A (0 <= x <= 1)."""
    if not (0.0 <= float(x) <= 1.0):
        raise ValueError("x must be between 0 and 1")
    return float(x) * a_A + (1.0 - float(x)) * a_B - float(bowing) * float(x) * (1.0 - float(x))


def composition_from_lattice(a: float, a_A: float, a_B: float, bowing: float = 0.0) -> dict:
    """Fraction x of end member A from a measured lattice parameter a.

    Returns dict(status, x) or dict(status, reason). status 'ok' only when the root lies in [0, 1]; a value outside the
    end-member range is reported as 'outside range' and no composition is returned."""
    a, a_A, a_B, b = float(a), float(a_A), float(a_B), float(bowing)
    if a_A == a_B and b == 0.0:
        return {"status": "no composition range", "reason": "both end members have the same lattice parameter"}
    # b x^2 + (a_A - a_B - b) x + (a_B - a) = 0
    if b == 0.0:
        roots = [(a - a_B) / (a_A - a_B)]
    else:
        qa, qb, qc = b, a_A - a_B - b, a_B - a
        disc = qb * qb - 4 * qa * qc
        if disc < 0:
            return {"status": "outside range", "reason": "no real composition for this lattice parameter"}
        sq = math.sqrt(disc)
        roots = [(-qb + sq) / (2 * qa), (-qb - sq) / (2 * qa)]
    inside = [r for r in roots if -1e-9 <= r <= 1 + 1e-9]
    if not inside:
        return {"status": "outside range", "reason": "the lattice parameter lies outside the two end members"}
    x = min(max(inside[0], 0.0), 1.0)
    return {"status": "ok", "x": x, "check_a": lattice_from_composition(x, a_A, a_B, b)}


# ---------------------------------------------------------------- v2.23: shift of reflections for the PXRD page (cubic end members)
def two_theta_for_lattice(a: float, d0_hkl_norm: float, wavelength: float) -> float:
    """2θ (degrees) of a cubic reflection with h²+k²+l² = d0_hkl_norm when the cubic cell parameter is a (Å).
    d = a / sqrt(h²+k²+l²); Bragg: λ = 2 d sin θ. Returns nan when the reflection cannot occur at this wavelength."""
    if a <= 0 or d0_hkl_norm <= 0:
        return float("nan")
    d = float(a) / math.sqrt(float(d0_hkl_norm))
    s = float(wavelength) / (2.0 * d)
    if s > 1.0:
        return float("nan")
    return math.degrees(2.0 * math.asin(s))


def shifted_two_theta(two_theta_ref: float, a_ref: float, a_new: float) -> float:
    """Shift a reflection from the reference cell (a_ref) to the new cell (a_new), same hkl: d scales with a.
    Pure scaling of the cubic lattice; no wavelength needed (d_new = d_ref · a_new / a_ref)."""
    if a_ref <= 0 or a_new <= 0 or not math.isfinite(two_theta_ref):
        return float("nan")
    th = math.radians(float(two_theta_ref) / 2.0)
    d_ref = 1.0 / (2.0 * math.sin(th))            # λ = 1 unit: only the ratio matters below
    d_new = d_ref * float(a_new) / float(a_ref)
    s = 1.0 / (2.0 * d_new)
    if s > 1.0:
        return float("nan")
    return math.degrees(2.0 * math.asin(s))


def x_from_two_theta(two_theta_obs: float, hkl_norm: float, wavelength: float, a_A: float, a_B: float, bowing: float = 0.0) -> dict:
    """Estimate the fraction x of end member A from one observed reflection (cubic, h²+k²+l² = hkl_norm).
    The cell is a_obs = wavelength · sqrt(h²+k²+l²) / (2 sin θ); then composition_from_lattice. Outside the end-member
    range nothing is estimated (status 'outside range', as in composition_from_lattice)."""
    if not math.isfinite(two_theta_obs) or two_theta_obs <= 0:
        return {"status": "no peak", "reason": "2θ must be a positive angle"}
    s = math.sin(math.radians(two_theta_obs) / 2.0)
    if s <= 0:
        return {"status": "no peak", "reason": "2θ must be a positive angle"}
    a_obs = float(wavelength) * math.sqrt(float(hkl_norm)) / (2.0 * s)
    out = composition_from_lattice(a_obs, a_A, a_B, bowing)
    out["a_obs"] = a_obs
    return out


def x_from_peaks(peaks, wavelength: float, a_A: float, a_B: float, bowing: float = 0.0) -> dict:
    """v2.23.1 (S-Structuur-2): composition x of end member A from SEVERAL observed cubic reflections.
    peaks: list of (two_theta_obs, hkl_norm). Each peak gives one x (x_from_two_theta). Returns the mean, the sample
    standard deviation and the individual values. Any peak outside the end-member range makes the result 'outside range'
    (no mean is returned): a partial answer would hide a bad peak. Fewer than two valid peaks gives 'too few peaks'."""
    vals, bad = [], []
    for th, n in peaks:
        r = x_from_two_theta(float(th), float(n), wavelength, a_A, a_B, bowing)
        if r["status"] == "ok":
            vals.append(float(r["x"]))
        else:
            bad.append({"two_theta": float(th), "status": r["status"]})
    if bad:
        return {"status": "outside range", "bad_peaks": bad, "n_ok": len(vals)}
    if len(vals) < 2:
        return {"status": "too few peaks", "n_ok": len(vals)}
    mean = sum(vals) / len(vals)
    sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / (len(vals) - 1))
    return {"status": "ok", "x": mean, "sd": sd, "n": len(vals), "values": vals}


# ---------------------------------------------------------------- v2.23.1 (S-Structuur-3): cubic, tetragonal, hexagonal end members
def system_of(cell) -> str:
    """'cubic', 'tetragonal', 'hexagonal' or '' for a cell (a, b, c, alpha, beta, gamma), tolerance 1e-6 (degrees and Å)."""
    a, b, c, al, be, ga = (float(v) for v in cell)
    if all(abs(v - 90.0) < 1e-6 for v in (al, be, ga)) and abs(a - b) < 1e-6 and abs(b - c) < 1e-6:
        return "cubic"
    if all(abs(v - 90.0) < 1e-6 for v in (al, be, ga)) and abs(a - b) < 1e-6:
        return "tetragonal"
    if abs(a - b) < 1e-6 and abs(al - 90.0) < 1e-6 and abs(be - 90.0) < 1e-6 and abs(ga - 120.0) < 1e-6:
        return "hexagonal"
    return ""


def d_spacing(cell, hkl) -> float:
    """d (Å) of reflection hkl in a cubic, tetragonal or hexagonal cell (axes as in system_of). Other systems: nan."""
    s = system_of(cell)
    h, k, l = (float(v) for v in hkl)
    if s == "cubic":
        a = float(cell[0])
        inv2 = (h * h + k * k + l * l) / (a * a)
    elif s == "tetragonal":
        a, c = float(cell[0]), float(cell[2])
        inv2 = (h * h + k * k) / (a * a) + (l * l) / (c * c)
    elif s == "hexagonal":
        a, c = float(cell[0]), float(cell[2])
        inv2 = (4.0 / 3.0) * (h * h + h * k + k * k) / (a * a) + (l * l) / (c * c)
    else:
        return float("nan")
    return math.sqrt(1.0 / inv2) if inv2 > 0 else float("nan")


def two_theta_for_cell(cell, hkl, wavelength: float) -> float:
    """2θ (degrees) of hkl in the given cell; nan when the reflection cannot occur at this wavelength."""
    d = d_spacing(cell, hkl)
    if not math.isfinite(d):
        return float("nan")
    s = float(wavelength) / (2.0 * d)
    return math.degrees(2.0 * math.asin(s)) if s <= 1.0 else float("nan")


def cell_from_composition(x: float, cell_A, cell_B):
    """Linear Vegard per axis for two end members of the SAME system: a(x) and c(x) (no bowing). Angles are taken from A.
    Raises ValueError for different systems."""
    sA, sB = system_of(cell_A), system_of(cell_B)
    if not sA or sA != sB:
        raise ValueError(f"end members must be of the same cubic, tetragonal or hexagonal system ({sA!r} vs {sB!r})")
    if not (0.0 <= float(x) <= 1.0):
        raise ValueError("x must be between 0 and 1")
    out = []
    for i in (0, 2):  # a and c (b = a for these systems)
        out.append(float(x) * float(cell_A[i]) + (1.0 - float(x)) * float(cell_B[i]))
    a_x, c_x = out
    return (a_x, a_x, c_x, float(cell_A[3]), float(cell_A[4]), float(cell_A[5]))


# ---------------------------------------------------------------- v2.24.0: x from several peaks, any system (grid search)
def estimate_x_grid(peaks, wavelength: float, cell_A, cell_B, tol_deg: float = 0.05, steps: int = 1001) -> dict:
    """Composition x (fraction of A) from observed reflections of a cubic, tetragonal or hexagonal solid solution.
    peaks: list of (two_theta_obs, (h, k, l)). For each x on a grid the cell a(x), c(x) is linear (cell_from_composition);
    the 2θ of every reflection is predicted; the x with the smallest sum of squares wins. rms = root-mean-square residual (°).
    Status: 'ok' (interior minimum with rms <= tol_deg); 'end member' (minimum at x = 0 or 1, rms <= tol_deg);
    'outside range' (minimum at an end and rms > tol_deg: the peaks are not explained by A and B); 'too few peaks'."""
    obs = [(float(t), tuple(float(v) for v in hkl)) for t, hkl in peaks if math.isfinite(float(t))]
    if len(obs) < 2:
        return {"status": "too few peaks", "n": len(obs)}
    best = None
    xs = [i / (steps - 1) for i in range(steps)]
    for x in xs:
        cell = cell_from_composition(x, cell_A, cell_B)
        pred = [two_theta_for_cell(cell, hkl, wavelength) for _, hkl in obs]
        if any(not math.isfinite(p) for p in pred):
            continue
        ss = sum((p - t) ** 2 for p, (t, _) in zip(pred, obs))
        if best is None or ss < best[0]:
            best = (ss, x)
    if best is None:
        return {"status": "outside range", "reason": "no composition predicts these reflections"}
    ss, x = best
    rms = math.sqrt(ss / len(obs))
    if x <= xs[0] or x >= xs[-1]:
        status = "end member" if rms <= tol_deg else "outside range"
    else:
        status = "ok" if rms <= tol_deg else "outside range"
    return {"status": status, "x": x, "rms_deg": rms, "n": len(obs)}


# ---------------------------------------------------------------- v2.24.0: end members from numbers or from a pattern
def symmetric_cell(system: str, a: float, c: float | None = None) -> tuple:
    """Exact cell of a system from its free parameters (removes the small asymmetry a refinement leaves):
    cubic (a), tetragonal (a, c), hexagonal (a, c; gamma = 120)."""
    a = float(a)
    if system == "cubic":
        return (a, a, a, 90.0, 90.0, 90.0)
    if system == "tetragonal":
        return (a, a, float(c), 90.0, 90.0, 90.0)
    if system == "hexagonal":
        return (a, a, float(c), 90.0, 90.0, 120.0)
    raise ValueError(f"unsupported system {system!r} (cubic, tetragonal, hexagonal)")


def cell_from_pattern(x, y, lo: float, hi: float, wavelength: float = 1.540593) -> dict:
    """Start cell by indexing the peaks of a pattern (cubic, tetragonal, hexagonal), then a Le Bail fit with the cell-only
    route (Phase.from_cell + quick_refine, mode 'lebail'). Returns status, the symmetric refined cell, its system, Rwp (%),
    and the indexing solution used. Status 'no cell' when indexing finds nothing."""
    import numpy as np
    from mfs import phaseid as PID
    from mfs import indexing as IX
    from mfs import rietveld as RV
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = (x >= lo) & (x <= hi)
    x, y = x[m], y[m]
    pk, _, _ = PID.exp_peaks(x, y)
    peaks = [float(v) for v in np.asarray(pk, float)]
    if len(peaks) < 4:
        return {"status": "no cell", "reason": "too few peaks in the range"}
    sols = IX.index(peaks, wavelength, systems=("cubic", "tetragonal", "hexagonal"), max_solutions=3)
    if not sols:
        return {"status": "no cell", "reason": "indexing found no cell"}
    top = sols[0]
    system = top["system"]
    if system not in ("cubic", "tetragonal", "hexagonal"):
        return {"status": "no cell", "reason": f"indexing gave a {system} cell (not supported by the panel)"}
    cell0 = tuple(float(v) for v in top["cell"])
    ph = RV.Phase.from_cell(cell0, "end member")
    res = RV.quick_refine(x, y, [ph], radiation="Cu Kα1", mode="lebail", tt_range=(lo, hi))
    refined = res["phases"][0]["cell"]
    if system == "cubic":
        a = (refined[0] + refined[1] + refined[2]) / 3.0
        cell = symmetric_cell("cubic", a)
    else:
        a = (refined[0] + refined[1]) / 2.0
        cell = symmetric_cell(system, a, refined[2])
    return {"status": "ok", "cell": cell, "system": system, "rwp": 100.0 * float(res["stats"]["Rwp"]), "index": top}
