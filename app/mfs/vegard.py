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
