"""v2.19.1: Vegard's law for solid solutions: composition from the lattice parameter, and the lattice parameter from the composition.

Vegard (1921): a(x) = x a_A + (1 - x) a_B, where x is the fraction of end member A and a_A, a_B are the lattice parameters of the
end members. Many solid solutions bend away from the straight line; the bowing parameter b gives
a(x) = x a_A + (1 - x) a_B - b x (1 - x) (Denton and Ashcroft 1991, Phys. Rev. A 43, 3161).
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
