r"""v2.18.4: small crystallographic helpers shared by the structure-solution workflow: X-ray scattering factors, structure factors of an atom list,
the Lorentz-polarisation factor of a powder pattern."""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np


@lru_cache(maxsize=None)
def _coeff(sym):
    from pymatgen.analysis.diffraction.xrd import ATOMIC_SCATTERING_PARAMS
    from pymatgen.core import Element
    return np.asarray(ATOMIC_SCATTERING_PARAMS[sym], float), Element(sym).Z


def f_xray(sym, s2):
    """X-ray scattering factor of an element (as in pymatgen: Z - 41.78214 s^2 sum a_i exp(-b_i s^2)), s2 = (sin theta / lambda)^2."""
    co, Z = _coeff(sym)
    s2 = np.asarray(s2, float)
    return Z - 41.78214 * s2 * np.sum(co[:, 0] * np.exp(-co[:, 1] * s2[:, None]), axis=1)


def atomic_number(sym):
    return _coeff(sym)[1]


def structure_factors(frac, elems, occ, hkl, Gs, B=0.0):
    """Complex structure factors of an atom list. frac (n,3); elems (n,); occ (n,); hkl (m,3); Gs = inverse metric tensor of the cell."""
    hkl = np.asarray(hkl, float)
    inv_d2 = np.einsum("ki,ij,kj->k", hkl, Gs, hkl)
    s2 = inv_d2 / 4.0
    F = np.zeros(len(hkl), complex)
    frac = np.asarray(frac, float)
    occ = np.asarray(occ, float)
    elems = np.asarray(elems)
    for sym in sorted(set(elems.tolist())):
        sel = elems == sym
        fs = f_xray(sym, s2) * np.exp(-B * s2)
        ph = np.exp(2j * np.pi * (hkl @ frac[sel].T))
        F += fs * (ph @ occ[sel])
    return F


def lorentz_polarisation(two_theta_deg, mono_2theta=0.0):
    th = np.radians(np.asarray(two_theta_deg, float) / 2.0)
    c2m = math.cos(math.radians(mono_2theta)) ** 2
    return (1.0 + c2m * np.cos(2 * th) ** 2) / ((1.0 + c2m) * np.sin(th) ** 2 * np.cos(th))
