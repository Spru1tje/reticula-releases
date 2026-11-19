r"""v2.18.4: space groups for the structure-solution workflow of a powder pattern.

  * the symmetry operations of the 230 space groups (pymatgen tables, standard settings); orthorhombic groups in the six axis settings;
  * the systematic absences (reflection conditions) of a group, computed from its operations: a reflection h is absent when the sum of exp(2*pi*i h.t)
    over the operations (R, t) that leave h unchanged (R^T h = h) vanishes - this gives the centring, glide-plane and screw-axis conditions in one rule;
  * the ranking of the groups of a crystal system against the extracted (Le Bail) intensities: a group is consistent when the reflections it forbids are
    really absent in the pattern; among the consistent groups the one that forbids the most observed reflections (the most informative extinction symbol)
    comes first, centrosymmetric and higher symmetry before lower.
Powder patterns overlap, so only ISOLATED forbidden reflections count as evidence."""
from __future__ import annotations

import itertools
import warnings
from functools import lru_cache

import numpy as np

SYSTEM_NUMBERS = {"triclinic": (1, 2), "monoclinic": (3, 15), "orthorhombic": (16, 74), "tetragonal": (75, 142),
                  "hexagonal": (143, 194), "cubic": (195, 230)}
PERMS = [(0, 1, 2), (1, 0, 2), (2, 0, 1), (2, 1, 0), (1, 2, 0), (0, 2, 1)]       # new axis m = standard axis perm[m]
PERM_LABEL = {(0, 1, 2): "abc", (1, 0, 2): "ba-c", (2, 0, 1): "cab", (2, 1, 0): "-cba", (1, 2, 0): "bca", (0, 2, 1): "a-cb"}


@lru_cache(maxsize=None)
def _group(number):
    warnings.filterwarnings("ignore")
    from pymatgen.symmetry.groups import SpaceGroup
    sg = SpaceGroup.from_int_number(int(number))
    ops = sorted(sg.symmetry_ops, key=lambda o: (tuple(np.round(o.rotation_matrix.ravel(), 2)), tuple(np.round(np.mod(o.translation_vector, 1.0), 4))))
    R = np.array([np.rint(o.rotation_matrix) for o in ops], float)
    t = np.array([np.mod(o.translation_vector, 1.0) for o in ops], float)
    return sg.symbol, R, t, sg.crystal_system


def group_symbol(number):
    return _group(number)[0]


def operations(number, perm=(0, 1, 2)):
    """(R (n, 3, 3), t (n, 3)) in the axis setting `perm` (identity = standard setting)."""
    _, R, t, _ = _group(number)
    if tuple(perm) == (0, 1, 2):
        return R.copy(), t.copy()
    P = np.zeros((3, 3))
    for m, k in enumerate(perm):
        P[m, k] = 1.0
    return np.einsum("ij,njk,lk->nil", P, R, P), np.mod(t @ P.T, 1.0)


def is_centrosymmetric(number):
    _, R, t, _ = _group(number)
    return bool(any(np.allclose(r, -np.eye(3)) for r in R))


def order(number):
    return len(_group(number)[1])


def candidates(system):
    """All space groups (with the axis settings that matter) of a crystal system as a list of dicts {number, symbol, perm, label, centro, n_ops}."""
    s = {"rhombohedral": "hexagonal", "trigonal": "hexagonal"}.get(system, system)
    lo, hi = SYSTEM_NUMBERS[s]
    out = []
    for n in range(lo, hi + 1):
        sym = group_symbol(n)
        perms = PERMS if s == "orthorhombic" else [(0, 1, 2)]
        seen = set()
        for pm in perms:
            R, t = operations(n, pm)
            key = _conditions_key(R, t)
            if key in seen:                                          # settings with the same conditions are the same group for a powder
                continue
            seen.add(key)
            out.append({"number": n, "symbol": sym, "perm": pm, "label": f"{sym}" + (f" ({PERM_LABEL[pm]})" if pm != (0, 1, 2) else ""),
                        "centro": is_centrosymmetric(n), "n_ops": len(R)})
    return out


def _conditions_key(R, t, hmax=3):
    """A hashable summary of the absences of a setting (the forbidden set for all |h|,|k|,|l| <= hmax)."""
    g = np.arange(-hmax, hmax + 1)
    H = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3).astype(float)
    return forbidden(R, t, H).tobytes()


def forbidden(R, t, hkl, tol=1e-6):
    """Boolean mask: reflection systematically absent. R (n,3,3), t (n,3), hkl (m,3) integers (as float)."""
    hkl = np.asarray(hkl, float)
    S = np.zeros(len(hkl), complex)
    for r, tt in zip(R, t):
        stab = np.all(np.abs(hkl @ r - hkl) < 1e-9, axis=1)
        if stab.any():
            S[stab] += np.exp(2j * np.pi * (hkl[stab] @ tt))
    return np.abs(S) < tol


def point_operations(number, perm=(0, 1, 2)):
    """One operation per distinct rotation part (the centring and other pure translations are dropped): R (k,3,3), t (k,3). With the extinction mask
    these are everything a structure-factor symmetrisation needs."""
    R, t = operations(number, perm)
    seen, Ro, to = [], [], []
    for r, tt in zip(R, t):
        key = tuple(np.rint(r.ravel()).astype(int))
        if key in seen:
            continue
        seen.append(key)
        Ro.append(r)
        to.append(tt)
    return np.array(Ro), np.array(to)


def rank(system, groups, min_abs=1):
    """Ranks the space groups of `system` against the pattern.
    groups: list of dicts {hkl (m,3) all members of one reflection group, I (extracted area), thr (detection threshold, same unit), iso (bool: isolated)}.
    Returns a list of dicts (best first): the candidate dict plus n_abs (isolated forbidden groups that are absent), n_viol (isolated forbidden groups that are
    present), n_forb (all fully forbidden groups), consistent (bool), score."""
    out = []
    iso = [g for g in groups if g.get("iso", True)]
    allh = np.concatenate([g["hkl"] for g in iso]) if iso else np.zeros((0, 3))
    owner = np.concatenate([np.full(len(g["hkl"]), i) for i, g in enumerate(iso)]) if iso else np.zeros(0, int)
    for c in candidates(system):
        R, t = operations(c["number"], c["perm"])
        f = forbidden(R, t, allh) if len(allh) else np.zeros(0, bool)
        n_abs = n_viol = 0
        viol_weight = 0.0
        for i, g in enumerate(iso):
            sel = owner == i
            if sel.any() and f[sel].all():
                if g["I"] > g["thr"]:
                    n_viol += 1
                    viol_weight += min(g["I"] / max(g["thr"], 1e-12), 20.0)
                else:
                    n_abs += 1
        consistent = n_viol == 0 or (n_viol <= 1 and n_viol <= 0.1 * (n_abs + n_viol))
        c2 = {**c, "n_abs": n_abs, "n_viol": n_viol, "viol_weight": viol_weight, "consistent": bool(consistent)}
        c2["score"] = (1000.0 if consistent else 0.0) + n_abs - 3.0 * viol_weight + (0.5 if c["centro"] else 0.0) + 0.01 * c["n_ops"]
        out.append(c2)
    out.sort(key=lambda d: -d["score"])
    return out
