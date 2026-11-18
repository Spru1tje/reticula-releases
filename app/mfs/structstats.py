"""v2.16: statistics of a refined crystal structure - bond lengths, bond angles and coordination with their estimated standard deviations (propagated
from the covariance matrix of the refinement), bond-valence sums and a plausibility verdict per atom.

All quantities are functions of the refined parameters (position shifts, cell); the esd of every distance / angle is sqrt(J C J^T) with J the numerical
derivative with respect to the refined parameters and C their covariance."""
from __future__ import annotations

import math

import numpy as np

from mfs.structmodel import _covalent, is_bond

# Brese & O'Keeffe (1991) bond-valence parameters R0 (Å) for cation–O bonds, b = 0.37 Å; formal oxidation state of the cation
BV_R0 = {"Zr": (1.937, 4), "Hf": (1.923, 4), "Ti": (1.815, 4), "Ce": (2.028, 4), "Zn": (1.704, 2), "Cu": (1.679, 2), "Al": (1.651, 3), "Fe": (1.759, 3),
         "Ca": (1.967, 2), "Mg": (1.693, 2), "P": (1.617, 5), "Si": (1.624, 4), "S": (1.624, 6), "Cr": (1.724, 3), "Mn": (1.790, 2), "Co": (1.692, 2),
         "Ni": (1.654, 2), "V": (1.743, 4), "La": (2.172, 3), "Y": (2.019, 3), "Zr_": (1.937, 4)}
BV_B = 0.37


def _images():
    return np.array([(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)], float)


def neighbours(model, p, G, extra: float = 0.45, dmax: float = 4.2):
    """For every representative atom the bonded neighbours: dict rep_atom_index -> list of (j atom, image n, distance). `G` = cell metric (Å²)."""
    X = model.atom_positions(p)
    rad = np.array([_covalent(e) for e in model.atom_el])
    imgs = _images()
    out = {}
    reps_atoms = [int(np.nonzero(model.atom_site == r)[0][0]) for r in model.reps]
    for i in reps_atoms:
        lst = []
        for n in imgs:
            d = (X + n) - X[i]
            dd = np.sqrt(np.maximum(np.einsum("ni,ij,nj->n", d, G, d), 0.0))
            for j in np.nonzero((dd < dmax) & (dd > 0.3))[0]:
                if (int(j) == i and not n.any()) or not is_bond(model.atom_el[i], model.atom_el[int(j)], float(dd[j])):
                    continue
                lst.append((int(j), n.copy(), float(dd[j])))
        out[i] = sorted(lst, key=lambda t: t[2])
    return out


def _angle(a, b, c, G):
    """Angle at b between the vectors b->a and b->c (fractional vectors a, c relative to b), degrees."""
    u, v = a, c
    cu = u @ G @ v / math.sqrt(max((u @ G @ u) * (v @ G @ v), 1e-30))
    return math.degrees(math.acos(max(-1.0, min(1.0, cu))))


def _propagate(f, names, base, cov, steps=None):
    """esd of the vector-valued f(values dict) from the covariance of the named parameters (numerical Jacobian)."""
    f0 = np.asarray(f(base), float)
    if not names:
        return np.zeros_like(f0)
    J = np.zeros((len(f0), len(names)))
    for k, n in enumerate(names):
        h = 1e-6 * max(1.0, abs(base[n]))
        b2 = dict(base)
        b2[n] = base[n] + h
        J[:, k] = (np.asarray(f(b2), float) - f0) / h
    var = np.einsum("ik,kl,il->i", J, cov, J)
    return np.sqrt(np.maximum(var, 0.0))


def report(model, p, cell_metric_fn, cov=None, names=None, base=None, labels=None, extra: float = 0.45):
    """Bond lengths, angles, coordination and bond-valence sums of the representative atoms.
    cell_metric_fn(values dict) -> 3x3 metric; `base` = values dict of ALL refined parameters; `cov`/`names` = covariance matrix and its parameter names;
    the position parameters are the entries 'p0'..'p{n-1}' of `base` (the shifts of model.pos_params)."""
    base = dict(base or {})
    for k, v in enumerate(p):
        base.setdefault(f"p{k}", float(v))
    names = [n for n in (names or []) if n in base]
    C = None
    if cov is not None and names:
        C = np.asarray(cov, float)
    labs = labels or model.labels()
    pvec = lambda vals: np.array([vals[f"p{k}"] for k in range(len(p))])
    G0 = cell_metric_fn(base)
    nb = neighbours(model, pvec(base), G0, extra)
    bonds, angles, coord, bvs = [], [], [], []
    seen_b, seen_a = {}, {}
    for i, lst in nb.items():
        li = labs[int(model.atom_rep_site[i])]
        if not lst:
            coord.append({"atom": li, "CN": 0, "mean_d": float("nan")})
            continue
        for (j, n, d) in lst:
            lj = labs[int(model.atom_rep_site[j])]
            key = (li, lj, round(d, 3))
            if key in seen_b:
                seen_b[key]["count"] += 1
                continue

            def fd(vals, i=i, j=j, n=n):
                X = model.atom_positions(pvec(vals))
                dv = (X[j] + n) - X[i]
                return [math.sqrt(max(dv @ cell_metric_fn(vals) @ dv, 0.0))]
            e = float(_propagate(fd, [x for x in names], base, C)[0]) if C is not None else float("nan")
            rec = {"atom1": li, "atom2": lj, "distance": d, "esd": e, "count": 1}
            seen_b[key] = rec
            bonds.append(rec)
        # angles at the central atom between its neighbours (within the first shell: bonded set)
        X = model.atom_positions(pvec(base))
        vecs = [(j, (X[j] + n) - X[i], d) for (j, n, d) in lst]
        for a in range(len(vecs)):
            for b in range(a + 1, len(vecs)):
                ja, va, da = vecs[a]
                jb, vb, db = vecs[b]
                ang = _angle(va, 0, vb, G0)
                la, lb = labs[int(model.atom_rep_site[ja])], labs[int(model.atom_rep_site[jb])]
                key = (la, li, lb, round(ang, 1))
                if key in seen_a or (lb, li, la, round(ang, 1)) in seen_a:
                    continue

                def fa(vals, i=i, ja=ja, jb=jb, na=lst[a][1], nb_=lst[b][1]):
                    Xv = model.atom_positions(pvec(vals))
                    return [_angle((Xv[ja] + na) - Xv[i], 0, (Xv[jb] + nb_) - Xv[i], cell_metric_fn(vals))]
                e = float(_propagate(fa, [x for x in names], base, C)[0]) if C is not None else float("nan")
                rec = {"atom1": la, "centre": li, "atom2": lb, "angle": ang, "esd": e}
                seen_a[key] = rec
                angles.append(rec)
        wts = [float(model.atom_occ[j]) for (j, n, d) in lst]               # a partially occupied neighbour counts with its occupancy
        cn = round(sum(wts), 2)
        cn = int(cn) if abs(cn - round(cn)) < 1e-6 else cn
        coord.append({"atom": li, "CN": cn, "mean_d": float(np.average([d for _, _, d in lst], weights=np.maximum(wts, 1e-9)))})
        el = model.atom_el[i]
        if el in BV_R0:
            R0, ox = BV_R0[el]
            s = sum(float(model.atom_occ[j]) * math.exp((R0 - d) / BV_B) for (j, n, d) in lst if model.atom_el[j] == "O")
            if s > 0:
                bvs.append({"atom": li, "BVS": s, "formal": ox, "deviation_pct": 100.0 * (s - ox) / ox})
    return {"bonds": bonds, "angles": angles, "coordination": coord, "bvs": bvs}


def verdict(rep: dict) -> list[dict]:
    """Plain-language plausibility notes: bond-valence sums far from the formal charge, very short / long bonds, angles far from typical."""
    out = []
    for b in rep.get("bvs", []):
        if abs(b["deviation_pct"]) > 25:
            out.append({"level": "warn", "text": f"{b['atom']}: bond-valence sum {b['BVS']:.2f} v.u. vs formal {b['formal']} ({b['deviation_pct']:+.0f} %) – the local geometry is questionable "
                                                 "(or the oxidation state / coordination differs from the assumed one)."})
        elif abs(b["deviation_pct"]) > 12:
            out.append({"level": "info", "text": f"{b['atom']}: bond-valence sum {b['BVS']:.2f} v.u. vs formal {b['formal']} ({b['deviation_pct']:+.0f} %)."})
    short = [b for b in rep.get("bonds", []) if b["distance"] < 0.9 * (_covalent(b["atom1"].rstrip("0123456789/")) + _covalent(b["atom2"].rstrip("0123456789/")))]
    for b in short[:3]:
        out.append({"level": "warn", "text": f"{b['atom1']}–{b['atom2']} = {b['distance']:.3f} Å is shorter than expected for a covalent bond."})
    if not out:
        out.append({"level": "ok", "text": "Bond lengths and bond-valence sums are within normal ranges."})
    return out
