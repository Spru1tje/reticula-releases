"""1H peak ↔ atom assignment and impurity check (v2.13).

assign(): predicted proton groups (symmetry-equivalent H-bearing heavy atoms, mfs.nmrpredict) are matched to observed
multiplets (mfs.nmrpeaks) by a small integer programme: shift error / prediction uncertainty, integral capacity (nH),
multiplicity vs. vicinal-H count; groups may stay unassigned (OH/NH often missing). Manual corrections (drag) are hard
constraints. impurity_check(): unassigned multiplets against the Fulmer 2010 / Gottlieb 1997 trace-impurity table.
"""
from __future__ import annotations

import json
import os

import numpy as np

from mfs import nmrpredict as HP

IMP_PATH = os.path.join(os.path.dirname(__file__), "data", "impurities_1h.json")
_IMP = None
SIGMA = {"high": 0.2, "medium": 0.35, "low": 0.5, "none": 2.0}
# impurity SMILES (for molar mass → wt %)
IMP_SMILES = {"water": "O", "acetic acid": "CC(=O)O", "acetone": "CC(C)=O", "acetonitrile": "CC#N", "benzene": "c1ccccc1",
              "tert-butyl alcohol": "CC(C)(C)O", "chloroform": "ClC(Cl)Cl", "18-crown-6": "C1COCCOCCOCCOCCOCCO1",
              "cyclohexane": "C1CCCCC1", "1,2-dichloroethane": "ClCCCl", "dichloromethane": "ClCCl",
              "diethyl ether": "CCOCC", "diglyme": "COCCOCCOC", "dimethylformamide": "CN(C)C=O", "1,4-dioxane": "C1COCCO1",
              "DME": "COCCOC", "ethane": "CC", "ethanol": "CCO", "ethyl acetate": "CCOC(C)=O", "ethylene": "C=C",
              "ethylene glycol": "OCCO", "hexamethylbenzene": "Cc1c(C)c(C)c(C)c(C)c1C", "n-hexane": "CCCCCC",
              "HMDSO": "C[Si](C)(C)O[Si](C)(C)C", "HMPA": "CN(C)P(=O)(N(C)C)N(C)C", "hydrogen": "[HH]", "imidazole": "c1c[nH]cn1",
              "methane": "C", "methanol": "CO", "nitromethane": "C[N+](=O)[O-]", "n-pentane": "CCCCC", "propane": "CCC",
              "2-propanol": "CC(C)O", "propylene": "C=CC", "pyridine": "c1ccncc1", "pyrrole": "c1cc[nH]c1",
              "pyrrolidine": "C1CCNC1", "tetrahydrofuran": "C1CCOC1", "toluene": "Cc1ccccc1", "triethylamine": "CCN(CC)CC",
              "dimethyl sulfoxide": "CS(C)=O", "tetramethylsilane": "C[Si](C)(C)C"}
_MULT_NEIGH = {"s": 0, "br s": 0, "d": 1, "t": 2, "q": 3, "p": 4, "dd": 2, "dt": 3, "td": 3, "ddd": 3, "dq": 4, "qd": 4,
               "tt": 4, "sext": 5, "sept": 6}


def impurity_table():
    global _IMP
    if _IMP is None:
        with open(IMP_PATH, encoding="utf-8") as f:
            _IMP = json.load(f)
    return _IMP


def groups(mol, preds=None):
    """Symmetry-equivalent proton groups: [{gid, atoms, labels, nH, shift, conf, sigma, exch, vicinal}]."""
    preds = preds if preds is not None else HP.predict(mol)
    cls = HP.equivalence_classes(mol)
    byc = {}
    for p in preds:
        byc.setdefault(cls[p["atom"]], []).append(p)
    out = []
    for gi, (c, ps) in enumerate(sorted(byc.items(), key=lambda kv: min(p["atom"] for p in kv[1]))):
        atoms = [p["atom"] for p in ps]
        sh = [p["shift"] for p in ps if p["shift"] is not None]
        conf = min((p["conf"] for p in ps), key=lambda c: ["none", "low", "medium", "high"].index(c))
        exch = any(p["exchangeable"] for p in ps)
        vic = 0
        a0 = mol.GetAtomWithIdx(atoms[0])
        for nb in a0.GetNeighbors():
            if nb.GetSymbol() in HP.EXCHANGEABLE:
                continue
            vic += nb.GetTotalNumHs() if nb.GetIdx() not in atoms else 0
        out.append(dict(gid=gi, atoms=atoms, labels=[a + 1 for a in atoms], nH=sum(p["nH"] for p in ps),
                        shift=float(np.mean(sh)) if sh else None, conf=conf,
                        sigma=1.5 if exch else max(SIGMA[conf], 1.3 * max([p["spread"] or 0 for p in ps] or [0])),
                        exch=exch, vicinal=0 if exch else vic,
                        spread=max([p["spread"] or 0 for p in ps] or [0])))
    # v2.14: proton labels are LETTERS (Ha, Hb ... from the most downfield group), so they cannot be mistaken for the atom numbers
    for k, gr in enumerate(sorted(out, key=lambda g: -(g["shift"] if g["shift"] is not None else -99.0))):
        gr["tag"] = _letters(k)
    return out


def _letters(k):
    """0 -> a, 1 -> b ... 25 -> z, 26 -> aa ..."""
    s = ""
    k += 1
    while k:
        k, r = divmod(k - 1, 26)
        s = chr(97 + r) + s
    return s


def proton_label(g):
    return "H" + g["tag"]


def _usable(mults):
    return [i for i, m in enumerate(mults) if not m.note or m.overlap]


def _n_large(m, jmin=5.0):
    """Number of coupling partners with J ≥ jmin from the multiplet name and J list (d=1, t=2, q=3 …)."""
    name = m.mult
    if not m.J or name in ("m", "s", "br s"):
        return None
    codes = {"sept": 6, "sext": 5, "p": 4}
    if name in codes:
        return codes[name] if m.J[0] >= jmin else 0
    per = {"d": 1, "t": 2, "q": 3}
    if len(name) != len(m.J) or any(ch not in per for ch in name):
        return None
    return sum(per[ch] for ch, J in zip(name, m.J) if J >= jmin)


def assign(mol_or_groups, mults, fixed=None, banned=None, lam=3.0, lam_over=6.0):
    """Assign proton groups to multiplets. fixed: {gid: mult index} (user drag), banned: {(gid, mi)}; returns
    {gid: mult index or None}. Sets m.label on the multiplets ('H-3, H-5')."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    G = mol_or_groups if isinstance(mol_or_groups, list) else groups(mol_or_groups)
    fixed = dict(fixed or {}); banned = set(banned or ())
    M = _usable(mults)
    pairs = []
    for g in G:
        for mi in M:
            m = mults[mi]
            if (g["gid"], mi) in banned:
                continue
            if fixed.get(g["gid"]) not in (None, mi):
                continue
            if g["shift"] is None:
                d = 0.0 if g["exch"] else 9.0
            else:
                lo, hi = min(m.lo, m.hi), max(m.lo, m.hi)
                d = 0.0 if lo <= g["shift"] <= hi and m.mult == "m" else abs(g["shift"] - m.centre)
            s = g["sigma"]
            if g["exch"] and m.mult == "br s":
                d = min(d, 1.0)                                   # exchangeable H: any broad singlet is plausible
            if fixed.get(g["gid"]) != mi and d > 5 * s + 0.5:
                continue
            c = min((d / s) ** 2, 6.0) + 0.3 * max(0.0, d / s - 3)
            if m.note:                                             # coupled lines on the solvent residual
                c += 1.0
            k = _MULT_NEIGH.get(m.mult)
            if k is not None and not g["exch"]:
                if (k == 0) != (g["vicinal"] == 0):
                    c += 0.8
                elif k:
                    nl = _n_large(m)
                    c += 0.6 * min(abs((nl if nl is not None else k) - min(g["vicinal"], 6)), 2.5)
            pairs.append((g["gid"], mi, c))
    nG, nM, nP = len(G), len(M), len(pairs)
    if nG == 0:
        return {}
    gidx = {g["gid"]: i for i, g in enumerate(G)}; midx = {mi: j for j, mi in enumerate(M)}
    # variables: x_p (nP), u_g (nG), s+_m (nM), s-_m (nM)
    nv = nP + nG + 2 * nM
    c = np.zeros(nv)
    for p, (_, _, cc) in enumerate(pairs):
        c[p] = cc
    for i, g in enumerate(G):
        c[nP + i] = 0.0 if g["gid"] in fixed and fixed[g["gid"]] is None else (1.0 if g["exch"] else 6.0 + 2.0 * g["nH"])
    c[nP + nG:nP + nG + nM] = lam_over                              # s+ : more H put into a multiplet than integrated
    c[nP + nG + nM:] = lam                                         # s− : multiplet not filled
    A, lb, ub = [], [], []
    for i, g in enumerate(G):                                     # each group once (or unassigned)
        row = np.zeros(nv)
        for p, (gid, _, _) in enumerate(pairs):
            if gid == g["gid"]:
                row[p] = 1
        row[nP + i] = 1
        A.append(row); lb.append(1); ub.append(1)
        if g["gid"] in fixed and fixed[g["gid"]] is not None:
            r2 = np.zeros(nv); r2[nP + i] = 1; A.append(r2); lb.append(0); ub.append(0)
        if g["gid"] in fixed and fixed[g["gid"]] is None:
            r2 = np.zeros(nv); r2[nP + i] = 1; A.append(r2); lb.append(1); ub.append(1)
    for j, mi in enumerate(M):                                    # capacity: Σ nH·x − s+ + s− = nH_obs
        m = mults[mi]
        row = np.zeros(nv)
        for p, (gid, mm, _) in enumerate(pairs):
            if mm == mi:
                row[p] = G[gidx[gid]]["nH"]
        obs = float(m.nH or 0.0)
        has = row[:nP].any()
        if not has:
            continue
        row[nP + nG + j] = -1; row[nP + nG + nM + j] = 1
        A.append(row); lb.append(obs); ub.append(obs)
    integ = np.zeros(nv); integ[:nP + nG] = 1
    res = milp(c, constraints=LinearConstraint(np.array(A), lb, ub), integrality=integ,
               bounds=Bounds(np.zeros(nv), np.r_[np.ones(nP + nG), np.full(2 * nM, np.inf)]),
               options={"time_limit": 10})
    out = {g["gid"]: None for g in G}
    if res.x is None:
        return out
    for p, (gid, mi, _) in enumerate(pairs):
        if res.x[p] > 0.5:
            out[gid] = mi
    for mi in M:
        tags = sorted((g["tag"] for g in G if out[g["gid"]] == mi), key=lambda t: (len(t), t))
        mults[mi].label = ", ".join("H" + t for t in tags)
        mults[mi].nH_assigned = float(sum(g["nH"] for g in G if out[g["gid"]] == mi)) or None
    return out


def unassigned(mults, assignment, min_nH=0.0):
    used = {mi for mi in assignment.values() if mi is not None}
    return [i for i in _usable(mults) if i not in used and not mults[i].note and (mults[i].nH or 0) >= min_nH]


def _shift(v):
    return float(np.mean(v)) if isinstance(v, list) else float(v)


def impurity_check(mults, idx, solvent, tol=0.04, unit_nH=1.0, mw_compound=None):
    """For multiplets mults[idx] (unassigned): candidate impurities. Returns
    {"per_peak": {i: [ {name, proton, delta, table_shift, score}, …]}, "best": [ {name, peaks, mol_ratio, wt_pct, score} ]}.
    mol_ratio = mol impurity per mol compound (integrals normalised to the compound's 1 H = unit_nH)."""
    from mfs import nmrproc as NP
    key = NP.solvent_key(solvent) or solvent
    T = impurity_table()["signals"]
    comps = {}
    for s in T:
        v = s["shifts"].get(key)
        if v is None or (s["name"] == "water" and False):
            continue
        comps.setdefault(s["name"], []).append((s, v))
    per = {i: [] for i in idx}
    best = []
    for name, sigs in comps.items():
        matched, tot_w, err = [], 0.0, []
        for s, v in sigs:
            t = tol * (3.5 if s["proton"] in ("OH", "NH") or name == "water" else 1.0)
            lo, hi = (v if isinstance(v, list) else (v, v))
            cand = []
            for i in idx:
                m = mults[i]
                d = 0.0 if lo - 0.005 <= m.centre <= hi + 0.005 else min(abs(m.centre - lo), abs(m.centre - hi))
                if d <= t:
                    cand.append((d, i))
            w = 0.5 if s["proton"] in ("OH", "NH") else 1.0
            tot_w += w
            if cand:
                d, i = min(cand)
                sc = 1.0 - d / (t + 1e-9) * 0.5
                k = _MULT_NEIGH.get(mults[i].mult); kt = _MULT_NEIGH.get(s["mult"].split(",")[0].strip())
                if k is not None and kt is not None and k != kt and mults[i].mult not in ("m",):
                    sc -= 0.3
                per[i].append(dict(name=name, proton=s["proton"], delta=round(mults[i].centre - _shift(v), 3),
                                   table_shift=v, mult=s["mult"], score=round(sc, 2), src=s.get("src", "")))
                matched.append((s, i, w, sc))
                err.append(d)
        if not matched:
            continue
        frac = sum(w for _, _, w, _ in matched) / tot_w
        # integral consistency among matched signals of the same impurity
        r = [(mults[i].nH or 0) / s["nH"] for s, i, _, _ in matched if (mults[i].nH or 0) > 0]
        cons = 1.0
        if len(r) >= 2 and min(r) > 0:
            cons = float(min(r) / max(r))
        score = frac * np.mean([sc for *_, sc in matched]) * (0.5 + 0.5 * cons)
        mol = float(np.median(r) / unit_nH) if r else None
        wt = None
        if mol is not None and mw_compound and name in IMP_SMILES:
            from rdkit import Chem
            from rdkit.Chem.Descriptors import MolWt
            mw = MolWt(Chem.MolFromSmiles(IMP_SMILES[name]))
            wt = round(100 * mol * mw / (mw_compound + mol * mw), 2)
        best.append(dict(name=name, peaks=sorted({i for _, i, _, _ in matched}), matched=f"{len(matched)}/{len(sigs)}",
                         score=round(float(score), 3), mol_ratio=None if mol is None else round(mol, 4), wt_pct=wt))
    for i in per:
        per[i].sort(key=lambda d: -d["score"])
    best.sort(key=lambda d: -d["score"])
    # greedy: each peak explained by its best compound
    taken, guess = set(), []
    for b in best:
        if b["score"] < 0.35:
            continue
        free = [p for p in b["peaks"] if p not in taken]
        if not free:
            continue
        taken |= set(free)
        guess.append(dict(b, peaks=free))
    return {"per_peak": per, "candidates": best[:15], "best": guess,
            "unexplained": [i for i in idx if i not in taken]}
