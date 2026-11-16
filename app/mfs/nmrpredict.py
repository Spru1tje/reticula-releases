"""1H chemical shift prediction from HOSE-style circular atom environments (v2.13).

Table built from experimental 1H assignments in nmrshiftdb2 (CC-BY-SA; https://nmrshiftdb.nmr.uni-koeln.de) by
`build_table` (calculated/predicted spectra excluded). Prediction: for each H-bearing heavy atom, the longest environment
(sphere 4 → 1) with enough examples gives the median shift; the spread (MAD) and sphere reached give the confidence.
Typical accuracy on held-out nmrshiftdb2 molecules: see NOTES / tests (MAE ≈ 0.2 ppm).
"""
from __future__ import annotations

import gzip
import json
import os
import re
from collections import defaultdict

import numpy as np

MAX_SPHERE = 4
_TABLE = None
TABLE_PATH = os.path.join(os.path.dirname(__file__), "data", "h1_hose.json.gz")
EXCHANGEABLE = {"O", "N", "S"}


def _rd():
    from rdkit import Chem
    return Chem


def _bsym(b):
    if b.GetIsAromatic():
        return ":"
    return {1.0: "-", 2.0: "=", 3.0: "#"}.get(b.GetBondTypeAsDouble(), "~")


def _atok(a):
    s = a.GetSymbol() + ("a" if a.GetIsAromatic() else "") + str(a.GetTotalNumHs())
    if a.GetFormalCharge():
        s += "+" if a.GetFormalCharge() > 0 else "-"
    return s


def env_codes(mol, idx, max_sphere=MAX_SPHERE):
    """Codes for spheres 0..max_sphere around heavy atom idx (each code includes all inner spheres)."""
    a0 = mol.GetAtomWithIdx(idx)
    ring = mol.GetRingInfo()
    rs = min([len(r) for r in ring.AtomRings() if idx in r] or [0])
    codes = [f"{_atok(a0)}r{rs}"]
    seen = {idx}
    frontier = [(idx, "")]
    cur = codes[0]
    for _sph in range(1, max_sphere + 1):
        toks, nxt = [], []
        for ai, ptok in frontier:
            a = mol.GetAtomWithIdx(ai)
            for b in a.GetBonds():
                j = b.GetOtherAtomIdx(ai)
                if j in seen:
                    continue
                nb = mol.GetAtomWithIdx(j)
                t = _bsym(b) + _atok(nb) + ("r" if nb.IsInRing() else "")
                toks.append(ptok + t)
                nxt.append((j, t[:3]))
        for j, _ in nxt:
            seen.add(j)
        cur = cur + "/" + ",".join(sorted(toks))
        codes.append(cur)
        frontier = nxt
        if not nxt:
            # pad: no more atoms → same code for all deeper spheres
            codes += [cur] * (max_sphere - _sph)
            break
    return codes


def _parse_spectra(props, nucleus="1H"):
    """Experimental spectra of one SD record → list of {atom: [shifts]}."""
    calc = set()
    for k in ("NMRMethod", "NMRProgram", "NMRBasisSet", "GeomMethod", "Program"):
        v = props.get(k, "")
        calc |= {int(m) for m in re.findall(r"(\d+):", v)}
    out = []
    for k, v in props.items():
        m = re.match(rf"Spectrum {nucleus} (\d+)$", k)
        if not m or int(m.group(1)) in calc:
            continue
        d = defaultdict(list)
        for part in v.strip().split("|"):
            f = part.split(";")
            if len(f) >= 3 and f[0] and f[2].strip().isdigit():
                try:
                    d[int(f[2])].append(float(f[0]))
                except ValueError:
                    pass
        if d:
            sol = props.get("Solvent", "")
            out.append((dict(d), sol))
    return out


def iter_training(sd_path, holdout=None):
    """Yield (mol, {atom: shift}, inchikey, solvent) for experimental 1H spectra."""
    Chem = _rd()
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    sup = Chem.ForwardSDMolSupplier(sd_path, sanitize=True, removeHs=True)
    for mol in sup:
        if mol is None:
            continue
        props = mol.GetPropsAsDict()
        props = {k: str(v) for k, v in props.items()}
        specs = _parse_spectra(props)
        if not specs:
            continue
        key = props.get("INChI key", "")
        for d, sol in specs:
            yield mol, {a: float(np.mean(v)) for a, v in d.items() if a < mol.GetNumAtoms()}, key, sol


def build_table(sd_path, out_path=TABLE_PATH, test_frac=0.1, seed=7, min_n=1):
    """Build the environment → shift table; molecules whose InChIKey hash falls in the test fraction are held out and
    written alongside (for validation). Returns (n_train_atoms, n_test_atoms)."""
    import hashlib
    acc = defaultdict(list)
    test = []
    ntr = 0
    for mol, shifts, key, sol in iter_training(sd_path):
        h = int(hashlib.md5((key or str(id(mol))).encode()).hexdigest()[:8], 16) / 0xFFFFFFFF
        is_test = h < test_frac
        items = []
        for a, s in shifts.items():
            at = mol.GetAtomWithIdx(a)
            if at.GetTotalNumHs() == 0 or not (-2 < s < 18):
                continue
            items.append((a, s))
        if is_test:
            from rdkit import Chem
            test.append({"smiles": Chem.MolToSmiles(mol, canonical=False), "molblock": Chem.MolToMolBlock(mol),
                         "shifts": items, "solvent": sol})
            continue
        for a, s in items:
            for c in env_codes(mol, a):
                acc[c].append(s)
            ntr += 1
    table = {}
    for c, v in acc.items():
        if len(v) >= min_n:
            v = np.asarray(v)
            med = float(np.median(v))
            table[c] = [round(med, 2), len(v), round(float(np.median(np.abs(v - med))), 2)]
    with gzip.open(out_path, "wt") as f:
        json.dump({"source": "nmrshiftdb2 (CC-BY-SA), experimental 1H only", "max_sphere": MAX_SPHERE, "table": table}, f,
                  separators=(",", ":"))
    with gzip.open(out_path.replace(".json.gz", "_test.json.gz"), "wt") as f:
        json.dump(test, f)
    return ntr, sum(len(t["shifts"]) for t in test)


def table():
    global _TABLE
    if _TABLE is None:
        with gzip.open(TABLE_PATH, "rt") as f:
            _TABLE = json.load(f)["table"]
    return _TABLE


def predict_atom(mol, idx, min_n=1):
    """(shift, sphere_used, n, spread) for heavy atom idx (None if no match)."""
    T = table()
    codes = env_codes(mol, idx)
    for sph in range(len(codes) - 1, -1, -1):
        r = T.get(codes[sph])
        if r and (r[1] >= min_n or (sph >= 3 and r[1] >= 1)):
            return r[0], sph, r[1], r[2]
    return None


def predict(mol):
    """1H prediction for an RDKit mol (implicit H). Returns list of dicts, one per H-bearing heavy atom:
    {atom, idx(1-based label), nH, shift, sphere, n, spread, exchangeable, conf}."""
    out = []
    for a in mol.GetAtoms():
        nh = a.GetTotalNumHs()
        if nh == 0:
            continue
        r = predict_atom(mol, a.GetIdx())
        exch = a.GetSymbol() in EXCHANGEABLE
        if r is None:
            out.append(dict(atom=a.GetIdx(), label=a.GetIdx() + 1, element=a.GetSymbol(), nH=nh, shift=None, sphere=0,
                            n=0, spread=None, exchangeable=exch, conf="none"))
            continue
        s, sph, n, spread = r
        conf = "high" if sph >= 3 and n >= 5 and spread <= 0.15 else ("medium" if sph >= 2 else "low")
        if exch:
            conf = "low"
        out.append(dict(atom=a.GetIdx(), label=a.GetIdx() + 1, element=a.GetSymbol(), nH=nh, shift=s, sphere=sph, n=n,
                        spread=spread, exchangeable=exch, conf=conf))
    return out


def equivalence_classes(mol):
    """Topologically equivalent atoms (canonical ranks, symmetry kept) → {atom: class id}."""
    Chem = _rd()
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    return {i: r for i, r in enumerate(ranks)}


def mol_from_any(text):
    """SMILES or molfile text → RDKit mol (None on failure)."""
    Chem = _rd()
    if not text or not text.strip():
        return None
    t = text.strip()
    mol = Chem.MolFromMolBlock(text, sanitize=True) if ("M  END" in t or "V2000" in t or "V3000" in t) else Chem.MolFromSmiles(t)
    if mol is None:
        return None
    return Chem.RemoveHs(mol)


def validate(test_path=None, limit=None):
    """MAE / RMSE / coverage on the held-out molecules (non-exchangeable H only)."""
    Chem = _rd()
    test_path = test_path or os.path.join(os.path.dirname(__file__), "..", "tests", "data", "nmr", "h1_hose_test.json.gz")
    with gzip.open(test_path, "rt") as f:
        test = json.load(f)
    err, miss, by_sph = [], 0, defaultdict(list)
    for t in test[:limit]:
        mol = Chem.MolFromMolBlock(t["molblock"])
        if mol is None:
            continue
        for a, s in t["shifts"]:
            if mol.GetAtomWithIdx(a).GetSymbol() in EXCHANGEABLE:
                continue
            r = predict_atom(mol, a)
            if r is None:
                miss += 1
                continue
            err.append(r[0] - s); by_sph[r[1]].append(abs(r[0] - s))
    e = np.abs(err)
    return {"n": len(err), "missing": miss, "MAE": round(float(e.mean()), 3), "median_AE": round(float(np.median(e)), 3),
            "RMSE": round(float(np.sqrt(np.mean(np.square(err)))), 3), "within_0.3": round(float(np.mean(e < 0.3)), 3),
            "by_sphere": {k: [len(v), round(float(np.mean(v)), 3)] for k, v in sorted(by_sph.items())}}
