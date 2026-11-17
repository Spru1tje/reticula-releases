"""v2.15: library of linkers / modulators / common species for the ¹H NMR of digested MOFs, and matching of found multiplets against it.

The shifts are TYPICAL literature values (digestion in NaOD/D₂O, DMSO-d₆ with acid, CDCl₃ …) and depend on pH, concentration and the
counter-ion: they are meant to suggest, not to prove. Every row can be corrected / extended by the user (own rows win and are kept in the project).
Format: name → role, formula and list of signals (δ ppm, multiplicity, nH) per solvent."""
from __future__ import annotations

# solvent keys: 'D2O' (also NaOD/D2O digests), 'DMSO' (DMSO-d6, with or without HF/acid), 'CDCl3'
LIB = {
    "terephthalic acid (BDC)": ("linker", "C8H6O4", {"D2O": [(7.75, "s", 4)], "DMSO": [(8.04, "s", 4)]}),
    "trimesic acid (BTC)": ("linker", "C9H6O6", {"D2O": [(8.35, "s", 3)], "DMSO": [(8.64, "s", 3)]}),
    "2-aminoterephthalic acid (BDC-NH2)": ("linker", "C8H7NO4", {"D2O": [(7.55, "d", 1), (7.15, "s", 1), (6.95, "d", 1)],
                                                              "DMSO": [(7.80, "d", 1), (7.30, "s", 1), (7.07, "d", 1)]}),
    "fumaric acid": ("linker", "C4H4O4", {"D2O": [(6.45, "s", 2)], "DMSO": [(6.63, "s", 2)]}),
    "succinic acid": ("linker", "C4H6O4", {"D2O": [(2.40, "s", 4)], "DMSO": [(2.42, "s", 4)]}),
    "formic acid / formate": ("modulator", "CH2O2", {"D2O": [(8.45, "s", 1)], "DMSO": [(8.25, "s", 1)], "CDCl3": [(8.05, "s", 1)]}),
    "acetic acid / acetate": ("modulator", "C2H4O2", {"D2O": [(1.90, "s", 3)], "DMSO": [(1.91, "s", 3)], "CDCl3": [(2.10, "s", 3)]}),
    "benzoic acid / benzoate": ("modulator", "C7H6O2", {"D2O": [(7.80, "d", 2), (7.50, "m", 3)],
                                                      "DMSO": [(7.95, "d", 2), (7.62, "t", 1), (7.50, "t", 2)]}),
    "DMF": ("solvent / residue", "C3H7NO", {"D2O": [(7.92, "s", 1), (2.96, "s", 3), (2.80, "s", 3)],
                                         "DMSO": [(7.95, "s", 1), (2.89, "s", 3), (2.73, "s", 3)], "CDCl3": [(8.02, "s", 1), (2.95, "s", 3), (2.88, "s", 3)]}),
    "dimethylamine (DMF decomposition)": ("residue", "C2H7N", {"D2O": [(2.72, "s", 6)]}),
    "methanol": ("solvent / residue", "CH4O", {"D2O": [(3.34, "s", 3)], "DMSO": [(3.17, "s", 3)], "CDCl3": [(3.49, "s", 3)]}),
}


def _rows(user_rows=None):
    rows = []                                   # (name, role, solvent, delta, mult, nH, source)
    for r in user_rows or []:
        try:
            rows.append((str(r["compound"]).strip(), str(r.get("role", "own")), str(r["solvent"]).strip(), float(r["δ (ppm)"]),
                         str(r.get("mult", "s")), float(r["nH"]), "own"))
        except (KeyError, TypeError, ValueError):
            continue
    for name, (role, _f, per) in LIB.items():
        for sol, sig in per.items():
            for d, m, n in sig:
                rows.append((name, role, sol, d, m, float(n), "typical"))
    return rows


def solvent_key(solvent: str | None) -> str:
    s = (solvent or "").lower()
    if "d2o" in s or "naod" in s or "water" in s:
        return "D2O"
    if "dmso" in s:
        return "DMSO"
    if "cdcl" in s or "chloroform" in s:
        return "CDCl3"
    return "D2O"


def match(mults, solvent=None, user_rows=None, tol=0.12):
    """mults: iterable of (centre ppm, integral). Returns one dict per library compound with ≥ 1 signal found:
    compound, role, found 'x/y', mean |Δδ|, amount (mean integral per H over the found signals), complete (all signals found), signals [(δlib, δfound, nH, integral)]."""
    sk = solvent_key(solvent)
    mults = [(float(c), float(i)) for c, i in mults]
    by = {}
    for name, role, sol, d, m, n, src in _rows(user_rows):
        if sol != sk:
            continue
        by.setdefault((name, role, src), []).append((d, m, n))
    out = []
    used = set()
    for (name, role, src), sig in by.items():
        found = []
        for d, m, n in sig:
            cand = [(abs(c - d), k, c, i) for k, (c, i) in enumerate(mults) if abs(c - d) <= tol]
            if cand:
                dd, k, c, i = min(cand)
                found.append((d, c, n, i, k, dd))
        if not found:
            continue
        amt = sum(f[3] / f[2] for f in found) / len(found)
        out.append({"compound": name, "role": role, "source": src, "found": f"{len(found)}/{len(sig)}", "complete": len(found) == len(sig),
                    "mean_dd": sum(f[5] for f in found) / len(found), "amount": amt, "signals": [(f[0], f[1], f[2], f[3]) for f in found]})
        used.update(f[4] for f in found)
    out.sort(key=lambda r: (not r["complete"], r["mean_dd"]))
    return out


def ratios(matches, reference):
    """Molar ratio of every found compound relative to `reference` (compound name), from the integral per H (only complete matches)."""
    ref = next((m for m in matches if m["compound"] == reference), None)
    if not ref or ref["amount"] <= 0:
        return {}
    return {m["compound"]: m["amount"] / ref["amount"] for m in matches if m["complete"] or m["compound"] == reference}
