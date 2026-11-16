"""Element molar ratios from ICP-OES / XRF wt-% with propagated (Gaussian, uncorrelated) uncertainties, and their
translation into per-node quantities for the defect formula (Hf fraction, S per M6 node)."""
from __future__ import annotations
import math

M = {"Zr": 91.224, "Hf": 178.49, "S": 32.06, "Ti": 47.867, "Ce": 140.116, "Na": 22.990, "Cl": 35.45, "P": 30.974}


def moles(w, sw, el):
    return w / M[el], sw / M[el]


def ratio(a, sa, b, sb):
    """a/b with σ = r √((σa/a)² + (σb/b)²)."""
    if not b:
        return float("nan"), float("nan")
    r = a / b
    return r, math.hypot(sa / b, a * sb / b ** 2)       # = |r|·√((σa/a)² + (σb/b)²), also valid for a = 0


def zr_hf_s(wt: dict, sd: dict):
    """wt, sd: {'Zr': wt%, 'Hf': wt%, 'S': wt%} (any subset). Returns molar ratios, Hf fraction, S per M6 node ± σ."""
    n = {el: moles(wt[el], sd.get(el, 0.0) or 0.0, el) for el in wt if el in M and wt[el] is not None}
    out = {}
    if "Zr" in n and "Hf" in n:
        nz, sz = n["Zr"]; nh, sh = n["Hf"]
        tot, stot = nz + nh, math.hypot(sz, sh)
        f = nh / tot
        sf = math.sqrt((nz * sh) ** 2 + (nh * sz) ** 2) / tot ** 2       # ∂f/∂nh = nz/tot², ∂f/∂nz = −nh/tot²
        out.update({"Zr:Hf": ratio(nz, sz, nh, sh), "Hf/(Zr+Hf)": (f, sf)})
    else:
        tot, stot = (n.get("Zr") or n.get("Hf") or (None, None))
    if "S" in n and tot:
        ns, ss = n["S"]
        r, sr = ratio(ns, ss, tot, stot)
        out.update({"S/(Zr+Hf)": (r, sr), "S per M6 node": (6 * r, 6 * sr)})
    if all(k in n for k in ("Zr", "Hf", "S")):
        nz, nh, ns = n["Zr"][0], n["Hf"][0], n["S"][0]
        out["Zr:Hf:S (per M6)"] = tuple(6 * v / (nz + nh) for v in (nz, nh, ns))
    return out
