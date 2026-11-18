"""v2.16: plausibility checks BETWEEN techniques. The numbers that different pages put into the sample store are compared with each other; a result
that is physically inconsistent (or suspiciously far apart) gets a plain-language finding with what it can mean and what to check.

Each finding: level ('ok' | 'info' | 'warn'), sample, rule, title, text, advice, values."""
from __future__ import annotations

import math
import re

from mfs.summary import TECH_NAMES


def _num(rec: dict, section: str, *keys):
    """First numeric value among `keys` (exact names, then names that contain the key) of one section of one sample, else None."""
    vals = ((rec.get("sections") or {}).get(section) or {}).get("values") or {}
    for k in keys:
        v = vals.get(k)
        if isinstance(v, (int, float)) and math.isfinite(v):
            return float(v)
    for k in keys:
        for name, v in vals.items():
            if k.lower() in str(name).lower() and isinstance(v, (int, float)) and math.isfinite(v):
                return float(v)
    return None


def _f(level, name, rule, title, text, advice="", **values):
    return {"level": level, "sample": name, "rule": rule, "title": title, "text": text, "advice": advice, "values": values}


def check_sample(rec: dict) -> list[dict]:
    name = rec.get("name", "?")
    out: list[dict] = []
    Dc = _num(rec, "pxrd", "crystallite_size_WH_nm", "crystallite_size_scherrer_nm")
    dp = _num(rec, "sem", "mean (nm)")
    rho = _num(rec, "pxrd", "density_g_cm3", "density")
    S = _num(rec, "n2", "S_BET")
    z = _num(rec, "dls", "z_avg_nm")
    # 1 crystallite size (PXRD) vs particle size (SEM / TEM)
    if Dc and dp and Dc > 0:
        q = dp / Dc
        v = dict(crystallite_nm=round(Dc, 1), particle_nm=round(dp, 1), ratio=round(q, 2))
        if q < 0.7:
            out.append(_f("warn", name, "size", "Particles smaller than the crystallites",
                          f"The mean particle size from the images ({dp:.0f} nm) is smaller than the crystallite size from PXRD ({Dc:.0f} nm). A particle cannot be smaller "
                          "than the crystal it consists of.", "Check the scale-bar calibration (nm/px) of the images, whether the segmentation splits particles, and the "
                          "instrument broadening used for the size (a too small instrumental width gives a too large size).", **v))
        elif q <= 1.5:
            out.append(_f("ok", name, "size", "Crystallite and particle size agree",
                          f"Particle {dp:.0f} nm vs crystallite {Dc:.0f} nm (ratio {q:.1f}): consistent with single-crystalline particles.", **v))
        elif q <= 4:
            out.append(_f("info", name, "size", "Particles somewhat larger than the crystallites",
                          f"Particle {dp:.0f} nm is {q:.1f} × the crystallite size {Dc:.0f} nm: mosaic / defective crystals, or peak broadening by strain and the instrument "
                          "that makes the crystallites look smaller.", "Compare with the Williamson–Hall size and the strain.", **v))
        else:
            out.append(_f("info", name, "size", "Particles are probably polycrystalline",
                          f"Particle {dp:.0f} nm is {q:.1f} × the crystallite size {Dc:.0f} nm: about {q ** 3:.0f} crystallites per particle if they are compact aggregates.",
                          "Look for aggregates in the images; report both sizes and say which is which.", **v))
    # 2 BET area vs the geometric (external) area of the particles
    if S and dp and rho and rho > 0:
        Sg = 6000.0 / (rho * dp)
        q = S / Sg
        v = dict(S_BET=round(S), S_geometric=round(Sg, 1), ratio=round(q, 1), density=round(rho, 2))
        if q < 0.5:
            out.append(_f("warn", name, "bet", "BET area below the external area of the particles",
                          f"S_BET = {S:.0f} m²/g is lower than the geometric surface of {dp:.0f} nm spheres of density {rho:.2f} g/cm³ ({Sg:.0f} m²/g).",
                          "Pores blocked or not activated, strongly aggregated particles, or a wrong particle size / density. Check the activation and the degassing.", **v))
        elif q < 20:
            out.append(_f("ok", name, "bet", "BET area plausible for the particle size",
                          f"S_BET = {S:.0f} m²/g vs {Sg:.0f} m²/g external: ratio {q:.1f}.", **v))
        else:
            out.append(_f("info", name, "bet", "BET area dominated by the internal pores",
                          f"S_BET = {S:.0f} m²/g is {q:.0f} × the external area ({Sg:.0f} m²/g): the surface is inside the framework (expected for a microporous solid).", **v))
    # 3 defect methods disagree
    rng_txt = ((rec.get("sections") or {}).get("defects") or {}).get("values", {}).get("L_range")
    if isinstance(rng_txt, str):
        m = re.match(r"\s*([\d.]+)\s*[–-]\s*([\d.]+)", rng_txt)
        if m:
            lo, hi = float(m.group(1)), float(m.group(2))
            if hi - lo > 1.0:
                out.append(_f("warn", name, "defects", "The defect methods disagree",
                              f"Linkers per node range from {lo:.2f} to {hi:.2f} between the methods (TGA / NMR / elemental). The spread is more than one linker per node.",
                              "Check the sample preparation (residual solvent, modulator, extra-framework oxide) and which methods were used; do not average methods that disagree.",
                              L_min=lo, L_max=hi))
            else:
                out.append(_f("ok", name, "defects", "The defect methods agree",
                              f"Linkers per node {lo:.2f}–{hi:.2f} across the methods used.", L_min=lo, L_max=hi))
    # 4 surface (XPS) vs bulk (ICP / XRF) metal ratio
    r_bulk = _num(rec, "icp", "Hf/(Zr+Hf)", "Hf_fraction") or _num(rec, "xrf", "Hf/(Zr+Hf)")
    r_surf = _num(rec, "xps", "Hf/(Zr+Hf)")
    if r_bulk and r_surf and r_bulk > 0:
        d = (r_surf - r_bulk) / r_bulk
        v = dict(bulk=round(r_bulk, 3), surface=round(r_surf, 3), relative_difference=round(d, 2))
        if abs(d) > 0.25:
            out.append(_f("info", name, "surface", "Surface and bulk metal ratio differ",
                          f"XPS (surface) {r_surf:.2f} vs ICP / XRF (bulk) {r_bulk:.2f} ({100 * d:+.0f} %).",
                          "Surface enrichment (coating, segregation) or a bulk inhomogeneity; XPS probes only the top few nm.", **v))
        else:
            out.append(_f("ok", name, "surface", "Surface and bulk metal ratio agree", f"XPS {r_surf:.2f} vs bulk {r_bulk:.2f}.", **v))
    # 5 hydrodynamic size (DLS) vs particle size (SEM)
    if z and dp:
        q = z / dp
        v = dict(z_avg_nm=round(z), particle_nm=round(dp), ratio=round(q, 2))
        pdi = _num(rec, "dls", "pdi")
        if q < 0.7:
            out.append(_f("warn", name, "dls", "DLS size smaller than the image size",
                          f"Z-average {z:.0f} nm vs particle size {dp:.0f} nm from the images. The hydrodynamic diameter is normally not smaller than the dry size.",
                          "Check the calibration of the images and whether the sample is the same; large particles settle and drop out of the DLS measurement.", **v))
        elif q <= 3:
            out.append(_f("ok", name, "dls", "DLS and image size are consistent", f"Z-average {z:.0f} nm vs {dp:.0f} nm from the images (ratio {q:.1f}).", **v))
        else:
            out.append(_f("info", name, "dls", "DLS size much larger than the image size",
                          f"Z-average {z:.0f} nm is {q:.1f} × the particle size from the images" + (f" (PDI {pdi:.2f})" if pdi else "") + ": aggregation in the dispersion.",
                          "Sonicate / change the medium or pH, and report the distribution, not only the Z-average.", **v))
    return out


def check(store: dict) -> list[dict]:
    """All findings for all samples of the sample store, warnings first."""
    res = []
    for sid, rec in (store or {}).items():
        res += check_sample(rec)
    order = {"warn": 0, "info": 1, "ok": 2}
    return sorted(res, key=lambda f: (order[f["level"]], f["sample"], f["rule"]))


def as_markdown(findings: list[dict]) -> str:
    if not findings:
        return "No cross-technique checks possible yet (need results of at least two techniques for the same sample).\n"
    icon = {"warn": "⚠️", "info": "ℹ️", "ok": "✅"}
    return "\n".join(f"- {icon[f['level']]} **{f['sample']} – {f['title']}.** {f['text']}" + (f" _{f['advice']}_" if f["advice"] else "") for f in findings) + "\n"
