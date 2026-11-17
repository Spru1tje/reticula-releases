"""v2.14: ready-to-edit Methods / SI sentences and data-quality notes generated from the data on a page.

The text states only what the data show (counts, ranges, mean ± SD); instrument, sample preparation and temperature are placeholders in
[square brackets] for the author to fill in - nothing about the experiment is invented."""
from __future__ import annotations

import numpy as np


def _rng(v, nd):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if v.size == 0:
        return None
    lo, hi = v.min(), v.max()
    return f"{lo:.{nd}f}" if round(lo, nd) == round(hi, nd) else f"{lo:.{nd}f}–{hi:.{nd}f}"


def dls_quality_notes(rows):
    """rows: iterable of dicts with sample, n, z_avg, pdi, zeta, zeta_sd. Returns a list of (level, text), level 'warn' or 'info'."""
    notes, low = [], []
    for r in rows:
        s = r["sample"]
        n = int(r.get("n") or 0)
        pdi, zeta = r.get("pdi"), r.get("zeta")
        if n == 1:
            notes.append(("info", f"{s}: one record only – no standard deviation; repeat the measurement (≥ 3 records) for a paper."))
        if pdi is not None and np.isfinite(pdi):
            if pdi > 0.7:
                notes.append(("warn", f"{s}: PDI {pdi:.2f} > 0.7 – very broad / aggregated; the Z-average is not a reliable size (report the distribution)."))
            elif pdi > 0.4:
                notes.append(("info", f"{s}: PDI {pdi:.2f} – broad distribution; interpret the Z-average with care."))
        if zeta is not None and np.isfinite(zeta) and abs(zeta) < 30:
            low.append(f"{s} ({zeta:.0f} mV)")
    if low:
        notes.append(("info", "|ζ| < 30 mV (colloidally not very stable by the usual rule of thumb): " + ", ".join(low) + "."))
    return notes


def dls_methods(rows, with_sd=True):
    """A Methods sentence + a results sentence for the DLS / zeta data of a page. rows as in dls_quality_notes."""
    rows = list(rows)
    if not rows:
        return ""
    n_all = [int(r.get("n") or 0) for r in rows]
    za = _rng([r["z_avg"] for r in rows], 0)
    pd_ = _rng([r["pdi"] for r in rows], 2)
    zp = _rng([r["zeta"] for r in rows], 1)
    meth = ("Hydrodynamic diameters (Z-average, intensity-weighted) and polydispersity indices (PDI) were measured by dynamic light scattering, "
            "and ζ-potentials by electrophoretic light scattering, on a Malvern Zetasizer [model] at [temperature] °C in [medium, pH] "
            f"(sample concentration [c]; {min(n_all)}–{max(n_all)} records per sample" if min(n_all) != max(n_all) else
            "Hydrodynamic diameters (Z-average, intensity-weighted) and polydispersity indices (PDI) were measured by dynamic light scattering, "
            "and ζ-potentials by electrophoretic light scattering, on a Malvern Zetasizer [model] at [temperature] °C in [medium, pH] "
            f"(sample concentration [c]; {n_all[0]} records per sample")
    meth += (", values are mean ± standard deviation of the records)." if with_sd else ", values are means of the records).")
    res = []
    if za:
        res.append(f"Z-average diameters were {za} nm")
    if pd_:
        res.append(f"PDI {pd_}")
    if zp:
        res.append(f"ζ-potentials {zp} mV")
    out = meth
    if res:
        out += " " + ", ".join(res) + f" ({len(rows)} samples)."
    return out
