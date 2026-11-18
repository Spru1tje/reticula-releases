"""v2.16: the refined structure as a CIF file ready for publication / deposition: cell and atom parameters with their estimated standard deviations in the
usual 1.2345(6) notation, the symmetry operations, the displacement parameters, bonds and angles, and the powder-refinement figures of merit."""
from __future__ import annotations

import datetime
import math
from fractions import Fraction

import numpy as np


def fmt_esd(value: float, esd: float | None, min_decimals: int = 0, max_decimals: int = 6) -> str:
    """1.2345(6): the esd is written as the last digits (two significant digits when its first digit is 1)."""
    if esd is None or not np.isfinite(esd) or esd <= 0:
        return f"{value:.{max(min_decimals, 5) if abs(value) < 100 else 3}f}".rstrip("0").rstrip(".") if value != int(value) else f"{value:.0f}"
    e = math.floor(math.log10(esd))
    dec = int(min(max(1 - e if str(f"{esd / 10 ** e:.1f}")[0] == "1" else -e, 0), max_decimals))
    digs = int(round(esd * 10 ** dec))
    if digs >= 20 and dec > 0:                                  # keep 1-19 as the esd digits
        dec -= 1
        digs = int(round(esd * 10 ** dec))
    return f"{value:.{dec}f}({max(digs, 1)})"


def _op_string(R, t) -> str:
    axes = "xyz"
    parts = []
    for i in range(3):
        s = ""
        for j in range(3):
            c = int(round(R[i][j]))
            if c:
                s += ("+" if c > 0 and s else "-" if c < 0 else "") + (axes[j] if abs(c) == 1 else f"{abs(c)}{axes[j]}")
        fr = Fraction(float(t[i])).limit_denominator(12) % 1
        if fr:
            s += ("+" if s else "") + f"{fr.numerator}/{fr.denominator}"
        parts.append(s or "0")
    return ",".join(parts)


def write_cif(phase_res: dict, model, name: str = "refined", radiation_wavelength: float | None = None, stats: dict | None = None, version: str = "") -> str:
    """CIF text of one refined phase. `phase_res` = the phase dict of a Rietveld result with a 'structure' entry; `model` = the phase's StructureModel."""
    st = phase_res["structure"]
    cell = phase_res["cell"]
    ce = phase_res.get("cell_esd") or {}
    names = ["a", "b", "c", "alpha", "beta", "gamma"]
    esd_of = {"a": ce.get("a"), "b": ce.get("b") if "b" in ce else ce.get("a") if phase_res["system"] in ("cubic", "tetragonal", "hexagonal", "rhombohedral") else None,
              "c": ce.get("c") if "c" in ce else (ce.get("a") if phase_res["system"] in ("cubic", "rhombohedral") else None), "alpha": ce.get("alpha"),
              "beta": ce.get("beta"), "gamma": ce.get("gamma")}
    safe = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in name)
    L = [f"data_{safe}", "#" + "-" * 70, f"# Refined with Reticula {version}".rstrip() + f" on {datetime.date.today().isoformat()} (Rietveld refinement of powder data)", "#" + "-" * 70,
         "_audit_creation_method 'Reticula Rietveld refinement'", f"_chemical_name_common '{name}'", f"_symmetry_space_group_name_H-M '{st['spacegroup']}'",
         f"_symmetry_Int_Tables_number {st['sg_number']}"]
    for k, n in enumerate(names):
        key = {"a": "length_a", "b": "length_b", "c": "length_c", "alpha": "angle_alpha", "beta": "angle_beta", "gamma": "angle_gamma"}[n]
        L.append(f"_cell_{key} {fmt_esd(float(cell[k]), esd_of.get(n))}")
    L.append(f"_cell_volume {fmt_esd(phase_res['volume'], phase_res.get('volume_esd'))}")
    if radiation_wavelength:
        L.append(f"_diffrn_radiation_wavelength {radiation_wavelength:.5f}")
    if stats:
        L += [f"_pd_proc_ls_prof_wR_factor {100 * stats['Rwp']:.2f}", f"_pd_proc_ls_prof_R_factor {100 * stats['Rp']:.2f}",
              f"_refine_ls_goodness_of_fit_all {math.sqrt(stats['chi2']):.3f}"]
    if phase_res.get("Rb") is not None:
        L.append(f"_refine_ls_R_Fsqd_factor {100 * phase_res['Rb']:.2f}")
    L += ["", "loop_", "_symmetry_equiv_pos_as_xyz"]
    seen = set()
    for R, t in zip(model.rotations, model.translations):
        op = _op_string(R, t)
        if op not in seen:
            seen.add(op)
            L.append(f"'{op}'")
    L += ["", "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_U_iso_or_equiv",
          "_atom_site_occupancy", "_atom_site_symmetry_multiplicity"]
    for a in st["atoms"]:
        el = a["element"].split("/")[0]
        u = a["B"] / (8 * math.pi ** 2)
        ue = (a["B_esd"] / (8 * math.pi ** 2)) if a.get("B_esd") else None
        L.append(f"{a['label']} {el} {fmt_esd(a['x'], a['x_esd'])} {fmt_esd(a['y'], a['y_esd'])} {fmt_esd(a['z'], a['z_esd'])} {fmt_esd(u, ue)} {a['occupancy']:.4f} {a['n_equivalent']}")
    if st.get("bonds"):
        L += ["", "loop_", "_geom_bond_atom_site_label_1", "_geom_bond_atom_site_label_2", "_geom_bond_distance"]
        for b in sorted(st["bonds"], key=lambda b: (b["atom1"], b["distance"])):
            L.append(f"{b['atom1']} {b['atom2']} {fmt_esd(b['distance'], b['esd'])}")
    if st.get("angles"):
        L += ["", "loop_", "_geom_angle_atom_site_label_1", "_geom_angle_atom_site_label_2", "_geom_angle_atom_site_label_3", "_geom_angle"]
        for g in sorted(st["angles"], key=lambda g: (g["centre"], g["atom1"], g["atom2"])):
            L.append(f"{g['atom1']} {g['centre']} {g['atom2']} {fmt_esd(g['angle'], g['esd'])}")
    return "\n".join(L) + "\n"
