"""v2.18.4: outputs of ANY powder refinement result (Rietveld with a known structure, Le Bail, structure solution): a CIF with the refined cell (and the atoms of the model),
a 3D view of the structure in the refined cell, the pattern as .xy files, the reflection list, and one ZIP with everything. Before this, the CIF, the 3D viewer and the tables
only existed after 'Rietveld + refine the atoms'.

Atoms are those of the model (positions not refined unless the structure refinement was used); the cell, the overall displacement parameter, the scale and the profile are refined.
The CIF says so in a comment, so nobody mistakes the starting positions for refined ones."""
from __future__ import annotations

import datetime
import io
import json
import math
import zipfile

import numpy as np
import pandas as pd

from mfs import structcif as SC


# ------------------------------------------------------------------------------------------------------------------------------ model of the phase in the refined cell
def model_for_phase(ph_obj, phase_res):
    """(StructureModel, atoms dict like Refinement._structure_table without esd) of a CIF-based phase in its REFINED cell; None when the phase has no structure."""
    s0 = getattr(ph_obj, "_structure", None)
    if s0 is None:
        return None
    from pymatgen.core import Lattice, Structure
    from mfs.structmodel import StructureModel
    cell = [float(v) for v in phase_res["cell"]]
    s = Structure(Lattice.from_parameters(*cell), s0.species_and_occu, s0.frac_coords)
    mdl = StructureModel(s)
    B = float(phase_res.get("B") or 0.0)
    labs = mdl.labels()
    atoms = []
    for r in mdl.reps:
        idx = np.nonzero(mdl.atom_site == r)[0]
        els = [mdl.atom_el[j] for j in idx]
        atoms.append({"label": labs[r], "element": "/".join(sorted(set(els))), "x": float(mdl.site_frac[r][0]), "y": float(mdl.site_frac[r][1]), "z": float(mdl.site_frac[r][2]),
                      "x_esd": 0.0, "y_esd": 0.0, "z_esd": 0.0, "B": B, "B_esd": phase_res.get("B_esd"), "occupancy": float(mdl.atom_occ[idx[0]]),
                      "n_equivalent": int(np.sum(mdl.site_rep == r)), "site_symmetry_order": int(mdl.site_sym_order[r]), "free_coordinates": 0})
    sd = {"spacegroup": mdl.spacegroup, "sg_number": mdl.sg_number, "n_atoms_cell": mdl.n_atoms, "atoms": atoms, "bonds": [], "angles": []}
    return mdl, sd


def cif_for_phase(res, i, ph_obj, name, version=""):
    """CIF text of phase i: refined cell (esd where available), the atoms of the model, symmetry; a comment states which parameters were refined. Cell-only phases (Le Bail
    without a structure) get a cell CIF."""
    ph = res["phases"][i]
    wl = _wavelength(res)
    mm = model_for_phase(ph_obj, ph) if ph_obj is not None else None
    if mm is None:
        return cell_cif(ph, name, res, version)
    mdl, sd = mm
    pr = dict(ph)
    pr["structure"] = sd
    txt = SC.write_cif(pr, mdl, name, wl, res["stats"], version)
    mode = "Rietveld refinement of the powder pattern" if res["mode"] == "rietveld" else "Le Bail fit"
    note = ("# Refined: cell, overall displacement parameter (B), scale, background and peak profile.\n"
            "# NOT refined: the atom positions (taken from the starting model). Use 'Rietveld + refine the atoms' for refined positions with esd.\n")
    return txt.replace("_audit_creation_method 'Reticula Rietveld refinement'", f"_audit_creation_method 'Reticula {mode}'") .replace("#" + "-" * 70 + "\n_audit", "#" + "-" * 70 + "\n" + note + "_audit", 1)


def cell_cif(ph, name, res, version=""):
    """A CIF with the cell only (Le Bail of an unknown phase: no structure yet)."""
    ce = ph.get("cell_esd") or {}
    cell = ph["cell"]
    sysm = ph.get("system", "triclinic")
    esd = {"a": ce.get("a"), "b": ce.get("b", ce.get("a") if sysm in ("cubic", "tetragonal", "hexagonal") else None),
           "c": ce.get("c", ce.get("a") if sysm == "cubic" else None)}
    safe = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in name)
    L = [f"data_{safe}", "#" + "-" * 70, f"# Cell from a Le Bail fit with Reticula {version}".rstrip() + f" on {datetime.date.today().isoformat()}",
         "# No atoms yet: the space group is not part of this file. Solve the structure with the PXRD structure-solution assistant.", "#" + "-" * 70,
         "_audit_creation_method 'Reticula Le Bail fit (cell only)'", f"_chemical_name_common '{name}'", "_symmetry_space_group_name_H-M 'P 1'", "_symmetry_Int_Tables_number 1"]
    keys = [("a", "length_a"), ("b", "length_b"), ("c", "length_c"), ("alpha", "angle_alpha"), ("beta", "angle_beta"), ("gamma", "angle_gamma")]
    for k, (n, key) in enumerate(keys):
        L.append(f"_cell_{key} {SC.fmt_esd(float(cell[k]), esd.get(n))}")
    L.append(f"_cell_volume {SC.fmt_esd(float(ph['volume']), ph.get('volume_esd'))}")
    L += [f"_pd_proc_ls_prof_wR_factor {100 * res['stats']['Rwp']:.2f}", f"_pd_proc_ls_prof_R_factor {100 * res['stats']['Rp']:.2f}"]
    return "\n".join(L) + "\n"


# ------------------------------------------------------------------------------------------------------------------------------ 3D
def payload_for_phase(res, i, ph_obj, repeat=(1, 1, 1)):
    """Data for the 3D viewer (mfs.structview.payload) of phase i in its refined cell, or None."""
    from mfs import structview as SV
    ph = res["phases"][i]
    mm = model_for_phase(ph_obj, ph) if ph_obj is not None else None
    if mm is None:
        return None
    mdl, sd = mm
    B = None if not ph.get("B") else np.full(mdl.n_atoms, float(ph["B"]))
    return SV.payload(mdl, np.zeros(mdl.n_pos), [float(v) for v in ph["cell"]], B, mdl.labels(), repeat)


# ------------------------------------------------------------------------------------------------------------------------------ pattern and reflections
def xy_text(x, y, header=""):
    buf = io.StringIO()
    if header:
        buf.write("# " + header + "\n")
    for a, b in zip(x, y):
        buf.write(f"{a:.5f} {b:.5f}\n")
    return buf.getvalue()


def pattern_files(res, name="pattern"):
    """{filename: bytes}: observed / calculated / difference / background as two-column .xy files, one CSV with all columns, the reflection tick positions per phase."""
    x = np.asarray(res["x"], float)
    files = {}
    for key, lab in (("y", "observed"), ("calc", "calculated"), ("resid", "difference (obs - calc)"), ("bg", "background")):
        short = {"y": "observed", "calc": "calculated", "resid": "difference", "bg": "background"}[key]
        files[f"{name}_{short}.xy"] = xy_text(x, np.asarray(res[key], float), f"{lab}; 2theta (deg), intensity; Reticula").encode("utf-8")
    df = pd.DataFrame({"two_theta": x, "observed": res["y"], "calculated": res["calc"], "background": res["bg"], "difference": res["resid"]})
    files[f"{name}_pattern.csv"] = df.to_csv(index=False).encode("utf-8")
    for k, p in enumerate(res["phases"]):
        t = np.asarray(p.get("ticks", []), float)
        if len(t):
            files[f"{name}_ticks_{k + 1}_{_safe(p['name'])}.txt"] = ("\n".join(f"{v:.4f}" for v in t) + "\n").encode("utf-8")
    return files


def reflection_table(res, i, ph_obj):
    """Reflection list of phase i: representative h k l of every d-spacing group, d (Å), 2θ (°), and the Le Bail intensity (relative) or multiplicity·|F|² (model)."""
    ph = res["phases"][i]
    if ph_obj is None or ph_obj.hkl is None:
        return pd.DataFrame(columns=["h", "k", "l", "d (Å)", "2theta (°)", "intensity"])
    cell = [float(v) for v in ph["cell"]]
    from mfs.rietveld import _cell_metric
    Gs = np.linalg.inv(_cell_metric(cell))
    hkl = np.asarray(ph_obj.hkl, float)
    d = 1.0 / np.sqrt(np.einsum("ki,ij,kj->k", hkl, Gs, hkl))
    lam = _wavelength(res)
    tt = 2 * np.degrees(np.arcsin(np.clip(lam / (2 * d), 0, 1))) if lam else np.full(len(d), np.nan)
    li = res.get("lebail_I")
    inten = np.asarray(li[i], float) if li and i in li else np.asarray(ph_obj.I0, float)
    df = pd.DataFrame({"h": hkl[:, 0].astype(int), "k": hkl[:, 1].astype(int), "l": hkl[:, 2].astype(int), "d (Å)": d, "2theta (°)": tt,
                       "intensity": inten / max(float(np.max(inten)), 1e-12) * 100.0})
    return df.sort_values("d (Å)", ascending=False).reset_index(drop=True)


def hkl_file(table: pd.DataFrame) -> str:
    """SHELX-style .hkl (h k l I sigma, 3I4 + 2F8.2) from the reflection list (relative intensities; sigma = 5 % + 1)."""
    L = []
    for r in table.itertuples(index=False):
        I = float(r[5])
        L.append(f"{int(r[0]):4d}{int(r[1]):4d}{int(r[2]):4d}{I:8.2f}{max(0.05 * abs(I), 1.0):8.2f}")
    L.append(f"{0:4d}{0:4d}{0:4d}{0.0:8.2f}{0.0:8.2f}")
    return "\n".join(L) + "\n"


def _safe(s):
    return "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in str(s))


def _wavelength(res):
    """Main wavelength (Å) of the radiation of a result, or None."""
    try:
        from mfs.rietveld import RADIATION
        r = res.get("radiation")
        return float(RADIATION[r][0][0]) if isinstance(r, str) and r in RADIATION else None
    except Exception:
        return None


# ------------------------------------------------------------------------------------------------------------------------------ Streamlit
def render(res, i, ph_obj, nm, g, version=""):
    """3D view + output downloads of phase i for a result WITHOUT atom refinement (Rietveld with a known structure, Le Bail)."""
    import streamlit as st
    from mfs import structview as SV
    dl = g["dl"]
    ph = res["phases"][i]
    has_struct = getattr(ph_obj, "_structure", None) is not None
    st.markdown(f"#### 🧱 Structure and files – {ph['name']}")
    tabs = st.tabs(["🧊 3D structure", "⬇️ Files"]) if has_struct else st.tabs(["⬇️ Files"])
    k = 0
    if has_struct:
        with tabs[0]:
            c1, c2 = st.columns([1, 3])
            rep = c1.selectbox("Show", ["1 × 1 × 1 cell", "2 × 2 × 2 cells"], key=f"ro_rep_{i}")
            c2.caption("The structure of the model in the refined cell (positions of the starting model unless the atoms were refined). Drag = rotate · wheel = zoom · click an atom for its bonds.")
            data = payload_for_phase(res, i, ph_obj, (2, 2, 2) if rep.startswith("2") else (1, 1, 1))
            if data:
                SV.show(st, data, 600)
        k = 1
    with tabs[k]:
        base = _safe(f"{nm}_{ph['name']}")
        f1, f2, f3 = st.columns(3)
        dl(f1, "⬇️ CIF (refined cell)" if has_struct else "⬇️ CIF (cell only)", cif_for_phase(res, i, ph_obj, f"{nm} {ph['name']}", version).encode("utf-8"),
           f"{base}.cif", f"ro_cif_{i}", mime="chemical/x-cif",
           help="The refined cell with esd" + (", the symmetry and the atoms of the model." if has_struct else ". A cell-only file: no atoms yet."))
        rt = reflection_table(res, i, ph_obj)
        if len(rt):
            dl(f2, "⬇️ Reflection list (CSV)", rt.to_csv(index=False).encode("utf-8"), f"{base}_reflections.csv", f"ro_ref_{i}",
               help="h k l, d, 2θ and the intensity (Le Bail: extracted from the pattern).")
            dl(f3, "⬇️ Reflections (.hkl)", hkl_file(rt).encode("utf-8"), f"{base}.hkl", f"ro_hkl_{i}", help="h k l I σ in the fixed-width format of structure-solution programs (relative intensities).")
        if len(rt):
            with st.expander("Reflection list", expanded=False):
                st.dataframe(rt, hide_index=True, width="stretch", height=300)


def bundle(res, ph_objs, nm, version="", extra_files=None):
    """ZIP of everything to continue working with: CIF per phase, pattern .xy files, reflection lists, summary JSON (+ extra files)."""
    base = _safe(nm)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for k, p in enumerate(res["phases"]):
            po = ph_objs[k] if ph_objs and k < len(ph_objs) else None
            tag = f"{base}_{_safe(p['name'])}"
            z.writestr(f"{tag}.cif", cif_for_phase(res, k, po, f"{nm} {p['name']}", version))
            rt = reflection_table(res, k, po)
            if len(rt):
                z.writestr(f"{tag}_reflections.csv", rt.to_csv(index=False))
                z.writestr(f"{tag}.hkl", hkl_file(rt))
        for fn, b in pattern_files(res, base).items():
            z.writestr(fn, b)
        summ = {"program": "Reticula " + version, "mode": res["mode"], "radiation": res.get("radiation"),
                "Rwp_pct": 100 * res["stats"]["Rwp"], "Rp_pct": 100 * res["stats"]["Rp"],
                "phases": [{"name": p["name"], "cell": [float(v) for v in p["cell"]], "cell_esd": {k: (None if v is None else float(v)) for k, v in (p.get("cell_esd") or {}).items()},
                            "volume": float(p["volume"]), "weight_fraction": None if p.get("weight_fraction") != p.get("weight_fraction") else p.get("weight_fraction")} for p in res["phases"]],
                "notes": list(res.get("notes", []))}
        z.writestr(f"{base}_summary.json", json.dumps(summ, indent=1, ensure_ascii=False, default=str))
        for fn, b in (extra_files or {}).items():
            z.writestr(fn, b)
    return buf.getvalue()
