"""v2.16: the results of a structure refinement on the PXRD page - 3D viewer, tables with esd, publication figure, CIF / Excel output and a Methods paragraph.
Rendered from mfs/pages/spectra.py (_rietveld_ui) with the page helpers passed in (show, newfig, dl)."""
from __future__ import annotations

import io
import math
import re

import numpy as np
import pandas as pd
import streamlit as st

from mfs import structcif as SC
from mfs import structview as SV


def _fmt(v, e, dec=5):
    return SC.fmt_esd(float(v), e if (e and e > 0) else None)


def atoms_table(sd: dict) -> pd.DataFrame:
    rows = []
    for a in sd["atoms"]:
        rows.append({"atom": a["label"], "element": a["element"], "x": _fmt(a["x"], a["x_esd"]), "y": _fmt(a["y"], a["y_esd"]), "z": _fmt(a["z"], a["z_esd"]),
                     "B (Å²)": _fmt(a["B"], a.get("B_esd")), "U_iso (Å²)": _fmt(a["B"] / (8 * math.pi ** 2), (a["B_esd"] / (8 * math.pi ** 2)) if a.get("B_esd") else None),
                     "occupancy": round(a["occupancy"], 3), "atoms in cell": a["n_equivalent"], "free coordinates": a["free_coordinates"]})
    return pd.DataFrame(rows)


def bonds_table(sd: dict) -> pd.DataFrame:
    return pd.DataFrame([{"atom 1": b["atom1"], "atom 2": b["atom2"], "distance (Å)": _fmt(b["distance"], b["esd"]), "× in the cell": b["count"]}
                         for b in sorted(sd.get("bonds", []), key=lambda b: (b["atom1"], b["distance"]))])


def angles_table(sd: dict) -> pd.DataFrame:
    return pd.DataFrame([{"atom 1": g["atom1"], "central atom": g["centre"], "atom 2": g["atom2"], "angle (°)": _fmt(g["angle"], g["esd"])}
                         for g in sorted(sd.get("angles", []), key=lambda g: (g["centre"], g["atom1"], g["atom2"]))])


def methods_paragraph(res: dict, phase: dict, version: str, radiation: str) -> str:
    sd = phase["structure"]
    st_ = res["stats"]
    return (f"The crystal structure of {phase['name']} was refined from the powder diffraction pattern by the Rietveld method with Reticula {version} "
            f"({radiation}, {res['x'].min():.1f}–{res['x'].max():.1f}° 2θ). The starting model was taken from the crystallographic information file of the phase; "
            f"the space group {sd['spacegroup']} (No. {sd['sg_number']}) was kept and the atom positions were refined within the freedom allowed by the site symmetry, "
            f"together with one isotropic displacement parameter per element, the lattice parameters, the zero / displacement offset, a Chebyshev background "
            f"({res['bg_order']} terms) and a pseudo-Voigt (Thompson–Cox–Hastings) peak profile. Soft restraints on the bond lengths of the starting model were used. "
            f"The refinement converged with R_wp = {100 * st_['Rwp']:.2f} %, R_p = {100 * st_['Rp']:.2f} %"
            + (f", GoF = {math.sqrt(st_['chi2']):.2f}" if res.get("counts") else "") + (f" and R_Bragg = {100 * phase['Rb']:.2f} %" if phase.get("Rb") is not None else "")
            + f". The lattice parameters are a = {_fmt(phase['cell'][0], phase['cell_esd'].get('a'))} Å"
            + (f", c = {_fmt(phase['cell'][2], phase['cell_esd'].get('c'))} Å" if phase["system"] in ("tetragonal", "hexagonal", "orthorhombic", "monoclinic") else "") + ".")


def render(res: dict, i: int, nm: str, g: dict, version: str = ""):
    """Everything about the refined structure of phase i of the result `res`."""
    show, newfig, dl = g["show"], g["newfig"], g["dl"]
    ph = res["phases"][i]
    sd = ph.get("structure")
    if not sd:
        return
    stt = sd.get("_state")
    st.markdown(f"#### 🧱 Refined structure – {ph['name']}")
    m = st.columns(5)
    m[0].metric("Space group", sd["spacegroup"], f"No. {sd['sg_number']}", delta_color="off")
    m[1].metric("Atoms in the cell", sd["n_atoms_cell"], f"{len(sd['atoms'])} independent", delta_color="off")
    m[2].metric("Rwp", f"{100 * res['stats']['Rwp']:.1f} %")
    m[3].metric("Bragg R", f"{100 * ph['Rb']:.1f} %" if ph.get("Rb") is not None else "–")
    dev = sd.get("bond_rms_deviation")
    m[4].metric("Geometry vs start", f"{dev:.3f} Å rms" if dev is not None else "–", help="Root-mean-square change of the bond lengths with respect to the starting model.")
    for v in sd.get("verdict", []):
        (st.success if v["level"] == "ok" else st.warning if v["level"] == "warn" else st.info)(v["text"])
    if sd.get("report_error"):
        st.caption(f"(statistics partly unavailable: {sd['report_error']})")
    tabs = st.tabs(["🧊 3D structure", "📋 Tables", "🖼️ Figure", "⬇️ Output"])
    # ------------------------------------------------------------------------------ 3D
    with tabs[0]:
        c1, c2 = st.columns([1, 3])
        rep = c1.selectbox("Show", ["1 × 1 × 1 cell", "2 × 2 × 2 cells"], key=f"sv_rep_{i}")
        c2.caption("Drag = rotate · wheel = zoom · shift+drag = move · click an atom for its bonds · buttons: views along a / b / c, thermal spheres (displacement size), labels, PNG.")
        if stt:
            r = (2, 2, 2) if rep.startswith("2") else (1, 1, 1)
            data = SV.payload(stt["model"], stt["p"], stt["cell"], stt["B"], stt["labels"], r, cn=sd.get("cn"))
            SV.show(st, data, 600)
    # ------------------------------------------------------------------------------ tables
    with tabs[1]:
        at = atoms_table(sd)
        st.markdown("**Atoms** (fractional coordinates, esd in brackets)")
        st.dataframe(at, hide_index=True, width="stretch")
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**Bond lengths**")
            st.dataframe(bonds_table(sd), hide_index=True, width="stretch", height=300)
        with c2:
            st.markdown("**Bond angles**")
            st.dataframe(angles_table(sd), hide_index=True, width="stretch", height=300)
        d1, d2 = st.columns(2)
        d1.markdown("**Coordination** (number of bonded neighbours, mean distance)")
        d1.dataframe(pd.DataFrame([{"atom": c["atom"], "CN": c["CN"], "mean distance (Å)": round(c["mean_d"], 3) if c["mean_d"] == c["mean_d"] else None}
                                   for c in sd.get("coordination", [])]), hide_index=True, width="stretch")
        d2.markdown("**Bond-valence sums** (Brese & O'Keeffe; cation–O)")
        d2.dataframe(pd.DataFrame([{"atom": b["atom"], "BVS (v.u.)": round(b["BVS"], 2), "formal charge": b["formal"], "deviation (%)": round(b["deviation_pct"], 1)}
                                   for b in sd.get("bvs", [])]), hide_index=True, width="stretch")
    # ------------------------------------------------------------------------------ figure
    with tabs[2]:
        if not stt:
            st.info("The structure is not available for the figure.")
        else:
            f1, f2, f3, f4 = st.columns(4)
            view = f1.selectbox("View", ["perspective", "along a", "along b", "along c"], key=f"sf_view_{i}")
            az = f2.slider("Azimuth (°)", 0, 360, 35, 5, key=f"sf_az_{i}", disabled=view != "perspective")
            el = f3.slider("Elevation (°)", -90, 90, 25, 5, key=f"sf_el_{i}", disabled=view != "perspective")
            size = f4.slider("Atom size", 0.4, 2.0, 1.0, 0.1, key=f"sf_size_{i}")
            g1, g2, g3 = st.columns(3)
            cell_on = g1.checkbox("Unit cell", True, key=f"sf_cell_{i}")
            lab_on = g2.checkbox("Atom labels", False, key=f"sf_lab_{i}")
            h_off = g3.checkbox("Hide hydrogen", True, key=f"sf_h_{i}")
            data = SV.payload(stt["model"], stt["p"], stt["cell"], stt["B"], stt["labels"], (1, 1, 1), cn=sd.get("cn"))
            along = {"along a": "a", "along b": "b", "along c": "c"}.get(view)

            def draw(ax, data=data):
                SV.draw_structure(ax, data, azim=az, elev=el, along=along, show_cell=cell_on, size=size, labels=lab_on, hide_h=h_off)
            fig, ax = newfig(1.0)
            draw(ax)
            show(fig, f"Structure_{re.sub(r'[^A-Za-z0-9]+', '_', ph['name'])}", at if False else atoms_table(sd), draw, "Crystal structure")
    # ------------------------------------------------------------------------------ output
    with tabs[3]:
        o1, o2 = st.columns(2)
        base = re.sub(r"[^A-Za-z0-9_.-]+", "_", f"{nm}_{ph['name']}")
        if stt:
            cif = SC.write_cif(ph, stt["model"], f"{nm} {ph['name']}", None, res["stats"], version)
            dl(o1, "⬇️ Refined structure (CIF with esd)", cif.encode("utf-8"), f"{base}_refined.cif", f"sx_cif_{i}", mime="chemical/x-cif",
               help="Cell, atoms with esd, symmetry operations, displacement parameters, bonds, angles and the figures of merit.")

        def _xlsx():
            buf = io.BytesIO()
            with pd.ExcelWriter(buf, engine="openpyxl") as xw:
                atoms_table(sd).to_excel(xw, sheet_name="Atoms", index=False)
                bonds_table(sd).to_excel(xw, sheet_name="Bonds", index=False)
                angles_table(sd).to_excel(xw, sheet_name="Angles", index=False)
                pd.DataFrame(sd.get("coordination", [])).to_excel(xw, sheet_name="Coordination", index=False)
                pd.DataFrame(sd.get("bvs", [])).to_excel(xw, sheet_name="Bond valence", index=False)
            return buf.getvalue()
        dl(o2, "⬇️ Tables (Excel: atoms, bonds, angles, coordination, BVS)", _xlsx, f"{base}_structure_tables.xlsx", f"sx_xlsx_{i}")
        st.markdown("**Methods text** (edit freely)")
        st.text_area("Methods", methods_paragraph(res, ph, version, res.get("radiation", "")), height=190, key=f"sx_methods_{i}", label_visibility="collapsed")
