"""Reticula page module: Structure (CIF) page.

Moved verbatim out of app.py in v2.13. Names that the code takes from the app (helpers such as show(), dl(), files_input(),
the sidebar values WIDTH / ASPECT / FONT, the cached reader wrapper R …) are listed in INJECTED and filled in by
mfs.pages.bind(globals()) on every run of app.py, so the code behaves exactly as it did inside app.py."""
# ruff: noqa: F821
from __future__ import annotations
import io, os, re, json, glob, zipfile, warnings, datetime, hashlib
import numpy as np
import pandas as pd
import streamlit as st
from mfs import numentry as NE  # every slider takes a typed number too (2.19.6)
import matplotlib.pyplot as plt

from mfs import analysis as A, plots as P, style as S, defects as D, composer as C, report as REP
from mfs import figctl as F, figctl_ui as FU, ui as UI, __version__
from mfs.chem import SPECIES, FRAMEWORKS, molar_mass

INJECTED = (
    'REF',
    'dl',
    'show',
)


@st.cache_data(show_spinner=False, max_entries=8)
def _cif_text(raw: bytes):
    return raw.decode("utf-8", "replace")


@st.cache_data(show_spinner=False, max_entries=8)
def _cif_pxrd(text, lo, hi, fw):
    from mfs import cifpxrd as CPX
    return CPX.simulate(text.encode(), (lo, hi), fwhm=fw)


def _cif_struct(text, name):
    from pathlib import Path  # noqa
    from mfs import cif as CX
    key = ("cifS", hash(text))
    if key not in st.session_state:
        st.session_state[key] = CX.read(text.encode(), name)
    return st.session_state[key]


def page_structure():
    from pathlib import Path
    from mfs import cif as CX, cifplot as CP, cifsi as CS, cifpxrd as CPX
    st.header("Structure (CIF)")
    st.caption("Single-crystal or MOF CIF → ORTEP / ball-and-stick / stick / space-filling figures, H-bond tables with "
               "symmetry codes, SI tables (crystal data, bonds, angles, torsions with ESDs) as DOCX/XLSX, simulated PXRD "
               "with hkl labels (and overlay on a measured pattern) and, for MOFs, geometric pore analysis. Pure matplotlib.")
    refs = {p.name: p for p in sorted(Path(REF).glob("*.cif"))}
    c1, c2 = st.columns([2, 1])
    up = c1.file_uploader("CIF file", type=["cif"], key="cx_up")
    ex = c2.selectbox("…or a bundled example", ["—"] + list(refs), key="cx_ex")
    if up is not None:
        text, name = _cif_text(up.getvalue()), up.name
    elif ex != "—":
        text, name = refs[ex].read_text(errors="replace"), ex
    else:
        st.info("Upload a CIF (SHELXL/Olex2 output or a COD/CSD file) or pick an example."); return
    try:
        blocks = _cif_struct(text, Path(name).stem)
    except Exception as e:
        st.error(f"Could not read the CIF: {e}"); return
    if not blocks:
        st.error("No data block with atoms found."); return
    struct = blocks[st.selectbox("Data block", range(len(blocks)), format_func=lambda i: blocks[i].name, key="cx_blk")] if len(blocks) > 1 else blocks[0]
    kn = CS.key_numbers(struct)
    st.success(" · ".join(f"{k} {v}" for k, v in kn.items() if v not in ("", None)))
    is_mof = any(a["el"] in CX.METALS for a in struct.atoms)
    tabs = st.tabs(["Figure", "H-bonds", "SI tables", "Simulated PXRD"] + (["Pores (MOF)"] if is_mof else []) + ["Hirshfeld (experimental)"])
    # ---- figure
    with tabs[0]:
        a1, a2, a3, a4 = st.columns(4)
        style = a1.selectbox("Style", ["ortep", "ball", "stick", "space"], key="cx_style",
                             format_func={"ortep": "ORTEP ellipsoids", "ball": "Ball-and-stick", "stick": "Stick", "space": "Space-filling"}.get)
        prob = NE.select_slider("Ellipsoid probability", [0.3, 0.5, 0.7], 0.5, key="cx_prob", format_func=lambda p: f"{p:.0%}")
        mode = a3.selectbox("Grow", ["asu", "molecules", "packing", "cell", "supercell"], index=1 if not is_mof else 3, key="cx_mode",
                            format_func={"asu": "Asymmetric unit", "molecules": "Complete molecules", "packing": "Packing (unit cell)",
                                         "cell": "Atoms in the cell", "supercell": "Supercell"}.get)
        labels = a4.selectbox("Labels", ["non-H", "all", "none"], key="cx_lab")
        b1, b2, b3, b4 = st.columns(4)
        hmode = b1.selectbox("Hydrogen atoms", ["all", "polar", "none"], key="cx_h",
                             format_func={"all": "All H", "polar": "Polar H only (O–H, N–H)", "none": "No H"}.get)
        showh = hmode != "none"
        has_dis = any(str(a["dis_group"]).strip(".?-") for a in struct.atoms)
        part = b2.selectbox("Disorder", ["major", "minor", "all"] if has_dis else ["all"], key="cx_part",
                            format_func={"major": "Major part (1)", "minor": "Minor part (2)", "all": "All parts"}.get)
        cell = b3.checkbox("Cell edges", True, key="cx_cell"); triad = b3.checkbox("Axis triad", True, key="cx_triad")
        hbo = b4.checkbox("H-bonds (dashed)", True, key="cx_hb")
        poly = b4.checkbox("Metal polyhedra", is_mof, key="cx_poly", disabled=not is_mof)
        e1, e2, e3 = st.columns([1, 1, 2])
        hsc = NE.slider("H atom size", 0.3, 1.5, 0.7, 0.1, key="cx_hsc", help="Scale of the H spheres (ORTEP / ball-and-stick)")
        ttl = e2.checkbox("Title on figure", False, key="cx_ttl", help="Off for publication exports (caption goes in the manuscript)")
        offs_txt = e3.text_input("Manual label offsets (Å, screen x y)", "", key="cx_offs", placeholder="C10A 0 -0.8; O5A 0.5 0",
                                 help="Shift single labels after the automatic placement, e.g. 'C10A 0 -0.8; O5A 0.5 0'")
        offs = {}
        for chunk in offs_txt.split(";"):
            tok = chunk.split()
            if len(tok) == 3:
                try:
                    offs[tok[0]] = (float(tok[1]), float(tok[2]))
                except ValueError:
                    st.warning(f"Offset '{chunk.strip()}' not understood (label dx dy)")
        sc = (1, 1, 1)
        if mode == "supercell":
            s1, s2, s3 = st.columns(3)
            sc = (s1.number_input("na", 1, 4, 2, key="cx_na"), s2.number_input("nb", 1, 4, 1, key="cx_nb"), s3.number_input("nc", 1, 4, 1, key="cx_nc"))
        v1, v2, v3 = st.columns([1, 1, 2])
        view = v1.selectbox("View along", ["a", "b", "c", "uvw", "hkl", "best"], key="cx_view",
                            format_func={"a": "a axis", "b": "b axis", "c": "c axis", "uvw": "vector [uvw]", "hkl": "plane normal (hkl)", "best": "best view (mean plane)"}.get)
        vec = v2.text_input("u v w / h k l", "1 1 0", key="cx_vec", disabled=view not in ("uvw", "hkl"))
        with v3:
            r1, r2, r3 = st.columns(3)
            rot = (NE.slider("rot x (°)", -180, 180, 0, 5, key="cx_rx"), NE.slider("rot y (°)", -180, 180, 0, 5, key="cx_ry"),
                   NE.slider("rot z (°)", -180, 180, 0, 5, key="cx_rz"))
        try:
            vv = [float(t) for t in vec.replace(",", " ").split()][:3]
        except ValueError:
            vv = [1, 1, 0]
        gmode = "packing" if (mode == "supercell" and not is_mof) else mode   # molecular crystals: whole molecules per cell
        sites, bonds = struct.grow(gmode, part=part, supercell=sc, include_h=showh)
        if len(sites) >= 6000:
            st.warning("Atom limit (6000) reached; drawing a partial structure.")
        if hmode == "polar":
            sites, bonds = CP.polar_h_only(sites, bonds)
        V = CP.view_matrix(struct, view, uvw=vv if view == "uvw" else None, hkl=vv if view == "hkl" else None, rot=rot, sites=sites)
        hb = CX.hbond_pairs(struct, sites) if (hbo and showh) else []
        n_mol = len(CX._components(sites, bonds)) if sites else 0
        st.caption(f"{len(sites)} atoms, {len(bonds)} bonds, {n_mol} fragment(s), {len(hb)} H-bond(s) drawn"
                   + (" · polymeric framework (molecules = connected network)" if getattr(struct, "polymeric", False) else ""))
        fig, ax = plt.subplots(figsize=(4.5, 4.5))
        CP.draw(ax, struct, sites, bonds, style=style, prob=prob, labels=labels, hbonds=hb, show_cell=cell, triad=triad, view=V,
                polyhedra=poly and is_mof, cells=sc, h_scale=hsc, label_offsets=offs)
        if ttl:
            ax.set_title(f"{struct.name}: {dict(asu='asymmetric unit', molecules='molecules', packing='packing', cell='unit cell', supercell='supercell')[mode]}"
                         f", view {view}", fontsize=7)
        if style == "ortep":
            ax.text(0.01, 0.01, f"{prob:.0%} probability ellipsoids", transform=ax.transAxes, fontsize=6, color="#555")
        xyz = pd.DataFrame([{"label": s["label"], "el": s["el"], "x": s["xyz"][0], "y": s["xyz"][1], "z": s["xyz"][2],
                             "symmetry": struct.symcode(s["op"], s["T"])} for s in sites])
        show(fig, f"{struct.name}_{style}_{mode}_{view}", data=xyz, controls=False)
    # ---- H-bonds
    hbr = CX.hbonds(struct)
    with tabs[1]:
        if hbr:
            d = pd.DataFrame(hbr).drop(columns=["_A"])
            d["D···A"] = [CX.fmt_esd(v, e) for v, e in zip(d["D···A (Å)"], d["σ D···A"])]
            st.dataframe(d.drop(columns=["σ D···A"]).round(3), hide_index=True, width="stretch")
            st.caption("D–H···A with H···A < Σr_vdW (Bondi) and ∠DHA ≥ 120°; symmetry code of the acceptor "
                       "(SHELX style: operator + lattice translation). D···A ESD from the coordinate s.u. only.")
            dl(st, "H-bonds (CSV)", lambda: d.to_csv(index=False).encode(), f"{struct.name}_hbonds.csv", "cx_hb_csv")
        else:
            st.info("No D–H···A hydrogen bonds found (no H on N/O donors, or no contacts within the vdW criterion).")
    # ---- SI tables
    with tabs[2]:
        ih = st.checkbox("Include bonds to H", False, key="cx_si_h")
        po = st.checkbox("Only rows flagged for publication (_geom_*_publ_flag yes)", False, key="cx_si_publ")
        geom = CX.geom_tables(struct, include_h=ih, publ_only=po)
        n0 = st.number_input("First table number (S#)", 1, 99, 1, key="cx_si_n")
        secs = CS.si_sections(struct, geom, hbr, int(n0))
        for title, df, cols, note in secs:
            with st.expander(f"{title} ({len(df)} rows)", expanded=title.startswith(f"Table S{int(n0)}.")):
                st.dataframe((df[cols] if cols else df).astype(str), hide_index=True, width="stretch")
                if note:
                    st.caption(note)
        c1, c2 = st.columns(2)
        dl(c1, "SI tables (DOCX)", lambda: CS.docx(secs), f"{struct.name}_SI_tables.docx", "cx_si_docx")
        dl(c2, "SI tables (XLSX)", lambda: CS.xlsx(secs), f"{struct.name}_SI_tables.xlsx", "cx_si_xlsx")
        _si_put = st.session_state.setdefault("si_tables_extra", {}); _si_put[struct.name] = secs
    # ---- PXRD
    with tabs[3]:
        p1, p2, p3 = st.columns(3)
        lo, hi = NE.slider("2θ range (°)", 2.0, 90.0, (3.0, 50.0), 0.5, key="cx_px_rng")
        fw = p2.number_input("FWHM (°)", 0.02, 1.0, 0.10, 0.01, key="cx_px_fw")
        nl = p3.number_input("Label strongest n", 0, 40, 12, key="cx_px_nl")
        with st.spinner("Simulating…"):
            try:
                sim, sticks = _cif_pxrd(text, lo, hi, fw)
            except Exception as e:
                st.error(f"Simulation failed: {e}"); sim = None
        if sim is not None:
            mf = st.file_uploader("Measured pattern for overlay (.xy / .txt / .csv)", key="cx_px_meas")
            meas = fitr = None
            if mf is not None:
                try:
                    from mfs import readers as RD
                    m = RD.read_xy(mf); mx, my = np.asarray(m.iloc[:, 0], float), np.asarray(m.iloc[:, 1], float)
                    sel = (mx >= lo) & (mx <= hi); meas = (mx[sel], my[sel]); fitr = CPX.fit_to(meas[0], meas[1], sim)
                    st.caption(f"Least-squares scale {fitr['scale']:.4g} + linear background; Rwp = {fitr['Rwp']:.3f} "
                               "(profile agreement only — not a Rietveld refinement).")
                except Exception as e:
                    st.error(f"Measured pattern not read: {e}")
            if meas is not None:
                fig, (ax, rax) = plt.subplots(2, 1, figsize=(5.5, 3.6), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
            else:
                fig, ax = plt.subplots(figsize=(5.5, 3.0)); rax = None
            CPX.draw(ax, sim, sticks, label_top=int(nl), meas=meas, fitres=fitr, rax=rax, name=struct.name)
            show(fig, f"{struct.name}_simulated_PXRD", data=pd.DataFrame({"2theta": sim.x, "I_sim": sim.y}))
            st.dataframe(CPX.table(sticks, 30), hide_index=True, width="stretch")
    # ---- pores
    if is_mof:
        with tabs[4]:
            from mfs import poretools as PT
            q1, q2, q3 = st.columns(3)
            spc = NE.select_slider("Grid spacing (Å)", [0.25, 0.3, 0.4, 0.5, 0.7], 0.4, key="cx_pt_sp")
            pr = q2.number_input("Probe radius (Å)", 0.5, 3.0, 1.82, 0.01, key="cx_pt_pr", help="1.82 Å ≈ N₂ (kinetic diameter 3.64 Å)")
            if q3.button("Analyse pores", key="cx_pt_go", type="primary"):
                with st.spinner("Grid + percolation + Monte Carlo…"):
                    st.session_state["cx_pt_res"] = (struct.name, PT.analyse(struct, spc, pr, part="major"))
            r = st.session_state.get("cx_pt_res")
            if r and r[0] == struct.name:
                res = r[1]
                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Di (largest included sphere)", f"{res['Di (Å)']:.2f} Å"); m2.metric("Df (largest free sphere)", f"{res['Df (Å)']:.2f} Å")
                m3.metric("Void fraction", f"{res['void fraction (probe 0)']:.3f}"); m4.metric("Accessible SA", f"{res['ASA (m² g⁻¹)']:.0f} m² g⁻¹")
                tbl = pd.DataFrame([(kk, (f"{v:.4g}" if isinstance(v, (float, np.floating)) else v)) for kk, v in res.items() if kk != "method"],
                                   columns=["quantity", "value"]).astype(str)
                st.dataframe(tbl, hide_index=True, width="stretch")
                st.caption("Method: " + res["method"] + " Compare with Zeo++/RASPA for publication; values depend on radii and grid.")
                dl(st, "Pore analysis (CSV)", lambda: pd.concat([tbl, pd.DataFrame([["method", res["method"]]], columns=tbl.columns)]).to_csv(index=False).encode(),
                   f"{struct.name}_pores.csv", "cx_pt_csv")
    # ---- Hirshfeld
    with tabs[-1]:
        st.warning("EXPERIMENTAL: approximate promolecule densities (one exponential per atom, not the Spackman tables of "
                   "CrystalExplorer). Shapes and contact percentages are indicative; use CrystalExplorer for publication.")
        if is_mof or getattr(struct, "polymeric", False):
            st.info("Hirshfeld surfaces need a discrete molecule; not available for frameworks.")
        elif st.button("Compute Hirshfeld surface (first molecule)", key="cx_hs_go"):
            from mfs import hirshfeld as HS
            with st.spinner("Ray search for w = 0.5…"):
                mol, mb = struct.grow("molecules", part="major")
                comps = CX._components(mol, mb); m = [mol[i] for i in comps[0]]
                env, _ = struct.grow("packing", part="major", supercell=(3, 3, 3))
                env = [dict(e, frac=e["frac"] - 1, xyz=struct.cart(e["frac"] - 1)) for e in env]
                st.session_state["cx_hs"] = (struct.name, HS.surface(struct, m, env, n_dir=160))
        r = st.session_state.get("cx_hs")
        if r and r[0] == struct.name:
            from mfs import hirshfeld as HS
            hs = r[1]
            fig, (a1_, a2_) = plt.subplots(1, 2, figsize=(6.5, 3.2))
            sc_ = HS.draw_dnorm(a1_, hs, CP.view_matrix(struct, "a")); fig.colorbar(sc_, ax=a1_, shrink=0.7, label=r"$d_\mathrm{norm}$")
            HS.draw_fingerprint(a2_, hs, highlight=["O···H", "H···O"])
            show(fig, f"{struct.name}_Hirshfeld_experimental", data=pd.DataFrame({"di": hs["di"], "de": hs["de"], "d_norm": hs["d_norm"], "pair": hs["pair"]}), controls=False)
            st.dataframe(HS.contacts(hs).round(1), hide_index=True)
