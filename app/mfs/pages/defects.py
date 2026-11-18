"""Reticula page module: Defect calculator page.

Moved verbatim out of app.py in v2.13. Names that the code takes from the app (helpers such as show(), dl(), files_input(),
the sidebar values WIDTH / ASPECT / FONT, the cached reader wrapper R …) are listed in INJECTED and filled in by
mfs.pages.bind(globals()) on every run of app.py, so the code behaves exactly as it did inside app.py."""
# ruff: noqa: F821
from __future__ import annotations
import io, os, re, json, glob, zipfile, warnings, datetime, hashlib
import numpy as np
import pandas as pd
import streamlit as st
from mfs import levels as LV
import matplotlib.pyplot as plt

from mfs import analysis as A, plots as P, style as S, defects as D, composer as C, report as REP
from mfs import figctl as F, figctl_ui as FU, ui as UI, __version__
from mfs.chem import SPECIES, FRAMEWORKS, molar_mass

INJECTED = (
    'EX',
    'EXPORT_DIR',
    'R',
    'REF',
    'STORE_PATH',
    'WIDTH',
    '_RENDER_LOCK',
    'dl',
    'newfig',
    'safe_name',
    'show',
    'store_put',
)


def ND_sites(nd):
    from mfs import nodes as N
    return N.sites(nd)


def _res_lab(metal, nd):
    from mfs import nodes as N
    if metal in ("Zr", "Hf", "mixed"):
        return "MO₂"
    r = nd.get("residue")
    return r[0] if r else N.residue_label(metal)


def guess_framework(name):
    """v2.11: framework + metal from a sample name ('MOF-801(Hf)', 'UiO-66-NH2 Zr', 'HKUST-1') -> (framework, metal|None)."""
    import re as _re
    n = (name or "").lower().replace("_", "-").replace(" ", "")
    fw = next((k for k in sorted(FRAMEWORKS, key=len, reverse=True) if k != "custom framework"
               and k.lower().replace(" ", "").replace("(", "").replace(")", "") in n.replace("(", "").replace(")", "")), None)
    met = None
    for m in ("Hf", "Zr", "Ce", "Th", "Cu", "Zn", "Al", "Fe", "Cr", "Mg", "Co", "Ni", "Ti"):
        if _re.search(r"(^|[^a-z])" + m.lower() + r"($|[^a-z])", (name or "").lower()):
            met = m; break
    return fw, met


def default_defect_input(sample="my sample", fw=None):
    """Neutral defaults (UiO-66 / Zr); a sample name that names a framework/metal pre-selects them."""
    gfw, gmet = guess_framework(sample)
    fw = fw or gfw or "UiO-66"
    d = FRAMEWORKS[fw]; lk = d["linker"]; f_, z_, _r, nH_, _c = SPECIES[lk]
    from mfs import nodes as N
    nd = N.get(d.get("node"))
    metal = gmet if gmet in nd["metals"] else d.get("metal", nd["metals"][0])
    sp = [D.Species(lk, f_, z_, "linker", nH_, None, False, True)]
    if z_ >= 2 and nd["core"] != "{M}":
        sp.append(D.Species("Formate", "CHO2", 1, "cap", 1, None, False, False))
    return REP.input_to_dict(D.DefectInput(sample=sample, framework=fw, metal=metal, species=sp))


def node_picker(where, fw, key, node_default="", metal_default=None, custom_default=None):
    """v2.11 node library picker shared by the defect calculator and the TGA page. Returns
    {"node": key, "metal": str, "node_custom": dict|None, "def": node definition}."""
    from mfs import nodes as N
    from mfs.chem import AW as _AW
    fwd = FRAMEWORKS.get(fw, {})
    keys = list(N.NODES) + [N.CUSTOM]
    nd0 = node_default or fwd.get("node") or N.DEFAULT_NODE
    c1, c2 = where.columns([2, 1])
    node = c1.selectbox("Node (SBU)", keys, index=keys.index(nd0) if nd0 in keys else 0, key=f"{key}_node|{fw}",
                        help="Metal-oxo/metal node per formula unit. All counts (linkers, caps, Cl, open sites) are per node.")
    custom = None
    if node == N.CUSTOM:
        cd = dict(custom_default or {})
        e = where.container(border=True)
        e.caption("Custom node: metal atoms per node, oxidation state, core O²⁻ / OH⁻ and the TGA residue. "
                  "Charge for ligands = n·val − 2·O − OH.")
        k1, k2, k3, k4, k5, k6 = e.columns(6)
        cm = k1.text_input("metal", cd.get("metal", "Zn"), key=f"{key}_cm|{fw}").strip() or "Zn"
        if cm not in _AW:
            e.error(f"Unknown element '{cm}'"); cm = "Zn"
        cn_ = k2.number_input("metals/node", 1, 24, int(cd.get("nM", 2)), key=f"{key}_cn|{fw}")
        cv = k3.number_input("oxidation state", 1, 6, int(cd.get("val", 2)), key=f"{key}_cv|{fw}")
        co = k4.number_input("core O²⁻", 0.0, 48.0, float(cd.get("core_O", 0)), 1.0, key=f"{key}_co|{fw}")
        coh = k5.number_input("core OH⁻", 0.0, 48.0, float(cd.get("core_OH", 0)), 1.0, key=f"{key}_coh|{fw}")
        rs = (cd.get("residue") or [N.residue_label(cm), N.RESIDUE.get(cm, (None, 1))[1]])
        rf = k6.text_input("TGA residue", rs[0], key=f"{key}_crf|{fw}|{cm}", help="oxide formula, e.g. Co3O4")
        try:
            from mfs.chem import molar_mass as _mm
            _mm(rf); rn = max(1, int(round(_mm_count(rf, cm))))
        except Exception:
            e.error(f"Cannot parse residue '{rf}'"); rf, rn = N.residue_label(cm), N.RESIDUE.get(cm, (None, 1))[1]
        custom = dict(nM=int(cn_), val=int(cv), core_O=co, core_OH=coh, metal=cm, residue=[rf, rn])
        nd = N.custom_node(**custom)
        e.caption(f"→ core {nd['core'].format(M=cm)}, {N.sites(nd):g} charges for ligands, residue {rf} ({_mm(rf)/rn:.2f} g mol⁻¹ per {cm})")
        metal = cm
    else:
        nd = N.NODES[node]
        ms = nd["metals"]
        m0 = metal_default if metal_default in ms else fwd.get("metal") if fwd.get("metal") in ms else ms[0]
        metal = c2.selectbox("Metal", ms, index=ms.index(m0), key=f"{key}_met|{fw}|{node}")
    return {"node": node, "metal": metal, "node_custom": custom, "def": nd}


def _mm_count(formula, el):
    """Number of atoms of element `el` in a simple formula (residue oxides)."""
    import re as _re
    tot = 0.0
    for sym, n in _re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", formula):
        if sym == el:
            tot += float(n) if n else 1.0
    return tot or 1.0


def page_defects():
    st.header("Defect calculator – shotgun every method, compare, get the minimal formula")
    st.caption("Per node (M₆O₄(OH)₄, paddlewheel, Zn₄O, M₃O, Al-OH chain, MOF-74, Ti₈, ZIF or your own – see the node library). Enter data once; every applicable method runs in parallel. Equations, symbols, assumptions, your own substitution and "
               "pitfalls are shown for each method. Unit tests reproduce the published worked examples (see README).")
    with LV.gate("Advanced", globals()):
        with st.expander("Load / save inputs (JSON)", expanded=False):
            exs = sorted(glob.glob(os.path.join(EX, "Defects", "defect_input_*.json")))
            exn = st.selectbox("Example input (synthetic EXAMPLE data)", ["(none)"] + [os.path.basename(p) for p in exs])
            upj = st.file_uploader("…or load your own JSON", type=["json"], key="dj")
            if st.button("Load"):
                src = json.load(upj) if upj else (json.load(open(os.path.join(EX, "Defects", exn), encoding="utf-8")) if exn != "(none)" else None)
                if src:
                    st.session_state["def_in"] = src; st.rerun()
    base = st.session_state.get("def_in") or default_defect_input("my sample")
    from mfs import samplestore as SS, si_export as SI
    with st.expander("🗂️ Sample store – pull TGA / ¹H NMR / ICP / XRF / N₂ results of one sample from the other pages", expanded=True):
        ids_ = SS.ids(st.session_state)
        if not ids_:
            st.caption("Empty. Open the TGA, NMR digestion, ICP/XRF and N₂ pages with your files (same sample names): every page "
                       f"stores its per-sample results here and in `{STORE_PATH}`.")
        if ids_:
            st.dataframe(pd.DataFrame(SS.summary_rows(st.session_state)), hide_index=True, width="stretch")
        pick = st.selectbox("Pull data for sample", ["(none – type the inputs)"] + list(ids_), format_func=lambda k: ids_.get(k, k), key="def_pull")
        k1, k2, k3 = st.columns(3)
        k1.download_button("⬇️ Store (JSON)", json.dumps(st.session_state.get("sample_store", {}), indent=1, ensure_ascii=False).encode(),
                           "sample_store.json", key="store_dl")
        upst = k2.file_uploader("Load a store (JSON)", type=["json"], key="store_up")
        if upst is not None and k2.button("Merge", key="store_merge"):
            st.session_state.setdefault("sample_store", {}).update(json.load(upst)); st.rerun()
        if k3.button("🗑️ Clear the store", key="store_clear"):
            SS.clear(st.session_state); SS.save(st.session_state, STORE_PATH); st.rerun()
    pulled, pull_note = {}, []
    if pick in ids_:
        rec = SS.get(st.session_state, ids_[pick]) or {"sections": {}}
        secs = {k: v["values"] for k, v in rec["sections"].items()}
        t_ = secs.get("tga") or {}
        if t_:
            pulled.update(tga_W_initial=t_.get("W_initial"), tga_W_200=t_.get("W200"), tga_W_plateau=t_.get("W_plateau"), tga_W_end=t_.get("W_end"))
            pull_note.append("TGA")
        i_ = secs.get("icp") or {}
        if i_:
            pulled.update(icp_Zr_wt=i_.get("Zr_wt"), icp_Hf_wt=i_.get("Hf_wt")); pull_note.append("ICP")
        x_ = secs.get("xrf") or {}
        if x_.get("xrf_S_per_M") is None and i_.get("xrf_S_per_M") is not None:
            x_ = {**x_, "xrf_S_per_M": i_["xrf_S_per_M"], "xrf_S_per_M_sd": i_.get("xrf_S_per_M_sd")}
        if x_.get("xrf_S_per_M") is not None:
            pulled["xrf_S_per_M"] = x_["xrf_S_per_M"]
            pull_note.append("S/M" + (f" = {x_['xrf_S_per_M']:.3f} ± {x_['xrf_S_per_M_sd']:.3f}" if x_.get("xrf_S_per_M_sd") is not None else ""))
        if (x_.get("Hf_frac") or i_.get("Hf_frac")) is not None:
            hf_ = x_.get("Hf_frac") if x_.get("Hf_frac") is not None else i_.get("Hf_frac")
            sd_ = x_.get("Hf_frac_sd") if x_.get("Hf_frac") is not None else i_.get("Hf_frac_sd")
            pull_note.append(f"Hf/(Zr+Hf) = {hf_:.3f}" + (f" ± {sd_:.3f}" if sd_ is not None else ""))
        n2_ = secs.get("n2") or {}
        if n2_.get("V_total") is not None:
            pulled["V_meas"] = n2_["V_total"]; pull_note.append("N₂ V_total")
        nm_ = secs.get("nmr") or {}
        if nm_:
            pulled["_nmr"] = nm_; pull_note.append("¹H NMR integrals")
        pulled["sample"] = rec.get("name", pick)
        if not st.session_state.get("def_in"):          # v2.11: framework/metal/linker defaults follow the sample name
            base = default_defect_input(pulled["sample"])
        st.success(f"Pulled for **{pulled['sample']}**: " + (", ".join(pull_note) or "nothing (no results stored yet)")
                   + ". Fields below are pre-filled and can still be edited.")
    base = dict(base); base.update({k: v for k, v in pulled.items() if not k.startswith("_") and v is not None})
    if pulled.get("_nmr"):      # integrals into the species table, matched by species name
        sp_ = [dict(x) for x in base.get("species", [])]
        for x in sp_:
            key = next((k for k in pulled["_nmr"] if k.startswith("I ") and x["name"].lower() in k.lower()), None)
            if key and pulled["_nmr"][key] is not None:
                x["integral"] = pulled["_nmr"][key]
        if not any(x.get("name", "").lower().startswith("tmsp") for x in sp_):
            pass
        base["species"] = sp_
        std = next((v for k, v in pulled["_nmr"].items() if k.startswith("I TMSP")), None)
        if std:
            base["nmr_std_integral"] = std
    sfx = "" if pick not in ids_ else "_" + pick[:20]
    g = lambda k, d=None: base.get(k, d)
    st.subheader("1 · Framework, node and metal")
    c1, c2 = st.columns([1, 1])
    sample = c1.text_input("Sample", g("sample"), key="def_sample" + sfx)
    fw0 = g("framework", "UiO-66"); fw0 = fw0 if fw0 in FRAMEWORKS else "custom framework"
    fw = c2.selectbox("Framework", list(FRAMEWORKS), index=list(FRAMEWORKS).index(fw0), key="def_fw" + sfx,
                      help="Picking a framework sets its node, metal and main linker (all editable below).")
    nsel = node_picker(st, fw, "def", g("node", "") if fw == fw0 else "", g("metal") if fw == fw0 else None,
                       g("node_custom") if fw == fw0 else None)
    metal, ND = nsel["metal"], nsel["def"]
    c1, c2, c3 = st.columns(3)
    cn_def = int(g("CN") or FRAMEWORKS[fw]["CN"]) if fw == fw0 else int(FRAMEWORKS[fw]["CN"] if not nsel["node_custom"] else ND["CN"])
    CN = c1.number_input("Ideal connectivity CN (charge units)", 1, 24, cn_def, key=f"def_cn|{fw}|{nsel['node']}" + sfx,
                         help="Linker donor groups per node in the defect-free framework, counted in charge units "
                              f"(12 for UiO-66, 8 NU-1000, 6 MOF-808 / Zn4O / M3O, 4 paddlewheel / MOF-74, 2 Al-OH chain / ZIF). "
                              f"This node offers {ND_sites(ND):g} charges for ligands.")
    f_hf = c2.number_input("Hf/(Hf+Zr) (mixed; blank → from ICP)", 0.0, 1.0, float(g("f_hf") or 0.0)) if metal == "mixed" else None
    extra = c3.checkbox("All Zr is extra-framework ZrO₂ (ZrO₂@Hf-MOF composites) – correct TGA residue, ICP, porosity",
                        bool(g("extraframework_Zr"))) if metal in ("Hf", "mixed") else False
    if fw != fw0 and not st.session_state.get("def_in"):
        base = dict(base); base["species"] = default_defect_input(sample, fw)["species"]
    with LV.gate("Advanced", globals()):
        with st.expander("🧱 Node library – metals, charges, TGA residues (pick above or build a custom node)", expanded=False):
            from mfs import nodes as _N
            st.dataframe(pd.DataFrame(_N.table_rows()), hide_index=True, width="stretch")
            st.caption("Charge for ligands = metals × oxidation state − 2·core O²⁻ − core OH⁻. TGA residue in air per metal "
                       "(Co → Co₃O₄, Fe/Cr/Al → M₂O₃, Mn → Mn₂O₃). Al-OH chains keep their µ-OH at the plateau; M₆O₄(OH)₄ and Ti₈O₈(OH)₄ "
                       "dehydroxylate first. Linker/cap molar masses are computed from the formula you type.")
    st.subheader("2 · Species (main linker, caps, guests) and ¹H NMR integrals")
    st.caption("Bound (deprotonated) formula; z = carboxylates bound (= node sites occupied, charge); role cap = monocarboxylate modulator. "
               "'at plateau' = still present at the ~400 °C TGA plateau (e.g. benzoate, sulfonated linkers). Integrals on any common scale "
               "(set the TMSP-d4 standard to 9 for the qNMR method).")
    sp = pd.DataFrame(g("species") if fw == fw0 else base["species"])
    sp = st.data_editor(sp, num_rows="dynamic", key="def_sp" + sfx + ("" if fw == fw0 else "|" + fw), column_config={
        "role": st.column_config.SelectboxColumn(options=["linker", "cap", "neutral"]),
        "name": st.column_config.TextColumn(help="Pick from the library below or type your own")})
    _mm = []
    for r_ in sp.to_dict("records"):
        try:
            _mm.append(f"{r_['name']}: {molar_mass(str(r_['formula'])):.2f}" if r_.get("name") and r_.get("formula") else None)
        except Exception as e_:
            _mm.append(f"⚠️ {r_.get('name')}: {e_}")
    st.caption("Molar masses (bound form, g mol⁻¹, computed from the formula – type any custom linker/cap): " + " · ".join(x for x in _mm if x))
    with LV.gate("Advanced", globals()):
        with st.expander("Species library"):
            lib = pd.DataFrame([(k, *v) for k, v in SPECIES.items()], columns=["name", "formula", "z", "role", "nH", "¹H signal"])
            lib["M (g/mol)"] = [round(molar_mass(f), 2) for f in lib.formula]
            st.dataframe(lib, width="stretch", hide_index=True)
    nmr = st.session_state.get("nmr_integrals")
    if nmr is not None:
        st.caption("Integrals from the NMR tab are available:"); st.dataframe(nmr.filter(regex="sample|^I "), width="stretch")
    st.subheader("3 · TGA (in air, % of initial mass)")
    up = st.file_uploader("Optional: TGA file → auto-read W(200 °C), plateau, residue", key="def_tga")
    auto = {}
    if up:
        t = R.read_tga(up, up.name); a = A.tga_analysis(t["T"], t["w"])
        auto = {"tga_W_initial": a["W_initial"], "tga_W_200": a["W200"], "tga_W_plateau": a["W_plateau"], "tga_W_end": a["W_end"]}
        st.info(f"Auto: W_i = {a['W_initial']:.2f}, W(200) = {a['W200']:.2f}, plateau {a['W_plateau']:.2f} % at {a['T_plateau']:.0f} °C (DTG minimum 250–480 °C), residue {a['W_end']:.2f} %")
    gv = lambda k: auto.get(k, g(k))
    c1, c2, c3, c4 = st.columns(4)
    num = lambda col, lab, key, v, **kw: col.number_input(lab, value=None if v is None else float(v), **kw, key="dn_" + key + sfx)
    W_i = num(c1, "W initial (RT)", "wi", gv("tga_W_initial") or 100.0)
    W200 = num(c2, "W at 200 °C (M-200)", "w200", gv("tga_W_200"))
    Wpl = num(c3, "W at plateau (M-400)", "wpl", gv("tga_W_plateau"))
    Wend = num(c4, f"W residue ({_res_lab(metal, ND)})", "wend", gv("tga_W_end"))
    st.subheader("4 · qNMR with internal standard (De Roo protocol)")
    c1, c2, c3, c4 = st.columns(4)
    std_c = num(c1, "[TMSP-d4] (mM)", "stdc", g("nmr_std_conc_mM"))
    std_I = num(c2, "TMSP integral", "stdi", g("nmr_std_integral") or 9.0)
    mdig = num(c3, "MOF mass digested (mg)", "mdig", g("nmr_sample_mg"), format="%.4f")
    vol = num(c4, "NMR volume (mL)", "vol", g("nmr_volume_mL"))
    st.subheader("5 · Chloride, ICP-OES, XRF")
    c1, c2, c3, c4 = st.columns(4)
    cl_node = num(c1, f"Cl per node (XRF/EDX: {ND['nM']}·Cl/M)", "cln", g("cl_per_node"))
    cl_mM = num(c2, "[Cl⁻] in digest (mM, UV-vis)", "clm", g("cl_conc_mM"))
    cl_wt = num(c3, "Cl wt%", "clw", g("cl_wt_pct"))
    s_m = num(c4, "XRF S/M atomic ratio (check)", "sm", g("xrf_S_per_M"))
    c1, c2, c3, c4 = st.columns(4)
    if metal in ("Zr", "Hf", "mixed"):
        zr = num(c1, "ICP Zr (wt%)", "zr", g("icp_Zr_wt"))
        hf = num(c2, "ICP Hf (wt%)", "hf", g("icp_Hf_wt")); icpM = None
    else:
        zr = hf = None
        icpM = num(c1, f"ICP {metal} (wt%)", "icpm_" + metal, g("icp_M_wt"))
        c2.caption(f"framework metal {metal}: {ND['nM']} per node")
    state = c3.selectbox("ICP sample weighed", ["as-measured (RT)", "dried 200 °C"], index=0 if "RT" in g("icp_state", "RT") else 1)
    lum = num(c4, "qNMR linker content (µmol mg⁻¹)", "lum", g("linker_umol_per_mg"), format="%.5f")
    st.subheader("6 · Porosity (semi-quantitative) and titration")
    c1, c2, c3, c4 = st.columns(4)
    Vm = num(c1, "V measured (cm³ g⁻¹)", "vm", g("V_meas"), format="%.4f")
    Vi = num(c2, "V ideal (cm³ g⁻¹)", "vi", g("V_ideal"), format="%.4f")
    vL = num(c3, "Volume freed per missing linker (Å³)", "vl", g("v_linker_A3"))
    if metal in ("Zr", "Hf") and fw in ("UiO-66", "MOF-801") and c4.button("Compute ideal V from bundled CIF"):
        from mfs.porosity import probe_pore_volume
        cif = os.path.join(REF, "UiO-66_Zr_COD4512074.cif" if fw.startswith("UiO") else "MOF-801_Zr_COD4121459.cif")
        r = probe_pore_volume(cif, substitute={"Zr": "Hf"} if metal == "Hf" else None)
        st.success(f"Probe-occupiable (N₂, 1.82 Å) pore volume: {r['V_cm3_g']:.3f} cm³ g⁻¹ (void fraction {r['void_fraction']:.2f})")
    c1, c2, c3, c4 = st.columns(4)
    tm = num(c1, "Titrated mass (mg)", "tm", g("tit_mass_mg"))
    tc = num(c2, "NaOH (M)", "tc", g("tit_c_NaOH") or 0.1)
    veq = c3.text_input("Equivalence volumes (mL, comma)", ", ".join(f"{v:.4f}" for v in (g("tit_V_eq") or [])))
    hcl = num(c4, "HCl added (mmol)", "hcl", g("tit_HCl_mmol"), format="%.4f")
    tup = st.file_uploader("Optional: titration curve (V, pH) → equivalence points from dpH/dV", key="tit")
    if tup:
        df = R.read_table(tup); eq = A.equivalence_points(df.iloc[:, 0].values, df.iloc[:, 1].values, 3)
        st.info("Equivalence volumes: " + ", ".join(f"{v:.3f}" for v in eq) + " mL – copy into the field above if correct.")
        fig, ax = newfig(0.6)
        ax.plot(df.iloc[:, 0], df.iloc[:, 1], "k-", lw=0.9); ax.set_xlabel("V(NaOH) (mL)"); ax.set_ylabel("pH")
        tw = ax.twinx(); tw.plot(df.iloc[:, 0], np.gradient(A.smooth(df.iloc[:, 1].values, "Savitzky-Golay", 15, 2), df.iloc[:, 0].values), color="#0072B2", lw=0.8)
        tw.set_ylabel("dpH/dV (mL$^{-1}$)", color="#0072B2"); ax.tick_params(right=False, which="both")
        for v in eq: ax.axvline(v, color="#D55E00", lw=0.6, ls="--")
        show(fig, "titration")
    reo = st.selectbox("Low-angle reo reflections (missing-cluster check) observed in PXRD?", ["not checked", "yes", "no"])
    try:
        species = [D.Species(**{k: (None if (isinstance(v, float) and np.isnan(v)) else v) for k, v in r.items()}) for r in sp.to_dict("records") if r.get("name")]
        for s in species:
            s.z = int(s.z); s.nH = int(s.nH or 0); s.at_plateau = bool(s.at_plateau); s.main = bool(s.main)
        inp = D.DefectInput(sample=sample, framework=fw, CN=int(CN), metal=metal, f_hf=f_hf, extraframework_Zr=extra, species=species,
                            node=nsel["node"] if nsel["node"] != (FRAMEWORKS[fw].get("node") or "") else "", node_custom=nsel["node_custom"], icp_M_wt=icpM,
                            tga_W_initial=W_i, tga_W_200=W200, tga_W_plateau=Wpl, tga_W_end=Wend, nmr_std_conc_mM=std_c, nmr_std_integral=std_I,
                            nmr_sample_mg=mdig, nmr_volume_mL=vol, cl_per_node=cl_node, cl_conc_mM=cl_mM, cl_wt_pct=cl_wt, icp_Zr_wt=zr, icp_Hf_wt=hf,
                            icp_state=state, linker_umol_per_mg=lum, xrf_S_per_M=s_m, V_meas=Vm, V_ideal=Vi, v_linker_A3=vL, tit_mass_mg=tm,
                            tit_c_NaOH=tc, tit_V_eq=[float(v) for v in veq.split(",") if v.strip()] or None, tit_HCl_mmol=hcl,
                            reo_peaks_seen=None if reo == "not checked" else reo == "yes")
        res = D.run_all(inp)
    except Exception as e:
        st.error(f"Calculation error: {e}"); return
    st.markdown("---")
    st.subheader("Results")
    cons = res["consensus"]
    if cons:
        c = st.columns(5)
        c[0].metric("Linkers / node (median)", f"{cons['L']:.2f}", f"range {cons['L_min']:.2f}–{cons['L_max']:.2f}", delta_color="off")
        c[1].metric("Missing linkers / node", f"{cons['x']:.2f}")
        c[2].metric("Caps / node", f"{cons['m']:.2f}")
        c[3].metric("Open sites / node", f"{cons['OS']:.2f}")
        c[4].metric("M per node (g mol⁻¹)", f"{cons['M']:.1f}")
        if not inp.is_zrhf() or ND["core"] != "{M}6O4(OH)4":
            st.caption(f"Node {ND['core'].format(M=inp.M_symbol())}: {inp.nM()} {inp.M_symbol()} per node, {inp.Q():g} charges for ligands, "
                       f"TGA residue {inp.residue_label()} ({inp.MMO2():.2f} g mol⁻¹ per metal). The equations below are written for "
                       "M₆O₄(OH)₄; the substitutions use these node values.")
        st.markdown(f"**Minimal formula** (one H₂O per terminal OH⁻/Cl⁻): `{cons['formula']}`  \n**Anhydrous**: `{cons['formula_anhydrous']}` "
                    f"({cons['M_anhydrous']:.1f} g mol⁻¹) · metal {cons['metal_wt_pct']:.1f} wt%")
        if s_m:
            st.caption(f"XRF cross-check: S per node = 6·S/M = {6*s_m:.2f} (compare with the sulfonated species per node in the formula).")
    tab = D.results_table(res)
    st.dataframe(tab, width="stretch", hide_index=True)
    if cons:
        used = [m.title.split(" ·")[0].split(" –")[0] for m in res["methods"] if m.status == "ok" and m.L is not None and not m.semi_quant]
        store_put(sample, "defects", {"framework": fw, "metal": metal, "node": ND["core"].format(M=inp.M_symbol()), "L": cons["L"], "x": cons["x"], "m": cons["m"], "OS": cons["OS"],
                                      "M": cons["M"], "formula": cons["formula"], "L_range": f"{cons['L_min']:.2f}–{cons['L_max']:.2f}",
                                      "methods": "; ".join(used), "metal_wt_pct": cons.get("metal_wt_pct")}, "Defect calculator")
    with st.expander("📑 SI table – one consistent defect report per sample (all samples in the store)", expanded=False):
        sit = SI.defect_rows(st.session_state.get("sample_store", {}))
        if len(sit):
            st.dataframe(SI._round(sit), width="stretch", hide_index=True)
            notes = ("Defect analysis per node (Figuryte " + __version__ + "). L = linkers per node (median of the methods listed), "
                     "x = missing linkers, caps = monocarboxylate modulators, open sites = 2x - caps. TGA in air; S_BET from the "
                     "Rouquerol-consistent range (ISO 9277); V_total at p/p0 = 0.95.")
            q1, q2 = st.columns(2)
            dl(q1, "⬇️ SI table (XLSX)", SI.to_xlsx(sit, tab.assign(sample=sample), notes), "SI_defects.xlsx", "si_xlsx")
            try:
                dl(q2, "⬇️ SI table (DOCX)", SI.to_docx(sit, notes=notes), "SI_defects.docx", "si_docx")
            except Exception as e:
                q2.caption(f"DOCX export unavailable ({e})")
        else:
            st.caption("No defect results stored yet – every computed sample is added automatically.")
    fig = P.fig_defect_comparison(res, max(WIDTH, 7.0))
    show(fig, f"defects_{re.sub(r'[^A-Za-z0-9]+', '_', sample)}")
    c1, c2, c3 = st.columns(3)
    sn = safe_name(sample)

    def _pdf(res=res):
        import tempfile
        with _RENDER_LOCK, tempfile.TemporaryDirectory() as td:
            pth = os.path.join(td, "defect_report.pdf")
            REP.pdf_report(res, pth)
            return open(pth, "rb").read()
    dl(c1, "📄 PDF report (full calculation)", _pdf, f"defect_report_{sn}.pdf", "def_pdf", help="generated when you click (takes a few s)")
    dl(c2, "CSV (table + steps)", REP.csv_report(res).encode("utf-8-sig"), f"defect_report_{sn}.csv", "def_csv")
    dl(c3, "JSON (inputs, reload later)", json.dumps(REP.input_to_dict(inp), indent=1).encode(), f"defect_input_{sn}.json", "def_json")
    if st.button("💾 Save report to folder", key="def_save"):
        os.makedirs(EXPORT_DIR, exist_ok=True)
        pth = os.path.join(EXPORT_DIR, f"defect_report_{sn}.pdf")
        open(pth, "wb").write(_pdf())
        open(os.path.join(EXPORT_DIR, f"defect_report_{sn}.csv"), "wb").write(REP.csv_report(res).encode("utf-8-sig"))
        st.success(f"Saved to `{EXPORT_DIR}`")
    st.subheader("Method by method")
    h = res["h"]
    for m in res["methods"]:
        icon = {"ok": "✅", "warn": "⚠️", "n/a": "▫️"}[m.status]
        with st.expander(f"{icon} {m.title}" + (f" — L = {m.L:.2f}, x = {m.x:.2f}" if m.L is not None else f" — {m.message}")):
            if m.semi_quant:
                st.warning("Semi-quantitative – excluded from the consensus.")
            if m.message:
                st.caption(m.message)
            st.markdown("**Equations**")
            for e in m.equations:
                st.latex(e)
            st.markdown("**Symbols**")
            st.table(pd.DataFrame(list(m.symbols.items()), columns=["symbol", "meaning"]))
            st.markdown("**Assumptions**\n" + "\n".join(f"- {a}" for a in m.assumptions))
            if m.steps:
                st.markdown("**Your numbers, step by step**")
                for s_ in m.steps:
                    st.latex(s_)
            if m.formula:
                st.markdown(f"**Formula**: `{m.formula}` (M = {m.M_formula:.1f} g mol⁻¹)")
            st.markdown("**Pitfalls**\n" + "\n".join(f"- {p}" for p in m.pitfalls))
            st.markdown("**References**\n" + "\n".join(f"- {c} [doi:{d}](https://doi.org/{d})" for c, d in m.refs))
    with st.expander(h["title"], expanded=True):
        for e in h["equations"]:
            st.latex(e)
        st.markdown("\n".join(f"- {a}" for a in h["assumptions"]))
        st.markdown("References: " + "; ".join(f"[{c.split('.')[0]}](https://doi.org/{d})" for c, d in h["refs"]))
    with LV.gate("Developer", globals()):
        with st.expander("When and why the methods disagree"):
            st.table(pd.DataFrame(D.SUMMARY_PITFALLS, columns=["pattern", "likely cause / what to do"]))
