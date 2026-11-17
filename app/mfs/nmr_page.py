"""¹H NMR page (v2.13): characterisation (structure ↔ spectrum) and digestion/composition mode.
Rendered from app.py with the app's helpers passed in (show, newfig, files_input, store_put, dl …)."""
from __future__ import annotations

import io
import json
import os
import zipfile

import numpy as np
import pandas as pd
import streamlit as st

from mfs import nmrproc as NP, nmrpeaks as PK

SOLVENT_CHOICES = ["auto (from file)", "CDCl3", "DMSO-d6", "Acetone-d6", "CD3OD", "D2O", "C6D6", "CD3CN", "CD2Cl2", "THF-d8",
                   "Toluene-d8"]
# qNMR internal standards: (name, M g/mol, nH of the integrated signal, typical δ, note)
QNMR_STD = [("Maleic acid", 116.07, 2, "6.0–6.4 (s)", "D2O, DMSO-d6, CD3OD"),
            ("1,3,5-Trimethoxybenzene (ArH)", 168.19, 3, "6.1 (s)", "CDCl3, DMSO-d6"),
            ("1,3,5-Trimethoxybenzene (OCH3)", 168.19, 9, "3.8 (s)", "CDCl3, DMSO-d6"),
            ("Dimethyl sulfone", 94.13, 6, "3.0–3.1 (s)", "D2O, DMSO-d6, CDCl3"),
            ("1,4-Bis(trimethylsilyl)benzene (SiCH3)", 222.47, 18, "0.25 (s)", "CDCl3, DMSO-d6"),
            ("Dimethyl terephthalate (ArH)", 194.18, 4, "8.1 (s)", "CDCl3, DMSO-d6"),
            ("1,2,4,5-Tetrachloro-3-nitrobenzene", 260.89, 1, "7.7–8.5 (s)", "CDCl3, DMSO-d6"),
            ("TMSP-d4 (sodium salt)", 172.27, 9, "0.00 (s)", "D2O"),
            ("DSS (sodium salt)", 218.32, 9, "0.00 (s)", "D2O"),
            ("Benzoic acid (o-H)", 122.12, 2, "8.0 (d)", "CDCl3, DMSO-d6"),
            ("Potassium hydrogen phthalate", 204.22, 4, "7.5–7.7 (m)", "D2O"),
            ("custom", 100.0, 1, "", "")]


@st.cache_data(show_spinner=False, max_entries=24)
def _load(data: bytes, name: str):
    low = name.lower()
    if low.endswith((".jdx", ".dx", ".jcamp", ".jdx.gz")):
        return [], NP.load_jcamp(data, name)
    return NP.load(data, name), None


@st.cache_data(show_spinner=False, max_entries=24)
def _process(data: bytes, name: str, exp_name: str, lb, p0, p1, reference, use_topspin):
    exps, jsp = _load(data, name)
    if jsp is not None:
        sp = jsp
        if p0 is not None:
            sp = NP.rephase(sp, p0, p1 or 0.0)
        return sp
    e = next(x for x in exps if x.name == exp_name)
    return NP.process(e, lb=lb, p0=p0, p1=p1, reference=reference, use_topspin=use_topspin)


@st.cache_data(show_spinner=False, max_entries=24)
def _analyse(_sp_key, _sp, k, solvent):
    ms = PK.analyse(_sp, k=k, solvent=solvent)
    return ms


def _mol(text):
    from mfs import nmrpredict as HP
    try:
        return HP.mol_from_any(text)
    except Exception:
        return None


def _canon(mol):
    """Canonical SMILES (None for no / unreadable structure) - to tell a redrawn molecule from a new one."""
    if mol is None:
        return None
    try:
        from rdkit import Chem
        return Chem.MolToSmiles(Chem.RemoveHs(mol))
    except Exception:
        return None


def _example_structures(ex_dir):
    p = os.path.join(ex_dir, "NMR_1H", "structures.json")
    try:
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def render(g):
    st.header("¹H NMR")
    mode = st.radio("Mode", ["🧬 Characterisation (structure ↔ spectrum)", "⚗️ Digestion / composition (ratios, qNMR)"],
                    horizontal=True, key="nmr_mode", label_visibility="collapsed")
    if mode.startswith("⚗️"):
        return g["page_nmr_digest"]()
    _characterisation(g)


def _characterisation(g):
    show, files_input, store_put, dl = g["show"], g["files_input"], g["store_put"], g["dl"]
    ss = st.session_state
    st.caption("Bruker raw data: the experiment folder (acqus + fid, optional pdata/1) as .zip — Avance II/III (XWIN-NMR, "
               "TopSpin 1–3) and Avance NEO (TopSpin 4) — processed automatically (phase, baseline, reference). "
               "JCAMP-DX spectra too. Add a .mol/.sdf with the same name for the structure, or draw it below.")
    files = files_input("Spectra (+ optional .mol)", ["zip", "jdx", "dx", "mol", "sdf"], "NMR_1H/*.zip", "nmr1h")
    if not files:
        return
    mols = {os.path.splitext(f.name)[0]: f.getvalue().decode("utf-8", "ignore") for f in files
            if f.name.lower().endswith((".mol", ".sdf"))}
    specs = [f for f in files if not f.name.lower().endswith((".mol", ".sdf"))]
    exst = _example_structures(g["EX"])
    entries = []
    for f in specs:
        try:
            exps, jsp = _load(f.getvalue(), f.name)
        except Exception as e:
            st.error(f"{f.name}: {e}"); continue
        if jsp is not None:
            entries.append((f, jsp.name or f.name, None, "1H" if jsp.nucleus == "1H" else jsp.nucleus))
            continue
        for e in exps:
            entries.append((f, e.name, e, e.nucleus))
    h1 = [x for x in entries if x[3] == "1H"]
    other = [x for x in entries if x[3] != "1H"]
    if other:
        st.caption(f"ℹ️ {len(other)} non-¹H experiment(s) in the upload ({', '.join(sorted({x[3] for x in other}))}) are "
                   "skipped: this version handles ¹H only.")
    if not h1:
        st.warning("No ¹H experiment found."); return
    labels = [f"{x[0].name} · {x[1]}" if x[2] is not None else x[0].name for x in h1]
    pick = st.selectbox("Spectrum", range(len(h1)), format_func=lambda i: labels[i], key="nmr1h_pick")
    f, ename, exp, _ = h1[pick]
    skey = f"{f.name}|{ename}"
    info = exp.info() if exp is not None else {}
    # ---------------- processing -----------------------------------------------------------------------------------------
    with st.expander("⚙️ Processing (automatic; override here)", expanded=False):
        c1, c2, c3, c4 = st.columns(4)
        lb = c1.number_input("Line broadening LB (Hz)", 0.0, 5.0, 0.3, 0.1, key="nmr1h_lb")
        ref_mode = c2.selectbox("Referencing", ["auto (TMS/DSS or solvent residual)", "as recorded (TopSpin SR)", "manual shift"],
                                key="nmr1h_ref")
        ref_shift = c2.number_input("Manual shift (ppm, added)", -2.0, 2.0, 0.0, 0.001, format="%.3f", key="nmr1h_refv",
                                    disabled=not ref_mode.startswith("manual"))
        man = c3.checkbox("Manual phase", key="nmr1h_manph")
        p0 = c3.slider("p0 (°)", -180.0, 180.0, 0.0, 0.5, key="nmr1h_p0", disabled=not man)
        p1 = c3.slider("p1 (°)", -360.0, 360.0, 0.0, 1.0, key="nmr1h_p1", disabled=not man)
        use_ts = c4.checkbox("Use TopSpin's processed 1r instead", key="nmr1h_ts",
                             disabled=exp is None or exp.proc_1r is None)
        sol_pick = c4.selectbox("Solvent", SOLVENT_CHOICES, key="nmr1h_sol")
        kpk = c4.slider("Peak threshold (× noise)", 3.0, 30.0, 6.0, 0.5, key="nmr1h_k")
    reference = "auto" if ref_mode.startswith("auto") else ("none" if ref_mode.startswith("as") else float(ref_shift))
    try:
        sp = _process(f.getvalue(), f.name, ename, lb if exp is not None else None, p0 if man else None,
                      p1 if man else None, reference, bool(use_ts))
    except Exception as e:
        st.error(f"Processing failed: {e}"); return
    solvent = sp.info.get("solvent") if sol_pick.startswith("auto") else sol_pick
    if not sol_pick.startswith("auto"):
        sp.info["solvent"] = sol_pick
    mhz = float(sp.info.get("spectrometer (MHz)") or sp.sf or 400)
    st.caption(f"**{ename}** · {mhz:.0f} MHz · {solvent or 'solvent ?'} · "
               + " · ".join(f"{k}: {v}" for k, v in info.items() if k in ("instrument", "software", "pulse program", "scans", "date"))
               + f" · phase p0 {sp.p0:.1f}°, p1 {sp.p1:.1f}° · {sp.ref_note}")
    mults = [PK.Multiplet(**m.as_dict()) for m in _analyse(f"{skey}|{lb}|{p0}|{p1}|{reference}|{use_ts}|{solvent}", sp, kpk, solvent)]
    # ---------------- structure ------------------------------------------------------------------------------------------
    mk = "nmr1h_mol|" + skey
    if mk not in ss:
        stem = os.path.splitext(f.name)[0]
        ss[mk] = mols.get(stem) or exst.get(f.name, "")
    g["PS"].track(mk)                    # after the default is in place: opening the page is not an "unsaved change"
    st.subheader("1 · Structure")
    t1, t2 = st.tabs(["✏️ Draw (editor)", "⌨️ SMILES / molfile"])
    with t1:
        try:
            from streamlit_ketcher import st_ketcher
            drawn = st_ketcher(ss[mk], height=420, molecule_format="MOLFILE", key="ketch_" + skey)
            if drawn and drawn.strip() and drawn != ss[mk]:
                m_ = _mol(drawn)
                # the editor echoes the structure back as a molfile: only a different molecule is a change
                if m_ is not None and m_.GetNumAtoms() and _canon(m_) != _canon(_mol(ss[mk]) if ss[mk] else None):
                    ss[mk] = drawn
        except Exception as e:
            st.info(f"Structure editor unavailable ({e}); paste SMILES or a molfile in the other tab.")
    with t2:
        txt = st.text_area("SMILES or molfile", ss[mk], height=110, key="nmr1h_txt|" + skey)
        if st.button("Use this structure", key="nmr1h_use|" + skey):
            if _mol(txt) is None:
                st.error("Could not read that structure.")
            else:
                ss[mk] = txt; st.rerun()
    mol = _mol(ss[mk]) if ss[mk] else None
    if mol is not None:
        from rdkit.Chem import rdMolDescriptors as RD, Descriptors
        st.caption(f"{RD.CalcMolFormula(mol)} · M = {Descriptors.MolWt(mol):.2f} g mol⁻¹ · "
                   f"{sum(a.GetTotalNumHs() for a in mol.GetAtoms())} H")
    # ---------------- integrals ------------------------------------------------------------------------------------------
    st.subheader("2 · Multiplets and integrals")
    c1, c2 = st.columns([2, 3])
    nh_tot = sum(a.GetTotalNumHs() for a in mol.GetAtoms()) if mol is not None else None
    norm_opts = ["best integer fit"] + (["total H of the structure (non-exchangeable)"] if mol is not None else []) + \
        ["reference multiplet"]
    norm = c1.radio("Normalise integrals by", norm_opts, key="nmr1h_norm")
    refi = None
    cand = [i for i, m in enumerate(mults) if not m.note]
    if norm == "reference multiplet" and cand:
        refi = c2.selectbox("Reference multiplet", cand, format_func=lambda i: f"δ {mults[i].centre:.2f} ({mults[i].mult})",
                            key="nmr1h_refm")
        refn = c2.number_input("… contains nH", 0.5, 60.0, 1.0, 0.5, key="nmr1h_refn")
        PK.normalise(mults, ref_index=refi, ref_n=refn)
    elif norm.startswith("total") and mol is not None:
        from mfs import nmrpredict as HP
        n_ne = sum(a.GetTotalNumHs() for a in mol.GetAtoms() if a.GetSymbol() not in HP.EXCHANGEABLE)
        PK.normalise(mults)                                   # first pass: find the compound's multiplets
        main = [m for m in mults if not m.note and (m.nH or 0) >= 0.5]
        if main:
            unit = sum(m.integral for m in main) / max(n_ne, 1)
            for m in mults:
                m.nH = round(m.integral / unit, 2)
    else:
        PK.normalise(mults)
    # ---------------- assignment -----------------------------------------------------------------------------------------
    assignment, groups, imp = {}, [], {}
    fk = "nmr1h_fix|" + skey
    g["PS"].track(fk); g["PS"].track_undo(fk)
    fixed = {int(k): v for k, v in (ss.get(fk) or {}).items()}
    if mol is not None:
        from mfs import nmrassign as A, nmrui as UI
        groups = A.groups(mol)
        fixed = {k: v for k, v in fixed.items() if k < len(groups) and (v is None or v < len(mults))}
        assignment = A.assign(groups, mults, fixed=fixed)
        st.subheader("3 · Assignment (drag to correct)")
        ev = UI.widget(sp, mults, mol, groups, assignment, fixed, skey)
        if ev and ev.get("t") != ss.get("_nmra_last|" + skey):
            ss["_nmra_last|" + skey] = ev.get("t")
            ss[fk] = {str(k): v for k, v in UI.handle(fixed, ev).items()}
            st.rerun()
        rows = []
        for gr in groups:
            mi = assignment.get(gr["gid"])
            m = mults[mi] if mi is not None else None
            rows.append({"H on atom": ", ".join(map(str, gr["labels"])), "nH": gr["nH"],
                         "predicted δ": None if gr["shift"] is None else round(gr["shift"], 2), "confidence": gr["conf"],
                         "assigned δ": None if m is None else round(m.centre, 2), "multiplet": "" if m is None else m.mult,
                         "Δδ": None if (m is None or gr["shift"] is None) else round(m.centre - gr["shift"], 2),
                         "set by": "you" if gr["gid"] in fixed else "auto"})
        with st.expander("Assignment table", expanded=False):
            st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
            st.caption("Prediction: HOSE-code environments from experimental ¹H data in nmrshiftdb2 (held-out test: mean error "
                       "0.22 ppm, 79 % within 0.3 ppm). Exchangeable OH/NH are often not observed — they may stay unassigned.")
    # ---------------- impurities -----------------------------------------------------------------------------------------
    from mfs import nmrassign as A
    un = A.unassigned(mults, assignment, 0.0) if mol is not None else [i for i, m in enumerate(mults) if not m.note
                                                                         and (m.nH or 0) < 0.5]
    unit = 1.0
    res = A.impurity_check(mults, un, solvent or "CDCl3", unit_nH=unit,
                           mw_compound=(__import__("rdkit.Chem.Descriptors", fromlist=["MolWt"]).MolWt(mol) if mol is not None else None))
    imp = {i: b["name"] for b in res["best"] for i in b["peaks"]}
    st.subheader("4 · Unassigned peaks and impurity check")
    if not un:
        st.success("Every multiplet is assigned to the structure (solvent/TMS lines excluded).")
    else:
        rows = []
        for i in un:
            m = mults[i]; c = res["per_peak"].get(i) or []
            best = next((b for b in res["best"] if i in b["peaks"]), None)
            rows.append({"δ": round(m.centre, 3), "mult": m.mult, "nH (rel.)": m.nH,
                         "best guess": best["name"] if best else "unknown",
                         "evidence": (f"{best['matched']} signals" + (f", {best['mol_ratio']:.3f} mol/mol" if best and best["mol_ratio"] else "")
                                      + (f", {best['wt_pct']:.2f} wt %" if best and best.get("wt_pct") else "")) if best else "",
                         "other candidates": "; ".join(f"{x['name']} ({x['proton']}, Δ {x['delta']:+.2f})" for x in c[:3])})
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        st.caption("Table: Fulmer et al., Organometallics 2010, 29, 2176 (Table 1) + DMSO from Gottlieb et al., J. Org. Chem. "
                   "1997, 62, 7512. Multi-signal impurities count only if their other signals are present with matching "
                   "integral ratios. Shifts move with concentration/temperature — treat as a best guess.")
    # ---------------- ACS string ----------------------------------------------------------------------------------------
    st.subheader("5 · Characterisation string")
    c1, c2 = st.columns([1, 4])
    with_lab = c1.checkbox("Add assignments", value=False, key="nmr1h_acslab")
    acs = PK.acs_string(mults, mhz=mhz, solvent=solvent or "", labels=with_lab, skip=set(imp) | set(un if mol is not None else []),
                        min_nH=0.5)
    c2.code(acs, language=None, wrap_lines=True)
    # ---------------- figure ---------------------------------------------------------------------------------------------
    st.subheader("6 · Labelled figure")
    c1, c2, c3, c4 = st.columns(4)
    s_str = c1.checkbox("Structure inset", True, key="nmr1h_fstr", disabled=mol is None)
    s_int = c2.checkbox("Integrals", True, key="nmr1h_fint")
    s_lab = c3.checkbox("Assignment labels", True, key="nmr1h_flab")
    from mfs import nmrfig as F
    h0, l0 = F.auto_xlim(sp, mults)
    rng = c4.slider("δ window", float(np.floor(sp.ppm.min())), float(np.ceil(sp.ppm.max())), (round(l0, 1), round(h0, 1)), 0.1,
                    key="nmr1h_xr|" + skey)

    def draw(ax, sp=sp, mults=mults, mol=mol, groups=groups, imp=imp, rng=rng):
        F.draw(ax, sp, mults, mol if s_str else None, groups, imp, xlim=(rng[1], rng[0]), show_integrals=s_int,
               show_labels=s_lab)
    fig, ax = g["newfig"](0.45); draw(ax)
    df = NP.to_frame(sp)
    show(fig, "NMR_1H_" + g["safe_name"](ename), df, draw, "¹H NMR (labelled)")
    # ---------------- qNMR -----------------------------------------------------------------------------------------------
    pur = _qnmr(mults, mol, un, imp, skey, g)
    # ---------------- SI export ------------------------------------------------------------------------------------------
    st.subheader("8 · SI export")
    mt = pd.DataFrame([{"delta_ppm": round(m.centre, 3), "from_ppm": round(m.hi, 3), "to_ppm": round(m.lo, 3), "mult": m.mult,
                        "J_Hz": "; ".join(f"{j:.1f}" for j in m.J), "nH": m.nH, "integral_raw": m.integral,
                        "assignment": m.label, "note": m.note or ("impurity: " + imp[i] if i in imp else
                                                                    ("unassigned" if i in un else ""))}
                       for i, m in enumerate(mults)])
    st.dataframe(mt, hide_index=True, width="stretch", height=220)
    sample = os.path.splitext(f.name)[0]
    store_put(sample, "nmr_1h", {"acs": acs, "solvent": solvent, "MHz": round(mhz), "purity_pct": pur,
                                 "smiles": _smiles(mol), "impurities": ", ".join(sorted(set(imp.values())))},
              "¹H NMR characterisation")

    def _zip():
        import matplotlib.pyplot as plt
        b = io.BytesIO()
        with zipfile.ZipFile(b, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("characterisation.txt", acs + "\n")
            z.writestr("multiplets.csv", mt.to_csv(index=False))
            z.writestr("spectrum_processed.csv", df.to_csv(index=False))
            z.writestr("processing.txt", f"{ename}\n{sp.ref_note}\np0 {sp.p0:.2f} deg, p1 {sp.p1:.2f} deg\n"
                       + "\n".join(f"{k}: {v}" for k, v in sp.info.items()) + "\n")
            if mol is not None:
                from rdkit import Chem
                z.writestr("structure.mol", Chem.MolToMolBlock(mol))
                z.writestr("assignment.csv", pd.DataFrame(rows).to_csv(index=False) if groups else "")
            fig2, ax2 = g["newfig"](0.45); draw(ax2)
            for fmt in ("png", "pdf", "svg"):
                z.writestr(f"spectrum_labelled.{fmt}", g["export_bytes"](fig2, fmt))
            plt.close(fig2)
            z.writestr("SI_page.pdf", si_pdf(sp, mults, mol, groups, imp, acs, ename, solvent, mhz, pur))
        return b.getvalue()
    c1, c2, c3 = st.columns(3)
    dl(c1, "⬇️ SI package (.zip)", _zip, f"NMR_SI_{g['safe_name'](ename)}.zip", "nmr1h_zip|" + skey)
    dl(c2, "⬇️ SI page (PDF)", lambda: si_pdf(sp, mults, mol, groups, imp, acs, ename, solvent, mhz, pur),
       f"NMR_SI_{g['safe_name'](ename)}.pdf", "nmr1h_pdf|" + skey)
    dl(c3, "⬇️ Multiplets (CSV)", mt.to_csv(index=False).encode("utf-8-sig"), f"NMR_multiplets_{g['safe_name'](ename)}.csv",
       "nmr1h_csv|" + skey)


def _smiles(mol):
    if mol is None:
        return ""
    from rdkit import Chem
    return Chem.MolToSmiles(mol)


def qnmr_purity(I_a, N_a, M_a, m_a, I_s, N_s, M_s, m_s, P_s):
    """Internal-standard qNMR: P_a = (I_a/I_s)(N_s/N_a)(M_a/M_s)(m_s/m_a)·P_s (masses in the same unit, P in %)."""
    return (I_a / I_s) * (N_s / N_a) * (M_a / M_s) * (m_s / m_a) * P_s


def _broad(m):
    return m.mult == "br s"


def integral_consistency(mults, idx):
    """Coefficient of variation (%) of integral per H over the compound multiplets idx (nH rounded; broad exchangeable
    singlets and lines on the solvent residual are left out)."""
    v = [mults[i].integral / max(round(mults[i].nH or 0), 1) for i in idx
         if (mults[i].nH or 0) >= 0.5 and not _broad(mults[i]) and not mults[i].note]
    if len(v) < 2:
        return None
    return float(np.std(v) / np.mean(v) * 100)


def _qnmr(mults, mol, un, imp, skey, g):
    st.subheader("7 · qNMR purity")
    comp = [i for i, m in enumerate(mults) if m.label] if mol is not None else \
        [i for i, m in enumerate(mults) if not m.note and (m.nH or 0) >= 0.5]
    cv = integral_consistency(mults, comp)
    if cv is not None:
        (st.warning if cv > 5 else st.caption)(f"Integral consistency over the compound multiplets: CV {cv:.1f} % per H"
                                               + (" — check phase/baseline, overlap or relaxation delay (D1 ≥ 5·T1)." if cv > 5 else " ✓"))
    with st.expander("Internal-standard purity (qNMR)", expanded=False):
        c1, c2, c3 = st.columns(3)
        names = [s[0] for s in QNMR_STD]
        si = c1.selectbox("Internal standard", range(len(names)), format_func=lambda i: names[i], key="qn_std")
        M_s = c1.number_input("M standard (g/mol)", 1.0, 5000.0, float(QNMR_STD[si][1]), key=f"qn_ms{si}")
        N_s = c1.number_input("nH of the standard signal", 1, 60, int(QNMR_STD[si][2]), key=f"qn_ns{si}")
        st.caption(f"Typical δ {QNMR_STD[si][3]} · solvents {QNMR_STD[si][4]}")
        m_s = c2.number_input("Mass standard (mg)", 0.0, 1000.0, 5.0, 0.01, format="%.3f", key="qn_mst")
        P_s = c2.number_input("Purity standard (%)", 0.0, 100.0, 99.9, 0.1, key="qn_pst")
        m_a = c3.number_input("Mass analyte (mg)", 0.0, 1000.0, 10.0, 0.01, format="%.3f", key="qn_ma")
        M_a0 = None
        if mol is not None:
            from rdkit.Chem import Descriptors
            M_a0 = Descriptors.MolWt(mol)
        M_a = c3.number_input("M analyte (g/mol)", 1.0, 5000.0, float(M_a0 or 300.0), key="qn_Ma")
        opts = list(range(len(mults)))
        fmt = lambda i: f"δ {mults[i].centre:.2f} {mults[i].mult} ({mults[i].nH}H){' · ' + mults[i].label if mults[i].label else ''}"
        sidx = st.selectbox("Standard signal", opts, format_func=fmt, key="qn_si",
                            index=next((i for i in un if i not in imp), 0) if un else 0)
        aopts = [i for i in opts if i != sidx]
        aidx = st.multiselect("Analyte signal(s)", aopts, default=[i for i in comp if i in aopts and not _broad(mults[i])][:3],
                              format_func=fmt, key="qn_ai")
        if aidx and mults[sidx].integral > 0 and m_a > 0:
            vals = []
            for i in aidx:
                n_a = max(round(mults[i].nH or 1), 1)
                vals.append(qnmr_purity(mults[i].integral, n_a, M_a, m_a, mults[sidx].integral, N_s, M_s, m_s, P_s))
            P = float(np.mean(vals))
            st.metric("Purity (wt %)", f"{P:.1f} %", help="Mean over the selected analyte signals")
            if len(vals) > 1:
                st.caption("Per signal: " + ", ".join(f"{v:.1f} %" for v in vals) + f" · SD {np.std(vals, ddof=1):.1f} %")
            st.caption("P = (Iₐ/Iₛ)·(Nₛ/Nₐ)·(Mₐ/Mₛ)·(mₛ/mₐ)·Pₛ. Requires full relaxation (D1 ≥ 5·T1, 90° pulse) and a flat "
                       "baseline around both signals.")
            return round(P, 2)
    return None


def si_pdf(sp, mults, mol, groups, imp, acs, name, solvent, mhz, purity=None):
    """One SI page (A4): labelled spectrum, structure, characterisation string, multiplet table."""
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from mfs import nmrfig as F
    import textwrap
    b = io.BytesIO()
    with PdfPages(b) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.07, 0.96, f"¹H NMR — {name}", fontsize=13, weight="bold")
        fig.text(0.07, 0.94, f"{mhz:.0f} MHz, {solvent or ''}, 298 K; processed in Reticula ({sp.ref_note})", fontsize=8,
                 color="#46607F")
        ax = fig.add_axes([0.07, 0.6, 0.86, 0.32])
        F.draw(ax, sp, mults, mol, groups, imp)
        fig.text(0.07, 0.55, "\n".join(textwrap.wrap(acs, 120)), fontsize=8, va="top", family="DejaVu Sans")
        if purity is not None:
            fig.text(0.07, 0.48, f"qNMR purity: {purity:.1f} wt %", fontsize=8)
        rows = [[f"{m.centre:.2f}", m.mult, ", ".join(f"{j:.1f}" for j in m.J), "" if m.nH is None else f"{m.nH:.2f}", m.label,
                 m.note or (f"impurity: {imp[i]}" if i in imp else "")] for i, m in enumerate(mults) if (m.nH or 0) >= 0.02 or m.label]
        rows = rows[:40]
        if rows:
            axt = fig.add_axes([0.07, 0.05, 0.86, 0.4]); axt.axis("off")
            t = axt.table(cellText=rows, colLabels=["δ (ppm)", "mult.", "J (Hz)", "nH", "assignment", "note"], loc="upper center",
                          cellLoc="left", colLoc="left")
            t.auto_set_font_size(False); t.set_fontsize(7); t.scale(1, 1.15)
        pdf.savefig(fig); plt.close(fig)
    return b.getvalue()
