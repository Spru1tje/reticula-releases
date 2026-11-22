"""Reticula page module: 1H NMR digestion page (the new NMR page is mfs/nmr_page.py).

Moved verbatim out of app.py in v2.13. Names that the code takes from the app (helpers such as show(), dl(), files_input(),
the sidebar values WIDTH / ASPECT / FONT, the cached reader wrapper R …) are listed in INJECTED and filled in by
mfs.pages.bind(globals()) on every run of app.py, so the code behaves exactly as it did inside app.py."""
# ruff: noqa: F821
from __future__ import annotations
import numpy as np
import pandas as pd
import streamlit as st
from mfs import numentry as NE  # every slider takes a typed number too (2.19.6)
from mfs import levels as LV

from mfs import analysis as A, plots as P

INJECTED = (
    'R',
    'dl',
    'files_input',
    'newfig',
    'offset_controls',
    'sample_table',
    'show',
    'store_put',
)


def _read_nmr_any(f):
    """v2.13: TopSpin processed 1r (as before) or, for raw folders without pdata, automatic processing of the FID."""
    raw = f.getvalue()
    try:
        return R.read_nmr(f, f.name)
    except Exception as e0:
        if raw[:2] != b"PK":
            raise
        from mfs import nmrproc as NP
        exps = [x for x in NP.load(raw, f.name) if x.nucleus == "1H" and x.fid is not None]
        if not exps:
            raise e0
        sp = NP.process(exps[0])
        out = pd.DataFrame({"x": sp.ppm, "y": sp.y})
        out.attrs["source"] = "FID processed automatically (" + sp.ref_note + ")"
        return out


_DIGEST_PRESETS = {      # (species, from, to, nH) — typical windows after digestion (NaOD/D2O or D2SO4/DMSO-d6); edit to your data
    "MOF-801 / Zr-fumarate (NaOD/D2O)": [("TMSP-d4 (standard)", -0.15, 0.15, 9), ("Fumarate", 6.40, 6.65, 2), ("Formate", 8.35, 8.55, 1),
                                         ("2-Sulfosuccinate (dicarboxylate, SO3H kept)", 2.40, 2.90, 2), ("Acetate", 1.80, 2.00, 3)],
    "UiO-66 / terephthalate (NaOD/D2O)": [("TMSP-d4 (standard)", -0.15, 0.15, 9), ("Terephthalate (BDC)", 7.70, 7.95, 4),
                                          ("Formate", 8.35, 8.55, 1), ("Acetate", 1.80, 2.00, 3), ("Benzoate (o-H)", 7.80, 7.95, 2)],
    "UiO-66-NH2 / 2-aminoterephthalate (NaOD/D2O)": [("TMSP-d4 (standard)", -0.15, 0.15, 9), ("2-Aminoterephthalate", 6.95, 7.75, 3),
                                                     ("Formate", 8.35, 8.55, 1), ("Acetate", 1.80, 2.00, 3)],
    "MOF-808 / trimesate (NaOD/D2O)": [("TMSP-d4 (standard)", -0.15, 0.15, 9), ("Trimesate (BTC)", 8.20, 8.45, 3),
                                       ("Formate", 8.45, 8.60, 1), ("Acetate", 1.80, 2.00, 3)],
    "Acid digestion in DMSO-d6 (D2SO4/HF)": [("Dimethyl sulfone (standard)", 2.95, 3.05, 6), ("Main linker", 7.90, 8.10, 4),
                                             ("Formic acid", 8.05, 8.20, 1), ("Acetic acid", 1.85, 1.95, 3), ("DMF (CHO)", 7.90, 7.99, 1)],
    "Custom (empty)": [("Standard", -0.15, 0.15, 9), ("Component A", 7.0, 7.5, 1)],
}


def page_nmr_digest():
    st.subheader("Digestion / composition – species ratios and qNMR against an internal standard")
    st.caption("Digested framework (MOF, COF, coordination polymer …) or any mixture: integrate each species, ratios to the main "
               "component, internal standard. Bruker experiment folder as .zip (processed pdata/1/1r, or raw FID → automatic "
               "processing), TopSpin *totxt* export, or ppm/intensity text.")
    files = files_input("Spectra", ["zip", "txt", "csv", "dat"], "NMR/*_totxt.txt", "nmr")
    if not files:
        return
    samples = sample_table(files, "nmr_tab")
    c1, c2, c3 = st.columns(3)
    lo, hi = NE.slider("δ range (ppm)", -1.0, 14.0, (1.5, 9.0), 0.1)
    cut = c2.text_input("Cut regions (e.g. solvent HDO 4.55-5.05)", "4.55-5.05")
    off = offset_controls(c3, "nmr_off", 1.1, 0.2, 2.0, 0.05, "Offset")
    st.session_state["_even|nmr_off"] = 1.05
    ann = c3.text_input("Annotations (ppm:label; …)", "8.44:HCOO$^-$; 6.52:fumarate; 3.62:SS CH; 2.68:SS CH$_2$")
    cuts = []
    for p in cut.split(";"):
        try:
            a, b = [float(v) for v in p.split("-")]; cuts.append((a, b))
        except Exception:
            pass
    traces, specs = [], {}
    for nm, f in samples:
        try:
            d = _read_nmr_any(f)
        except Exception as e:
            st.error(f"{f.name}: {e}"); continue
        specs[nm] = d
        m = np.ones(len(d), bool)
        for a, b in cuts:
            m &= ~((d.x > a) & (d.x < b))
        y = np.where(m, d.y, np.nan)
        traces.append((nm, d.x.values, y))
    annotations = []
    for p in ann.split(";"):
        if ":" in p:
            try:
                a, b = p.split(":", 1); annotations.append((float(a), b.strip()))
            except Exception:
                pass
    def draw(ax, traces=traces):
        P.draw_stacked(ax, traces, "nmr", offset=off, xlim=(hi, lo), annotations=annotations)
    fig, ax = newfig(0.9); draw(ax)
    show(fig, "NMR_digest", None, draw, "1H NMR")
    st.subheader("Integration → component ratios (feeds the defect calculator)")
    st.caption("Integrals are baseline-corrected (linear between region ends). Pick a preset (typical windows after digestion; "
               "edit them to your spectra), the internal standard row and the main component (= 1 in the ratios).")
    presets = _DIGEST_PRESETS
    pc1, pc2, pc3 = st.columns([2, 2, 2])
    pre = pc1.selectbox("Preset", list(presets), key="nmr_preset")
    if st.session_state.get("_nmr_preset_last") != pre:          # new preset → new editable table
        st.session_state["_nmr_preset_last"] = pre
        st.session_state.pop("nmr_regions", None)
    regions = pd.DataFrame(presets[pre], columns=["species", "from_ppm", "to_ppm", "nH"])
    regions = st.data_editor(regions, num_rows="dynamic", key="nmr_regions")
    names = [str(x) for x in regions["species"].dropna()]
    if not names:
        return
    std_name = pc2.selectbox("Internal standard", ["(none)"] + names, index=1 if names[0].lower().startswith(("tmsp", "dss", "maleic", "dimethyl sulfone", "std")) else 0,
                             key="nmr_std_pick|" + pre)
    main_opts = [n for n in names if n != std_name]
    main_name = pc3.selectbox("Main component (= 1)", main_opts, key="nmr_main_pick|" + pre) if main_opts else None
    if specs:
        with LV.gate("Advanced", globals()):
            with st.expander("🔎 Peaks found in the first spectrum (to set the windows)", expanded=False):
                try:
                    from mfs import nmrpeaks as _PK, nmrproc as _NP
                    d0 = next(iter(specs.values()))
                    sp0 = _NP.Spectrum(d0.x.values, d0.y.values, "1H", "1H", float(d0.attrs.get("SF") or 400.0), {})
                    ms0 = _PK.analyse(sp0, k=8.0); _PK.normalise(ms0)
                    st.dataframe(pd.DataFrame([{"δ": round(m.centre, 3), "from": round(m.lo, 2), "to": round(m.hi, 2), "mult": m.mult,
                                                "rel. integral": round(m.integral / max(x.integral for x in ms0), 4)} for m in ms0
                                               if m.integral > 0.003 * max(x.integral for x in ms0)]), hide_index=True, width="stretch")
                except Exception as e:
                    st.caption(f"(peak list unavailable: {e})")
    std_nH = float(regions.loc[regions["species"] == std_name, "nH"].iloc[0]) if std_name in names else None
    out = []
    for nm, d in specs.items():
        row = {"sample": nm}
        for r in regions.itertuples():
            row[f"I {r.species}"] = A.integrate(d.x.values, d.y.values, r.from_ppm, r.to_ppm)
        std = row.get(f"I {std_name}") if std_nH else None
        std = std if std else np.nan
        tag = f"rel. to {std_nH:g}H std" if std_nH else "per H"
        for r in regions.itertuples():
            if r.species == std_name:
                continue
            if std == std and std_nH:
                row[f"I {r.species}"] = row[f"I {r.species}"] / std * std_nH
            row[f"n {r.species} ({tag})"] = row[f"I {r.species}"] / r.nH
        main = f"n {main_name} ({tag})" if main_name else None
        if main and main in row and row[main]:
            for r in regions.itertuples():
                k = f"n {r.species} ({tag})"
                if k in row and k != main:
                    row[f"ratio {r.species}:{main_name}"] = row[k] / row[main]
            parts = [(main_name.split(" (")[0], 1.0)] + [(r.species.split(" (")[0], row.get(f"ratio {r.species}:{main_name}"))
                                                         for r in regions.itertuples()
                                                         if r.species not in (std_name, main_name) and f"ratio {r.species}:{main_name}" in row]
            row["composition"] = " : ".join(f"{v:.2f}" for _, v in parts if v is not None) + "  (" + " : ".join(n_ for n_, v in parts if v is not None) + ")"
        out.append(row)
    if out:
        df = pd.DataFrame(out)
        st.dataframe(df, width="stretch")
        st.caption("Ratios per main component (= 1). The integrals go to the sample store and are pre-filled in the defect "
                   "calculator ('Pull data for sample') for methods (b), (c), (d) and (e).")
        st.session_state["nmr_integrals"] = df
        for row in out:
            store_put(row["sample"], "nmr", {k: v for k, v in row.items() if k != "sample"}, "1H NMR digest page")
        dl(st, "CSV integrals", df.to_csv(index=False).encode("utf-8-sig"), "NMR_integrals.csv", "nmr_int_csv")
