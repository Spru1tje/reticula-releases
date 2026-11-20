"""Automatic SI tables from a CIF: crystal data and structure refinement, selected bonds / angles / torsions (ESDs),
hydrogen bonds with symmetry codes. DOCX (plain WordprocessingML, no python-docx) and XLSX."""
from __future__ import annotations
import io, zipfile
import numpy as np
import pandas as pd
from .cif import num, fmt_esd
from .si_export import _run


def _g(S, *tags):
    for t in tags:
        v = S.items.get(t)
        if v not in (None, "?", "."):
            return v
    return None


def _e(S, tag, unit=""):
    v = _g(S, tag)
    if v is None:
        return "n.r."
    x, e = num(v)
    return (fmt_esd(x, e) if np.isfinite(x) else str(v)) + (f" {unit}" if unit else "")


def _plain(v, d=None):
    if v is None:
        return "n.r."
    x, _ = num(v)
    if not np.isfinite(x):
        return str(v)
    return f"{x:.{d}f}" if d is not None else (f"{x:g}")


def crystal_table(S):
    g = lambda *t: _g(S, *t)
    lam = g("_diffrn_radiation_wavelength")
    rad = g("_diffrn_radiation_type") or ""
    radn = {"CuK\\a": "Cu Kα", "MoK\\a": "Mo Kα", "CuKa": "Cu Kα", "MoKa": "Mo Kα"}.get(rad.replace(" ", ""), rad)
    tmin, tmax = g("_diffrn_reflns_theta_min"), g("_diffrn_reflns_theta_max")
    hkl = [(g(f"_diffrn_reflns_limit_{c}_min"), g(f"_diffrn_reflns_limit_{c}_max")) for c in "hkl"]
    rint = g("_diffrn_reflns_av_r_equivalents"); rsig = g("_diffrn_reflns_av_uneti/neti", "_diffrn_reflns_av_sigmai/neti")
    twin = "twin" in (g("_refine_special_details") or "").lower()
    rint_s = _plain(rint, 4) if rint else ("n.r." + (" (twinned data)" if twin else ""))
    nt, ngt = g("_reflns_number_total"), g("_reflns_number_gt")
    thf = g("_diffrn_reflns_theta_full")
    comp = g("_diffrn_reflns_laue_measured_fraction_full", "_diffrn_measured_fraction_theta_full")
    ccdc = g("_database_code_depnum_ccdc_archive", "_database_code_ccdc")
    cell = [fmt_esd(*c) if np.isfinite(c[1]) else f"{c[0]:g}" for c in S.cell]
    rows = [
        ("Identification code", S.name), ("Empirical formula", g("_chemical_formula_sum") or S.formula()),
        ("Formula weight", _plain(g("_chemical_formula_weight"), 2)), ("Temperature/K", _e(S, "_diffrn_ambient_temperature") if g("_diffrn_ambient_temperature") else _e(S, "_cell_measurement_temperature")),
        ("Crystal system", g("_space_group_crystal_system", "_symmetry_cell_setting") or "n.r."), ("Space group", S.sg),
        ("a/Å", cell[0]), ("b/Å", cell[1]), ("c/Å", cell[2]), ("α/°", cell[3]), ("β/°", cell[4]), ("γ/°", cell[5]),
        ("Volume/Å³", _e(S, "_cell_volume")), ("Z", _plain(g("_cell_formula_units_z"))),
        ("ρcalc/g cm⁻³", _plain(g("_exptl_crystal_density_diffrn"), 3) if g("_exptl_crystal_density_diffrn") else f"{S.density():.3f}"),
        ("μ/mm⁻¹", _plain(g("_exptl_absorpt_coefficient_mu"), 3)), ("F(000)", _plain(g("_exptl_crystal_f_000"))),
        ("Crystal size/mm³", " × ".join(_plain(g(f"_exptl_crystal_size_{d}"), 3) for d in ("max", "mid", "min")) if g("_exptl_crystal_size_max") else "n.r."),
        ("Crystal colour, shape", ", ".join(x for x in (g("_exptl_crystal_colour"), g("_exptl_crystal_description")) if x) or "n.r."),
        ("Radiation", f"{radn} (λ = {_plain(lam, 5)} Å)"),
        ("2Θ range for data collection/°", f"{2 * num(tmin)[0]:.2f} to {2 * num(tmax)[0]:.2f}" if tmin and tmax else "n.r."),
        ("Index ranges", ", ".join(f"{lo} ≤ {c} ≤ {hi}" for c, (lo, hi) in zip("hkl", hkl)) if hkl[0][0] else "n.r."),
        ("Reflections collected", _plain(g("_diffrn_reflns_number"))),
        ("Independent reflections", f"{_plain(nt)} [Rint = {rint_s}, Rsigma = {_plain(rsig, 4)}]"),
        ("Reflections with I > 2σ(I)", _plain(ngt)),
        ("Completeness to θ = " + (f"{num(thf)[0]:.3f}°" if thf else "θfull"), f"{100 * num(comp)[0]:.1f} %" if comp else "n.r."),
        ("Absorption correction", (g("_exptl_absorpt_correction_type") or "n.r.") + (f"; Tmin/Tmax {_plain(g('_exptl_absorpt_correction_t_min'), 3)}/{_plain(g('_exptl_absorpt_correction_t_max'), 3)}" if g("_exptl_absorpt_correction_t_min") else "")),
        ("Data/restraints/parameters", f"{_plain(g('_refine_ls_number_reflns'))}/{_plain(g('_refine_ls_number_restraints'))}/{_plain(g('_refine_ls_number_parameters'))}"),
        ("Goodness-of-fit on F²", _plain(g("_refine_ls_goodness_of_fit_ref"), 3)),
        ("Final R indexes [I > 2σ(I)]", f"R1 = {_plain(g('_refine_ls_r_factor_gt'), 4)}, wR2 = {_plain(g('_refine_ls_wr_factor_gt'), 4)}"),
        ("Final R indexes [all data]", f"R1 = {_plain(g('_refine_ls_r_factor_all'), 4)}, wR2 = {_plain(g('_refine_ls_wr_factor_ref'), 4)}"),
        ("Largest diff. peak/hole/e Å⁻³", f"{_plain(g('_refine_diff_density_max'), 3)}/{_plain(g('_refine_diff_density_min'), 3)}"),
        ("Flack parameter", _e(S, "_refine_ls_abs_structure_flack") if g("_refine_ls_abs_structure_flack") else "n.r."),
        ("CCDC deposition number", ccdc or "not deposited"),
    ]
    if twin:
        bas = g("_twin_individual_mass_fraction_refined")
        rows.append(("Note", "Twinned structure (refinement against twinned data)"))
    return pd.DataFrame(rows, columns=["Parameter", "Value"])


def key_numbers(S):
    t = dict(crystal_table(S).values)
    return {k: t[k] for k in ("Empirical formula", "Space group", "a/Å", "b/Å", "c/Å", "β/°", "Volume/Å³", "Z", "Radiation",
                              "Independent reflections", "Goodness-of-fit on F²", "Final R indexes [I > 2σ(I)]",
                              "Final R indexes [all data]", "Largest diff. peak/hole/e Å⁻³", "Flack parameter", "CCDC deposition number")}


def _tbl(df, cols=None):
    cols = cols or list(df.columns)
    cell = lambda t, b=False: f'<w:tc><w:p>{_run(t, b, 16)}</w:p></w:tc>'
    bd = "".join(f'<w:{e} w:val="single" w:sz="4" w:space="0" w:color="000000"/>' for e in ("top", "bottom", "insideH"))
    rows = ['<w:tr>' + "".join(cell(c, True) for c in cols) + '</w:tr>']
    rows += ['<w:tr>' + "".join(cell("" if pd.isna(v) else str(v)) for v in r) + '</w:tr>' for r in df[cols].itertuples(index=False)]
    return (f'<w:tbl><w:tblPr><w:tblW w:w="0" w:type="auto"/><w:tblBorders>{bd}</w:tblBorders></w:tblPr>'
            f'<w:tblGrid>{"".join("<w:gridCol/>" for _ in cols)}</w:tblGrid>{"".join(rows)}</w:tbl>')


def docx(sections):
    """sections: list of (title, DataFrame, columns|None, note)."""
    body = ""
    for title, df, cols, note in sections:
        body += f'<w:p>{_run(title, True, 20)}</w:p>' + (_tbl(df, cols) if df is not None and len(df) else f'<w:p>{_run("(none)", False, 16)}</w:p>')
        if note:
            body += f'<w:p>{_run(note, False, 14)}</w:p>'
    body += ('<w:sectPr><w:pgSz w:w="11906" w:h="16838"/><w:pgMar w:top="1134" w:right="1134" w:bottom="1134" w:left="1134" '
             'w:header="708" w:footer="708" w:gutter="0"/></w:sectPr>')
    W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
    doc = f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document xmlns:w="{W}"><w:body>{body}</w:body></w:document>'
    ct = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
          '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
          '<Default Extension="xml" ContentType="application/xml"/>'
          '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/></Types>')
    rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/></Relationships>')
    b = io.BytesIO()
    with zipfile.ZipFile(b, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", ct); z.writestr("_rels/.rels", rels); z.writestr("word/document.xml", doc)
    return b.getvalue()


def si_sections(S, geom, hb, n=1):
    hbd = pd.DataFrame(hb).drop(columns=["_A"], errors="ignore")
    if len(hbd):
        hbd = pd.DataFrame({"D–H···A": [f"{r['D']}–{r['H']}···{r['A']}" for _, r in hbd.iterrows()],
                            "d(D–H)/Å": hbd["D–H (Å)"].map(lambda v: f"{v:.2f}"), "d(H···A)/Å": hbd["H···A (Å)"].map(lambda v: f"{v:.2f}"),
                            "d(D···A)/Å": [fmt_esd(v, e) for v, e in zip(hbd["D···A (Å)"], hbd["σ D···A"])],
                            "∠DHA/°": hbd["D–H···A (°)"].map(lambda v: f"{v:.1f}"), "Symmetry (A)": hbd["symmetry (A)"]})
    c = lambda k, col: (geom[k], ["atoms", col, "symmetry"] if "symmetry" in geom[k] else ["atoms", col])
    b_, bc = c("bonds", "Distance (Å)"); a_, ac = c("angles", "Angle (°)"); t_, tc = c("torsions", "Torsion (°)")
    return [(f"Table S{n}. Crystal data and structure refinement for {S.name}.", crystal_table(S), None, ""),
            (f"Table S{n + 1}. Bond lengths (Å) for {S.name}.", b_, bc, "ESDs from the CIF (refinement); * computed from coordinates."),
            (f"Table S{n + 2}. Bond angles (°) for {S.name}.", a_, ac, ""),
            (f"Table S{n + 3}. Torsion angles (°) for {S.name}.", t_, tc, ""),
            (f"Table S{n + 4}. Hydrogen bonds for {S.name}.", hbd, None,
             "Criteria: H···A < Σr_vdW, ∠D–H···A ≥ 120°. D···A ESDs approximate (coordinate su only). Symmetry codes of the acceptor.")]


def xlsx(sections):
    b = io.BytesIO()
    with pd.ExcelWriter(b, engine="openpyxl") as w:
        for k, (title, df, cols, note) in enumerate(sections):
            nm = ["Crystal data", "Bonds", "Angles", "Torsions", "H-bonds"][k] if k < 5 else f"Table{k}"
            d = (df[cols] if cols else df) if df is not None and len(df) else pd.DataFrame({"": ["(none)"]})
            pd.DataFrame({title: [note or ""]}).to_excel(w, sheet_name=nm, index=False)
            d.to_excel(w, sheet_name=nm, index=False, startrow=3)
    return b.getvalue()
