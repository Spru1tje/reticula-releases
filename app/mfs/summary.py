"""v2.15: a Methods paragraph for every technique, and the sample summary table (xlsx / csv / LaTeX).

Only what the program does is stated; everything about the experiment is a [placeholder] for the author to fill in."""
from __future__ import annotations


TECHNIQUE_METHODS = {
    "PXRD": "Powder X-ray diffraction patterns were recorded on a [diffractometer] with [radiation, λ = … Å] in [Bragg–Brentano / transmission] geometry over "
            "[2θ range] with a step of [step] °. Patterns were plotted with [background treatment, normalisation] and compared with patterns simulated from "
            "the crystallographic information files [CIF references].",
    "Structure (CIF)": "Reference patterns and structure figures were calculated from the crystallographic information files [CIF references] "
                       "(radiation [λ]; peak shape [pseudo-Voigt, FWHM]).",
    "N2 sorption (BET/PSD)": "N2 physisorption isotherms were measured at 77 K on a [instrument] after degassing at [T] °C for [t] h. The BET surface area was "
                             "calculated in the pressure range [p/p0 range] selected with the Rouquerol criteria (positive C constant, increasing n(1−p/p0)); "
                             "pore volumes were taken at p/p0 = [value] and pore size distributions were calculated with [method, kernel].",
    "FTIR / Raman": "Infrared / Raman spectra were recorded on a [spectrometer] in [ATR / transmission] mode from [range] cm−1 at [resolution] cm−1 ([n] scans). "
                    "Spectra were [baseline-corrected, normalised to …] before plotting; band positions were determined by peak picking.",
    "UV/Vis": "UV/Vis [absorption / diffuse-reflectance] spectra were recorded on a [spectrometer] from [range] nm. Diffuse-reflectance data were converted with "
              "the Kubelka–Munk function; optical band gaps were determined from Tauc plots [direct / indirect transition, linear range].",
    "Photoluminescence": "Photoluminescence spectra were recorded on a [spectrofluorometer] with an excitation wavelength of [λex] nm, slits of [slit] nm and "
                         "[temperature range / cryostat]. Band areas were obtained by [integration / peak fitting]; thermal quenching and luminescence-"
                         "thermometry parameters were obtained by non-linear least-squares fitting.",
    "1H NMR (digest)": "1H NMR spectra of the digested samples were recorded on a [spectrometer, MHz] in [solvent] after digestion in [digestion medium]. "
                       "Spectra were phased and baseline-corrected, referenced to [reference], and the signals were integrated to determine [linker / "
                       "modulator ratios].",
    "XPS": "X-ray photoelectron spectra were recorded on a [instrument] with [Al Kα / Mg Kα] radiation and a pass energy of [E] eV. Binding energies were "
           "referenced to the C 1s signal at [284.8] eV; a [Shirley / linear] background was subtracted and peaks were fitted with [line shape].",
    "ICP-OES / XRF": "Elemental composition was determined by [ICP-OES / XRF] on a [instrument] after [digestion / pellet preparation]. Values are the "
                     "mean of [n] measurements.",
    "TGA": "Thermogravimetric analysis was carried out on a [instrument] from [T1] to [T2] °C at [rate] K min−1 under [gas, flow]. Mass losses were "
           "determined between the temperatures given in the figure; the residual mass was taken at [T] °C.",
    "DSC": "Differential scanning calorimetry was carried out on a [instrument] from [T1] to [T2] °C at [rate] K min−1 under [gas]. Transition "
           "temperatures were taken as [onset / peak] values; enthalpies were obtained by integration with a [baseline type] baseline.",
    "SEM / TEM / EDX": "Bright-field TEM images were recorded on a JEOL JEM-1010 at 100 kV (tungsten filament, CCD camera, point resolution 0.4 nm); "
                       "high-resolution TEM images on a JEOL JEM-2200FS field-emission microscope at 200 kV (image-lens Cs corrector, EELS and EDX). "
                       "SEM images were recorded on a [SEM model] at [kV]. Particle sizes were determined from [n] particles per sample "
                       "by automatic segmentation of calibrated images (scale bar [value]) and are given as mean ± standard deviation. Elemental "
                       "maps / spectra were acquired by EDX.",
    "Kinetics & adsorption": "Kinetic data were fitted with the [pseudo-first / pseudo-second order] model by non-linear least squares; uptake capacities "
                             "were obtained from [isotherm model]. Conditions: [concentration, pH, temperature, dose].",
    "SDS-PAGE (Image Lab)": "Proteins were separated by SDS-PAGE on [gel %] gels and stained with [stain]; bands were quantified by densitometry "
                            "in Image Lab [version] relative to [reference].",
    "Zeta / DLS": "Hydrodynamic diameters (Z-average, intensity-weighted) and polydispersity indices (PDI) were measured by dynamic light scattering, and "
                  "ζ-potentials by electrophoretic light scattering, on a Malvern Zetasizer [model] at [temperature] °C in [medium, pH].",
    "Defect calculator": "Linker and modulator defects were calculated from the TGA mass loss and the NMR digestion ratios using the formulas given in "
                         "the Supporting Information [reference].",
}


# Journal wording templates. They set only the heading and the software sentence; the experimental text above stays the same.
# These are starting points: check the current author guidelines of the journal before submitting.
SI_STYLES = {
    "ACS": {"name": "ACS (J. Am. Chem. Soc., Chem. Mater.)", "heading": "Experimental Section", "label": "Data processing",
            "software": "Data were processed and plotted with Reticula{v} (open source; figure settings are stored in the figure files and in the Supporting Information)."},
    "Wiley-VCH": {"name": "Wiley-VCH (Adv. Mater., Angew. Chem.)", "heading": "Experimental Section", "label": "Data evaluation",
                  "software": "Data evaluation and figures were made with Reticula{v} (open-source software; settings are given in the Supporting Information)."},
    "RSC": {"name": "RSC (J. Mater. Chem. A, Chem. Sci.)", "heading": "Experimental", "label": "Data analysis",
            "software": "Data were analysed and plotted using Reticula{v} (open source). Figure settings are provided in the ESI."},
    "Elsevier": {"name": "Elsevier (Microporous Mesoporous Mater., Chem. Eng. J.)", "heading": "Materials and methods", "label": "Data processing",
                 "software": "Data were processed and plotted with Reticula{v} (open source); the settings are given in the Supplementary material."},
}
DEFAULT_SI_STYLE = "ACS"


def methods_text(pages, version: str = "", style: str = DEFAULT_SI_STYLE) -> str:
    """A Methods paragraph per technique for the given page names (unknown pages are skipped) + the software sentence, worded for the chosen journal style."""
    st_ = SI_STYLES.get(style, SI_STYLES[DEFAULT_SI_STYLE])
    out = [f"**{p}.** {TECHNIQUE_METHODS[p]}" for p in dict.fromkeys(pages) if p in TECHNIQUE_METHODS]
    if out:
        out.append(f"**{st_['label']}.** " + st_["software"].format(v=(" " + version.strip()) if version.strip() else ""))
    if out:
        out.insert(0, f"### {st_['heading']}")
    return "\n\n".join(out)


def sample_summary(store: dict, sections=None):
    """One row per sample, one column per 'technique: quantity' (scalar values only)."""
    import pandas as pd
    from mfs.si_bundle import _flat
    rows = []
    for sid, rec in (store or {}).items():
        r = {"Sample": rec.get("name", sid)}
        for sec, d in rec.get("sections", {}).items():
            if sections and sec not in sections:
                continue
            for k, v in _flat(d.get("values") or {}).items():
                if isinstance(v, (int, float, str, bool)) and not (isinstance(v, str) and len(v) > 80):
                    r[f"{sec}: {k}"] = round(v, 4) if isinstance(v, float) else v
        rows.append(r)
    return pd.DataFrame(rows)


def summary_latex(df, caption: str = "Summary of the characterisation of the samples.", label: str = "tab:summary") -> str:
    """A booktabs LaTeX table of the sample summary (empty cells '--', special characters escaped)."""
    def esc(x):
        s = "--" if x is None or str(x) == "nan" else str(x)
        for a, b in (("\\", "\\textbackslash{}"), ("&", "\\&"), ("%", "\\%"), ("_", "\\_"), ("#", "\\#"), ("$", "\\$")):
            s = s.replace(a, b)
        return s
    cols = list(df.columns)
    lines = ["\\begin{table}[ht]", "\\centering", f"\\caption{{{esc(caption)}}}", f"\\label{{{label}}}", "\\small",
             "\\begin{tabular}{" + "l" + "r" * (len(cols) - 1) + "}", "\\toprule", " & ".join(esc(c) for c in cols) + " \\\\", "\\midrule"]
    for _, r in df.iterrows():
        lines.append(" & ".join(esc(r[c]) for c in cols) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}", "\\end{table}"]
    return "\n".join(lines)


TECH_NAMES = {"tga": "TGA", "nmr": "¹H NMR (digest)", "nmr_1h": "¹H NMR", "icp": "ICP-OES", "xrf": "XRF", "xps": "XPS", "sem": "SEM / TEM sizing", "n2": "N₂ sorption",
              "pxrd": "PXRD", "kinetics": "Kinetics", "dsc": "DSC", "defects": "Defect calculator", "pl": "Photoluminescence", "uvvis": "UV/Vis"}


def summary_long(store: dict, sections=None):
    """One row per number: Sample, Technique, Quantity, Value, Source (page), Data files, Time, Version - the provenance of every value of the summary."""
    import pandas as pd
    from mfs.si_bundle import _flat
    rows = []
    for sid, rec in (store or {}).items():
        for sec, d in rec.get("sections", {}).items():
            if sections and sec not in sections:
                continue
            for k, v in _flat(d.get("values") or {}).items():
                if isinstance(v, (int, float, str, bool)) and not (isinstance(v, str) and len(v) > 80):
                    rows.append({"Sample": rec.get("name", sid), "Technique": TECH_NAMES.get(sec, sec), "Quantity": k, "Value": round(v, 4) if isinstance(v, float) else v,
                                 "Source": d.get("source", ""), "Data files": "; ".join(d.get("files") or []), "Time": d.get("updated", ""), "Version": d.get("version", "")})
    return pd.DataFrame(rows, columns=["Sample", "Technique", "Quantity", "Value", "Source", "Data files", "Time", "Version"])


def summary_notes(store: dict, sections=None):
    """The footnotes of the summary: one per (sample, technique) with page, files, time. Returns (list of (label, text), {(sample, technique): label})."""
    notes, ref = [], {}
    for sid, rec in (store or {}).items():
        for sec, d in rec.get("sections", {}).items():
            if sections and sec not in sections:
                continue
            lab = chr(ord("a") + len(notes)) if len(notes) < 26 else f"a{len(notes)}"
            tx = (f"{rec.get('name', sid)}, {TECH_NAMES.get(sec, sec)}: {d.get('source') or 'Reticula'}"
                  + (f", file(s) {'; '.join(d.get('files') or [])}" if d.get("files") else "") + (f", {d.get('updated', '')[:16].replace('T', ' ')}" if d.get("updated") else "")
                  + (f", Reticula {d.get('version')}" if d.get("version") else ""))
            notes.append((lab, tx))
            ref[(rec.get("name", sid), sec)] = lab
    return notes, ref


def summary_latex_sourced(store: dict, caption: str = "Summary of the characterisation of the samples.", label: str = "tab:summary") -> str:
    """The summary as a LaTeX table whose every value carries a footnote letter that names its source (page, files, time, program version)."""
    wide = sample_summary(store)
    if not len(wide):
        return ""
    notes, ref = summary_notes(store)
    sec_of = {TECH_NAMES.get(k, k): k for k in TECH_NAMES}

    def esc(x):
        s = "--" if x is None or str(x) == "nan" else str(x)
        for a_, b_ in (("\\", "\\textbackslash{}"), ("&", "\\&"), ("%", "\\%"), ("_", "\\_"), ("#", "\\#"), ("$", "\\$")):
            s = s.replace(a_, b_)
        return s
    cols = list(wide.columns)
    lines = ["\\begin{table}[ht]", "\\centering", f"\\caption{{{esc(caption)}}}", f"\\label{{{label}}}", "\\small",
             "\\begin{tabular}{" + "l" + "r" * (len(cols) - 1) + "}", "\\toprule", " & ".join(esc(c) for c in cols) + " \\\\", "\\midrule"]
    for _, r in wide.iterrows():
        cells = [esc(r["Sample"])]
        for c in cols[1:]:
            tech_key = c.split(":")[0]
            lab = ref.get((r["Sample"], tech_key))
            v = r[c]
            cells.append(esc(v) + (f"$^{{{lab}}}$" if lab and str(v) != "nan" else ""))
        lines.append(" & ".join(cells) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}", "\\\\[2pt]", "\\begin{minipage}{\\linewidth}\\footnotesize"]
    lines += [f"$^{{{lab}}}$ {esc(tx)}\\\\" for lab, tx in notes]
    lines += ["\\end{minipage}", "\\end{table}"]
    return "\n".join(lines)
