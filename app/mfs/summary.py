"""v2.15: a Methods paragraph for every technique, and the sample summary table (xlsx / csv / LaTeX).

Only what the program does is stated; everything about the experiment is a [placeholder] for the author to fill in."""
from __future__ import annotations

import numpy as np

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
    "SEM / TEM / EDX": "Electron micrographs were recorded on a [microscope] at [kV]. Particle sizes were determined from [n] particles per sample "
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


def methods_text(pages, version: str = "") -> str:
    """A Methods paragraph per technique for the given page names (unknown pages are skipped) + the software sentence."""
    out = [f"**{p}.** {TECHNIQUE_METHODS[p]}" for p in dict.fromkeys(pages) if p in TECHNIQUE_METHODS]
    if out:
        out.append(f"**Data processing.** Data were processed and plotted with Reticula {version}".rstrip() +
                   " (open source; the settings of each figure are stored in the figure files and in the Supporting Information).")
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
