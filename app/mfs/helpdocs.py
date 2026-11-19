"""v2.18.3: the info documents of the Help and info page (one per module). English. Every abbreviation is written out where it is first used.
Sources: the user guide (app/README.md), the program's own code and tests. Advice that is general practice (not something the program checks)
is marked 'General practice'. A document = {"title", "subtitle", "sections": [(kind, text), ...]}, kinds: h p b n f k (see helppdf.py)."""
from __future__ import annotations

GP = "General practice (not checked by the program): "


def _doc(title, sub, *sections):
    return {"title": title, "subtitle": sub, "sections": [s for s in sections]}


def H(t): return ("h", t)
def P(t): return ("p", t)
def B(t): return ("b", t)
def N(t): return ("n", t)
def F(t): return ("f", t)
def K(t): return ("k", t)


def terms(*items):
    return [H("Terms used in this document")] + [B(i) for i in items]


DOCS: dict[str, dict] = {}

DOCS["Getting started"] = _doc(
    "Getting started with Reticula", "Read this first: it applies to every module.",
    H("What Reticula does"),
    P("Reticula turns raw instrument files into journal-style figures and does the usual analysis on the way. It runs on your computer; your data are never uploaded. "
      "Every technique has its own module in the sidebar. They all work the same way."),
    *terms("Supporting information (SI): the extra material that goes with a paper.",
           "Crystallographic information file (CIF): the standard text file with a crystal structure.",
           "Arbitrary units (a.u.): intensities that have no absolute scale.",
           "Dots per inch (dpi): the resolution of a picture."),
    H("Tutorial: from file to figure"),
    N("Pick a technique in the sidebar (grouped by Structure, Porosity, Spectroscopy, Thermal, Microscopy, Kinetics / Bio, Tools)."),
    N("Upload one or more files, or tick 'Use the bundled EXAMPLE data' to try the module first. The example files are synthetic."),
    N("Rename the samples and put them in the order you want in the sample list. Colours follow that order (colour-blind-safe Okabe-Ito palette by default). "
      "A colour rule in the sidebar keeps a colour for a sample in every figure."),
    N("Choose the processing (baseline, smoothing, normalisation, offset) and the figure width (ACS single column 3.25 in, double column 7 in, Wiley 8.5 or 17.5 cm)."),
    N("Open 'Figure controls' under the figure to change title, axes, units, ticks, legend and colours. Presets can be saved as JSON; the sidebar holds a global house style."),
    N("Download the figure as SVG, PDF, TIFF (600 dpi), PNG (600 dpi), the processed data as CSV, or everything as one ZIP. "
      "If your browser blocks downloads, use 'Save to folder'."),
    N("Use 'Add to combined figure' to put panels of different techniques into one figure, and tick 'Collect figures for the SI bundle' to get the SI export."),
    H("Detail levels"),
    P("The sidebar switch 'Detail level' has three settings. Standard shows only what a publishable figure needs. Advanced adds fit options, corrections, ticks and presets. "
      "Developer adds diagnostics. You can set it for all modules or for one module."),
    H("Simple mode"),
    P("Under a figure, Figure controls → Style & colours, 'Simple mode' makes the figure black only and tells the samples apart by shape: one black symbol per sample for data "
      "points, one line style per sample for continuous curves, grey bars with a hatch per series. Legends follow. Switch it off to get the colours back."),
    H("Style rules"),
    B("Arial or Helvetica 7–9 pt, inward ticks on four sides, no grid, no title by default (each can be changed in Figure controls)."),
    B("Units in parentheses, for example 2θ (°) or Wavenumber (cm⁻¹). Unit conversion only relabels the ticks; the CSV keeps the original units."),
    H("Limits"),
    B("The file readers were built from published format descriptions and tested on synthetic and some real files. If an import fails, export plain text as described in the Export guide; two-column text always works."),
    B("If the Arial font is missing, Liberation Sans or DejaVu Sans is used: check fonts in the SVG/PDF before submission."),
)

DOCS["PXRD"] = _doc(
    "PXRD – powder X-ray diffraction", "Patterns, simulated references, crystallite size, Rietveld and structure refinement.",
    H("What it does"),
    P("Plots one or more patterns (stacked or overlaid), simulates a reference pattern from a crystal structure, and can fit the pattern: crystallite size from peak broadening, "
      "Rietveld / Le Bail refinement, and a one-click identification and structure refinement."),
    *terms("Powder X-ray diffraction (PXRD): intensity versus the diffraction angle 2θ of a powder.",
           "Full width at half maximum (FWHM): the width of a peak at half its height.",
           "Crystallographic information file (CIF): the text file that holds a crystal structure.",
           "Rietveld refinement: fitting the whole pattern with a crystal structure. Le Bail fit: the same without a structure; the intensities are free, only cell and peak shape are fitted.",
           "Residual weighted profile (Rwp) and goodness of fit (GoF): numbers that tell how well the calculated pattern matches the measured one.",
           "Estimated standard deviation (esd): the uncertainty of a refined value.",
           "Instrument profile: the peak width the diffractometer itself adds. Caglioti parameters U, V, W (Gaussian part) and X, Y (Lorentzian part) describe it as a function of angle."),
    H("Tutorial: a figure"),
    N("Upload the patterns (STOE WinXPOW .xy/.dat/.txt/.csv; Rigaku .rasx/.ras). Order and rename them in the sample list."),
    N("Choose the processing: background removal (asymmetric least squares, rolling ball or polynomial), smoothing, normalisation, offset."),
    N("Optional: add a simulated pattern from a CIF (Cu Kα1, pseudo-Voigt peaks) as a grey trace or sticks, and highlight bands."),
    N("Adjust the figure in Figure controls and download."),
    H("Tutorial: crystallite size and microstrain (opt-in)"),
    N("Switch on 'Crystallite size and microstrain'. Peaks are fitted automatically."),
    N("Choose the instrumental broadening: typical laboratory value (0.05°), a standard you measured (LaB₆ or Si), or a value you know."),
    N("Read the Scherrer, Williamson–Hall and Halder–Wagner results with their uncertainties."),
    F(r"$D = \dfrac{K\,\lambda}{\beta\,\cos\theta}$"),
    P("D is the crystallite size, K the shape factor (about 0.9), λ the wavelength, β the sample broadening in radians and θ the Bragg angle. "
      "The sample broadening is the measured width with the instrument's contribution removed (Gaussian: β² = β_obs² − β_inst²; Lorentzian: β = β_obs − β_inst)."),
    K("A too small instrumental width gives a too large size. The program says so when the width was only assumed (0.05°) and when the sample width is as narrow as the instrument (size above about 100 nm: a lower limit)."),
    H("Tutorial: Rietveld / Le Bail (opt-in)"),
    N("Switch on 'Rietveld / Le Bail refinement'. Choose the phases from the library or upload CIF files. Or press 'Identify the phases'; or the 'Do everything automatically' button."),
    N("Choose the method: Rietveld (structure known), Rietveld + refine the atoms, or Le Bail. Choose the 2θ range."),
    N("Press Refine. A pattern takes about a second. Read the figure, Rwp, GoF, cell ± esd and weight fractions."),
    F(r"$R_{wp} = \sqrt{\dfrac{\sum_i w_i\,(y_i^{obs} - y_i^{calc})^2}{\sum_i w_i\,(y_i^{obs})^2}}$"),
    P("Advanced level adds: background polynomial order, monochromator angle, instrument offset (sample displacement, zero shift or both), preferred orientation, and bond-restraint strength."),
    H("Background, glass holder and standard (Advanced)"),
    B("The background is a polynomial (default order 6, range 2–14). On top of it, broad humps between 15° and 35° 2θ model the amorphous signal of a glass sample holder; this is automatic and can be switched off."),
    B("You can load a measured background (an empty holder pattern); it is scaled and subtracted in the refinement."),
    B("Instrument profile: refined from the pattern by default. Or fixed from a Si standard (SRM 640) measured on the same instrument: the standard is fitted with Le Bail at its certificate lattice parameter (5.4312 Å, not refined) and U, V, W, X, Y are then kept fixed for the sample."),
    H("What is checked"),
    B("The sample displacement and the zero shift are not refined together by default (they are strongly correlated and bias the cell); 'both' needs a standard."),
    B("The program warns when the instrument broadening was not removed from the apparent size, and when the model misses something (large difference curve)."),
    B("Checked on synthetic data: the fixed Si profile gives the right cell for a second material. The fitted profile is less reliable above about 80° 2θ: check it against your own standard."),
    H("Limits"),
    B("Binary STOE .raw files need the optional xylib package; export .xy instead."),
    B("Indexing covers cubic, tetragonal, hexagonal and orthorhombic cells, not monoclinic or triclinic."),
    B("Only Si is built in as a standard profile; LaB₆ is not yet."),
)

DOCS["Structure (CIF)"] = _doc(
    "Structure (CIF) – crystal structure drawings and tables", "From a CIF to figures, SI tables and a simulated pattern.",
    H("What it does"),
    P("Reads a CIF and draws the structure, makes the tables for the SI, simulates a powder pattern and gives some pore and surface-contact quantities."),
    *terms("Crystallographic information file (CIF): the text file that holds a crystal structure (from SHELXL, Olex2, COD or CSD).",
           "Oak Ridge thermal-ellipsoid plot (ORTEP): atoms drawn as ellipsoids of their thermal motion (30, 50 or 70 % probability).",
           "Asymmetric unit (ASU): the smallest part of the cell from which symmetry builds the rest.",
           "Estimated standard deviation (esd): the uncertainty of a bond length or angle.",
           "Hirshfeld surface: a surface around a molecule that shows how it touches its neighbours."),
    H("Tutorial"),
    N("Upload a CIF."),
    N("Choose the drawing style (ORTEP, ball-and-stick, stick, space-filling), the view (along an axis, a [uvw] direction, a (hkl) plane or best view), and what to show (ASU, molecules, packing, supercell, hydrogen bonds)."),
    N("Download the figure, or the SI tables (crystal data, bond lengths, angles, torsions with esd) as DOCX or XLSX."),
    N("Optional: simulated powder pattern with hkl labels and overlay on your measured pattern; MOF pore quantities (Di, Df, void fraction, accessible surface area); Hirshfeld fingerprint (experimental)."),
    H("Limits"),
    B("The drawings here have no colour controls, because the colours carry information; Simple mode does not apply."),
    B("A NU-1000 CIF is not bundled; load your own."),
)

DOCS["N2 sorption (BET/PSD)"] = _doc(
    "N₂ sorption – isotherms, BET and pore size distribution", "From the instrument export to surface area, micropore volume and pore sizes.",
    H("What it does"),
    P("Plots nitrogen isotherms and calculates the BET surface area with the Rouquerol criteria, the t-plot micropore volume and external area, and pore size distributions."),
    *terms("Brunauer–Emmett–Teller (BET) theory: the monolayer capacity from the adsorption isotherm gives the surface area.",
           "Relative pressure (p/p₀): the pressure divided by the saturation pressure.",
           "Pore size distribution (PSD): how the pore volume is spread over pore widths.",
           "Density functional theory (DFT), Horváth–Kawazoe (HK), Saito–Foley (SF), Barrett–Joyner–Halenda (BJH): methods to calculate a PSD from the isotherm.",
           "t-plot (with the Harkins–Jura thickness): separates micropore volume from external area.",
           "Standard temperature and pressure (STP): the reference for gas volumes (cm³ g⁻¹ STP)."),
    H("Tutorial"),
    N("Upload the MicroActive Excel export (.xls/.xlsx, full report or the Isotherm Tabular table) or a CSV. The 'File structure' expander shows which tables were found."),
    N("Plot the isotherms (adsorption filled, desorption open; linear and/or logarithmic p/p₀)."),
    N("Read the BET table: the app's Rouquerol BET next to the value in the file, and the app's value on the file's p/p₀ range. Check the Rouquerol criteria panel."),
    N("Choose a PSD method (DFT kernel fit, HK, SF or BJH) and read the t-plot values."),
    F(r"$\dfrac{1}{v\,[(p_0/p)-1]} = \dfrac{c-1}{v_m\,c}\,\dfrac{p}{p_0} + \dfrac{1}{v_m\,c}$"),
    P("v is the adsorbed quantity, v_m the monolayer capacity and c a constant; the straight-line fit over a p/p₀ window gives v_m, from which the area follows."),
    H("What is checked"),
    P("The Rouquerol criteria (ISO 9277) choose the fitting window: (1) the constant c must be positive; (2) the quantity v·(1 − p/p₀) must keep increasing with p/p₀ in the window; "
      "(3) the relative pressure of the monolayer must lie inside the window; (4) it must agree with the value predicted by the BET equation. The app also shows the range of areas over all windows."),
    B("Ultramicroporous MOFs (MOF-801, UiO-66) often cannot meet criterion 3 or 4. The app flags this; then report a 'BET-equivalent area' and the p/p₀ range used."),
    B("The instrument usually fits a fixed window (0.05–0.30). For microporous MOFs the Rouquerol window lies at much lower p/p₀, so the two values differ by design."),
    H("Limits"),
    B("The DFT kernel is a carbon-slit N₂ 77 K kernel; it is not an oxide NLDFT/QSDFT kernel, so PSDs of oxides are approximate."),
    B("MicroActive .SMP/.QPS files cannot be read; use the Excel export."),
)

DOCS["FTIR / Raman"] = _doc(
    "FTIR and Raman – vibrational spectra", "Stacked spectra, band labels, band fitting and difference spectra.",
    H("What it does"),
    P("Plots infrared and Raman spectra, processes them, labels bands and fits overlapping bands."),
    *terms("Fourier-transform infrared spectroscopy (FTIR), here with attenuated total reflection (ATR): infrared absorption of a solid pressed on a crystal.",
           "Fourier-transform Raman (FT-Raman): scattered light with a 1064 nm laser.",
           "Transmittance (T) and absorbance (A): A = −log₁₀ T.",
           "Wavenumber (cm⁻¹): the horizontal axis, the number of waves per centimetre.",
           "Asymmetric least squares (ALS): a baseline-removal method."),
    H("Tutorial"),
    N("Choose the spectrum type: ATR-FTIR (transmittance), ATR-FTIR (absorbance) or FT-Raman."),
    N("Upload the files: Bruker OPUS binaries (.0, .1, …), OMNIC .spa (one spectrum per file), .dpt, or CSV/TXT. Header-less OMNIC CSV exports are read as one spectrum."),
    N("If a text file holds several spectra (a wavenumber column plus one column per spectrum), tick 'The file(s) contain several spectra'. Each column then becomes a spectrum, named by its header."),
    N("Set the range, the baseline (none for FTIR, ALS for Raman), normalisation, and highlight bands ('lo-hi:label')."),
    N("Advanced: normalise to a reference band, automatic peak labels, a zoomed inset, difference spectra, and a band-fit (assignment table and simultaneous multi-band fit)."),
    F(r"$A = -\log_{10} T$"),
    H("What is checked"),
    B("OMNIC .spa files were checked against their own CSV export: same points and intensities (wavenumber within 0.001 cm⁻¹); the .spa holds the full measured range."),
    B("An OMNIC .spg (a series) is not read."),
)

DOCS["UV/Vis"] = _doc(
    "UV/Vis – absorbance, assays, band gaps", "Spectra, Beer–Lambert calibration, assay presets, kinetics and Tauc band gaps.",
    H("What it does"),
    P("Plots ultraviolet–visible spectra, and calculates concentrations, assay results, single-wavelength kinetics and optical band gaps."),
    *terms("Ultraviolet–visible spectroscopy (UV/Vis): light absorption between about 200 and 800 nm.",
           "Absorbance (A), transmittance (%T) and reflectance (%R).",
           "Limit of detection (LOD) and limit of quantification (LOQ): the smallest amounts that can be detected or measured reliably.",
           "Kubelka–Munk function: turns diffuse reflectance into a quantity proportional to absorption.",
           "Tauc plot: a plot to read the optical band gap E_g from the absorption edge.",
           "Bicinchoninic acid (BCA), o-phthalaldehyde (OPA): reagents of protein and amine assays."),
    H("Tutorial"),
    N("Upload spectra (Implen NanoPhotometer Excel, Cary or generic CSV, JCAMP-DX)."),
    N("Choose the quantity (A, %T or %R), blank, smoothing and the maximum wavelength (λmax)."),
    N("Beer–Lambert: enter standards, fit the calibration (with LOD and LOQ) and read unknowns. Pick an assay preset (BCA, Bradford, Lowry, ninhydrin, OPA)."),
    N("Band gap: use Kubelka–Munk and a Tauc plot on reflectance data."),
    F(r"$A = \varepsilon\, c\, l$"),
    F(r"$F(R) = \dfrac{(1-R)^2}{2R}\qquad (F(R)\,h\nu)^{n}\ \mathrm{versus}\ h\nu$"),
    P("ε is the molar absorption coefficient, c the concentration, l the path length; R the reflectance, hν the photon energy. The exponent n depends on the type of transition (2 for direct, 1/2 for indirect)."),
    B(GP + "read the band gap from the linear part of the Tauc plot and state which n you used."),
)

DOCS["Photoluminescence"] = _doc(
    "Info implementatie ratiometric thermometry in Reticula", "Photoluminescence module: spectra, thermal quenching and ratiometric thermometry.",
    H("1. The idea in plain words"),
    P("Some luminescent materials (for example Er³⁺-doped upconversion particles or phosphors) emit light in two bands whose energy levels lie close together. "
      "When the material is heated, the upper level is populated more strongly, so the band at the shorter wavelength grows relative to the other one. "
      "We therefore do not use the absolute brightness (it depends on laser power, sample position and concentration) but the RATIO of the two bands. "
      "That ratio is the thermometer: ratiometric luminescence thermometry."),
    H("2. Terms, defined before they are abbreviated"),
    B("Photoluminescence (PL): light emitted by a material after it was excited with light."),
    B("Intensity of band 1 and band 2 (I1, I2): the area under the shorter- and the longer-wavelength emission band."),
    B("Ratio (Δ): I1 divided by I2."),
    B("Energy gap (ΔE): the energy difference between the two emitting levels, in cm⁻¹ or eV."),
    B("Boltzmann constant (kB): 0.695 cm⁻¹ per kelvin; the temperature T is in kelvin (K)."),
    B("Relative sensitivity (Sr): how many percent the ratio changes per kelvin."),
    B("Temperature uncertainty (σT): the smallest temperature change that can be told from the noise, in K."),
    B("Full width at half maximum (FWHM): the width of a band at half of its height."),
    B("Bayesian information criterion (BIC): a score that compares two models and punishes extra parameters; the lower score wins."),
    B("Coefficient of determination (R²) and reduced chi-squared (χ²red): numbers that tell how well a fit describes the data."),
    B("Standard deviation (SD), arbitrary units (a.u.), supporting information (SI): used in the usual way."),
    H("3. The formulas"),
    F(r"$\Delta = \dfrac{I_1}{I_2} = \alpha\,\exp\!\left(-\dfrac{\Delta E}{k_B T}\right)\ (+\,C)$"),
    B("α is a scale factor of the material; C is an optional constant (tick 'Add a constant C')."),
    P("Relative sensitivity (percent per kelvin):"),
    F(r"$S_r = 100\%\cdot\left|\dfrac{1}{\Delta}\dfrac{\partial\Delta}{\partial T}\right| = 100\%\cdot\dfrac{\Delta E}{k_B T^{2}}$"),
    P("Temperature uncertainty with photon counts:"),
    F(r"$\sigma_T = \dfrac{1}{S_r}\sqrt{\dfrac{1}{I_1}+\dfrac{1}{I_2}}$"),
    P("If the error of the ratio itself is known (from the global fit), the general form is σT = (σΔ / Δ) / Sr. The absolute sensitivity Sa = |dΔ/dT| = (Δ − C)·ΔE / (kB T²) is calculated as well."),
    B("Consequence: Sr falls with 1/T². A thermometer is more sensitive at low temperature, and a larger ΔE gives a higher Sr."),
    H("4. The steps (tab 'Thermometry (band ratio)')"),
    N("Load the spectra: one table (wavelength plus one column per temperature) or separate files. The temperature is read from the column header or file name ('293.15 K', '25C', 'T=300K'); numbers above 150 are read as K, smaller ones as °C. At least 3 temperatures are needed."),
    N("Find the bands: automatically (the two strongest separate bands, also hidden shoulders) or choose two integration windows yourself."),
    N("Determine the ratio at every temperature (section 6)."),
    N("Fit the Boltzmann law to Δ versus T. Result: ΔE (cm⁻¹ and eV, with error), α, R²."),
    N("Calculate Sr and σT per temperature and draw the figures: emission map, Δ versus T with fit, Sr(T), σT(T). They can go into the combined figure."),
    H("5. The background polynomial"),
    P("Under the emission bands there is always a background: detector offset, stray light, or the tail of a neighbouring emission. It is described by a polynomial in wavelength x:"),
    F(r"$\mathrm{background}(x) = b_0 + b_1 x + b_2 x^2 + b_3 x^3$"),
    P("The order is the highest power used. Order 0 is a flat offset, order 1 a straight sloping line, order 2 a curve, order 3 an S-shaped curve."),
    B("Default: order 1 (a straight line). In the Advanced level the slider 'Background polynomial order' allows 0 to 3."),
    B("Every spectrum gets its own background over the fit window (the bands plus 3 FWHM on each side), built from Chebyshev polynomials (a numerically stable form of the same curve)."),
    B("The coefficients are solved exactly together with the band areas in every step of the fit (the problem is linear in them)."),
    B("The integration-window method uses no polynomial: it uses the area above a straight line through the two ends of each window."),
    P(GP + "keep order 1. Use 0 for a flat baseline. Go to 2 only if the residual shows a smooth curve. A high order can absorb band tails and shift the ratio: check that the ratio does not change when the order changes by one."),
    H("6. Two ways to determine the ratio"),
    P("A. Global fit (recommended). All spectra of the series are described together by the SAME bands (pseudo-Voigt shapes with unit area). Position and width may drift linearly with temperature; areas and background of each spectrum are solved exactly. "
      "An extra 'nuisance' band (a shoulder or third peak) can be added so that it does not contaminate the ratio. The simple model (fixed shape) and the drifting model are compared with the BIC; the better one is used."),
    P("B. Integration windows. The classic method: the area in two windows above a straight line through the window ends. Drawback: overlapping bands are cut at an arbitrary valley and the result depends on your windows."),
    P("'Compare both' shows the two next to each other. If they differ strongly, the bands overlap and the global fit is preferred."),
    H("7. What is checked"),
    B("Error of the ratio = statistical error (including the covariance of the two areas) + systematic error. The systematic part is the spread of refits with narrower or wider windows and with the other model; total = √(statistical² + systematic²)."),
    B("Quality flags per spectrum: R², χ²red and the autocorrelation of the residuals. Badly fitting spectra are marked: do not trust those points."),
    B("The program looks for structure left in the residuals and proposes an extra band only if the BIC improves (at most 3 extra bands)."),
    B("The ratio must be positive for every spectrum; otherwise a warning appears and the Boltzmann fit stops."),
    B("σT from I1 and I2 is calculated only for photon counts. For intensities in a.u. leave 'counts' unticked; σT then follows from the error of the global fit."),
    B("Saturation: the extras panel measures the power dependence (slope n < 1 at high power means saturation)."),
    H("8. The other tabs of the module"),
    B("Spectra: coloured by temperature, waterfall, emission map; nm, eV or cm⁻¹ axes; peak position and FWHM."),
    B("Temperature dependence: thermal quenching (activation energy E_a, T₅₀, Arrhenius plot)."),
    H("9. Practical tips"),
    B(GP + "use as many temperatures as possible over a wide range; a Boltzmann fit through 3 points is weak."),
    B(GP + "correct the spectrometer response (detector and grating); otherwise the ratio is wrong even when the fit looks good."),
    B(GP + "laser heating falsifies the temperature: check that the ratio does not change with power."),
    B(GP + "report ΔE, Sr at a reference temperature and σT, not only the curve."),
)

DOCS["1H NMR (digest)"] = _doc(
    "¹H NMR – spectra, assignment and digestion analysis", "From a Bruker folder to a labelled figure, a peak list, and linker : modulator ratios.",
    H("What it does"),
    P("Processes a one-dimensional proton NMR spectrum, finds the multiplets, assigns peaks to atoms of a structure you draw, and (in the digestion mode) turns integrals into ratios and concentrations."),
    *terms("Nuclear magnetic resonance (NMR): the radio-frequency response of nuclei in a magnetic field; ¹H is the proton.",
           "Chemical shift (δ, in ppm): the position of a signal relative to a reference.",
           "Coupling constant (J, in Hz): the splitting of a signal by neighbouring protons.",
           "Quantitative NMR (qNMR): concentrations from integrals compared with an internal standard.",
           "Trimethylsilylpropanoic acid (TMSP): a common internal standard (here TMSP-d4).",
           "Metal–organic framework (MOF): a porous crystalline material that is digested (dissolved) before measurement."),
    H("Tutorial: characterisation"),
    N("Upload a Bruker experiment folder as .zip (fid or processed pdata/1/1r), a JCAMP-DX file, a TopSpin 'totxt' export, or x,y text."),
    N("Phase, baseline and referencing are automatic; adjust if needed. Multiplets and J couplings are found."),
    N("Draw the structure in the built-in editor (or paste a molfile or SMILES). The program predicts shifts (nmrshiftdb2 table) and assigns peaks to atoms; drag to correct."),
    N("Check impurities against the Fulmer 2010 and Gottlieb solvent tables. Make the labelled figure and the ACS-style peak list for the SI."),
    P("'Simple mode (black)' draws the whole figure in black."),
    H("Tutorial: digestion and qNMR"),
    N("Switch to the digestion mode. Choose the integration windows for linker, modulator and the internal standard."),
    N("Read linker : modulator ratios and concentrations versus the standard. They feed the defect calculator."),
    F(r"$n_i \propto \dfrac{A_i}{N_{H,i}}\qquad \dfrac{n_i}{n_{std}} = \dfrac{A_i\,/\,N_{H,i}}{A_{std}\,/\,N_{H,std}}$"),
    P("A is the integral of a signal and N_H the number of protons it represents."),
    H("What is checked"),
    B("Phasing, baseline and referencing are automatic; the shift prediction is a database lookup, not a quantum calculation: check assignments by eye."),
    B(GP + "use a long enough relaxation delay (D1) for quantitative integrals."),
)

DOCS["XPS"] = _doc(
    "XPS – X-ray photoelectron spectroscopy", "Region spectra, background, doublet fitting and atomic percentages.",
    H("What it does"),
    P("Plots XPS regions (binding energy axis reversed), subtracts a background, fits doublets and gives atomic percentages."),
    *terms("X-ray photoelectron spectroscopy (XPS): the energy of electrons emitted by X-rays tells the elements and their chemical state.",
           "Binding energy (BE): the energy axis, in eV.",
           "Shirley and Tougaard backgrounds: two common models of the step background under a peak.",
           "Spin–orbit doublet: a pair of peaks (for example Zr 3d, Hf 4f, S 2p) with a fixed splitting and area ratio.",
           "Relative sensitivity factor (RSF): the correction that turns peak areas into amounts (Scofield values by default)."),
    H("Tutorial"),
    N("Upload CSV/TXT region files (CasaXPS, Avantage) or simple VAMAS (.vms)."),
    N("Choose the region, the background (Shirley or Tougaard) and the components. For a doublet the splitting and the area ratio are fixed."),
    N("Read peak positions and areas; give or edit RSFs to get atomic percent."),
    F(r"$x_i = \dfrac{I_i / \mathrm{RSF}_i}{\sum_j I_j / \mathrm{RSF}_j}\times 100\%$"),
    B("VAMAS support is minimal (single-block regions). Use the multi-region EXAMPLE to try the module."),
)

DOCS["ICP-OES / XRF"] = _doc(
    "ICP-OES and XRF – elemental analysis", "From concentrations to weight percent and molar ratios.",
    H("What it does"),
    P("Reads elemental analysis tables and calculates weight percent, mean ± standard deviation bars and molar ratios (also Zr : Hf : S with propagated uncertainty)."),
    *terms("Inductively coupled plasma optical emission spectrometry (ICP-OES): elements in a digested sample, as mg/L.",
           "X-ray fluorescence (XRF): elements in a solid from their characteristic X-rays.",
           "Standard deviation (SD): the spread over replicates.",
           "Metal node: the metal cluster of a MOF (for example M₆)."),
    H("Tutorial"),
    N("Upload a Syngistix CSV (long or wide) or any table with sample, element and value."),
    N("Enter the digested mass and the volume to convert mg/L to weight percent."),
    N("Choose the bars (one panel per element or grouped) and the molar ratios (pick the reference metal, per n atoms)."),
    N("The Zr : Hf : S case gives Hf/(Zr + Hf) and S per M₆ node with propagated uncertainty. The values go to the sample store and the defect calculator."),
    F(r"$w\,[\mathrm{wt\%}] = \dfrac{c\,[\mathrm{mg/L}]\cdot V\,[\mathrm{L}]}{m\,[\mathrm{mg}]}\times 100\%$"),
)

DOCS["TGA"] = _doc(
    "TGA – thermogravimetric analysis", "Mass loss curves, derivative, steps and the oxide residue.",
    H("What it does"),
    P("Plots mass versus temperature with the derivative, labels steps, normalises to the oxide residue and passes values to the defect calculator."),
    *terms("Thermogravimetric analysis (TGA): the mass of a sample while it is heated.",
           "Derivative thermogravimetry (DTG): the rate of mass change; its minima mark the steps.",
           "Oxide residue: the metal oxide (MO₂) that remains after combustion in air.",
           "Quadrupole mass spectrometry (QMS): a mass spectrometer coupled to the TGA to see what leaves."),
    H("Tutorial"),
    N("Upload Netzsch Proteus ASCII, TA Universal Analysis / TRIOS text, Mettler STARe export, or generic CSV (mg or %)."),
    N("Choose the DTG display (second panel or twin axis) and optional renormalisation to the oxide residue = 100 %."),
    N("Read the automatic step labels; set the plateau if needed. The values W(200 °C), plateau and end go to the defect calculator. A TGA–MS overlay is available."),
    B("The plateau is picked automatically at the DTG minimum between 250 and 480 °C. Check it, and override it in the defect tab if dehydroxylation and combustion overlap."),
    B(GP + "measure in air so that the end product is the oxide MO₂."),
)

DOCS["DSC"] = _doc(
    "DSC – differential scanning calorimetry", "Peaks, glass transitions, cycles and heat capacity.",
    H("What it does"),
    P("Plots heat flow versus temperature or time and calculates transition temperatures, enthalpies, glass transitions, cycle comparisons and heat capacity."),
    *terms("Differential scanning calorimetry (DSC): the heat flow into a sample while it is heated or cooled.",
           "Onset, peak and endset: the start, top and end of a thermal event.",
           "Enthalpy (ΔH): the heat of a transition, from the area of the peak (J/g).",
           "Glass transition (Tg) and heat-capacity step (Δcp).",
           "Heat capacity (cp): heat needed per gram and kelvin; measured against a sapphire reference (ASTM E1269).",
           "Simultaneous thermal analysis (STA): TGA and DSC together."),
    H("Tutorial"),
    N("Upload Mettler STARe, TA TRIOS / Universal Analysis, Netzsch Proteus, PerkinElmer Pyris text exports, or generic CSV/XLSX."),
    N("Choose the units, the exo direction (up or down) and the segment colours."),
    N("Pick a peak and a baseline (linear or sigmoidal) to get onset, peak, endset and ΔH; or Tg (midpoint, inflection, Δcp)."),
    N("Full cycle: heating and cooling segments are found from dT/dt; compare Heat 1 with Heat 2; see the supercooling ΔT = Tm − Tc."),
    N("Heat capacity: load the blank, the sapphire and the sample; the ratio method gives cp."),
    F(r"$\Delta H = \dfrac{1}{\beta\, m}\int (\Phi - \Phi_{base})\,dT$"),
    P("Φ is the heat flow, β the heating rate and m the mass."),
)

DOCS["SEM / TEM / EDX"] = _doc(
    "SEM, TEM and EDX – images, scale bars and particle sizes", "Calibrated images, scale bars, particle statistics and element maps.",
    H("What it does"),
    P("Shows electron-microscope images with a correct scale bar, measures particle sizes and draws element maps."),
    *terms("Scanning electron microscopy (SEM) and transmission electron microscopy (TEM).",
           "High-resolution TEM (HRTEM) and bright-field TEM: image modes of a TEM.",
           "Energy-dispersive X-ray spectroscopy (EDX): element maps and spectra from X-rays.",
           "Pixel size (nm/px): the calibration that turns pixels into nanometres.",
           "Instruments named in the documentation: JEOL JEM-1010 (bright-field TEM, 100 kV, tungsten filament, CCD camera, 0.4 nm resolution, up to 500 000×) and JEOL JEM-2200FS (HRTEM, 200 kV field-emission gun, image-lens spherical-aberration corrector, electron energy loss spectrometer and EDX)."),
    H("Tutorial"),
    N("Upload the ORIGINAL image: TIFF (with the databar), TEM .dm3/.dm4/.ser/.emd/.mrc, or PNG/JPG with a manual pixel size. EDX maps are greyscale images, one per element."),
    N("Check the calibration panel: it compares the file metadata, the scale bar burned into the info strip (measured, with text recognition) and the magnification. If all fail, click the two bar ends and type its value."),
    N("Set the scale bar (value or length, position, colour, thickness, font, box), build the (a)–(d) panel grid, and use matched bars or one calibration for a batch."),
    N("Particle sizing: automatic recognition with background flattening; click to exclude or to measure; log-normal histogram. EDX maps can be overlaid on the SE image; line scans are available."),
    B("JEM-1010: save the camera TIFF at full resolution with the scale bar on (Olympus SIS / analySIS / MegaView III files are read). JEM-2200FS: save the raw Digital Micrograph or Velox file, not a screenshot."),
    B("The 127 mm magnification estimate is only valid for Philips/FEI XL SEMs."),
)

DOCS["Kinetics & adsorption"] = _doc(
    "Kinetics and adsorption", "Rate constants, activation parameters, turnover frequencies and isotherms.",
    H("What it does"),
    P("Fits concentration-versus-time data and adsorption isotherms."),
    *terms("Observed rate constant (k_obs) and half-life (t½).",
           "Activation energy (E_a), activation enthalpy (ΔH‡) and entropy (ΔS‡).",
           "Turnover frequency (TOF): conversions per open site per hour.",
           "Confidence interval (CI): the range that holds the true value with 95 % probability.",
           "Langmuir and Freundlich: two adsorption isotherm models."),
    H("Tutorial"),
    N("Upload a table: sample, replicate, time, concentration; or raw NMR integrals with an internal standard."),
    N("Read the first-order fits (k_obs ± SD, t½, R²), the C/C₀ plots and the k_obs bars. Blank subtraction is available."),
    N("With several temperatures: Arrhenius and Eyring analysis with confidence intervals; initial rates; TOF per open site (open sites from the defect calculator); pD–rate profile (pKa)."),
    N("Adsorption table (sample, Ce, qe): Langmuir and Freundlich fits."),
    F(r"$C = C_0\,e^{-k t}\qquad t_{1/2} = \dfrac{\ln 2}{k}$"),
    F(r"$k = A\,e^{-E_a/RT}\qquad k = \dfrac{k_B T}{h}\,e^{\Delta S^\ddagger/R}\,e^{-\Delta H^\ddagger/RT}$"),
    F(r"$q = \dfrac{q_{max} K C_e}{1 + K C_e}\ \ (\mathrm{Langmuir})\qquad q = K_F C_e^{1/n}\ \ (\mathrm{Freundlich})$"),
)

DOCS["SDS-PAGE (Image Lab)"] = _doc(
    "SDS-PAGE – protein cleavage from gel bands", "Band percentages versus time.",
    H("What it does"),
    P("Reads the band table of a gel and plots how much of each band (intact protein and fragments) there is over time."),
    *terms("Sodium dodecyl sulfate polyacrylamide gel electrophoresis (SDS-PAGE): proteins separated by size in a gel.",
           "Image Lab: the Bio-Rad software that quantifies the bands (Gel Doc EZ)."),
    H("Tutorial"),
    N("Upload the Image Lab 'Lane and Bands' table (CSV or XLSX)."),
    N("Map each lane to a time point."),
    N("Read the band percentage versus time for each band."),
)

DOCS["Zeta / DLS"] = _doc(
    "Zeta potential and dynamic light scattering", "Size, polydispersity and surface charge from Malvern Zetasizer data.",
    H("What it does"),
    P("Reads Zetasizer data (CSV, Excel or the .dts file), and plots zeta potentials, size distributions and the extra quantities of the instrument report."),
    *terms("Dynamic light scattering (DLS): particle size from the flicker of scattered light.",
           "Z-average: the intensity-weighted mean hydrodynamic diameter (nm), from a cumulant fit of the correlation function.",
           "Polydispersity index (PDI): the width of the size distribution (a small number means uniform).",
           "Zeta potential (ζ, in mV): the electrical potential at the particle's slipping plane, from the electrophoretic mobility.",
           "Zeta deviation: the width (standard deviation) of the zeta distribution of one record.",
           "Conductivity (mS/cm) and count rate (kcps, thousands of counts per second)."),
    H("Tutorial"),
    N("Upload the export or the .dts file. For .dts files of software 8.x every record is read; size records give names, the intensity distribution, Z-average and PDI as stored by the software; zeta records give ζ, deviation, conductivity, count rate, the three zeta peaks (mean, area %, SD) and the zeta distribution."),
    N("Untick runs (records) you want to leave out of the averages and graphs."),
    N("Read the table (also saved for the SI), choose what to show in the figures (Z-average, PDI, ζ, ± SD), and download the figures."),
    N("The Methods text is generated from the data; fill in the bracketed items yourself."),
    B("Checked on real files against the Zetasizer reports: Z-average, PDI, ζ, deviation, conductivity and the peaks are equal to the reported values."),
    B("The zeta distribution is shown only for records whose stored block reproduces the software's ζ."),
    B("Software 7.11 files: Z-average and PDI are calculated from the stored correlation function; ζ needs the Excel/CSV export."),
    F(r"$d_H = \dfrac{k_B T}{3\pi\,\eta\,D}$"),
    P("The Stokes–Einstein relation behind the Z-average: d_H is the hydrodynamic diameter, η the viscosity and D the diffusion coefficient."),
)

DOCS["Defect calculator"] = _doc(
    "Defect calculator – missing linkers in MOFs", "Every method the data allows, side by side, with a consensus.",
    H("What it does"),
    P("Takes the results from TGA, NMR, ICP, N₂ sorption and titration (or typed values) and calculates the number of linkers per metal node by every possible method."),
    *terms("Linkers per node (L), missing linkers (x), modulator or cap per node (m).",
           "Open sites (OS): metal positions free for catalysis.",
           "Thermogravimetric analysis (TGA), nuclear magnetic resonance (NMR), inductively coupled plasma (ICP), Brunauer–Emmett–Teller (BET).",
           "Consensus: the median of the quantitative methods, with the minimum, maximum and standard deviation."),
    H("Tutorial"),
    N("Fill the data (or pull it from the sample store, which holds results from the other modules). Choose the framework and node type."),
    N("Read each method: equation, symbols, assumptions, a step-by-step substitution with your own numbers, pitfalls and the DOI."),
    N("Read the consensus; download the PDF report, CSV table or JSON input."),
    H("Methods"),
    B("a: TGA plateau (ideal-normalised). b: M-200 (mass at 200 °C + NMR ratios + Cl with charge balance) and b′ the missing-cluster model. c: TGA plateau corrected with the NMR modulator ratio. d: experimental minimal formula (oxide residue, qNMR, Cl, water by difference). "
      "e: ICP metal : linker ratio and mass balance. f: porosity versus the ideal pore volume (semi-quantitative). g: potentiometric titration. h: open-site rule."),
    F(r"$\mathrm{OS} = (12 - CN) + t\,x - \sum_i z_i\,n_i - n_{Cl}$"),
    P("For 12-connected nodes this reduces to OS = 2x − m − Cl."),
    H("What is checked"),
    B("The code reproduces published values (Valenzano 2011, Shearer 2014/2016, De Roo 2025) and recovers the compositions used to build the synthetic examples."),
    B("Methods f and b′ are semi-quantitative and are excluded from the consensus. The Klet titration is only checked for its conventions, not against raw data."),
    B("Extra-framework ZrO₂ in composites is removed from the TGA masses before methods a–d."),
)

DOCS["Figure composer"] = _doc(
    "Combined figure – several panels in one figure", "Panels from different techniques with journal templates.",
    H("What it does"),
    P("Puts graphs of different techniques into one multi-panel figure with panel letters."),
    *terms("Panel: one graph of the combined figure.", "Journal template: preset column width, font size and label style (ACS, RSC, Wiley)."),
    H("Tutorial"),
    N("Under any graph click 'Add to combined figure'."),
    N("Open 'Combined figure' and arrange the panels (1–4 columns). Choose the panel letters (A/B/C, a/b/c or a)) and their place."),
    N("Choose a journal template; the program checks the width and the minimum text size."),
    N("Export as vector graphics."),
    B("Panels use the default style of each figure; the composer has its own controls."),
)

DOCS["Reticon"] = _doc(
    "Reticon – schemes with lab and materials icons", "A drawing editor for graphical abstracts and method schemes.",
    H("What it does"),
    P("A scheme editor with about 800 icons (CC0, bundled) for lab equipment, materials, instruments and symbols. Everything is offline."),
    *terms("Scalable vector graphics (SVG), portable document format (PDF): vector export formats.",
           "Portable network graphics (PNG), tagged image file format (TIFF): picture export formats with a chosen resolution (dpi)."),
    H("Tutorial"),
    N("The editor opens with the module. Press the prominent 'Full screen' button (or Ctrl+Shift+M) to work in a large window."),
    N("Drag icons from the library onto the canvas. Widen the library by dragging its edge (double-click resets it)."),
    N("Add arrows (connect icons), text and boxes; align, group and distribute items; undo with Ctrl+Z."),
    N("Start from a template or from the content tab (methods schematic, your data, molecule, structure, combined figure)."),
    N("Export to SVG, PDF, PNG or TIFF at the width and dpi of your journal. Schemes are saved with the project file."),
)

DOCS["Batch processing"] = _doc(
    "Batch processing – one preset for many files", "A whole folder to a set of figures.",
    H("What it does"),
    P("Applies one set of figure settings to every file of a folder."),
    *terms("Preset: a saved set of figure settings (JSON) or the house style."),
    H("Tutorial"),
    N("Choose the technique (PXRD, FTIR, Raman, XPS, TGA, DSC, N₂ sorption or SEM/TEM particles) and a folder (sub-folders included), an upload or the example folder."),
    N("Choose the preset, the journal width and the layout (one figure per file, all overlaid or stacked)."),
    N("Run. Download the ZIP (PNG, TIFF, PDF, SVG and CSV) or write to an output folder. A log lists every file."),
)

DOCS["SI export bundle"] = _doc(
    "SI export bundle – everything for the supporting information", "Figures, data, parameters and tables in one ZIP.",
    H("What it does"),
    P("Collects every figure shown in the session into one ZIP."),
    *terms("Supporting information (SI): the extra material that goes with a paper."),
    H("Tutorial"),
    N("Keep 'Collect figures for the SI bundle' ticked in the sidebar while you work."),
    N("Open the SI export page and choose what to include."),
    N("Download the ZIP: figures (PDF, PNG, SVG), data CSVs, SI_parameters.xlsx (figure list, analysis and figure settings, results per sample) and the defect SI table (DOCX/XLSX)."),
    N("'Dossier per material' makes one ZIP per material (figures, data, methods, results). The cross-checks between techniques (for example crystallite size versus particle size) are listed here too."),
    B("Every number in the sample summary carries its page, data files, time and program version."),
)

DOCS["Sample database"] = _doc(
    "Sample database – your samples and their results", "A local list of sample IDs linked to protocols, results and data files.",
    H("What it does"),
    P("A local database (SQLite, in the Reticula data folder) of your samples."),
    *terms("Sample ID: the code of a sample; data files are linked when it appears in their names (<ID>-B<batch>_…).",
           "Digital object identifier (DOI): the permanent link of a source paper."),
    H("Tutorial"),
    N("Load your Excel sample list (a 'Sample ID' column; optional Name, Paper/Project, Protocol family, Parent / starting material, Protocol basis, Main source (DOI), protocol files)."),
    N("Point it at your data folders. Browse the samples with their protocol files, family, parent and children, source DOI, results from the sample store and linked data files; add notes."),
    B("A synthetic example list is bundled."),
)

DOCS["Export guide"] = _doc(
    "Export guide – which file to export from each instrument", "Open the 'Export guide' module for the full table.",
    H("What it does"),
    P("Lists for every instrument which file to export, which formats Reticula reads, and the menu path in the instrument software. It includes the JEOL JEM-1010 and JEM-2200FS electron microscopes and the Zetasizer."),
    *terms("Comma-separated values (CSV): a plain text table; the safest export format.",
           "Transmission electron microscope (TEM) and its high-resolution mode (HRTEM).",
           "Dynamic light scattering (DLS): the Zetasizer size measurement."),
    H("Tip"),
    B("If an import fails, export plain text as described there; two-column text always works."),
)
