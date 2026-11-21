# Changelog

## 2.18.8 - 2026-11-21

### Fixes (photoluminescence tab)
- Power dependence (`fit_power`): equal power values for all spectra gave `numpy.linalg.LinAlgError: Singular matrix` and blocked the page. The page now says it needs at least two different powers.
- Panel composer (`composer._group_boxes`), and the label and plot helpers that measure text: `fig.canvas.get_renderer()` does not exist on a plain canvas (`AttributeError`). A new helper `mfs/render.py` (`renderer_for`) falls back to an Agg renderer of the same size and resolution.
- Tests: `tests/test_v2187_hotfix.py` reproduces both crashes (they fail on 2.18.7 and pass now).

## 2.18.7 - 2026-11-20

### Code cleanup (no change to computed results)
- Duplicate, identical helper `_grid_for` removed from `mfs/indexing.py`; the indexing path keeps using the one remaining copy.
- Dead code removed: a no-op assignment in `mfs/pages/thermal.py`, unused imports (F401) and empty f-strings (F541) across `mfs/`.
- Local variable `S` renamed to `struct` in `mfs/pages/structure.py` so it no longer hides the module-level `S`. Text in strings is unchanged.
- `Image` imported only under `TYPE_CHECKING` in `mfs/reticon/qa.py` (and its runtime copy), so the linter reports no false error.
- `tests/test_v215_plextras.py`: `assert False` replaced by `raise AssertionError(...)`, because `assert` is removed under `python -O`.
- Release workflow: the example tag in the header comment of `release.yml` is now `vX.Y.Z`; `macos-build.yml` accepts the same version spacing as `release.yml`.

### Checks
- Unit tests: 652 passed, 50 skipped, 1 failed (`test_calibrate_synthetic_with_ocr[zeiss]`, the SEM scale-bar OCR test excluded from the gate by the user).
- Self-test (`check_selftest.py`): 41/41 ok, exit 0.
- Known-answer coverage of the indexing path: `tests/test_v216_indexing.py::test_known_structures_are_indexed` (known cells are recovered).

## 2.18.6 - 2026-11-20

### Photoluminescence: thermometry (ratiometric)
- **Poisson maximum-likelihood fit on raw counts** (new `mfs/plmle.py`): N1 ~ Poisson(B1 + D·S2), N2 ~ Poisson(B2 + S2), D = α·exp(−ΔE/k_B T). No background subtraction before the ratio, no logarithm, no dropped points. Uncertainty from the observed information (numerical Hessian, inverted), checked against the analytic Poisson Fisher bound.
- **Wired to the thermometry tab** (integration windows, counts on): the maximum-likelihood ΔE with its uncertainty is compared with the subtraction result in a new check "Poisson maximum likelihood (counts)".
- **Emission-masked baseline** (`mfs/plbackground.masked_poly`): when no signal-free window exists, the baseline is a polynomial of order 1 at most, fitted outside the emission (mask derived from the data by iteration). Default instead of the rubber band; the rubber band is an alternative. A safety check (overshoot above the spectrum by more than three noise units) falls back to the rubber band when the emission covers the whole range, and the notes say so.
- **Changed default (behaviour change):** without a signal-free window the program now uses the emission-masked polynomial instead of the rubber band. Measured on synthetic spectra (area error): sharp peaks on a slope rubber +78 % / masked +1 %; noisy spectra rubber +89 % / masked +0 %; broad band rubber +6 % / masked -2 %; concave background: both fail (+715 % / +485 %). A test of release 2.18.4 was adapted accordingly.
- **Analysis report (PDF)**: a third download next to the CSV and the JSON record (`mfs/plreport.py`): verdict, settings, equations, calibration figure, table, uncertainty budget, checks with advice, limitations.
- Help texts (`helpdocs.py`) updated for the three items.

### Measured (synthetic, known answer)
| Test | Result |
|---|---|
| Ratio, 300 runs, true ΔE 700 cm⁻¹: subtraction + ln-fit (released method) | bias −3,05 %; 33 % within 1σ, 53 % within 2σ |
| Ratio, 150 runs: Poisson maximum likelihood | bias +0,09 %; 70 % within 1σ, 95 % within 2σ; z-sd 1,06 |
| Broad emission band (area error): SNIP / arPLS / rubber band / masked polynomial order 1 | +13 % / −26 % / +5 % / −1,5 % |
| Masked polynomial with order 2 (not used) | −30 % (extrapolation into the band) |

### Known issues and open items
- The maximum-likelihood fit has no real-data validation yet: no online dataset with raw counts was usable (Sm²⁺ gated spectra are counts already background-corrected by the instrument; Yb³⁺ is OriginLab .opju and not readable; Pr³⁺ is normalised, without counts).
- Mn⁴⁺ (Science Advances table, Fig-3c): the table gives ΔE = 196–202 cm⁻¹ over every temperature range tested; the paper states 214 (fitted) and 226.5 cm⁻¹ (calculated from the ν₆ phonon energy). The cause is not identified.
- Dy³⁺ and YOF: no reference ΔE in the sources; the Boltzmann check correctly refuses a reliable gap for the Dy³⁺ ratio.
- Test suite (local): 653 passed, 0 failed before release (see the release notes).
- The PDF methodology document is not updated; it waits for approval.

## 2.18.5 - 2026-11-20

### Photoluminescence: ratiometric thermometry
- **Empirical noise when the counting error is too small.** When the reduced chi-squared of the Boltzmann fit is above 10, the counting error is no longer used as a weight. The fit is repeated with a constant relative error equal to the measured scatter of the ratio, and the report shows a "Ratio noise" warning. Before this change the error bars were scaled up by a factor of about 10^8 on the anti-Stokes test data.
- **Birge scaling only when significant.** The errors are multiplied by the square root of chi-squared only when chi-squared exceeds 1 by more than two standard deviations of chi-squared (threshold 1 + 2 sqrt(2/dof)).
- **Boltzmann-law check (R-squared).** R-squared below 0.9 fails, 0.9 to 0.98 warns, 0.98 and above passes. A ratio that does not follow the law is reported as having no reliable energy gap.
- **Simulation check renamed.** The known-answer self-check is now a simulation consistency check. It tests that the analysis recovers a gap it was given; it does not test accuracy.
- **Ratio error includes shape parameters.** The covariance of the line-shape parameters is propagated into the ratio error.
- **Fluoracle export fixes.** Trailing separators are accepted and rows without numbers are skipped.

### Validation (documented, not a release claim of accuracy)
- Visible LIR of the YOF data set reproduces an independent regression (0.0562 eV against 0.0558 eV).
- Mn4+ ratio table from the authors' supplement: the module gives 195.8 cm-1 (statistical error 1.8 cm-1) against 214 and 226.5 cm-1 quoted in the paper. The difference is not resolved.
- Cr3+ MOF relative sensitivity at 160 K: 2.79 %/K against 2.84 %/K in the paper.
- Dy3+ ratios fail the Boltzmann-law check; the source manuscript reports sensitivity peaking at 193 K, which a single gap cannot produce.
- The PDF "Reticula PL mathematics and implementation" is updated to match the method.

### Other
- Example sample names are masked: MOF-801(Hf) becomes Sample A, Sulfo-MOF-801(Hf) Sample B, ZrO2@Sulfo-MOF-801(Hf) Sample C. The example files are renamed to match.
- DSC examples and the DSC generator are unchanged from 2.18.4.
- Known open item: two defect-calculator snapshot tests (test_pagestate) depend on which example file is linked to which sample; decision pending.

## 2.18.4 — 2026-11-19

### Photoluminescence (visible and near-infrared)
- **Edinburgh Instruments `.FS` files are read** (Fluoracle software): single scans and "Multiple Scans" files, one spectrum per scan, with sample name, detector, excitation source and the temperature taken from the sample name. Checked on 9 real files (900–1250 nm, 351 points, near-infrared photomultiplier tube).
- **Standard background subtraction** for the whole series (new setting *Background*, default **Automatic**): signal-free windows are found over all spectra, a polynomial (order 0–2, chosen by the Bayesian information criterion, a higher order must earn its parameter) is fitted in those windows only, and the same plan is subtracted from every spectrum. The block shows method, windows, mean background and what is left in the windows (should be 0 ± 1 noise level), and a candidate table ("Why this method?").
- Manual control: your own windows with polynomial, rubber band, ModPoly, arPLS or SNIP; a measured dark / background file (scaled in the windows); constant pedestal options; optional removal of single-point spikes. Without any signal-free window the cautious **rubber band** (lower convex hull) is used.
- **The ratio fit uses the background-subtracted spectra** (Poisson weights still use the raw counts), and **the choice of background is part of the ratio's systematic error**: the series is refitted with the other reasonable background plans (other polynomial orders, shortened windows; or the other flexible methods) and the spread is added, shown in its own column ("of which background choice"). On the real NIR test series this part is about 0.4 % of the ratio. A file with 10 000 counts subtracted gives exactly the same corrected spectra as the original (maximum difference 0.0 counts).
- **Confidence verdict and quality report** above the thermometry results (HIGH / MODERATE / LOW, every check with its reason and what to do): saturation, settings that differ, whether the ratio carries temperature information at all, monotonic ratio, band-fit quality, Boltzmann residuals and runs test, ΔE against the spacing of the bands and the usual coupled-level range, an independent method (windows against fit), heating against cooling, reproducibility between runs, the self-check.
- **Manifold fit** for crowded, Stark-split bands: every line its own shape, the lines in two dragged ranges are summed into I₁ and I₂ (full covariance in the error). Proposed automatically when the two-band fit fails on more than 30 % of the spectra. On a real lanthanide NIR test series the two-band fit flagged every spectrum; the manifold fit describes all seven (χ² ≈ 1).
- **Uncertainty budget of ΔE**: Monte Carlo over the ratio errors, temperature-sensor random error and offset (typed in), a Boltzmann refit on every processing alternative (fit range / margin, peak model, background, windows), reproducibility between runs, Birge scaling. Two temperature uncertainties: resolution (one reading) and calibration (accuracy of T(Δ)); the temperature read back from every ratio is in the table.
- Integration windows: counting error from the raw counts; integrated above zero after the standard background subtraction (the straight line through the window ends cut into neighbouring lines and could even give negative areas). The ratio is turned the right way round automatically when it falls with temperature.
- **Laser-heating check** (spectra at one temperature, several powers → heating per watt), **heating / cooling hysteresis**, **known-answer self-check** (series simulated from your own fit with a known ΔE, analysed again: bias and coverage of the error bars), **analysis record** (JSON with settings, file fingerprints, results, budget and checks).
- **Decay curves (lifetimes)**: tail fit, rise-and-decay, reconvolution with a measured instrument response; 1–3 exponentials chosen by the information criterion; amplitude- and intensity-weighted averages; τ(T) with Eₐ and Sᵣ. Checked on simulated decays with known lifetimes.
- Spectral-response correction from a calibration lamp (measured + certified spectrum), coverage warning; removal of the excitation line and its second order; settings read from Fluoracle text exports and Ocean Optics files (excitation, dwell time, repeats); normalisation per second; data-check panel with saturation (photon-counting detectors exempt from converter limits); labels of typical lanthanide / Cr³⁺ transitions; CIE dominant wavelength and colour purity; thermal-quenching Eₐ and T₅₀ with errors and warnings for few points on the drop or an extrapolated T₅₀.
- The Photoluminescence info document (Help and info and the PDF button) explains the file reading, the background steps and what else influences ratiometric thermometry.

### Defect calculator (thermogravimetric analysis)
- **Lázaro et al. 2020 multi-faceted method** (methods i1 and i2): per-metal residue ratio R_exp, charge balance with OH⁻ / O²⁻ / Cl⁻, partial decomposition residue fraction, and composites (inert component, component burning before the linker, component burning with the linker). Validated against the worked examples of the paper's supporting information. The info document is rewritten.

### PXRD
- **Structure solution of an unknown compound** (new expander *Solve the crystal structure of an unknown compound*): the route of the field as one guided, automatic decision tree with a progress bar and a log of every decision. Data check → unit cell (indexing, then a Le Bail fit of each candidate on the real pattern; lowest Rwp, smallest volume among equals) → space group from the systematic absences (matched-filter intensities of isolated reflections; the 230 groups, orthorhombic groups in all six axis settings) → **charge flipping** with the space group imposed (E-values, re-partitioning of overlapping reflections, many random starts in parallel threads) → atoms (symmetry orbits, elements from your formula, number of formula units and positions chosen by the R factor of the amplitudes, difference-Fourier completion, position refinement) → Rietveld refinement of the model. A silicon (or other standard) pattern gives the instrument zero shift and the fixed peak profile. Outputs: CIF, 3D structure, fit, `.xy` files, reflection lists, a decision report, all in one ZIP. Checked on synthetic powder patterns of public-domain COD structures (SrTiO₃, rutile, ZnO, anatase, tetragonal ZrO₂, SnO₂ recovered); harder cases are less reliable, see the limits in the info document.
- Indexing: the centring that leaves fewer observed lines unindexed now wins (a centring that forbids an observed line is wrong however good its M20 over the other lines looks); this fixed rutile being reported as a four-times larger body-centred cell.
- **Rietveld / Le Bail results can be exported by hand**: 3D structure viewer in the refined cell, CIF (the atoms are not refined, stated in the file), `.xy` files for observed / calculated / difference / background, reflection list and `.hkl`, and one ZIP.
- Rietveld: glass-holder background (automatic or from a loaded background file), and the instrument profile and zero shift can be fixed from a measured Si standard.
- **Indexing** rewritten for a higher M20 (de Wolff figure of merit): linear refinement of the reciprocal metric with zero-shift, centring-aware seed ranking, staged search per crystal system, robust M20 over indexed lines.

## 2.18.3 — 2026-11-19

### PXRD: no assumed material
- The reference pattern now starts at *none*, the metal in the simulation at *As in the CIF (no substitution)* (Hf replacement is an option), and no phase is preselected in Rietveld. The label no longer adds `(Zr)` by itself.

### Help and info
- New module **Help and info**: an info document for every module (what it does, the terms written out before they are abbreviated, a tutorial, the formulas, what the program checks, the limits). Readable in the app, searchable, and downloadable as a PDF in the Reticula colours (one module, or all as a ZIP).
- Every module has a **📖 Info about this module (PDF)** button at the top.
- The Photoluminescence document explains ratiometric thermometry, including the background polynomial (default order 1, range 0–3).

## 2.18.2 — 2026-11-19

### Zetasizer (.dts)
- **Software 8.x .dts files are read**: sample names, the intensity size distribution, and **Z-average and PDI as stored by the software** for every record. Checked on two real files (two DLS files, software 8.02): all six records equal the values the Zetasizer software shows (167.0 / 166.0 / 166.3 nm and PDI 0.139 / 0.133 / 0.126; 2238 / 2363 / 2219 nm and PDI 0.155 / 0.527 / 0.631).
- **Zeta potential records are read** (software 8.x): ζ-potential, zeta deviation, conductivity, count rate, the three zeta peaks (mean, area %, SD) and the zeta distribution. Checked on 11 records of a real file against the Zeta Potential Reports of the software (ζ, deviation, conductivity, peaks and areas equal). The Zeta / DLS page adds conductivity, count rate and deviation to its table (and so to the SI data) and draws the ζ distribution. A distribution is only shown for records whose stored block reproduces the software's ζ (8 of 11 in the test file).
- Software 7.11 files keep the earlier route (Z-average and PDI calculated from the stored correlation function).

### Reticon (schemes)
- **The scheme editor loads when the page opens** (no extra "open the editor" step).
- **Full screen** is the prominent button at the start of the top bar (accent colour).
- In full screen the top bar sits below the screen edge; the handles follow it.
- The icon library can be widened by dragging its edge, in the browser and in full screen; entering full screen starts it at least 340 px wide (a wider width you chose stays).

### FTIR / Raman
- **Several spectra in one text file**: tick *The file(s) contain several spectra* and every column after the wavenumber column becomes its own spectrum, named by its column header. Without the tick, the behaviour is unchanged.
- **OMNIC `.spa` files are read** (one spectrum per file; series files `.spg` are not). Verified on a real file against its OMNIC CSV export: same points, same intensities (x within 0.001 cm⁻¹). The full measured range is read, also beyond the 400–4000 cm⁻¹ of the CSV export.
- Header-less OMNIC CSV exports (scientific notation, e.g. `3.996813e+002`) are read as one spectrum.

### Methods and export guide
- The SEM / TEM / EDX Methods paragraph names the instruments: **JEOL JEM-1010** (bright-field TEM, 100 kV, tungsten filament, CCD camera, 0.4 nm point resolution, up to 500 000×, dual-sample holder) and **JEOL JEM-2200FS** (HRTEM, 200 kV FEG, image-lens Cs corrector, EELS and EDX). SEM stays a placeholder.
- `docs/EXPORT_GUIDE.md`: separate rows for the JEM-1010 (camera TIFF) and the JEM-2200FS (raw Digital Micrograph / Velox files).

### Tagline
- The tagline is now "From raw data to figures, the easy way." (app, launcher, installer, README).

### Simple mode (every technique)
- The *Simple mode – all black* switch in the Style tab of every figure is redrawn for a sleek journal look: data points get one black symbol per sample, continuous curves (spectra, patterns, fits) one black line style per sample and no markers, scatter groups filled / open, bars and histograms light grey with a hatch per series, filled bands light grey.
- Legends that already exist follow the same look (they used to keep the colours).
- Not covered: CIF structure views and the Hirshfeld map (their colours carry the information).

### Rietveld: instrument profile from a standard
- *Instrument profile → fixed from a standard* (Advanced): upload a Si SRM 640 pattern measured on the same instrument. The standard is fitted by Le Bail with its certificate cell (a = 5.4312 Å, not refined); U, V, W, X, Y are then kept fixed for the sample refinement.
- Checked on synthetic data: the Si cell stays at the certificate; a Ge-type sample (a = 5.658 Å) refines to 5.6580 Å with the fixed profile, the same Rwp as with a refined profile. The fitted profile reproduces the peak width up to about 65° 2θ; above 80° the Lorentzian term is less well determined. Check the standard's width curve against your own data.
- LaB₆ is not built in yet (its cell can be given as a Phase; see `instrument_profile_from_standard`).

### Tests
- The Word (DOCX) export test no longer needs python-docx (the app writes DOCX without it); the Reticon tests follow the editor loading on page open.

## 2.18.1 — 2026-11-19 (quick fix)

- **The Schemes & Icons module is now called Reticon** (menu, page title, figure composer label, gallery tile).
- **Icon library resizable**: drag the edge next to the library to widen it; double-click resets it. Wider library shows more
  icons per row (it used to show about two and a half).
- **Full screen** for the Reticon editor: the *Full screen* button (or Ctrl+Shift+M) makes the editor fill the window, with
  larger tile labels and buttons. The browser's full-screen mode is used where the host allows it.

## 2.18.0 — 2026-11-19

### Schemes & Icons
- **500 new vector icons** (1054 in total). Each new icon has the five standard variants, so 6 shapes per icon: 3000 shapes
  in all with the variants. Families: spectra and analysis traces (150), lab vessels (120), instrument bodies (96),
  material morphologies (96), process arrows (38). All are original vector drawings, CC0, no images, masks or text.
  Source: `tools/icon_sources/reticula_v218.py`.
- The starting template chosen in *Start from* is shown as a preview (SVG). Hover previews are not possible in Streamlit.
- Editor work since 2.17.0: floating contextual bar above the selection, status line, Enter / Shift+Enter to add the next step,
  Repeat →, Number, Line up, Page ↔ / Page ↕, Copy style / Paste style, Copy PNG, focus mode (F), template start points,
  locked items, undo toast, draft kept in the browser, resizable right panel, connection points on hover.
- Fixed in the editor: flip could not be undone; typing in the side panel lost focus; the align bar showed without a selection.

### Thermal analysis (TGA)
- New plateau rule *start of decomposition (EJIC 2020)*, next to the flattest-stretch rule. It takes the knee of the mass curve
  (1 % below the pre-decomposition line, fitted between 150 °C and 320 °C). Based on Lázaro et al., Eur. J. Inorg. Chem. 2020,
  doi:10.1002/ejic.202000656. A point picked by hand on the graph still overrides either rule.
- Example: UiO-66 TGA (Leuven, CC-BY 4.0, Zenodo 15013202) in `example_data/TGA`.
- Known limit: the knee rule has not been checked against hand-picked onsets on real data; the example shows about 5–10 °C
  from the visible bend on UiO-66.

### Nuclear magnetic resonance (1H NMR)
- *Simple mode (black)* checkbox on the figure tab: the spectrum, labels, integrals, marks, the molecule inset, and the SI page
  PDF are drawn in black only.

### Tests and checks
- The batch-run test counts the example TGA files from the folder instead of a fixed number.
- New tests: icon library (counts, five variants, vector-only), EJIC onset rule (synthetic and real UiO-66 data), NMR simple mode.

### Known limits
- Zetasizer 8.02 `.dts` files are still not read (waiting for a sample file).
- Icons from the toolkit's `variants.py` need `rsvg-convert`; the app uses its own `raster.py` version instead (see `tools/build_icons.py`).

## 2.17.0 — 2026-11-18

### Name and project files: Reticula
- The program is called **Reticula** in the window, the start page, the menus, the exports, the methods texts and the provenance.
- **Project files are `.reticula`** (zip format). Projects from v2.16 and earlier open as before.
- Double-click or *Open with* works for `.reticula` files (Windows installer and macOS bundle register the new extension too).
- **Start page: continue or save a project.** *Continue a project* opens a `.reticula` file directly from the home page; *Save as .reticula* saves the whole work from there as well.
- Not yet renamed (technical, to be done with the release): the installer file name `Reticula_Setup`, the macOS app name and bundle identifier, the data folder `Documents\Reticula`, the GitHub release names and the update check. These keep working; changing them would break existing installs and the update check.

### Across the program
- **Tabs more visible**: every technique's tab bar has the app colours, bold labels and a clear bar.
- **Simple mode** (Style tab): full black graph; samples are told apart by their own point symbol.
- **SI export – journal style**: choose ACS, Wiley-VCH, RSC or Elsevier for the Methods heading and the software sentence. These are templates; check the author guidelines of your journal.
- **Progress with time left**: the Rietveld refinement shows elapsed time and an estimate of the remaining time.

### NMR
- *Unassigned & impurities* is the first figure; the molecule picture no longer overlaps the spectrum.
- Protons are labelled with numbers (literature-style locants) instead of letters.

### Rietveld / structures
- **ICSD import**: manual only. A button opens the ICSD login page; the query is yours (queries are paid), Reticula sends nothing.
- **CIF picture**: centred on the molecule; labels have a white background, so they do not run through the molecule.

### Defect calculator (TGA)
- **Plateau point** for the defect calculation is chosen automatically at the flattest part of the mass-loss curve between the last low-temperature step and the linker combustion, or picked on the graph (double-click resets).
- The results table says how the point was chosen: *automatic (flat)*, *picked on graph* or *automatic – no clear plateau*.
- The plateau values of Sample A and the consensus changed by a few degrees; the baseline test was updated for these two keys only.

### Zetasizer `.dts` (known limits)
- Files of software version **8.02** are read: sample names and the intensity size distribution.
- **Z-average and PDI are not calculated for 8.02 yet.** On the 8.02 test files the calculation gives impossible values (negative PDI, Z-averages that disagree between runs). They are only calculated for 7.11. For 8.02 use the Excel / CSV export from the software.
- ζ-potential is still not read from `.dts`; use the export.

### Schemes & Icons (new module, 🎨 in the menu)
- **Scheme editor** with the Reticon icon library (about 800 CC0 icons, bundled, offline): search in EN/NL/FR/DE, icon families (6 variants per icon), drag and drop, snapping, arrows that stay attached to icon anchors (straight, elbow or curved, both heads, dashed, label), text, panel boxes, undo / redo.
- **Multi-select**: drag on empty canvas or Shift-click; move, delete, duplicate and nudge the whole selection.
- **Align** (left, centre, right, top, middle, bottom), **distribute** (horizontally or vertically) and **group / ungroup** (clicking one grouped item selects the group).
- **Copy / paste** (Ctrl+C / Ctrl+V, Ctrl+D duplicates): also between schemes.
- **Arrow ends**: drag an end dot onto another icon's anchor to re-attach the arrow.
- **Subscript and superscript** in text: `H_{2}O`, `CO^{2}`; the same in the export.
- **Recently used** icons at the top of the picker; **keyboard shortcut list** with the ? key.
- **Export** of the scheme to PDF (vector), SVG, PNG or TIFF at the width of a single column, 1.5 column, double column or TOC graphic, at 300, 600 or 1200 dpi; warnings for thin lines and non-vector icons.
- Schemes are saved with the project file. Start from a template or an empty canvas, or import a scheme as JSON.

### Schemes & Icons: arrows, content from Reticula and the combined figure
- **Arrow bends**: add a bend, drag a bend to move it, or **Route around icons** (an automatic route that avoids the other icons, stored as bends so the export matches the screen).
- **Embedded figures**: a plot of your own data (x, y text file), a molecule from a SMILES, or a crystal structure from a CIF become vector icons in a scheme.
- **Methods schematic**: one click makes a flow from a sample and the techniques it went through.
- **Scheme in the combined figure**: adds the scheme as a panel next to the graphs of the other techniques.
- **Style and theme per icon** (Advanced/Developer layers of the export): Style and Theme in the inspector.
- **Periodic table view** of the 118 element tiles in the picker.

### Known limits (summary)
- Accepted as they are for this version: the TGA plateau rule (SI tables not checked against it) and the journal SI wording (templates, not checked against current author guidelines).
- Zetasizer 8.02: Z-average and PDI not calculated (use the software's export); ζ not read from .dts.
- Browser tests of the scheme editor and of 8 older browser tests are not run here (Playwright missing). PyInstaller build not run with these changes.
- Full to-do list for the next version: TODO_v2.19.md.

- The browser preview shows each icon in its default style; the export (PDF, SVG, PNG, TIFF) uses the chosen style and theme.
- The scheme panel in the combined figure is an image at 600 dpi, so that one panel is not vector.
- The automatic route works on an 8 px grid around the icon boxes; check the result after routing.
- The editor was tested by hand in a browser; automated browser tests are not in the repository yet.
- The PyInstaller build has not been run with these changes.
- TGA plateau: the Supporting Information table was not checked against the rule (see the defect calculator section).
- Zetasizer 8.02: Z-average and PDI are not calculated for 8.02 yet (see the .dts section).

## 2.16.0 — 2026-11-18

### PXRD: from a pattern to a refined, checked structure (almost automatic)
- **✨ One click**: *Do everything automatically* searches a library of 18 reference phases (MOFs, oxides, standards, minerals; COD) with a small search over the cell size and the zero shift, proposes the phases with a confidence, and refines them – cell, peak shape, background, atom positions, displacement parameters.
- **Structure refinement** (*Rietveld + refine the atoms*): symmetry-aware (space group found from the structure, atoms stay on their special positions), one displacement parameter per element, soft bond restraints that keep linkers chemically sensible (strength adjustable), esd from the covariance. A reflection is computed once per symmetry class, so a 328-atom MOF refines in under 10 s. Validated on two **real laboratory patterns** (fluorapatite: atom positions within 0.0007 of the GSAS result, Rwp 9.1 → 8.6 %; PbSO₄: Rwp 9.1 → 8.7 %, three different starting CIFs reach the same result).
- **Statistics of the model**: bond lengths and angles with esd, coordination numbers, bond-valence sums, a plain-language verdict (BVS far from the formal charge, too short bonds).
- **3D viewer** (own, no external libraries, works offline): rotate, zoom, click an atom for its bonds and coordination, views along a / b / c, thermal spheres, labels, 1×1×1 or 2×2×2 cells, PNG.
- **Publication output**: refined **CIF with esd** (cell, atoms, symmetry operations, U_iso, bonds, angles, figures of merit), tables (Excel), a structure figure (view, angles, atom size, labels) that also goes into the combined figure, and a Methods paragraph.
- **Unknown phase**: *indexing* finds the unit cell from the peak positions (cubic, tetragonal, hexagonal, orthorhombic; de Wolff M20, physical limits) and a Le Bail fit of the cell checks it. Monoclinic / triclinic are not searched (stated in the page).

### Across the program
- **Home page redesigned**: dark hero, "what do you want to do?" cards, the five steps, a coloured tile per tool with what you get, the three detail levels (switch on the page), the gallery.
- **Sample summary with sources**: every number carries its page, data files, time and program version (Excel with a *Sources* sheet, LaTeX with footnotes); also in each material dossier.
- **Cross-checks between techniques** (SI export page): crystallite size vs particle size, BET vs the external area of the particles, defect methods that disagree, surface (XPS) vs bulk (ICP / XRF) metal ratio, DLS vs image size – with what it can mean and what to check.
- **Problem report** in the sidebar: a small ZIP with the version, settings and the last error – no data – to send to the developer.

## 2.15.1 — 2026-11-18

### PXRD / Rietveld
- **Preferred orientation** (March–Dollase): give the axis (h k l) per phase under *Advanced → Preferred orientation*; the factor r (r < 1 plates, r > 1 needles) is refined after a scan that finds the right basin (the objective has two). Recovered on a synthetic texture within 0.05.
- Zero shift and sample displacement are no longer refined together by default (they are strongly correlated and bias the cell, 0.06 % on a real pattern): displacement only (reflection geometry), option *Instrument offset*. Validated on a real laboratory pattern (fluorapatite, Cu Kα1+Kα2): Rwp 9 %, GoF 1.6, cell within 0.003 % of the GSAS result.

### Photoluminescence (Advanced)
- **Spectral-response correction** file (wavelength, factor) applied to every spectrum.
- **Temperature uncertainty from repeat measurements**: σT from the scatter between spectra at the same temperature (the experimental temperature resolution).
- **Power dependence**: power read from the name (`250mW`) or typed, slope n of log I vs log P ± error, saturation warning, log–log figure.
- **CIE 1931 colour**: x, y and colour temperature of every spectrum and a chromaticity diagram.

### Dossier per material
- SI export page: a *material* column (guessed when a graph shows one sample) and **one ZIP per material**: figures, plotted data in one Excel, settings, methods text, the results of that material from the sample store; or all materials in one ZIP with a folder each.

### Interface
- **Readable sidebar**: the expanders of the sidebar (figure style, colour rules, house style, About, project, saved work) were white with pale text on the dark sidebar; they now have a dark background with light text.
- **Detail level** is now a coloured box at the TOP of the sidebar with three buttons (🟢 Standard, 🔵 Advanced, 🟣 Developer) and the colour of the current level; the override for the current module sits right under it.
- Standard level of *Customise this figure*: title size, sub/superscript rule and font sizes moved to Advanced (they keep their defaults).

### Build and tests
- macOS: pyarrow's `libarrow_*` libraries are no longer removed from the app (the self-test of v2.15.0 failed on this); a manual `macos-build` workflow attaches the dmg to an existing release; `SHA256SUMS.txt` without the build path.
- Real-browser tests updated (composer button / page names, sample-list shift-click timing); the old NMR drag test is skipped until it is rewritten for the new layout.

## 2.15.0 — 2026-11-17

### New: detail levels (Standard / Advanced / Developer)
- The app starts in **Standard**: only what a publishable figure needs. **Advanced** adds the details (ticks and presets of *Customise this figure*, fit options, series styles, libraries, quality checks), **Developer** shows everything.
- One switch in the sidebar for all modules, and an override per module (e.g. Standard everywhere, Advanced for PXRD). The choice is saved in the project file.

### New: PXRD – crystallite size and Rietveld / Le Bail (both hidden until you switch them on)
- **Crystallite size and microstrain**: automatic peak detection and pseudo-Voigt fits of every sample, instrumental broadening removed (typical value, a measured LaB₆ / Si standard, or a value you know), Scherrer per peak (mean ± SD), **Williamson–Hall** and **Halder–Wagner** with uncertainties from weighted regression, W–H figure, plain-language warnings (too few peaks, negative strain, size above the limit of a laboratory diffractometer).
- Rietveld on a REAL laboratory pattern (fluorapatite, Cu Kα1+Kα2, fully automatic): Rwp 9 %, GoF 1.6, cell within 0.003 % of the GSAS result. Zero shift and sample displacement are strongly correlated, so by default only the displacement is refined (option *Instrument offset*).
- **Rietveld / Le Bail refinement** (own implementation): pick the phases (bundled CIFs or your own, any number) and press *Refine*. Automatic staged strategy (start-value search for cell scale and peak width, then cell, zero shift, peak shape, displacement, asymmetry), scale factors and background solved exactly in every step, about **one second per pattern**. Gives Rwp, Rp, GoF, Bragg R, cell parameters ± esd, weight fractions (Hill & Howard), the standard Rietveld figure (observed, calculated, background, difference, reflection ticks), an Excel report and results into the sample store. Validated on synthetic patterns (cell within the esd, weight fractions within 0.5 %, nanocrystalline broad peaks, Kα1 + Kα2) and against the structure factors of pymatgen. Not included yet: preferred orientation, atom positions, anisotropic broadening (stated in the report).

### New: reproducibility and journals
- **Provenance in every figure file**: program version, page, source files, figure settings and time are embedded in PNG, TIFF, PDF and SVG; the ZIP and *Save to folder* add `<figure>_provenance.json`.
- **SI export bundle**: all figure data in **one Excel** (a sheet per figure + index), a **journal check** (final width, smallest font and thinnest line at the printed size, colour-blind distinguishability, resolution) with **one-click fix** (width, fonts, lines, colour-blind safe palette), **Methods text** for every technique (placeholders in [brackets] for what only the author knows), and a **sample summary table** (one row per sample from the sample store; Excel, CSV, LaTeX).

### Photoluminescence: robust band ratio (thermometry)
- **Automatic global fit** (default): all spectra of the series are described by the same bands (pseudo-Voigt; position and width may drift with temperature); areas and background of every spectrum are solved exactly. Hidden shoulders of overlapping bands are found (second derivative), extra "nuisance" bands are added automatically when the residual shows a clear structure and the BIC improves, and the shared-shape and drifting-shape models are compared with the BIC. A single dirty spectrum does not bend the result (robust loss).
- **Honest uncertainties**: the ratio error includes the covariance of the two areas and a systematic part from refits with a narrower / wider fit range and the other model; the Boltzmann fit is weighted with these errors and gives ΔE ± error; σT follows from the ratio error. Quality flag per spectrum (good / check / poor from R², reduced χ², residual correlation) with a plain warning.
- **Select ranges on the graph**: drag over the mean spectrum and press *→ band 1 / band 2 / fit range*; the method *Integration windows* and *Compare both* (shows how much the windows deviate from the fit) remain available.
- On synthetic series with overlapping bands, a shoulder and a weak neighbour the fit recovers the ratio within 2 %, where valley windows are off by > 8 %.

### Modules
- **FTIR / Raman**: band assignment from an editable list of typical ranges (your own rows first), simultaneous multi-band pseudo-Voigt fit with uncertainties and a figure.
- **XPS**: batch fit – the same component model applied to every spectrum of a region, table + stacked-share figure. (UV/Vis already fits the band gap of all selected spectra.)
- **SEM / TEM particle sizer**: sizing **templates** (save the detection settings per sample type, apply with one click, kept between sessions and in the project) and a confidence bar per image (high / medium / low).
- **¹H NMR**: library of linkers, modulators and solvents (typical shifts in D₂O/NaOD, DMSO-d₆, CDCl₃; editable) that matches the found multiplets and gives molar ratios.

### Combined figure
- Panels no longer overlap: the gaps grow automatically until colour bars, axis titles and panel letters of neighbouring panels are clear of each other.
- New *Fine-tune the layout*: spacing between columns / rows, wide panels (span all columns), relative column widths and row heights, panel-letter size, automatic overlap check on / off.

### ¹H NMR
- **nmrML** files (`.nmrML`) can be uploaded: with an FID they go through the same automatic processing as Bruker data, a stored processed spectrum is used as is.
- Assignment structure: the proton chips (Hₐ, H_b …) sit **beside** their atom with a dotted leader line instead of on top of it, do not overlap each other, are larger and have readable text on every colour.

### Development
- `docs/DEVELOPMENT.md`: every new page must go through `show()`, be added to the figure-controls suite and have a methods text.

## 2.14.0 — 2026-11-18

### New: photoluminescence page
- Spectra from any spectrofluorometer (CSV / TXT / DAT / Excel; two columns or one table with a column per temperature; the temperature is read from the header or the file name).
- Temperature-coloured overlay with colour bar, waterfall, emission map; x axis in nm / eV / cm⁻¹ (with the Jacobian); baseline, smoothing, normalisation; peak position and FWHM.
- **Thermal quenching**: Eₐ and T₅₀ from I(T) = I₀/(1 + A·exp(−Eₐ/k_BT)), Arrhenius plot, peak shift and width against T (A–D figure).
- **Ratiometric thermometry**: Δ = I₁/I₂ = α·exp(−ΔE/k_BT), relative sensitivity Sᵣ, absolute sensitivity, temperature uncertainty σT (photon counts); emission map, Δ vs T with fit, Sᵣ(T), σT(T) as an A–D figure.

### New: combined figure
- Every graph has a **🧩 Add to combined figure** button; the sidebar shows the panels ("Open and arrange"). Pick graphs of different techniques, automatic layout, panel letters A/B/C (inside or outside),
  caption under the figure, long legends hidden, the per-figure settings (axis titles, ranges, units) travel with the panel.

### ¹H NMR
- Protons are labelled with **letters** (Hₐ, H_b …, most downfield first) in the structure, the chips, the spectrum and the tables - no longer confusable with carbon numbers.
- Figure: smaller labels with computed head-room (no clipping), integral numbers in their own band under the baseline (option to hide the numbers), thinner line.
- Page: the structure editor is tall enough for all its tools, with a side panel (SMILES, use / clear, normalisation); the other sections are tabs.
- Analysis: ¹³C satellites are recognised, nH ± uncertainty, a flag for integrals that are not integers (overlap / exchange), and a reliability column in the assignment table (ok / check / ambiguous / doubtful).

### Zeta / DLS
- **Runs can be switched off**: a table of all records (runs) with a tick per run; the averages (Z-average, PDI, ζ) and the graphs follow. `.dts` names like `Sample 3` become sample `Sample` + run 3.
- Excel (.xlsx / .xls) exports as well as CSV, and `.dts` files (sample names, intensity size distribution, and Z-average / PDI calculated from the stored correlation function – checked on Zetasizer software 7.11; ζ needs the Excel / CSV export); volume and number distributions; Z-average, PDI and ζ as movable texts in the DLS figure; quality notes; generated Methods / SI text.

### All figures
- **Font family** per figure and for the legend, **legend-title size**; legend font size, frame and title now work without ticking *Customise the legend* (the font size also sizes direct labels).
- A log axis without positive data no longer breaks the figure.
- **New test suite** (`tests/test_figure_controls.py`): every control of *Customise this figure* is applied to the real figures of all 16 pages and every widget is set in the running page.

### Particle sizer
- **Parallel**: images are segmented at the same time on all CPU cores (16 images: 60 s → 11 s, identical results), with real progress; falls back to serial.
- Edge of non-round particles on the same calibrated edge as the circle fits (mask / circle diameter ratio 1.12 → 1.01).
- Elongated particles (rods) are no longer cut into pieces or flagged as aggregates (synthetic rods found: 32/72 → 69/72); cubes and spheres are unaffected.
- Warning for low working resolutions (the sphere model switches off below ~1000 px).

### N₂ / BET
- The page starts with the BET surface area; the Rouquerol consistency, the comparison with the file's own value, the BET / Rouquerol plots and the t-plot are behind collapsed blocks.

### Code and build
- Windows: `run_windows.bat` picks Python 3.11-3.13; the frozen app's self-test checks the parallel particle sizer.
- Release workflow no longer waits for macOS; Mac build kit (`build_mac.command`) in every release.

## 2.13.0 — 2026-11-17

Built by Stijn Germonpré.

### New: ¹H NMR page
- Bruker raw import (Avance II + NEO, TopSpin 2-4; fid or processed 1r; JCAMP-DX), automatic phase / baseline / referencing.
- Offline structure editor (Ketcher); ¹H shift prediction from an nmrshiftdb2 HOSE-code table (held-out MAE 0.22 ppm, median 0.09 ppm).
- Automatic peak → atom assignment with drag-to-correct, impurity check (Fulmer 2010 + Gottlieb tables), labelled figure, qNMR, SI export (zip / PDF).
- The former digestion page is kept as the *digestion* sub-mode (linker : modulator ratios feed the defect calculator).
- ¹³C / DEPT import exists in the code but is hidden in the interface.

### N₂ / BET
- Redesigned graphs: isotherm, BET / Rouquerol, t-plot, PSD with a cumulative axis (the PSD inset was removed: no graphs in graphs).
- **BET range sensitivity**: the table shows the range of areas over all fitting windows that were evaluated (only the valid ones when any satisfies all
  Rouquerol criteria). `analysis.bet_range_spread()`.
- The Rouquerol window search is vectorised: ~11× faster (0.44 s → 0.04 s for 150 points) with identical results (tested against the previous
  implementation on the examples and 40 random isotherms).

### For people outside the original lab
- **No lab-specific defaults any more**: sample colours follow the palette in list order; name → colour rules are optional (sidebar *Sample → colour
  rules*, with an *Example-set rules* button for the bundled data). A file name is only turned into a display name for the bundled example files.
- **Sample database** reads any Excel sample list (**Sample ID** column; other columns optional). The bundled list is a synthetic example; the author's
  synthesis pipeline is no longer part of the program.
- ICP / XRF: molar ratios ± σ for **any elements** with a chosen reference metal (the Zr/Hf/S special case is unchanged).
- Requirements list completed (rdkit, streamlit-ketcher, pytesseract, streamlit-image-coordinates, h5py).

### Code and build
- `app.py` (3,800 lines) split: the page code now lives in `mfs/pages/*.py`, moved verbatim; tests check that every name a page takes from the app is declared.
- One cross-platform launcher; Windows installer and macOS dmg are built by GitHub Actions on a tag (tests → build → self-test of the built app → release).
- The scheduled 30-minute macOS retry workflow is gone. The release no longer waits for a macOS runner: the Windows installer is published first and the macOS `.dmg` is attached later; a Mac build kit (`build_mac.command`, `BUILD_ON_MAC.txt`) is part of every release.

### Fixes before release
- ¹H NMR page: just opening it no longer counts as an "unsaved change" (the structure editor echoed the structure back as a
  molfile), so switching to another page no longer stops at the *Unsaved work* dialog. Home tile and dialog use the label "¹H NMR".
- Windows: text files (impurity table, example structures, defect inputs, tool scripts) are read and written as UTF-8 instead of
  the system code page (cp1252 broke on characters such as ⁻).
- Windows build: PyInstaller with Python 3.12.0 breaks scipy.stats in the frozen app (15 of 40 self-tests failed); the spec now
  refuses 3.12.0. Inno Setup script saved with a UTF-8 BOM again (author name). Only the program files go into the build
  (no tests, tools, .venv or exports from a developer folder).
- Example-figure scripts find a sans-serif font on Windows and macOS too (were Linux-only paths).
- Tests: private-data skip also matches `SampleASGA2_TriStar_full_report_REAL.xls`; start-up smoke test talks to Streamlit 1.67
  with `websockets` (tornado is no longer installed by Streamlit).

## 2.12.0 — 2026-11-16

Built by Stijn Germonpré.

- **Particle sizer — sphere model** (auto-detected from the image's own isolated particles; Particle shape: auto / spheres / off):
  - outlines sit on the real edge (sub-pixel, calibrated on manually sized SiO2-NH2 TEM: O10 pooled 211 ± 23 nm vs manual 209 ± 19 nm;
    v2.10 overestimated by ~+30 nm at low magnification);
  - touching doubles / triples are split by circle fits (dotted outlines);
  - occlusion: an overlapping particle keeps its full diameter (darker overlap lens);
  - each fit gets a confidence score; below the threshold (slider, default 0.5) it is auto-excluded and drawn orange-dashed —
    click to put it back. Cubic / irregular particles (e.g. MOF-801) are unchanged.
- **Movable, editable labels on every graph**: drag sample labels / legends, double-click to edit text; positions follow each
  trace's offset; "Save this offset for every sample"; reset; exports use the moved positions; undo/redo and projects keep them.
- Update check on macOS prefers the download for the Mac's own architecture (macOS_arm64 / macOS_x86_64).
- O10-type images with tiny specks no longer break the automatic size estimate.

## 2.11.0 — 2026-11-16

Built by Stijn Germonpré.

- **Loading screen** instead of a white page on start (progress + stage), heavy libraries pre-loaded meanwhile.
- **Faster pages**: figure preview cache, decimated previews of large data, cached readers, fragment-based redraws.
- **Unified sample list** on every plotting page: show/hide (shift-click ranges), colour, rename, drag to reorder, "Stack evenly".
- **Shape model** clearly visible next to "Detect particles" (toggle + learned-template hint); also in Batch processing ("SEM/TEM particles").
- **Defect calculator for any node/metal**: node library (Zr6/Hf6/Ce6/Th6, Cu/Zn paddlewheel, Zn4O, Al-OH chains, Fe3O/Cr3O trimers,
  MOF-74, Ti8, ZIF) + custom nodes and linkers; molar masses computed automatically; TGA residue per metal. Zr/Hf results unchanged.
- **Project files (.reticula)**: save / open the whole session; Windows file association.
- **Global figure style** with journal presets (ACS, RSC, Nature, Elsevier, Wiley, …) applied on every page.
- **Undo / redo** for exclusions and edits on every page.
- **Update check** (opt-in, asked once): banner with Download / Later / Skip this version; never installs by itself; silent offline.
- Single-file installer Reticula_Setup_<version>.exe next to the split download.

## 2.10.0 — 2026-11-16

Built by Stijn Germonpré.

### Particle sizer: better recognition on real TEM images
- New default splitting **"adaptive v2.10"**: background flattening (smooth background from the support pixels, so
  uneven illumination does not matter) → particle contrast; hysteresis threshold from the background noise (seeds above
  max(Otsu, 8σ), grown down to the *Sensitivity* × σ, default 4) so light-grey particles next to dark ones are found;
  particle width estimated from the image (or typed: *Expected particle width*); one marker per particle from a blob
  filter at the particle scale plus the distance map of the locally darker cores, greedy non-maximum suppression;
  marker watershed on the contrast. Regions > 2× the particle area or with solidity < 0.80 are aggregates (excluded by
  default). Edge-cut particles are still removed (*Exclude particles cut by the image edge*, on by default).
- Optional **shape prior** (off by default, *Shape prior* checkbox, *Prior strength*, *Prior k*): learns the outline
  (radial profile, corner count from its Fourier harmonics, orientation spread, aspect ratio) and the size distribution
  from the confident isolated particles; weight w = n/(n+k) · consistency · strength. Rotated templates at the learned
  scale (normalised cross-correlation) seed a watershed inside aggregates and oversized particles; pieces are kept only
  if their size is inside the learned distribution and the outline matches the template. Panel *Learned shape prior*
  shows the template, example particles, the size histogram and the weight; prior splits have a dotted outline.
- User's TEM image (MOF-801 nanocubes, 120 hand-checked cubes outside the dense clumps): recall 58 % → 82 % (prior off)
  → 91 % (prior on); precision 77 % → 92 % / 91 %; ECD 29.6 ± 16.5 nm (v2.9, many fragments) → 38.3 ± 7.5 nm / 37.7 ± 6.8 nm.
- Cleaner overlay (picker and exported figure): thin teal outline = included, dotted = shape-prior split, red dashed
  outline + light fill = aggregate / excluded; numbers off by default (*# Numbers* button / *Show particle numbers*,
  hover for a tooltip). Working resolution default 1536 px.
- The v2.8 watershed and the v2.6 legacy split remain selectable.

## 2.9.0 — 2026-11-13

Built by Stijn Germonpré.

### Work is kept when you switch pages
- Every page keeps its uploads, settings and results when you open another page and come back (un-keyed widgets get a
  deterministic key from page + widget + parameters; uploaded files are kept per uploader and shown as "kept from
  earlier", ✖ forgets them; particle exclusions are kept).
- **Unsaved-work dialog**: leaving a page with changes asks *Do you want to save your current work?* –
  **Save** (kept in the session and written to `Documents\Reticula\workspace\<page>.json` + files, restored on the next
  start), **Don't save** (page reverts to its last saved state, or to the defaults), **Cancel** (stay). Home tiles go
  through the same guard. Sidebar: "● n unsaved change(s) · 💾 Save page" and "💾 Saved page work" (forget all).
- Defect calculator: methods a–h and the consensus are identical to v2.8 (regression test on the stored example inputs and
  on data pulled from the TGA / ICP / NMR / N₂ pages after navigating).

### SEM / TEM particle sizer: click to exclude
- Interactive picker on the image (offline HTML/SVG component, no full rerun per click): green = included,
  red = excluded by you, orange = aggregate (auto-excluded); one click toggles, drag a box to exclude (Shift / "Box:
  include" to include), undo (Ctrl+Z), reset to auto, include all, exclude all aggregates; live N, mean ± SD and
  histogram next to the image. Statistics, figures, CSV and the sample store use only included particles; the
  exclusion list is saved with the sample. The outline editor (merge / split / delete / add) is the secondary mode.
- Default polarity is now **dark particles on a lighter background (TEM)** everywhere (sizer UI, `segment()` API,
  helper script); "bright (SEM)" stays an option – the bundled SEM example images need "bright".

### UV/Vis (new page)
- Readers: Implen NanoPhotometer / NPOS Excel export (results sheet + per-sample spectra, layout-agnostic), Implen
  `.ids` (zip/XML/text containers; encrypted archives give export advice), Cary-style paired columns, generic CSV/TXT
  (decimal comma), wavelength-across tables, JCAMP-DX.
- Overlay (ACS style), absorbance ⇄ %T ⇄ %R, Kubelka–Munk, blank/baseline subtraction, zeroing, Savitzky–Golay
  smoothing, normalisation, peak picking with λmax table → sample store.
- Beer–Lambert calibration: linear / through origin / quadratic, R², s_y/x, LOD (3.3 s/m), LOQ (10 s/m), unknowns with
  SE and dilution, ε; assay presets BCA, Bradford, Lowry, Biuret, ninhydrin, OPA (+ degree of hydrolysis), A280, A205.
- Kinetics at one λ (times from the names), first-order fit, "send to the Kinetics page".
- Diffuse reflectance → Tauc (direct / indirect) with automatic or manual linear window, Eg table → sample store.
- 1st / 2nd derivative spectra with zero-crossing / minimum band positions. All figures go to CSV / SI bundle.

### PXRD
- Reads Rigaku SmartLab Studio II **.rasx** (Data*/Profile*.txt, intensity = counts × attenuator factor,
  MesurementConditions*.xml metadata) and legacy **.ras** (*RAS_INT_START blocks); several scans per file;
  scans whose axis is not 2θ/θ are flagged; target, Kα1/Kα2, kV/mA, range, step and speed go to the sample store.

### Credits
- "Built by Stijn Germonpré" under the version in the sidebar, the About box, the launcher window and the installer
  welcome page; installer publisher = Stijn Germonpré.

## 2.8.0 — 2026-11-13

### Structure (CIF) (new page)
- Reads SHELXL/Olex2, COD and CSD CIFs (multi-block, symmetry operators, Uij, disorder groups, `_geom_*` loops with ESDs).
- Figures in pure matplotlib: **ORTEP** ellipsoids from the Uij at 30/50/70 % probability with principal axes, ball-and-stick,
  stick, space-filling; labels for all / non-H atoms with **automatic placement** (away from the bonded neighbours, repel
  pass, leader lines) plus manual per-label offsets; H atoms all / polar only / none and an H-size slider; disorder parts;
  views along a, b, c, [uvw], (hkl) normal or a **best view** (mean plane, tilted to separate overlapping atoms) with
  rotation sliders; cell edges, axis triad; grow to the ASU, complete molecules, packing, cell or supercell; metal
  coordination polyhedra for MOFs; optional title (off by default); PNG/TIFF/SVG/PDF export.
- **H-bonds** D–H···A (H···A < Σr_vdW, angle ≥ 120°) dashed in the figure and as a table with symmetry codes.
- **SI tables**: crystal data and structure refinement (incl. Rint, completeness, GooF, R1/wR2 for I > 2σ(I) and all data,
  peak/hole, Flack, CCDC), bond lengths, angles, torsions with ESDs, H-bonds → DOCX and XLSX.
- **Simulated PXRD** with hkl labels (close reflections share a label) and an overlay on a measured pattern (scale + linear
  background, residual panel, Rwp).
- **MOF pore tools**: largest included (Di) and free (Df, periodic percolation) sphere, void fraction, probe-occupiable
  fraction and Monte-Carlo accessible surface area (m² g⁻¹, m² cm⁻³), method stated.
- **Hirshfeld surface (experimental)**: approximate d_norm point cloud and 2D fingerprint with contact percentages.
- py3Dmol was not added: it loads 3Dmol.js from a CDN and would not work offline in the exe.

### SEM / TEM particle sizing
- New **watershed** splitting (default): h-maxima markers of the distance map with a minimum marker distance tied to the
  automatically estimated particle radius (sliders), priority-flood watershed on the distance map, the image gradient
  (facet edges) or both; regions with solidity < 0.85 or > 2× the median area are re-split, otherwise classed
  **aggregate** (excluded from the statistics unless included). Bright round **support holes** (holey carbon) are
  excluded before the threshold; info strip and edge-touching particles excluded as before.
- **Manual corrections**: click to delete, merge, split along a cut line or add a particle; undo/clear; coloured,
  numbered labels per particle ('A' = aggregate). The v2.6 method stays available ("legacy").

### Multi-image SEM panels
- Scale-bar label, magnification and panel letters are sized to the panel (4.5–9 pt, one uniform font for all panels,
  user override), the bar thickness is set in points per panel, the label is aligned to the bar end, and every label is
  shrunk until it lies inside its own panel with padding (no more overlap with the neighbouring image).

### Windows launcher (white page)
- MIME types for .js/.mjs/.css/… forced before Streamlit starts (a registry mapping of .js to text/plain gave a white
  page), 127.0.0.1 everywhere, browser opened only after /_stcore/health is ok and the JS bundle is served correctly,
  old Reticula servers on 8501–8510 stopped, cache-busting URL, proxy variables ignored for the local server;
  status line, log path, "Reload app" and "Open log" in the launcher window.

## 2.7.0 — 2026-11-13

### Figure composer
- **Journal templates** ACS (3.25 / 7.0 in), RSC (8.3 / 17.1 cm), Wiley (8.5 / 17.5 cm): width, font size and panel-label
  style set from the template; a **check** flags a wrong width and any text below the journal minimum.

### Batch processing (new page)
- One technique (PXRD, FTIR, Raman, XPS, TGA, DSC, N₂ isotherm), all files of a folder (recursive), an upload or the
  example folder; one figure per file, overlay or stacked; **one preset** (saved figure-settings JSON or the house style);
  journal width; PNG/TIFF/PDF/SVG + CSV in a ZIP and optionally in an output folder; per-file log; preset saved in the ZIP.

### SI export bundle (new page)
- Every figure shown is collected (sidebar switch, on by default) with its data, its figure settings and the analysis
  settings at capture time. One ZIP: figures/ (PDF, PNG, SVG), data/, **SI_parameters.xlsx** (Figures, Parameters,
  Results per sample), SI_parameters.csv, the defect SI table (DOCX/XLSX), README.

### Sample database (new page)
- Local **SQLite** database `Documents\Reticula\reticula.db`. Imports the sample IDs of the synthesis pipeline
  (Master_sample_overview.xlsx; a copy is bundled, point it to your Synthesis_Pipeline folder to link the protocol
  PDF/DOCX files), shows paper, family, parent and children, basis and source DOI, opens/downloads the protocol; syncs every
  result of the sample store; links data files whose names contain a sample ID (batch from `-B<n>`, technique from the
  name/folder); notes per sample; download of the database file.

## 2.6.0 — 2026-11-13

### SEM / TEM
- **Particle sizing**: automatic (Otsu or manual threshold, smoothing, opening, hole filling, splitting of touching
  particles with distance-transform markers – over-split pieces that share a long boundary are merged back –, edge
  particles excluded, min/max size) or **click-to-measure** (two clicks per particle, editable table). ECD, Feret max/min,
  aspect ratio; outlines overlay; **histogram with log-normal fit and mean ± SD** (n, median, σg, mode, KS p); stored
  per sample (section "sem"). Several images of one sample give one combined histogram.
- **Batch options**: **matched scale bars** (one physical bar value on every panel) and **one calibration for the whole
  batch** (first image's nm/px, or a typed value).

### EDX
- Colour picker and gain per element, **overlay on the SE/BSE image** (additive blend; SE image auto-detected by
  "_SE" in the name, or any loaded SEM image), legend, scale bar; **line scan** (start/end, band width, normalised
  profiles; distance in nm/µm when the pixel size is set).

### ICP-OES / XRF
- **Zr : Hf : S molar ratios with propagated uncertainties** (Gaussian; replicate SD vs pooled instrument SD, the
  larger is used): Zr:Hf, Hf/(Zr+Hf), S/(Zr+Hf), S per M₆ node ± σ, Zr:Hf:S per M₆; CSV; bar chart of S per node.
  Values go into the sample store; the defect calculator's *Pull data* shows S/M and Hf/(Zr+Hf) with σ.
  The older "S per Hf6 node" column is now "S per M6 node" (Zr + Hf).

## 2.5.0 — 2026-11-13

### FTIR / Raman
- **Normalise to a reference band** (height or area above a local chord, e.g. νas(COO⁻) 1540–1620 cm⁻¹ or ν(C=C) for Raman);
  traces then share one intensity scale so relative intensities survive the stacking.
- **Peak labels** (prominence threshold, max. labels, chosen traces) with a peak table; **zoomed inset** of any range;
  **difference spectra** (sample − f × reference, f = 1, manual or least-squares in a window) as a separate figure.
- ALS baseline was already available (Processing → Baseline → ALS).

### ¹H NMR
- Integration table now has a "linker : modulators" ratio string; integrals go to the sample store and pre-fill the defect calculator.

### XPS
- Peak fitting with **Shirley** (iterative) or **Tougaard** (universal cross-section) background, GL components and
  **spin–orbit doublets with fixed splitting and area ratio**: Zr 3d (2.43 eV, 2:3), Hf 4f (1.71 eV, 3:4), S 2p (1.18 eV, 1:2),
  plus Hf 4d / Ti 2p; residual shown. **Atomic-% table** with editable RSFs (Scofield defaults), Hf/(Zr+Hf), S/(Zr+Hf);
  multi-region EXAMPLE (ZrO₂@Sample B: Zr 3d, Hf 4f, S 2p, O 1s, C 1s).

### TGA
- **Automatic step labels** (mass loss and DTG peak per step, steps < threshold merged), step table, **plateaus normalised
  to the oxide residue** with the ideal plateau and linkers per node (same equations as the defect calculator), and a
  **TGA–MS overlay** (ion currents vs temperature, m/z labels such as H₂O, CO₂, SO₂; time-based files mapped with T₀ + βt).

### DSC
- **Cycle-to-cycle comparison** (Heat 1 vs Heat 2 …): table with ΔT peak/onset and ΔH ratio vs cycle 1, overlay figure;
  heat–cool–heat EXAMPLE (REAL heating 1, SIMULATED cooling and heating 2).
- **cp with a sapphire reference** (ASTM E1269 ratio method, sapphire table 150–1000 K) + direct cp for comparison; EXAMPLE runs.

## 2.4.0 — 2026-11-13

### Shared sample store → one defect report per sample
- TGA, ¹H NMR digestion, ICP-OES / XRF, N₂ sorption and kinetics put their per-sample results into a shared **sample store**
  (session + `Documents\Reticula\sample_store.json`; names matched ignoring case, spaces, TeX and punctuation).
- Defect calculator: **Pull data for sample** pre-fills TGA weights, ICP metal wt-%, XRF S/M, NMR integrals and N₂ V_total;
  the consensus result is stored back. **SI table** (all samples: L, x, caps, open sites, formula, methods, S_BET, V_total)
  as XLSX or DOCX (landscape Word table; written without extra libraries). Store download / merge / clear.

### N₂ sorption
- Rouquerol criteria panel (all four criteria per sample, chosen range, C, n_m, R²); **S_BET in the legend** option.
- t-plot (Harkins–Jura) with adjustable t range: micropore volume, external surface area, plot.
- **Computed PSD** from the isotherm: DFT kernel fit (pyGAPS N₂/77 K carbon-slit kernel, MIT licence, NNLS + Tikhonov
  regularisation – *not* an NLDFT/QSDFT kernel for oxides), Horváth–Kawazoe (slit, carbon), Saito–Foley (cylinder,
  oxide-ion parameters) and BJH (Harkins–Jura thickness, Kelvin). The method and its limits are written on the figure caption.

### PXRD
- Scherrer table for several peaks with instrumental correction (constant or from a LaB₆/Si standard pattern; Gaussian,
  Lorentzian or Voigt width correction; K editable). Crystallinity index vs a reference pattern (integrated reflection
  windows). Peak-shift tracking: Δ2θ vs a reference trace, d-spacing, cubic a from hkl, with a plot.

### Kinetics
- k_obs with **95 % confidence intervals** (t-distribution × SE), CI or SE error bars on the bar chart.
- Optional temperature / pD columns: one fit per (sample, T, pD); **Arrhenius** (Ea, ln A) and **Eyring** (ΔH‡, ΔS‡, ΔG‡298)
  with CIs directly from the fits; temperature-series EXAMPLE.
- Initial-rate tool (linear fit up to a chosen conversion, r0 ± CI); **TOF per open site** with open sites per node and
  M(node) pulled from the defect calculator; **pD–rate profile** fit (sigmoidal or bell-shaped → pKa ± CI).

## 2.3.0 — 2026-11-13

### DSC – full heating / cooling cycle
- Segments are detected from the sign of dT/dt (heating, cooling, isothermal) and numbered Heat 1, Cool 1, Heat 2, Iso 1 …;
  ramp/hold boundaries are sharpened so switching transients do not end up in a ramp.
- Works for (a) one file with all segments, (b) separate heating and cooling exports of one sample – **Combine files of one
  sample into one cycle**, ordered by the file time stamp or the table order, time axis made continuous – and (c) STARe files
  with several 'Curve Name:' blocks (restarted clocks are joined).
- Views: **Full cycle vs temperature** (closed loop, direction arrows, colour + line style per sample and segment kind),
  **Full cycle vs time** (temperature programme Ts or Tr on a twin axis, shaded segments), **Stacked: heating above cooling**.
- Peaks on every segment (melting on heating, crystallisation on cooling) and a **cycle summary**: Tm / Tc onset and peak,
  ΔHm, ΔHc, supercooling ΔT = Tm − Tc (peak and onset) and |ΔHc/ΔHm|.
- Full-cycle EXAMPLE: the real DSC-cycle_1.1 heating + a simulated isothermal and cooling exotherm (clearly labelled).
- Peak limits: broad peak tops no longer stop the limit search early; end levels are taken robustly.

## 2.2.0 — 2026-11-13

### DSC (new page under Thermal, next to TGA)
- Readers: Mettler STARe (Curve Values blocks; footer mass, method and the STARe *Results* are shown next to Reticula's values),
  TA TRIOS and Universal Analysis, Netzsch Proteus (exo-down handled, `Segment` column), PerkinElmer Pyris (endo-up flipped),
  generic CSV/XLSX. Decimal commas, several curves/segments (also split by heating/cooling rate).
- Heat flow in mW, W/g or mW/mg (mass normalisation from the file or the table), vs temperature (Ts) or time; exo up/down with
  arrow; segment selection with colours; offset stacking.
- Peaks with linear or sigmoidal baseline: onset / peak / endset (tangent construction) and ΔH (J/g, integrated over time);
  Tg (ISO 11357-2 midpoint, inflection, Δcp) with an auto-suggested window; TGA–DSC (STA) dual-axis overlay; annotations.
- Validated on a real STARe export (DSC-cycle_1.1, 1 K/min): onset 48.30 °C / peak 48.85 °C / ΔH −105.6 J/g vs STARe 48.29 / 48.85 / −105.38.

### SEM / TEM — fixes after user feedback
- **⚡ Auto-calibrate** button, a preview of what was detected on the full image, and an "Applied to figure: X nm/px" line
  showing the new bar length in original pixels, the crop row and the magnification label.
- Burned-in (overlay) bars are cropped away, so the old bar no longer appears in the final figure; no silent fall-back to a
  default 1 µm when the label is not read (a warning asks for the value); widgets reset after re-detection; duplicate file
  names supported. End-to-end browser test that clicks Auto on real images and measures the exported bar.

## 2.1.0 — 2026-11-12

### SEM / TEM — calibration from the original image
- New `mfs/semcal.py` + a **📐 Calibration** panel per image. It shows the metadata pixel size, the detected scale bar (px),
  the label read by OCR, the resulting nm/px and the magnification (metadata and/or OCR), with cross-check notes.
- Sources, in order of reliability: metadata (FEI/Thermo PixelWidth, Zeiss *Image Pixel Size*, ImageJ, DM/EMD, and now
  **JEOL .txt / Hitachi .txt / TESCAN .hdr sidecars** uploaded with the image) › **scale bar detected in the info strip**
  (longest solid bar, measured to the pixel, label read by Tesseract OCR, bundled in the Windows app) › magnification ›
  **scale bar (mark ends)**: click the two bar ends on a zoomable strip ("Select the scale bar length"), then type the value.
- Every result can be overridden (bar length, value, unit) and ticked as checked; the new scale bar uses it.
- The magnification is read from the strip ('x10,000', '10.0kx', '50 000 x', 'Mag = 10.00 K X', '×5.00k' …), cross-checked
  against the metadata, and can be printed on the final figure (editable text, any corner).
- The info strip is detected and cropped automatically. A built-in example without metadata is included.
- Olympus SIS / analySIS TIFFs (e.g. JEOL JEM-1010 + MegaView III): palette images are mapped to grey, the SIS pixel size and
  magnification are read, bars burned into the image (no strip) are read with an ink-colour OCR pass, the old bar is cropped,
  and a warning appears when a burned-in bar disagrees with the metadata or matches another image's calibration (stale bar).

### N₂ sorption — TriStar/3Flex full reports with side-by-side report bands
- Sheets holding several report bands separated by `|` columns are split before parsing; decimal commas in
  metadata are read; BET (Rouquerol) falls back to 3 points for steep Type I isotherms and says so; isotherm quality warnings.

## 2.0.0 — 2026-11-12

### N₂ sorption (3Flex / MicroActive) — bug fix
- **Cause of the zig-zag isotherm.** The v1 reader appended *every* table in the Excel export that had a
  p/p₀ and a quantity column into one series. For a MicroActive *full report*, that includes Isotherm Tabular,
  the BET report points, the linear plot overlay (side-by-side ads/des columns) and the log plot overlay, which is often in **mmol/g**.
  The plot then sorted all of these points by p/p₀. The result was two interleaved curves (one about 22× smaller), a zig-zag adsorption
  branch, a squashed y range, and a BET computed on mixed data.
- **Fix.** There is a new table-aware reader (`mfs/readers.py: read_sorption`). It:
  - finds every header and numeric block and its title, and classifies each one (Isotherm Tabular, plot overlay, BET, Langmuir, t-plot, DFT);
  - uses **one** isotherm table by priority;
  - splits adsorption/desorption by labels, by column pairs or by measurement order (turning point at max p/p₀);
  - converts mmol/g → cm³/g STP (×22.414);
  - reads the reported BET (± err, C, Qm, range, number of points), Langmuir, single-point, V_total, t-plot, sample name, mass, adsorptive, bath temperature;
  - reads `.xls` with xlrd.
- **Plot (Type I MOF convention).**
  - closed symbols for adsorption and open for desorption;
  - lines in measurement order, with desorption starting from the last adsorption point;
  - optional log p/p₀, or linear and log side by side;
  - PSD as a separate panel;
  - **one** legend with adsorption/desorption key entries (v1 drew two legend boxes).
- **Per-sample style.** Adsorption colour, desorption colour, marker, desorption marker, marker size, line style and line width.
- **BET.** The app's own Rouquerol/ISO 9277 calculation is unchanged. The value reported in the file is shown **for reference only**, with the file's fitting range. For comparison the table also shows the app's value on the file's range (classical 0.05–0.30). The difference between the two values is normally the p/p₀ window, not an error.
- New file-structure inspector ("which tables were found and which one was used").

### SEM / TEM / EDX scale bar
- Mode **number**: presets 200 nm, 500 nm, 1 µm, 2 µm, or a custom value with unit nm, µm or mm. The bar length is recomputed from nm/px.
- Mode **length**: set the bar length (fraction slider ↔ px number input). The label updates.
- Mode **auto**: a 1-2-5 value at about 20 % of the width.
- Manual override of magnification (with the reference display width) or nm/px.
- Position (4 corners), colour (white/black), thickness, font size, background box, label on/off and custom label text.
- An option to give every image the same bar style.
- TEM: Gatan `.dm3/.dm4`, FEI/TIA `.ser`, Velox/Berkeley `.emd` and `.mrc` via **ncempy**, with the pixel size read from the file. ImageJ-calibrated TIFFs are read too. TIFF/PNG without calibration take a manual nm/px.

### Figure controls for every technique (`mfs/figctl.py`, `mfs/figctl_ui.py`)
- Title and axis titles (Unicode, or mathtext such as `cm^3 g^{-1}`).
- Axis min/max and display-unit conversion. The data and CSV are unchanged; the ticks are relabelled. Conversions:
  - cm³ g⁻¹ STP ↔ mmol g⁻¹ / mg g⁻¹
  - °C ↔ K
  - wt % ↔ mg (with m₀)
  - cm⁻¹ ↔ µm/nm
  - 2θ ↔ d (Å) / q (Å⁻¹)
  - h ↔ min/s/d
  - nm ↔ Å
  - mM ↔ µM
- Legend:
  - on/off, rename, reorder, show/hide entries;
  - every matplotlib location plus outside right, top, bottom and left;
  - columns, font size, frame, marker scale.
- Ticks:
  - sides: all 4, any 3, 2 (bottom+left, bottom+right, top+bottom), 1, or custom;
  - direction in/out/inout;
  - major increment per axis (**slider and number input, kept in sync**);
  - minor ticks on/off with an increment;
  - tick labels per side;
  - frame follows the ticks or stays a full box.
- Per-series colour pickers on every page.
- **House-style presets.** Save and load all settings as JSON, globally in the sidebar or per figure.

### Downloads — bug fix
- **Cause.** Every download click re-ran the whole script. On every rerun the v1 app eagerly regenerated PNG, TIFF, SVG, PDF and CSV for every figure on the page, which took several seconds on the N₂ and defect pages. Streamlit then invalidated the previous media files, so a click during or just after a rerun returned *"Download canceled / file wasn't available on site"*. Reproduced with headless Chromium on v1.
- **Fix.**
  - Exports are **deferred**: they are generated only when clicked, through a callable.
  - Download buttons use `on_click="ignore"`, so they cause no rerun.
  - Explicit MIME types.
  - Windows-safe file names.
  - Deterministic SVG/PDF.
  - A render lock.
  - A **ZIP (all formats)** button.
  - A **💾 Save to folder** button that writes to `exports/` next to the app and shows the full paths.
  - The defect PDF report is also generated on click.
- Verified with Playwright and Chromium: every download on six pages arrives with the right type and a valid size, including rapid clicks.

### Other
- Duplicate sample names are made unique.
- The colour cascade with duplicate families is fixed.
- New tests: `tests/test_v2_features.py` and `tests/test_downloads_browser.py` (Playwright, skipped if Playwright is not installed), plus app-interaction tests for the figure controls and the scale bar.
- v1 reader and plot code kept for reference in `tools/legacy/`.
- `tools/n2_before_after.py <file.xls>` writes v1-vs-v2 comparison PNGs for any 3Flex export.

## 1.0.0
- First release: all Paper 1 techniques, defect calculator, figure composer.
