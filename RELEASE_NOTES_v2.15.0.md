## Reticula 2.15.0

**Downloads**
- Windows 10/11: `Reticula_Setup_2.15.0.exe` (not code-signed - SmartScreen: *More info* → *Run anyway*)
- macOS: no free macOS runner on GitHub; build it yourself in a few minutes with `BUILD_ON_MAC.txt` (Dutch: `BUILD_ON_MAC_NL.txt`).
- `SHA256SUMS.txt` to verify the downloads. Everything else: run from source (see the README).

**New**
- **Detail levels**: Standard (default, only what a publishable figure needs), Advanced (more details) and Developer (everything) – one switch for all modules plus an override per module.
- **PXRD, opt-in (hidden until you switch them on)**:
  - **Crystallite size and microstrain**: Scherrer, Williamson–Hall, Halder–Wagner with uncertainties, instrumental broadening, plain-language warnings.
  - **Rietveld / Le Bail**: choose the phases (CIF), press *Refine*; automatic staged strategy, ≈ 1 s per pattern, Rwp / GoF / Bragg R, cell ± esd, weight fractions, standard Rietveld figure, Excel report. Validated on synthetic patterns with known answers.
- **Every figure file carries its provenance** (version, page, source files, settings, time) in PNG / TIFF / PDF / SVG; the ZIP adds `provenance.json`.
- **SI bundle**: all figure data in one Excel, **journal check with one-click fix** (width, fonts, lines, colour-blind safe), Methods text for every technique, sample summary table (Excel / CSV / LaTeX).
- **Photoluminescence thermometry**: robust automatic global fit of the band ratio (overlapping, dirty or doubled bands), errors and quality flags, ranges selectable on the graph.
- **FTIR / Raman**: band assignment from an editable list, multi-band fit. **XPS**: batch fit of all spectra. **Particle sizer**: sizing templates per sample type and a confidence bar. **¹H NMR**: linker / modulator library with molar ratios.

**New and not yet seen on many real files** – please report anything that looks wrong: nmrML import, the proton chips next to the structure in the ¹H NMR assignment, range selection on the PL graph, the automatic spacing of the combined figure.

**Next version (2.16)**: one project per material, more photoluminescence tools (calibration, repeat-based σT, power dependence, CIE colour), Rietveld preferred orientation.

Full list: `app/CHANGELOG.md`.
