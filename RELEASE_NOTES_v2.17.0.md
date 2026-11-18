## Reticula 2.17.0

**Downloads**
- Windows 10/11: `Reticula_Setup_2.17.0.exe` (not code-signed - SmartScreen: *More info* → *Run anyway*)
- macOS, Apple Silicon: `Reticula_2.17.0_macOS_arm64.dmg` (not notarised - first start: right-click → *Open*; see `README_macOS.txt`)
- `SHA256SUMS.txt` to verify the downloads. Everything else: run from source (see the README).

**New**
- **Schemes & Icons** (🎨 in the menu): a scheme editor with about 800 offline CC0 icons. Multi-select, align, distribute, group, copy and paste, arrows with bends and automatic routing around icons, subscript and superscript, a periodic table of the elements, and export to SVG, PDF, PNG or TIFF at journal widths and dpi.
- **Content from Reticula in a scheme**: a methods schematic from a sample, your own x/y data, a molecule from a SMILES, a crystal structure from a CIF, and a scheme as a panel in the combined figure.
- **Rietveld and structures**: ICSD import is manual (a button opens the ICSD login; Reticula sends no query); the CIF picture is centred on the molecule; the refinement shows elapsed time and time left.
- **Defect calculator**: the TGA plateau for the defect calculation is chosen automatically at the flattest part, or picked on the graph; the results table says how it was chosen.
- **SI export**: choose the journal style for the Methods heading and wording: ACS, Wiley-VCH, RSC or Elsevier.
- **Simple mode** (black graph style), more visible work tabs, NMR figure first with literature-style proton labels.
- **Zetasizer .dts**: files of software version 8.02 are read (sample names and size distribution).
- **Project files** are `.reticula` (zip format). The start page can open and save them.

**Limits, stated plainly**
- Zetasizer 8.02: Z-average and PDI are not calculated for 8.02; use the software's export. ζ is not read from `.dts`.
- TGA plateau rule: based on the main text of the source paper; the Supporting Information tables were not checked against it.
- Journal SI wording: templates; check against the current author guidelines of your journal.
- Schemes: the browser preview shows the default icon style (the export uses the chosen style); the panel in the combined figure is a 600 dpi image, not vector; own icons live in the scheme where they were added.
- Technical names (installer file name, macOS bundle id, data folder `Documents\Reticula`, update check) keep the old name in this version.

Full to-do list for the next version: `TODO_v2.19.md`.
