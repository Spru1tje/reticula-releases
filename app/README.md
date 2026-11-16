# Reticula

*Publication figures and analysis for framework materials (MOFs, COFs, zeolites, …).*

One local web app (Streamlit) that turns raw instrument exports into journal-style figures (ACS, RSC, Wiley, Nature, Elsevier, …) and
does the usual analysis on the way: PXRD, FTIR/Raman, ¹H NMR, N₂ sorption (BET/Rouquerol, t-plot, PSD), TGA, DSC, SEM/TEM particle sizing,
XPS, UV/Vis, ICP/XRF, kinetics, crystal structures from CIF, and a defect calculator for MOFs. It runs entirely on your computer; your
data are never uploaded.

> All files in `example_data/` are **EXAMPLE DATA** (synthetic, written in each instrument's export format) so that you can try every
> page. The two files with `_REAL` in their name are a real DSC heating run of a MOF-801 sample. None of this is a reference measurement.

## 1. Install and start

**Easiest - the installers** (no Python needed): download from the [Releases](../../../releases) page

| System | File | Notes |
|---|---|---|
| Windows 10/11 (64-bit) | `Reticula_Setup_<version>.exe` | Not code-signed: SmartScreen → *More info* → *Run anyway*. Installs for your user, no admin rights. |
| macOS (Apple Silicon) | `Reticula_<version>_macOS_arm64.dmg` | Not notarised: first start via right-click → *Open* (see `README_macOS.txt`). |

**From source** (any system; Python 3.11-3.13, 3.12 recommended):

| OS | Double-click | What it does |
|---|---|---|
| Windows | `run_windows.bat` | First run: creates `.venv` and installs `requirements.txt` (about 3-5 min, internet needed). Then it opens the app in your browser. |
| macOS | `run_mac.command` | The same. If macOS blocks it, right-click → Open once, or run `chmod +x run_mac.command`. |
| Linux | `./run_linux.sh` | The same. Needs the system package `tesseract-ocr` for the scale-bar OCR (optional). |

Later starts take a few seconds. The app opens at http://localhost:8501 (only reachable from your own computer). Close the terminal window to stop it.
By hand: `python -m venv .venv`, activate, `pip install -r requirements.txt`, `streamlit run app.py`.

**Tests**: `pip install pytest` (and `playwright` + `playwright install chromium` for the real-browser tests), then `python -m pytest -q` in this folder.
The suite covers the readers, analyses, plots, exports, the defect calculator against literature values, every page through Streamlit's AppTest
and the real-browser downloads. A few tests use private measurements that are not in the repository; they are skipped automatically.
`python tools/make_example_data.py` regenerates the synthetic examples.

**Releases** are built by GitHub Actions: push a tag `vX.Y.Z` (equal to `__version__` in `mfs/__init__.py`) and the workflow runs the tests, builds the
Windows installer and the macOS dmg, runs a self-test inside each built app and publishes the release (`.github/workflows/release.yml`, `pkg/`).

## 2. How to use it

1. Pick a technique in the sidebar.
2. Upload one or more files. You can also tick *Use the bundled EXAMPLE data*.
3. Rename the samples and change their order in the sample list (show/hide, colour, drag to reorder). Samples take the palette colours in list order
   (Okabe-Ito by default, colour-blind safe; other palettes under *Global figure style*). To keep a colour for a sample in every figure, add a **colour
   rule** under *Sample → colour rules* in the sidebar (a regular expression on the sample name and a colour, e.g. `uio-66` → `#0072B2`) or pick the colour
   in the list. *Example-set rules* loads the rules that match the bundled example data.
4. Choose processing (baseline, smoothing, normalisation, offset) and figure width:
   - 3.25 in, ACS single column
   - 7 in, ACS double column
   - 8.5 cm or 17.5 cm, Wiley/Angew
5. Optional: open **🎛 Figure controls** under the figure to change title, axes, units, ticks, legend and colours (see section 2a).
6. Download the figure as SVG, PDF (fonts embedded as TrueType), 600 dpi TIFF (LZW), 600 dpi PNG, CSV of the processed data, or **ZIP (all)**. Files are generated when you click. If your browser blocks downloads, click **💾 Save to folder**: the files are written to `exports/` next to `app.py`, and the full paths are shown.
7. Click **➕ Add this panel to the figure composer** under a figure. Then open *Figure composer* to build a multi-panel figure with (a)/(b)/(c) labels and export it as vector graphics.

**Style.** Every figure follows the same rules:
- Arial or Helvetica 7–9 pt. On Linux, Liberation Sans is used because it has identical metrics.
- Axis and tick linewidth 0.6–0.8 pt; inward ticks on all four sides; no grid; no title (all of these are defaults, and each can be changed in Figure controls).
- Units in parentheses, e.g. *2θ (°)*, *Wavenumber (cm⁻¹)*, *δ (ppm)*, *Quantity adsorbed (cm³ g⁻¹ STP)*.
- Offset spectra have no y tick labels ("a.u."). The legend or direct labels are placed at the trace ends.

## 2a. Figure controls and scale bars (details of every release in `CHANGELOG.md`)

**Figure controls (every technique).** The expander under each figure has these tabs:

- *Title & axes*
  - Title and axis titles. Unicode works; so does `^`/`_`, e.g. `cm^3 g^{-1}`.
  - Min/max and display unit per axis:
    - cm³ g⁻¹ STP ↔ mmol g⁻¹ / mg g⁻¹
    - °C ↔ K
    - wt % ↔ mg
    - cm⁻¹ ↔ µm
    - 2θ ↔ d / q
    - h ↔ min
    - nm ↔ Å
    - mM ↔ µM
  - Unit conversion relabels the ticks. The CSV stays in the original units.
- *Ticks*
  - Sides: all 4, 3, 2, 1 or custom.
  - Direction in/out/inout.
  - Major increment per axis: a slider and a number box, kept in sync; 0 = automatic.
  - Minor ticks and their increment.
  - Tick labels per side.
- *Legend*
  - On/off.
  - Rename, reorder and hide entries in the table.
  - Location: every matplotlib position, or outside right/top/bottom/left.
  - Columns, font size, frame and marker scale.
- *Colours*: a colour picker per series.
- *Preset*: save or load this figure's settings as JSON. The sidebar holds a global **house style** (JSON) that applies to every figure.

**SEM / TEM / EDX scale bar.** Each image has its own scale-bar block:

- Mode *number*: 200 nm, 500 nm, 1 µm, 2 µm, or a custom value and unit. The bar length is computed from nm/px.
- Mode *length*: drag the bar length; the number follows.
- Mode *auto*.
- Pixel size from metadata, from the magnification (reference width editable), or typed in manually.
- Position, white/black, thickness, font size, background box, custom label.
- *Same style for all images* applies one bar style to every image.
- The caption under each image shows the resulting bar in px and nm.
- TEM files: `.dm3/.dm4/.ser/.emd/.mrc` are read with ncempy.

**N₂ sorption.**

- Full MicroActive report exports (.xls/.xlsx) are read table by table. Only the Isotherm Tabular table is plotted.
- Branches follow the measurement order: filled symbols for adsorption, open for desorption.
- Optional log p/p₀, or linear and log side by side.
- One legend.
- Adsorption and desorption colours, markers and lines can be set per sample.
- The *File structure* expander shows which tables were found and which one was used.
- The BET table shows this app's Rouquerol BET next to the value reported in the file (reference only) and the app's value on the file's p/p₀ range.

To compare v1 and v2 on any export, run `python tools/n2_before_after.py your_file.xls out/prefix`. It writes `prefix_before_v1.png` and `prefix_after_v2.png`.

## 3. Modules

| Module | Input formats | What it produces |
|---|---|---|
| **PXRD** (STOE STADI P, Cu Kα1) | WinXPOW `.xy/.dat/.txt/.csv`. Binary STOE `.raw` only with optional `xylib-py`. | Stacked or overlay patterns. Simulated pattern from any **CIF** via pymatgen (Cu Kα1, pseudo-Voigt FWHM, optional Zr→Hf substitution), shown as grey trace and/or sticks. Background removal (ALS, rolling ball, polynomial), smoothing, normalisation, highlight bands. **Scherrer table** for several peaks with instrumental correction (constant or LaB₆/Si standard; Gaussian/Lorentzian/Voigt), **crystallinity index** vs a reference pattern, **peak-shift tracking** (Δ2θ, d, cubic a). Bundled CIFs: MOF-801(Zr) COD 4121459 and UiO-66(Zr) COD 4512074. |
| **ATR-FTIR / FT-Raman** (Vertex 70 + RAM II) | OPUS binary (`.0`, `.1` …) via brukeropusreader; `.dpt`; CSV | Transmittance or absorbance (with T→A conversion), Raman 1064 nm. Stacked plots with labelled bands such as νas(COO⁻) and ν(S=O) at about 1040/1180 cm⁻¹, ATR baseline, normalisation. **Normalise to a reference band**, automatic **peak labels**, **zoomed insets**, **difference spectra**. |
| **¹H NMR** (Bruker; any 1D ¹H) | Bruker experiment folder as **.zip** (fid or processed `pdata/1/1r`; Avance II and NEO, TopSpin 2-4), JCAMP-DX, TopSpin `totxt`, or x,y | **Characterisation**: automatic phase / baseline / referencing, multiplet analysis (J couplings), structure from the built-in editor (Ketcher) or a molfile/SMILES, ¹H shift prediction (nmrshiftdb2 HOSE table), automatic peak → atom assignment with drag-to-correct, impurity check (Fulmer 2010 / Gottlieb tables), qNMR, labelled figure, ACS-style peak list, SI export. **Digestion** sub-mode: integration windows give linker : modulator ratios and concentrations vs an internal standard (e.g. TMSP-d4), which feed the defect calculator. |
| **N₂ sorption, BET, PSD** (3Flex / MicroActive) | MicroActive **Excel export** (full report or Isotherm Tabular + DFT, .xls/.xlsx, single or multiple sheets), CSV (also semicolon/decimal comma) | Isotherms (ads filled, des open, measurement order, linear and/or log p/p₀, per-sample colours) and PSD (dV/dw, cumulative). **BET with the 4 Rouquerol criteria** (ISO 9277) and a Rouquerol plot, plus the **range of areas over all fitting windows** (how much the BET area depends on the chosen p/p₀ range). t-plot (Harkins–Jura) micropore volume and external area, V_total at p/p₀ = 0.95, BET table CSV. Rouquerol criteria panel, S_BET in the legend. **Computed PSD** from the isotherm: DFT kernel fit (pyGAPS N₂ 77 K carbon-slit kernel, MIT – not an oxide NLDFT/QSDFT kernel), Horváth–Kawazoe, Saito–Foley, BJH. |
| **TGA** | Netzsch Proteus ASCII, TA Universal Analysis/TRIOS text, Mettler STARe export, generic CSV (mg or %) | TGA with DTG (second panel or twin axis), optional renormalisation to the oxide residue = 100 %. **Automatic step labels**, plateaus normalised to the oxide residue with linkers per node, **TGA–MS** overlay. W(200 °C), plateau and end values are passed to the defect calculator. |
| **DSC** | Mettler STARe, TA TRIOS / Universal Analysis, Netzsch Proteus, PerkinElmer Pyris text exports, generic CSV/XLSX (decimal commas, several segments) | Heat flow (mW, W/g, mW/mg) vs T or time, exo arrow (up/down), segment colours, offsets. Peaks with linear/sigmoidal baseline: onset, peak, endset, ΔH (J/g); Tg (midpoint, inflection, Δcp); STARe-reported values shown for comparison; TGA–DSC overlay for STA files. **Cycle-to-cycle comparison** (Heat 1 vs Heat 2) and **cp with a sapphire reference** (ASTM E1269). |
|  | *Full cycle* | Heating / cooling segments from dT/dt (Heat 1, Cool 1 …), combine separate exports of one sample, full cycle vs T (loop with arrows) or vs time (temperature programme, shaded segments), stacked heating/cooling, cycle summary with supercooling ΔT = Tm − Tc. |
| **Defect calculator** | Values from the other tabs or typed in; JSON save/load | See section 4. **Sample store**: TGA, NMR, ICP/XRF, N₂ and kinetics results are stored per sample ID and pulled into the calculator; **SI table** for all samples as XLSX/DOCX. |
| **Kinetics and adsorption** | CSV/XLSX: sample, (replicate), time, concentration, or raw NMR integrals + TMSP | First-order fits with k_obs ± SD (from the covariance), t½ ± SD and R². C/C₀ plots with fits, k_obs bar chart, blank subtraction. Langmuir and Freundlich fits (q_max, K ± SD). **95 % CI** on k_obs; per-temperature fits → **Arrhenius** (Ea) and **Eyring** (ΔH‡, ΔS‡, ΔG‡) with CIs; initial rates; **TOF per open site** from the defect calculator; **pD–rate profile** (sigmoidal/bell → pKa). |
| **ICP-OES / XRF** (Syngistix CSV or any table) | Syngistix CSV (long or wide), any sample–element–value table | mg/L → wt% (digested mass and volume), mean ± SD bars per element (one panel per element or grouped), molar ratios, Zr/Hf fraction; **molar ratios ± σ for any elements** (choose the reference metal, per n atoms), and the Zr:Hf:S special case **± propagated σ** (Hf/(Zr+Hf), S per M₆ node). They feed the defect calculator. |
| **SEM / TEM / EDX** (Philips XL30 FEG; any TEM) | TIFF (FEI/Philips `[Scan] PixelWidth`, Zeiss, ImageJ calibration, magnification) with the databar cropped automatically; TEM `.dm3/.dm4/.ser/.emd/.mrc` (ncempy); PNG/JPG with manual pixel size; EDX element maps as images | Fully adjustable scale bar (value or length, position, colour, thickness, font, box), (a)–(d) panel grid, matched scale bars and one calibration for a batch; **particle sizing** (v2.10 adaptive recognition with background flattening, particle-scale splitting and an optional learned shape prior; click-to-exclude picker; or click-to-measure) with log-normal histogram; EDX maps in your colours, **overlay on the SE image** and **line scans**. |
| **SDS-PAGE** (Gel Doc EZ / Image Lab) | Image Lab *Lane and Bands* table (CSV/XLSX) | Band % vs time per band (intact protein and fragments). Lane → time mapping. |
| **Zeta / DLS** (Malvern) | Zetasizer CSV export (summary + size distribution) | ζ-potential bars ± SD, DLS intensity distributions (log x). |
| **XPS** (optional) | CSV/TXT (CasaXPS, Avantage), simple VAMAS | Region spectra with reversed BE axis, stacked samples. **Peak fitting** with Shirley or Tougaard background and spin–orbit doublets (Zr 3d, Hf 4f, S 2p: fixed splitting and area ratio); **atomic-%** with editable RSFs (Scofield defaults). |
| **Figure composer** | Panels added from any tab | Multi-panel figure (1–4 columns) with (a)/(b)/(c), a/b/c, A/B/C or a) labels, exported as vector graphics. **Journal templates** ACS / RSC / Wiley (column width, font size, label style) with a check of width and minimum text size. |
| **Batch processing** | A folder (sub-folders included), an upload or the example folder | PXRD, FTIR, Raman, XPS, TGA, DSC or N₂ files → one figure per file, all overlaid or stacked, ONE preset (figure-settings JSON or house style), journal width, PNG/TIFF/PDF/SVG + CSV; ZIP and optional output folder; log of every file. |
| **SI export bundle** | Every figure shown in the session (sidebar: *Collect figures for the SI bundle*) | ZIP with figures (PDF/PNG/SVG), data CSVs, **SI_parameters.xlsx** (figure list, analysis + figure settings per figure, results per sample) and the defect SI table (DOCX/XLSX). |
| **Sample database** | Your own Excel sample list (a **Sample ID** column; optional Name, Paper/Project, Protocol family, Parent / starting material, Protocol basis, Main source (DOI), Protocol files) + your data folders; a synthetic example list is bundled | Local SQLite `reticula.db` in your Reticula data folder: sample IDs with protocol files, family, parent/children, source DOI; all results from the sample store; data files linked by the sample ID in their names (`<ID>-B<batch>_…`); notes. |
| **Structure (CIF)** | A CIF (SHELXL/Olex2, COD, CSD) | ORTEP (30/50/70 %), ball-and-stick, stick, space-filling with auto-placed labels; views along axes/[uvw]/(hkl)/best; ASU, molecules, packing, supercell; H-bonds with symmetry codes; SI tables (crystal data, bonds, angles, torsions with ESDs) as DOCX/XLSX; simulated PXRD with hkl and overlay; MOF pores (Di, Df, void fraction, ASA); Hirshfeld fingerprint (experimental). |
| **UV/Vis** | Implen NanoPhotometer Excel (.ids only if not encrypted), Cary / generic CSV, JCAMP-DX | Overlay, A ⇄ %T ⇄ %R, blank, smoothing, λmax; Beer–Lambert (LOD/LOQ, unknowns), assay presets (BCA, Bradford, Lowry, ninhydrin, OPA + DH); single-λ kinetics; Kubelka–Munk + Tauc band gaps; derivatives. |
| **Export guide** | – | Which file to export from each instrument (also in `docs/EXPORT_GUIDE.md`). |

## 4. Defect calculator (Zr/Hf/Ce/Th₆ nodes, paddlewheels, Zn₄O, Al-OH chains, M₃O trimers, MOF-74, Ti₈, ZIFs, custom nodes)

You enter the data once. The calculator runs it through **every method that the data allows** and shows them side by side:
- linkers per node L
- missing linkers x
- caps/modulator per node m
- Cl⁻ per node
- open sites
- minimal formula and molar mass
- a **consensus** value: the median of the quantitative methods, with min–max and SD

For each method the app shows:
- the equation, rendered in LaTeX
- a symbol table and the assumptions
- a **step-by-step substitution using your own numbers**
- pitfalls and the reasons why methods disagree
- the DOI

A PDF report (summary, chart, table, one page per method), a CSV table and a JSON input file can be downloaded.

| Key | Method | Data needed | Main reference (DOI) |
|---|---|---|---|
| a | TGA plateau, "M-400", ideal-normalised weight loss | TGA in air: plateau and end mass | Valenzano *Chem. Mater.* 2011, 10.1021/cm1022882; Shearer *Chem. Mater.* 2014, 10.1021/cm501859p; 2016, 10.1021/acs.chemmater.6b00602 |
| b | Lillerud/Shearer/Sannes **M-200**: mass at 200 °C + NMR ratios + Cl, solved together with charge balance (b), plus the missing-cluster (reo) model (b′) | TGA W(200 °C), end, NMR ratios, Cl | Sannes *Chem. Mater.* 2023, 10.1021/acs.chemmater.2c03744; Øien *Cryst. Growth Des.* 2014, 10.1021/cg501386j; Cliffe *Nat. Commun.* 2014, 10.1038/ncomms5176 |
| c | TGA plateau corrected with the ¹H NMR modulator:linker ratio (species that survive to the plateau) | TGA + NMR | Shearer 2016; Atzori *J. Phys. Chem. C* 2017, 10.1021/acs.jpcc.7b00483 |
| d | **Pulparayil Mathew/De Roo** experimental minimal formula: M_exp from the oxide residue, qNMR vs internal standard, Cl, water by difference | TGA (initial, end), qNMR with TMSP/standard, Cl | *Adv. Sci.* 2025, 10.1002/advs.202504713 |
| e | ICP metal : qNMR linker ratio (e), and ICP wt% mass balance (e′) | ICP (+ qNMR, TGA for the dry basis) | De Roo 2025; Sannes 2023 |
| f | Porosity: measured vs ideal pore volume or BET (ideal V calculated from the CIF with a 1.82 Å N₂ probe). **Semi-quantitative**, not part of the consensus. | V_total or S_BET, CIF | Wu *JACS* 2013, 10.1021/ja404514r; Vermoortele *JACS* 2013, 10.1021/ja405078u; Shearer 2016 |
| g | **Klet potentiometric titration**: μ₃-OH, M-OH₂ and M-OH from three equivalence points (found automatically from dpH/dV) | NaOH titration curve, sample mass, HCl added | Klet *J. Mater. Chem. A* 2016, 10.1039/C5TA07687K; Planas *J. Phys. Chem. Lett.* 2014, 10.1021/jz501521x |
| h | Open-site rule applied to every method: OS = (12 − CN) + t·x − Σz_i·n_i(caps) − n_Cl. For 12-connected nodes this reduces to OS = 2x − m − Cl. | – | Shearer 2016; De Roo 2025 |

**Hf and mixed Zr/Hf** are handled:
- The ICP Hf/(Zr+Hf) fraction sets the node mass and the MO₂ residue mass.
- **Extra-framework ZrO₂** in ZrO₂@Hf-MOF composites is removed from the TGA masses before methods a–d. The framework fraction comes from ICP Zr/Hf.

**Literature checks** (`tests/test_defects_literature.py`, all pass):

| Check | Published | This code |
|---|---|---|
| Ideal dehydroxylated UiO-66(Zr), plateau vs ZrO₂ (Shearer 2014/16) | 220.2 % | 220.2 % ✔ |
| Valenzano 2011: M of Zr₆O₆(BDC)₆ / Zr₆O₄(OH)₄(BDC)₆ | 1628 / 1664 g mol⁻¹ | 1628.0 / 1664.1 ✔ |
| Shearer-type defective example: plateau 200 % → L | 4.99 | 4.99 ✔ |
| De Roo 2025 MOF-808 worked example: M_exp from 45.6 % ZrO₂ | 1619.5 g mol⁻¹ | 1621.3 (within 0.15 %; the paper value corresponds to 45.65 %, i.e. rounding) ✔ |
| De Roo 2025 MOF-808 minimal formula: BTC 2, acetate 2, Cl 2.8, w 1.2, M_anhy 1341, n(H₂O) 15.5 | as listed | all reproduced ✔ |
| De Roo 2025 Table 3 masses: 1447.9, 1664.1, 1620.2, 1852.6 g mol⁻¹; 33.8 % Zr | as listed | reproduced; UiO-66 %Zr 29.55 % vs 29.6 % in the paper (rounding) ✔ |
| Sannes 2023 M-200 equations | no raw numbers in the main text | exact round trip (L, m, β) ✔ |
| Benzoate that survives 400 °C: (a) overestimates L, (c) corrects it | qualitative statement in Shearer 2016 | ✔ |
| Klet 2016 titration | raw volumes are only in the ESI | **convention check only** (pKa set, L/x bookkeeping), not a raw-data reproduction |
| ZrO₂@MOF composite correction; consensus logic | – | internal consistency ✔ |

On the synthetic examples the calculator recovers the compositions used to generate the TGA curves (`example_data/TGA/example_design_truth.json`):

| Sample | Designed L | Consensus L (range) |
|---|---|---|
| MOF-801(Hf) | 5.2 | 5.20 (5.02–5.55, 7 methods) |
| Sulfo-MOF-801(Hf) | 4.7 | 4.67 (4.50–4.70) |
| ZrO₂@Sulfo-MOF-801(Hf) | 4.7 | 4.69 (4.65–4.81) |

## 5. Limitations and caveats

- **Parsers built from published format descriptions.** They were tested on synthetic files written in those formats, not on files from your instruments. If an import fails, export plain text as described in `docs/EXPORT_GUIDE.md`. Two-column text always works.
- **STOE binary `.raw`** needs the optional `xylib-py`, which needs a compiler. Export `.xy` from WinXPOW instead.
- **OPUS binary** is read through brukeropusreader. It has not been tested on a real Vertex 70 file. Use `.dpt` if it fails.
- **MicroActive `.SMP/.QPS`** files are not readable. Use the Excel export.
- **VAMAS** support is minimal (single-block regions).
- **SEM pixel size.** The calibration panel compares the file metadata with the scale bar burned into the info strip (measured + OCR of its label) and the magnification. Without metadata the detected scale bar is used; if it cannot be read, click the two bar ends ("scale bar (mark ends)") and type the value. The 127 mm magnification estimate is only valid for Philips/FEI XL instruments.
- **The real 3Flex file from the v2 bug report was not available while building 2.0.** The N₂ fix was developed and tested on a synthetic MicroActive full-report export that reproduces the reported symptoms (`example_data/Sorption/MOF801Hf_3Flex_full_report_EXAMPLE.xls`). If your file still looks wrong, open *File structure* on the N₂ page and send that table.
- **Reported vs computed BET.** MicroActive usually fits a fixed window (0.05–0.30). For Type I MOFs the app uses the Rouquerol window, which lies at much lower p/p₀. The two numbers differ by design. The app does not tune its value to match the report.
- **Unit conversion** in Figure controls relabels the ticks only. Exported CSVs keep the original units.
- **Figure composer panels** use the default style of each figure; the composer has its own controls.
- **TEM readers** were tested on the openNCEM sample files (`tests/data`), not on files from your microscope.
- **BET.** Ultramicroporous MOFs (MOF-801, UiO-66) often cannot meet Rouquerol criterion 3 or 4. The app flags this. In that case report a "BET-equivalent area" and give the p/p₀ range used.
- **Defect methods f (porosity) and b′ (missing-cluster)** are semi-quantitative and are excluded from the consensus.
- **Klet titration** is only checked for its conventions; see the table above.
- A **NU-1000 CIF** is not bundled (the download failed). Load your own CIF; the framework constants for NU-1000 (CN 8, 2 linkers per node) are built in.
- The **TGA plateau** is picked automatically at the DTG minimum between 250 and 480 °C. Check it, and override it in the defect tab if the dehydroxylation and combustion steps overlap.
- **Arial** is used if it is installed. Otherwise Liberation Sans or DejaVu Sans is used, so check fonts in the SVG/PDF before submission.
