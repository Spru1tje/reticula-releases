# Reticula v2.14 – plan and status

## 0. Check: is everything from the v2.13 next-steps in v2.13?
| Item (RETICULA_NEXT_STEPS.md) | Status |
|---|---|
| ¹H NMR page (Bruker import, auto processing, Ketcher, shift prediction, assignment, impurities, figure, qNMR, SI export, digestion sub-mode) | in the repo; checked in the app (aspirin example: 6/6 H groups assigned, structure can be redrawn / replaced) |
| BET/N₂ redesign, PSD inset removed | in the repo |
| ¹³C/DEPT hidden in the UI | yes |
| 1. Final pass NMR page against TopSpin on own Phobos/Atlas spectra | **open – needs your real spectra** |
| 2. Tests | done (v2.13: 395 pass locally; CI on 3 OS) |
| 3. Release notes + version bump | done (RELEASE_NOTES_v2.13.0.md, CHANGELOG, `__version__`) |
| 4. Windows installer | built by GitHub Actions on the tag (local PyInstaller build + self-test 39/40, OCR needs Tesseract) |
| 5. GitHub release v2.13.0 + SHA256SUMS | pending: tag after green CI |
| "Run from source needs Python 3.11–3.12" | now 3.11–3.13 (run_windows.bat picks it) |
| Mac: scheduled 30-min retry workflow | replaced by `release.yml` (tag → tests → Windows + macOS build → release) |
| Open inputs: Phobos/Atlas spectra, UV-Vis Excel, real .rasx, raw particle diameters, Hf-MOF-801 quantities | still open |
| Later ideas: ¹³C/DEPT, ³¹P/¹⁹F, 2D NMR, ICP/XRF Zr/Hf/S special case, data-editor tables in projects/undo | moved to the v2.15 list below |

## 1. v2.14 scope (your list)
| # | Request | Status |
|---|---|---|
| 1 | Zeta/DLS tab updated | **done (first version)** |
| 2 | Upload Excel exports of the Zetasizer Ultra | **done** (CSV, .xlsx, .xls; header row found automatically) |
| 2b | `.dts` files | **needs a sample file.** `.dts` is the *Classic* Zetasizer (Nano) software format; ZS Xplorer (Ultra) does not read it either (Malvern release notes). Reader to follow when a file is available; until then export to Excel/CSV |
| 3 | PDI and Z-average (and ζ) in the DLS figure with tick boxes, movable | **done** (coloured text per sample, ± SD option, drag with "Move / edit labels on the figure") |
| 4 | Legend movable | **done** (the v2.12 label mover moves legends; legend tick boxes added to the zeta and DLS figures) |
| 5 | Legend-title font size | **done** (Customise → Legend → "Legend-title size") |
| 6 | More font types | **done** (per figure and for the legend; installed extra families are offered, e.g. Verdana, Georgia, Palatino, Garamond, Consolas) |
| 7 | Particle sizer progress bar with time left | **done** (bar per batch: image k of n, seconds elapsed and left; the estimate learns from earlier runs in the session) |
| 8 | Ar and CO₂ sorption tabs | **designed below – needs example files** |

## 2. Ar / CO₂ sorption – proposal
Instruments: 3P (Ar), Quantachrome Autosorb 6100 and Autosorb iQ (physisorption).
Open-source readers show the layouts that exist: Quantachrome ASiQwin **text/Excel export** (header block with analysis gas, bath temperature,
sample weight, outgas, then p/p0 and volume @ STP columns); 3P **.xlsx report**. The native Quantachrome `.qps` is binary (no public spec found) – use the export.
*I will not guess columns: please send one real export of each instrument (Ar 87 K, CO₂ 273 K, N₂ 77 K).*

One "Gas sorption" page with a gas selector (N₂ 77 K / Ar 87 K / CO₂ 273 K / CO₂ 298 K) and the same flow for all gases:
1. Load files → sample list (existing).
2. Results table that is a paper table (copy as Word/Excel/LaTeX): S_BET (+ fitted p/p0 range, C, Rouquerol criteria), V_total (Gurvich, p/p0 = 0.95 / 0.99), V_micro (t-plot / DR), median pore width, uptake at chosen pressures.
3. Figures: isotherm (ads/des, log or linear), BET plot with range, t-plot/αs, PSD (differential + cumulative), CO₂: DR plot, uptake bar chart.
4. Gas-specific: **Ar 87 K** – BET with σ(Ar) = 0.142 nm², NLDFT/QSDFT pore size, HK/SF, t-plot, hysteresis type (IUPAC). **CO₂ 273 K** – Dubinin–Radushkevich ultramicropore volume, NLDFT CO₂ PSD (0.35–1.5 nm), uptake at 0.15 / 1 bar (mmol/g, wt %, cm³/g), Henry constant. **Several temperatures** – isosteric heat Q_st (Clausius–Clapeyron / virial), and with an N₂ isotherm IAST CO₂/N₂ selectivity.
5. "Methods paragraph" generated from the settings (instrument, degassing, model, range, kernel) + the figure/table in one click for the SI bundle.
Kernels: the MIT-licensed DFT kernel already bundled is N₂/77 K – Ar/87 K and CO₂/273 K kernels have to be added from a source with a compatible licence (to be checked).

## 3. Suggestions (my input)
1. **Gas sorption as ONE page** (above) rather than separate Ar and CO₂ tabs: same table/figure logic, fewer places to maintain.
2. **Methods + SI text generator** for every technique (what was done, with which settings) – saves the most time when writing.
3. **Everything-in-one-workbook export** (all figure data as Excel sheets, Origin-friendly) next to the figure files.
4. **Figure checks per journal** before export (minimum font size at final width, line width, colour-blind check) with a one-click fix.
5. **Zeta/DLS**: volume/number distributions, lognormal fit, correlation function, zeta-vs-pH curve with isoelectric point (MPT titration data).
6. **Data audit**: flag outlier records, suspicious PDI (> 0.7), negative isotherm points, non-monotonic adsorption branches before they end up in a figure.
7. **Project recipe**: one `.reticula` file that restores data, settings and figure layout (exists since v2.11) + "re-run with new data" for the same sample set.
8. v2.15 candidates: ¹³C/DEPT, ³¹P/¹⁹F, 2D NMR, generic ICP/XRF page, data-editor tables in projects/undo.
