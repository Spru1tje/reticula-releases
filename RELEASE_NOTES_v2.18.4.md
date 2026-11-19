# Reticula 2.18.4

**Release date:** 2026-11-19

## Downloads

- Windows 10/11: `Reticula_Setup_2.18.4.exe` (not code-signed - SmartScreen: *More info* → *Run anyway*)
- macOS, Apple Silicon: `Reticula_2.18.4_macOS_arm64.dmg` (not notarised - first start: right-click → *Open*; see `README_macOS.txt`)
- `SHA256SUMS.txt` to verify the downloads.

## What changed

- **Photoluminescence:** Edinburgh Instruments `.FS` files (single scans and "Multiple Scans") are read. A new standard **background subtraction** for visible and near-infrared spectra: automatic by default (signal-free windows found over the whole series, polynomial order chosen by the Bayesian information criterion, the same plan for every spectrum), with manual windows, several methods, a dark / background file and spike removal as options.
- **Photoluminescence – trustworthy thermometry:** a confidence verdict with a quality report, a manifold fit for crowded lanthanide bands, a full uncertainty budget (statistics, temperature sensor, every processing choice, runs) with resolution and calibration uncertainty, laser-heating and heating/cooling checks, a known-answer self-check and an analysis record (JSON) for the SI.
- **Photoluminescence – lifetimes:** decay curves with tail, rise-and-decay or reconvolution fits, average lifetimes and τ(T).
- **Photoluminescence – extras:** lamp-based spectral correction, excitation-line removal, settings from the files, saturation check, transition labels, dominant wavelength and colour purity.
- **Defect calculator:** the multi-faceted thermogravimetric method of Lázaro et al. 2020 (methods i1 and i2), including composites and charge balance.
- **PXRD – structure solution of an unknown compound:** a guided automatic route (cell → space group from the systematic absences → charge flipping → atoms from your formula → Rietveld refinement) with a progress bar, a log of every decision, an optional silicon reference for the zero shift and profile, and CIF, 3D view, pattern files and a report in one ZIP. A starting model for a careful refinement, not a finished structure.
- **PXRD:** Rietveld / Le Bail results give a 3D structure, CIF, `.xy` files, reflection list and a ZIP; glass-holder background and a measured Si standard for the instrument profile and zero shift; indexing rewritten for a higher M20.
- The info documents (Help and info, 📖 PDF button) of Photoluminescence, Defect calculator and PXRD are updated.

Full list: `app/CHANGELOG.md`.
