# Reticula 2.18.2

**Release date:** 2026-11-19

## Downloads

- Windows 10/11: `Reticula_Setup_2.18.2.exe` (not code-signed - SmartScreen: *More info* → *Run anyway*)
- macOS, Apple Silicon: `Reticula_2.18.2_macOS_arm64.dmg` (not notarised - first start: right-click → *Open*; see `README_macOS.txt`)
- `SHA256SUMS.txt` to verify the downloads.

## What changed

- **Zetasizer .dts (software 8.x)**: Z-average and PDI are read as the software stored them; checked against the Zetasizer values on two real files.
- **Reticon**: the editor loads when the page opens; the *Full screen* button is prominent; in full screen the top bar has a margin; the icon library can be widened by dragging, in the browser and in full screen.
- **FTIR / Raman**: tick *The file(s) contain several spectra* to read one spectrum per column of a CSV/TXT. OMNIC `.spa` files are now read.
- **Methods / export guide**: the SEM / TEM / EDX Methods paragraph names the JEOL JEM-1010 (bright-field TEM) and JEOL JEM-2200FS (HRTEM); the export guide has a row for each.
- **Simple mode, every technique**: the *Simple mode – all black* switch redraws every figure in a sleek style: one black symbol per sample for data points, one black line style per sample for curves (no markers), grey bars with a hatch per series, legends follow.
- **Rietveld**: optional instrument profile fixed from a Si SRM 640 standard (Advanced → *Instrument profile*). The standard's cell stays at the certificate value.

Full list: `app/CHANGELOG.md`.
