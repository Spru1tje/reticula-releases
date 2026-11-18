## Reticula 2.15.1

**Downloads**
- Windows 10/11: `Reticula_Setup_2.15.1.exe` (not code-signed - SmartScreen: *More info* → *Run anyway*)
- macOS, Apple Silicon: `Reticula_2.15.1_macOS_arm64.dmg` (not notarised - first start: right-click → *Open*; see `README_macOS.txt`); the Mac build kit is attached as well.
- `SHA256SUMS.txt` to verify the downloads. Everything else: run from source (see the README).

**New**
- **Rietveld**: preferred orientation (March–Dollase); zero shift and displacement are no longer refined together (biased the cell); validated on a real laboratory pattern against GSAS.
- **Photoluminescence**: spectral-response correction, temperature uncertainty from repeat measurements, power dependence (slope n), CIE 1931 colour and colour temperature.
- **Dossier per material**: one ZIP per material with its figures, data (one Excel), methods text and results.
- **Detail level** switch: a clear coloured box at the top of the sidebar (Standard / Advanced / Developer).
- **macOS**: the app builds and passes its self-test again (v2.15.0 needed a manual fix).

Full list: `app/CHANGELOG.md`.
