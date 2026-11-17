## Reticula 2.14.0

**Downloads**
- Windows 10/11: `Reticula_Setup_2.14.0.exe` (not code-signed - SmartScreen: *More info* → *Run anyway*)
- macOS, Apple Silicon: `Reticula_2.14.0_macOS_arm64.dmg` (not notarised - first start: right-click → *Open*; see `README_macOS.txt`). GitHub often has no free macOS runner: the `.dmg` is attached to this release as soon as one is available. **Not there yet? Build it yourself in ± 30 minutes:** download `Reticula_2.14.0_mac_build_kit.zip`, double-click `build_mac.command` (instructions: `BUILD_ON_MAC.txt`, Nederlands: `BUILD_ON_MAC_NL.txt`).
- `SHA256SUMS.txt` to verify the downloads. Everything else: run from source (see the README).

**New**
- **Photoluminescence**: spectra coloured by temperature, waterfall, emission map; thermal quenching (Eₐ, T₅₀, Arrhenius plot); ratiometric thermometry (Boltzmann fit, Sᵣ, σT) with ready-made A–D figures.
- **Combined figure**: put graphs of different techniques into ONE figure (A, B, C …) with a button under every graph.
- **¹H NMR**: proton labels as letters (Hₐ, H_b …), a cleaner figure (no clipped labels, integral numbers in their own band, optional), a bigger structure editor with all tools visible, ¹³C satellites, integral uncertainties and an assignment reliability column.
- **Zeta / DLS**: Excel import, volume / number distributions, Z-average / PDI / ζ in the figure, generated Methods text.
- **Particle sizer**: parallel on all CPU cores (≈ 5× faster for batches), calibrated edge for non-round particles, rods stay whole.
- **All figures**: font per figure and legend, legend-title size; the legend options now work without an extra tick; a test checks every option on every page.
- **N₂ / BET**: starts with the BET surface area; the Rouquerol details and plots are one click away.

Full list: `app/CHANGELOG.md`.
