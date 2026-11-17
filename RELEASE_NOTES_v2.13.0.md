## Reticula 2.13.0

**Downloads**
- Windows 10/11: `Reticula_Setup_2.13.0.exe` (not code-signed - SmartScreen: *More info* → *Run anyway*)
- macOS, Apple Silicon: `Reticula_2.13.0_macOS_arm64.dmg` (not notarised - first start: right-click → *Open*; see `README_macOS.txt`). GitHub often has no free macOS runner: the `.dmg` is attached to this release as soon as one is available. **Not there yet? Build it yourself in ± 30 minutes:** download `Reticula_2.13.0_mac_build_kit.zip`, double-click `build_mac.command` (instructions: `BUILD_ON_MAC.txt`, Nederlands: `BUILD_ON_MAC_NL.txt`).
- `SHA256SUMS.txt` to verify the downloads. Everything else: run from source (see the README).

**New**
- **¹H NMR page**: Bruker raw import, automatic processing, structure editor, shift prediction, automatic peak → atom assignment, impurity check, qNMR, labelled figure, SI export.
- **N₂ / BET**: redesigned graphs, BET area range over all fitting windows, ~11× faster Rouquerol search.
- **Made for everyone's data**: no lab-specific colour/name defaults, your own sample list in the Sample database, molar ratios for any elements in the ICP/XRF page.
- Code: page code split into `mfs/pages/`; releases are built and self-tested by GitHub Actions.

Full list: `app/CHANGELOG.md`.
