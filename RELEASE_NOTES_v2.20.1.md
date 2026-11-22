# Reticula 2.20.1

**Release date:** 2026-11-22

## Downloads

- Windows 10/11: `Reticula_Setup_2.20.1.exe` (not code-signed - SmartScreen: *More info* → *Run anyway*)
- macOS, Apple Silicon: `Reticula_2.20.1_macOS_arm64.dmg` (not notarised - first start: right-click → *Open*)
- `SHA256SUMS.txt` to verify the downloads.

## What changed

- **Rietveld refinement is faster, with identical results**: parallel finite-difference Jacobian and an exact cache of the structure factors.

## Known limits

- An analytic Jacobian (a further several-fold speed-up) is not in this release: it would change the numbers slightly.

Full list: `app/CHANGELOG.md`.
