# Reticula 2.18.0

**Release date:** 2026-11-19

## Downloads

- Windows 10/11: `Reticula_Setup_2.18.0.exe` (not code-signed - SmartScreen: *More info* → *Run anyway*)
- macOS, Apple Silicon: `Reticula_2.18.0_macOS_arm64.dmg` (not notarised - first start: right-click → *Open*; see `README_macOS.txt`)
- `SHA256SUMS.txt` to verify the downloads. Everything else: run from source (see the README).

## What is new

- **Schemes & Icons: 500 new vector icons** (1054 in the library). Every new icon has five variants, so six shapes per icon.
  All original, CC0, vector only. Families: spectra and analysis traces, lab vessels, instrument bodies, material morphologies,
  process arrows.
- **Scheme editor**: a floating bar above the selection, a status line that says what the next click or key does, Enter and
  Shift+Enter to add the next step, Repeat →, Number, Line up, Page centring, Copy/Paste style, Copy PNG, focus mode (F),
  template start points and a preview of the chosen template, locked items, an undo toast, a draft kept in the browser.
- **TGA**: a new plateau rule, *start of decomposition (EJIC 2020)*, based on Lázaro et al., Eur. J. Inorg. Chem. 2020
  (doi:10.1002/ejic.202000656). A point picked by hand on the graph still overrides the rule.
- **1H NMR**: *Simple mode (black)* for the final figure, the SI page PDF included.

## Fixed

- Flip in the scheme editor can be undone. Typing in the side panel keeps the focus. The align bar no longer shows without a selection.
- The batch-run test no longer assumes three example TGA files.

## Known limits

- Zetasizer 8.02 `.dts` files are not read yet (waiting for a sample file).
- The start-of-decomposition rule is not yet checked against hand-picked onsets on real data.
- Hover previews of templates are not possible in Streamlit; the preview shows after a template is chosen.

Full list: `app/CHANGELOG.md`.
