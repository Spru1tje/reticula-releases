# Reticula 2.19.4

## New methods and checks (ten additions)
- Hill plot for cooperative binding: theta = x^n / (K^n + x^n), with K, n and their errors.
- Michaelis–Menten fit of Vmax and Km, with the Lineweaver–Burk (double-reciprocal) estimate for comparison.
- Friedman differential isoconversional activation energy (no integral approximation).
- Heating-rate series check before Kissinger or KAS: warns on a narrow or non-monotone heating-rate spread.
- Kissinger against Friedman: warning when the two activation energies differ by more than 15 %.
- Gurvich total pore volume from N₂ uptake near P/P0 = 0.95.
- TGA: uncertainty of each mass-loss step from the balance resolution (√2 × resolution).
- Problem report: names the detected file family of each upload (no data).
- Missing Tesseract: the message names the install command and `$MFS_TESSERACT`.

## Sidebar
- Each technique button shows its icon at the left edge, recoloured for the selected row.

## Checks
- Known-answer tests for every new method: `tests/test_v2192_items.py`.
- Sidebar icon tests: `tests/test_v2192_sidebar_icons.py`.
- Excluded from the gate: the SEM scale-bar OCR check and its OCR tests.
