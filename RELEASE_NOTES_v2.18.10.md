# Reticula 2.18.10

## Navigation
- The technique list is now grouped into collapsible categories, with an icon for every page. The current page's category opens by itself, so the sidebar no longer needs long scrolling.

## Instrument files
- Automatic reading of instrument text exports: metadata lines are skipped, the delimiter and decimal comma are detected, and recognised units (°C, K, mW, mg, min, mbar, and more) are converted to base units. Unknown units are reported, not guessed.

## Methods (available from code; pages to follow)
- Kissinger and Ozawa-Flynn-Wall activation energies; Sips isotherm fit.

## Checks
- Unit tests: 674 passed, 50 skipped; one excluded SEM OCR test failed (not part of the gate).
- Self-test: 40/41 with the excluded SEM OCR check; exit 0 for the gate.
