# Reticula 2.18.8

## Fixes
- Photoluminescence tab: the power-dependence panel no longer crashes with "Singular matrix" when all spectra have the same power; it asks for at least two different powers instead.
- Panel composition (PL figure and similar) no longer crashes with "'FigureCanvasBase' object has no attribute 'get_renderer'".

No change to computed results for valid input. Details: `app/CHANGELOG.md`, section 2.18.8.

## Checks
- New regression tests: 4 passed (they fail on 2.18.7).
- Related existing tests (power fit, composer, label move, plots, PL): 242 passed, 16 skipped.
