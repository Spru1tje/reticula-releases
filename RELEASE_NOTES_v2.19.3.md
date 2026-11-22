# Reticula 2.19.3

## CIE chart
- sRGB colour gamut drawn; points inside it are shown in their own colour, points outside are marked and listed as "not showable".

## Speed
- Faster start-up (scipy loaded on first use in the analysis module).
- Pseudo-Voigt for all bands in one array operation (same numbers).
- The manifold fit is still slow (about 85 s on the composite sample); an analytic Jacobian is the next step.

## Checks
- Unit suite: 733 passed, 50 skipped; the excluded SEM OCR test failed (not part of the gate).
