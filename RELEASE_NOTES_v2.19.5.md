# Reticula 2.19.5

## PL manifold fit: 12–40× faster, and more exact
The whole manifold analysis of the thermometry tab (line search, both shape models, robustness refits with other windows, the other model and other backgrounds), measured on the same computer:

| Data set | v2.19.4 | v2.19.5 | Speed-up |
|---|---|---|---|
| Test series (7 spectra × 376 points, 6 lines) | 10.7 s | 0.28 s | 39× |
| Crowded series (13 spectra × 710 points, 10 lines) | 73.7 s | 1.9 s | 40× |
| Long series (41 spectra × 376 points, 6 lines) | 94.9 s | 3.1 s | 31× |
| Er³⁺ example | 1.34 s | 0.11 s | 12× |
| Yb/Ho composite | 1.96 s | 0.14 s | 14× |

What changed:
- Exact Jacobian (variable projection) instead of finite differences: one evaluation per optimiser step instead of one per parameter.
- Areas by QR; non-negative least squares only when an area would be negative, first as a batched solution that is checked with the Karush–Kuhn–Tucker conditions.
- All spectra in one batched array operation; one BLAS thread for the small matrices; robustness refits in parallel; both shape models from the same warm start, at the same time.
- Warm starts (mean spectrum, previous trial, main solution) with a safety net: a warm fit that does not converge is redone from the ordinary start, and the drifting model may never end above the shared model.
- A fixed shape model on the page reuses the fit of the automatic choice instead of a new fit.

## More exact results
- A stuck fit of the old code is fixed (offset start): on the test series with the window 500–575 nm the shared-shape fit ended at rss 131 916 instead of 2 603. The test series now gets the correct shared shape, and the systematic errors no longer come from stuck refits.
- Ratio errors: the old finite-difference step was relative and gave rounding noise for parameters at zero; parameters held on a bound are now left out of the error propagation. Monte Carlo over noisy copies: the reported errors match the true scatter (see the info document).

## Info document
- The PL info document has a new section: every optimisation and why, the errors that were found, the performance comparison with figures, and what was tried and rejected.

## Checks
- Known-answer tests: `tests/test_v2195_plfit.py` (exact Jacobian against finite differences, the NNLS shortcut against NNLS, noiseless data recovered exactly, the stuck minimum, Monte Carlo error bars, warm against cold, serial against parallel).
- Excluded from the gate: the SEM scale-bar OCR check and its OCR tests.
