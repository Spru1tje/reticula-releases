## 2.18.5 - 2026-11-20

### Photoluminescence: ratiometric thermometry
- **Empirical noise when the counting error is too small.** When the reduced chi-squared of the Boltzmann fit is above 10, the counting error is no longer used as a weight. The fit is repeated with a constant relative error equal to the measured scatter of the ratio, and the report shows a "Ratio noise" warning. Before this change the error bars were scaled up by a factor of about 10^8 on the anti-Stokes test data.
- **Birge scaling only when significant.** The errors are multiplied by the square root of chi-squared only when chi-squared exceeds 1 by more than two standard deviations of chi-squared (threshold 1 + 2 sqrt(2/dof)).
- **Boltzmann-law check (R-squared).** R-squared below 0.9 fails, 0.9 to 0.98 warns, 0.98 and above passes. A ratio that does not follow the law is reported as having no reliable energy gap.
- **Simulation check renamed.** The known-answer self-check is now a simulation consistency check. It tests that the analysis recovers a gap it was given; it does not test accuracy.
- **Ratio error includes shape parameters.** The covariance of the line-shape parameters is propagated into the ratio error.
- **Fluoracle export fixes.** Trailing separators are accepted and rows without numbers are skipped.

### Validation (documented, not a release claim of accuracy)
- Visible LIR of the YOF data set reproduces an independent regression (0.0562 eV against 0.0558 eV).
- Mn4+ ratio table from the authors' supplement: the module gives 195.8 cm-1 (statistical error 1.8 cm-1) against 214 and 226.5 cm-1 quoted in the paper. The difference is not resolved.
- Cr3+ MOF relative sensitivity at 160 K: 2.79 %/K against 2.84 %/K in the paper.
- Dy3+ ratios fail the Boltzmann-law check; the source manuscript reports sensitivity peaking at 193 K, which a single gap cannot produce.
- The PDF "Reticula PL mathematics and implementation" is updated to match the method.

### Other
- Example sample names are masked: MOF-801(Hf) becomes Sample A, Sulfo-MOF-801(Hf) Sample B, ZrO2@Sulfo-MOF-801(Hf) Sample C. The example files are renamed to match.
- DSC examples and the DSC generator are unchanged from 2.18.4.
- Known open item: two defect-calculator snapshot tests (test_pagestate) depend on which example file is linked to which sample; decision pending.