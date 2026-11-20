## 2.18.6 - 2026-11-20

### Photoluminescence: thermometry (ratiometric)
- **Poisson maximum-likelihood fit on raw counts** (new `mfs/plmle.py`): N1 ~ Poisson(B1 + D·S2), N2 ~ Poisson(B2 + S2), D = α·exp(−ΔE/k_B T). No background subtraction before the ratio, no logarithm, no dropped points. Uncertainty from the observed information (numerical Hessian, inverted), checked against the analytic Poisson Fisher bound.
- **Wired to the thermometry tab** (integration windows, counts on): the maximum-likelihood ΔE with its uncertainty is compared with the subtraction result in a new check "Poisson maximum likelihood (counts)".
- **Emission-masked baseline** (`mfs/plbackground.masked_poly`): when no signal-free window exists, the baseline is a polynomial of order 1 at most, fitted outside the emission (mask derived from the data by iteration). Default instead of the rubber band; the rubber band is an alternative. A safety check (overshoot above the spectrum by more than three noise units) falls back to the rubber band when the emission covers the whole range, and the notes say so.
- **Changed default (behaviour change):** without a signal-free window the program now uses the emission-masked polynomial instead of the rubber band. Measured on synthetic spectra (area error): sharp peaks on a slope rubber +78 % / masked +1 %; noisy spectra rubber +89 % / masked +0 %; broad band rubber +6 % / masked -2 %; concave background: both fail (+715 % / +485 %). A test of release 2.18.4 was adapted accordingly.
- **Analysis report (PDF)**: a third download next to the CSV and the JSON record (`mfs/plreport.py`): verdict, settings, equations, calibration figure, table, uncertainty budget, checks with advice, limitations.
- Help texts (`helpdocs.py`) updated for the three items.

### Measured (synthetic, known answer)
| Test | Result |
|---|---|
| Ratio, 300 runs, true ΔE 700 cm⁻¹: subtraction + ln-fit (released method) | bias −3,05 %; 33 % within 1σ, 53 % within 2σ |
| Ratio, 150 runs: Poisson maximum likelihood | bias +0,09 %; 70 % within 1σ, 95 % within 2σ; z-sd 1,06 |
| Broad emission band (area error): SNIP / arPLS / rubber band / masked polynomial order 1 | +13 % / −26 % / +5 % / −1,5 % |
| Masked polynomial with order 2 (not used) | −30 % (extrapolation into the band) |

### Known issues and open items
- The maximum-likelihood fit has no real-data validation yet: no online dataset with raw counts was usable (Sm²⁺ gated spectra are counts already background-corrected by the instrument; Yb³⁺ is OriginLab .opju and not readable; Pr³⁺ is normalised, without counts).
- Mn⁴⁺ (Science Advances table, Fig-3c): the table gives ΔE = 196–202 cm⁻¹ over every temperature range tested; the paper states 214 (fitted) and 226.5 cm⁻¹ (calculated from the ν₆ phonon energy). The cause is not identified.
- Dy³⁺ and YOF: no reference ΔE in the sources; the Boltzmann check correctly refuses a reliable gap for the Dy³⁺ ratio.
- Test suite (local): see the test summary before release. Known failures: `test_pagestate` (defect baseline, decision pending); browser end-to-end tests of batch A, B and C still use the pre-masking sample names; `test_pl_extras_panels_exist_and_run` failed once and passed on rerun.
- The PDF methodology document is not updated; it waits for approval.

### Known open items
- test_pagestate (2 cases): defect baseline decision pending (masked example names).
- Browser end-to-end tests of batch A and B still need the masked sample names (tests-extra workflow).
