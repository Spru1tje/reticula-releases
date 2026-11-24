# Reticula 2.21.0

Feature release: structure solution from powder data (experimental), PL publication checks and PDF.

- **Structure solution (experimental)**: a net matcher picks fcu, pcu, dia or tbo from the indexed reflections and the
  charge-flip density peaks. On simulated UiO-66 data it selects fcu and recovers 100 % of the carbons within 0.7 Å.
  An info document with six figures explains the method and its limits (Help page).
- **PL**: publication checks for ratiometric thermometry (model choice by AIC/BIC, validation against a reference
  temperature, response correction) and a PL report PDF with methods and figures.
- **Speed**: the Rietveld peak profile runs in float32 (lead-sulfate refinement 89 s → 48 s, same Rwp to 2e-6).

Known limits:
- The structure solver is validated on simulated data only. HKUST-1 and MIL-53 are not finished.
- The tbo net is taken from the HKUST-1 CIF, so the HKUST-1 net choice is not independent.
- A Zetasizer `.dts` test on a local file gives z-averages 1–2 nm above the value the test expects; under investigation.

Full list: app/CHANGELOG.md.
