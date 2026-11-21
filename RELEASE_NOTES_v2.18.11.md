# Reticula 2.18.11

## New methods
- Zeta potential from electrophoretic mobility (Smoluchowski and Hückel), with a validity warning.
- Avrami (JMAK) crystallisation fit: n, k and half-time.
- Toth and Dubinin-Radushkevich isotherm fits; Akaike information criterion for choosing between models.
These are in the code with tests; the pages do not use them yet.

## Instrument files
- Heating rate and concentration units are converted automatically; key/value metadata above the table is read.

## Checks
- Known-answer tests for every new method; navigation and icon checks.
- Details: `docs/suggestions/2026-11-21-r2.md`.
