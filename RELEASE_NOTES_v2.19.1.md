# Reticula 2.19.1

## New methods (each with a known-answer test)
- KAS (Kissinger–Akahira–Sunose) isoconversional activation energy, next to Ozawa–Flynn–Wall.
- Vegard's law: composition from the lattice parameter, with optional bowing.
- Sauerbrey QCM mass per area from the frequency shift (5 MHz quartz: about 17.7 ng per Hz and cm²).
- Birch–Murnaghan equation of state (V0, B0, B0′) from pressure–volume data.
- Randles–Sevcik peak current and diffusion coefficient; specific capacity and coulombic efficiency for battery cycling.

## Reading files
- Heat flow in mW is converted to W; wavelength in µm, Å or m is converted to nm.
- Tab-separated exports with decimal commas are read correctly.
- Format detection before parsing (JCAMP-DX is detected but not read yet).

## Limits
- The new methods are library functions; their page controls come in a later round.
- Vegard and Sauerbrey are approximations (see the info document, Kinetics and adsorption).

## Checks
- 17 known-answer tests (`app/tests/test_v2191_known_answers.py`).
