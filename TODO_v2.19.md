# To do for v2.19 (carried over from v2.18)

Accepted for v2.18 as they are (known limits, see CHANGELOG):
- TGA plateau (defect calculator): the rule is based on the main text of the source paper. The SI tables are not checked against it.
- Journal SI styles (ACS, Wiley-VCH, RSC, Elsevier): templates only. Not checked against the current author guidelines.

Open for the next version:
- Zetasizer .dts 8.02: Z-average and PDI are not calculated for 8.02 (the calculation gives impossible values). Needs an Excel/CSV export of the same measurements to verify. ζ is not read from .dts.
- Schemes & Icons:
  - Browser preview shows the default icon style; export uses the chosen style and theme.
  - Scheme panel in the combined figure is a 600 dpi bitmap, not vector.
  - Automatic routing works on an 8 px grid; check results.
  - Own icons are not yet a reusable library; they only live in the scheme where they were added.
  - Composer, methods schematic and element picker: first versions only.
  - Automated browser tests for the editor are missing.
- Test environment: Playwright is not installed on the development PC, so 8 browser tests (test_v212_e2e_browser, test_v29_e2e_browser) have not been run.
- PyInstaller build (Windows and macOS) has not been run with the v2.18 changes.
- Rename of technical names to Figuryte, with the release: installer file name (Reticula_Setup), macOS app name and bundle identifier, data folder Documents\Reticula, GitHub release names and update check. Changing them breaks existing installs and the update check, so they need a planned migration.
