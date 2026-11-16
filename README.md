# Reticula

**Publication figures and analysis for framework materials** (MOFs, COFs, zeolites, …) - one local app, your data never leave your computer.

PXRD · FTIR / Raman · ¹H NMR · N₂ sorption (BET / Rouquerol, PSD) · TGA · DSC · SEM / TEM particle sizing · XPS · UV/Vis · ICP / XRF · kinetics ·
crystal structures from CIF · defect calculator for MOFs · figure composer · batch processing · SI export bundle.
Journal presets for ACS, RSC, Wiley, Nature, Elsevier and more; colour-blind-safe palettes; SVG / PDF / TIFF / PNG export.

## Download

Get the latest installer from **[Releases](../../releases/latest)**:

| System | File |
|---|---|
| Windows 10/11 (64-bit) | `Reticula_Setup_<version>.exe` |
| macOS (Apple Silicon) | `Reticula_<version>_macOS_arm64.dmg` |

Both are unsigned: Windows SmartScreen → *More info* → *Run anyway*; macOS → right-click → *Open* the first time. Details: `pkg/README_Windows.txt`, `pkg/README_macOS.txt`.
Intel Macs and Linux: run from source (below).

## Run from source

```bash
pip install -r requirements.txt        # Python 3.11-3.13
cd app && streamlit run app.py
```
(or double-click `app/run_windows.bat`, `app/run_mac.command`, `app/run_linux.sh`). Tick *Use the bundled EXAMPLE data* on any page to try it.

## Repository layout

| Path | Content |
|---|---|
| `app/` | the program: `app.py` (start page, sidebar, helpers), `mfs/` (readers, analyses, plots; `mfs/pages/` = one module per technique page), `example_data/` (synthetic), `tests/`, `docs/`, `tools/` |
| `mfs_launcher.py` | desktop launcher of the installers (starts the local server and the browser) |
| `pkg/` | PyInstaller specs, Inno Setup script, dmg scripts, build helpers |
| `.github/workflows/` | `ci.yml` (tests on every push) and `release.yml` (tag `vX.Y.Z` → tests → Windows + macOS builds → GitHub release) |

User guide and technique reference: [`app/README.md`](app/README.md) · what's new: [`app/CHANGELOG.md`](app/CHANGELOG.md) · which file to export from each instrument: [`app/docs/EXPORT_GUIDE.md`](app/docs/EXPORT_GUIDE.md).

## Development

```bash
pip install -r requirements-dev.txt
cd app && python -m pytest -q
```
Release: bump `__version__` in `app/mfs/__init__.py`, add `app/CHANGELOG.md` and `RELEASE_NOTES_v<version>.md`, then `git tag v<version> && git push origin v<version>`.
The workflow refuses a tag that does not match `__version__`.

## License

CC0 1.0 (public domain dedication) for the code in this repository. Bundled third-party material keeps its own terms where it has any
(for example the DFT kernel in `app/references/kernels` is MIT-licensed, see the licence file there, and the CIFs come from the Crystallography Open Database).
Check the provenance of the data in `app/mfs/data` and `app/tests/data/nmr` before redistributing them separately.
