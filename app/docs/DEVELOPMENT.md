# Development notes

## Rules for every new page / module

1. **Every graph goes through `show()`** in `app.py` (never `st.pyplot` / `st.image` directly). That gives it the *Customise this figure*
   panel, the exports (PNG/TIFF/SVG/PDF/ZIP with provenance), the SI collection and the *Add to combined figure* button.
2. **Add the page to the figure-controls suite**: put its name in `PAGES` of `tests/test_figure_controls.py`. The suite runs every
   customise option (legend font size, line width, axis limits, …) on every page and fails when an option has no effect.
3. **Register the page** in `mfs/ui.py` (`NAV`), in `app.py` (page function) and give it a text in `mfs/summary.py`
   (`TECHNIQUE_METHODS`; a test fails when a technique has none).
4. **Put per-sample results in the sample store** (`samplestore.put`) so they appear in the sample summary table.
5. **Source files** are read through `files_input()` – their names then end up in the provenance record of each figure.
6. Optional / heavy tools (e.g. Rietveld, crystallite size) stay **hidden until the user clicks them open**.
7. Tests: `pytest tests/` – new behaviour needs a test on synthetic data with a known answer.
8. **Info documents**: after every big change (new module, method, default, file format or output) update the module's document in `mfs/helpdocs.py` (the *Help and info*
   module and the 📖 PDF button), in the same release. English; write abbreviations out first; mark unchecked advice as *General practice*. See also `AGENTS.md` at the repo root. Check the numbers with `pytest tests/test_helpdoc_numbers_sourced.py`: it counts the prose numbers that occur nowhere in the code, tests or notes, and fails above the ceiling `MAX_UNSOURCED` (a count only; the numbers are not shown).
