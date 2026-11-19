# Instructions for anyone (person, bot, other PC) working on Reticula

## Standing rule: keep the info documents current
After every **big change** – a new module, a new analysis method, a changed default, a new file format, a new output – update the info documents of the affected modules **in the same release**,
without being asked:

- `app/mfs/helpdocs.py`: one document per module (what it does, terms, tutorial, formulas, what is checked, limits). They are shown in the module *Help and info* and downloaded as a PDF
  (`app/mfs/helppdf.py`, Reticula colours); every module also has a 📖 info-PDF button.
- Language: English. Write every abbreviation out where it is first used, then abbreviate. Mark advice that the program does not check as *General practice*.
- A test (`tests/test_v2182.py::test_help_documents_cover_every_module_and_render_as_pdf`) fails when a module has no document or no *Terms* section.
- Also update `app/CHANGELOG.md` and `RELEASE_NOTES_v<version>.md`.

## Other rules
- See `app/docs/DEVELOPMENT.md` for the rules for a new page (figures through `show()`, registration, tests on synthetic data with a known answer).
- A release is a tag `vX.Y.Z` equal to `__version__` in `app/mfs/__init__.py`; `.github/workflows/release.yml` builds it. Do not commit the files that tests modify
  (`app/reticula.db`, `app/sample_store.json`, `app/run_windows.bat`).
