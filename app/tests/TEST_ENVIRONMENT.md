# Test environment for the TEM module

Documentation only. It describes what the CI totals depend on; it changes no test.

- Command (CI): from `app/`, with `MFS_NO_UPDATE_CHECK=1 PYTHONUTF8=0 PYTHONPATH=.`:
  `python -m pytest -q -rs -p no:cacheprovider --ignore-glob='tests/*_e2e_browser.py' --ignore=tests/test_downloads_browser.py`
- Tesseract: needed for real OCR tests. The code looks for it in this order: `$MFS_TESSERACT`, the copy bundled with the
  Windows app, `tesseract` on PATH, then the default install folders (Windows: `C:\Program Files\Tesseract-OCR`;
  macOS: Homebrew and MacPorts). Measured with Tesseract v5.4.0. The OCR tests in `test_v2222_*` replace OCR with a fake
  and run without it.
- `threadpoolctl`: optional. Without it `test_v2195_plfit` is skipped (`could not import 'threadpoolctl'`). Measured with
  Python 3.12.0 (Windows venv).
- Private measurement files: `MFS_TEST_DTS`, `MFS_TEST_DTS_A`, `MFS_TEST_DTS_B`, `MFS_TEST_DTS_ZETA`, `MFS_TEST_FS`.
  Without them the tests that need them are skipped (`-rs` lists the reason). The repository holds no such data.
- Totals: the number of passed tests depends on which branch is checked out, because each TEM branch adds tests.
  Record the branch, the commit hash, the totals and the skip count with every result.
- Skip count: 50 in the TEM-branch CI logs (`ci_*.log`, every TEM branch checked so far); the reasons are private data files,
  missing optional packages, and missing CIF/sample files. This count is for the TEM logs only, not for the whole suite.
