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

## Measured CI totals (full suite, CI command above)

Measured on this machine; each run uses the command above, with `-rs`. Totals differ because each branch adds tests.

| Tip | Passed | Skipped | Log |
|---|---|---|---|
| tem/regressie-tests fb8ab49 | 841 | 50 | ci_fb8ab49.log |
| tem/suggesties 7cb8768 | 850 | 50 | ci_7cb8768.log |
| tem/ruis-sweep 5ae2356 | 853 | 50 | ci_5ae2356.log |
| tem/ocr-meta-warn b70b288 | 855 | 50 | ci_b70b288.log |
| tem/raw-route 104c8bb | 860 | 50 | raw_route_full_ci.log |
| integration set (10 tips, local merge) | 888 | 50 | ci_final_set.log |

The totals and skip counts come from the final pytest line of each log. The logs are kept by the author, outside the repository.

## Skip reasons (from the same logs)

- Private measurement files not in the repository (the `MFS_TEST_DTS*` and `MFS_TEST_FS` tests, and the real-data tests).
- Optional packages not installed (`threadpoolctl` for `test_v2195_plfit`).
- User-supplied CIF files not on the machine.

The skip count is 50 in every run in the table. The reasons were read from the logs of the runs listed; the three groups above are what those logs show, not a per-tip check of every skip.
