# Agent guide: Reticula release workflow

Repo: the Reticula releases repo (`app/`, `pkg/`, `.github/workflows/`). Read this before touching a release.

## Rules
- Release only after the user says so.
- Never commit the test-modified files: `app/reticula.db`, `app/sample_store.json`, `app/run_windows.bat`.
- For every big change (new module, method, changed default, file format, output), update the info documents in `app/mfs/helpdocs.py` in the same release. Write in English, spell out each abbreviation at first use, and mark unchecked advice as *General practice*. A test checks this.

## Steps for a release
1. Run the tests: `cd app && python -m pytest -q`. Fix failures before continuing.
2. Set `__version__` in `app/mfs/__init__.py` (for example `"2.18.7"`).
3. Write `app/CHANGELOG.md` and `RELEASE_NOTES_v<version>.md` (in the repo root).
4. Update the info documents in `app/mfs/helpdocs.py` if the change is big.
5. Commit on top of the last tag. Leave out the three files listed above.
6. Create and push the tag: `git tag v<version> && git push origin v<version>`.
7. Check the `release` workflow run in GitHub Actions. Expect these jobs in order:
   - `prepare`: fails if the tag does not equal `v` + `__version__`, or if `RELEASE_NOTES_v<version>.md` is missing.
   - `tests`: runs `ci.yml` (pytest on Ubuntu for Python 3.11, 3.12, 3.13 and Windows for 3.12).
   - `windows`: PyInstaller build, bundled Tesseract, headless self-test, Inno Setup installer, SHA256 checksum.
   - `macos`: currently disabled (`if: false`, v2.18.6).
   - `publish`: creates or updates the GitHub release with the installer, `SHA256SUMS.txt`, `README_macOS.txt` and the macOS build kit zip.
8. Confirm the release page has the installer, `SHA256SUMS.txt` and the build kit.

## Other workflows
- `ci.yml`: blocking tests on pushes to `main` or `fix/**` and on pull requests. It also starts the app from source on Windows and Ubuntu.
- `ci-extra.yml`: macOS tests and Playwright browser tests. Informational; does not block a release.
- `macos-build.yml`: manual only (Actions, Run workflow). Builds the Apple Silicon `.dmg` from `main`. Set `attach_tag` to attach it to an existing release, for example `v2.18.6`.

## Known gaps
- No macOS `.dmg` is built in a release until the `macos` job is re-enabled. Use `macos-build.yml` with `attach_tag`, or build on a Mac with the build kit.
- The example tag in the header comment of `release.yml` is stale (`v2.13.0`). It is only a comment.
- `macos-build.yml` matches `__version__ = "` with exact spacing. `release.yml` allows any spacing. Keep the single-space format in `app/mfs/__init__.py`.

## If a step fails
- Report the exact failing job and its log line. Do not retry blindly.
- Do not re-tag a version that is already published. Use a new patch version instead.
