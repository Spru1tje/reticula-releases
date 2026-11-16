# Reticula macOS build source (v2.11.0)

This branch holds a **snapshot of the v2.11.0 source** plus the files needed
to build the unsigned macOS app on GitHub Actions (`macos-14`, Apple Silicon).

The default `main` branch stays a releases-only landing page. Do not develop
here; later versions will get their own snapshot / tag.

- Workflow: `.github/workflows/macos.yml`
- App source: `app/` (package `mfs`, version 2.11.0)
- Launcher: `mfs_launcher.py` (macOS data dir: `~/Library/Application Support/Reticula`)
- PyInstaller spec: `pkg/mfs_mac.spec`

Run the workflow with **Actions → macOS build → Run workflow**, or push to
`mac-build`.
