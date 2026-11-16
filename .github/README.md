The macOS GitHub Actions workflow lives on `main` (`.github/workflows/macos.yml`)
because scheduled workflows only run from the default branch. This branch holds
the v2.11.0 source snapshot that workflow checks out (`ref: mac-build`).
