"""Block a commit that stages files which must never be committed.

Forbidden: the generated sample store, test data, and core.py/nav_icons changes made by this icon work.
Exception: if a forbidden path must change on purpose, remove it from FORBIDDEN in this file in the same
commit and say why in the commit message; do not work around the guard any other way.
Scope: the guard checks only the files staged with `git add`. It does not look at unstaged or untracked files.
Usage (from the repo root, before committing):  python app/tools/check_staged.py
Exit code 1 if a forbidden file is staged; 0 otherwise.
"""
import fnmatch
import subprocess
import sys

FORBIDDEN = ["app/sample_store.json", "app/tests/data/*", "app/tests/data/**", "app/mfs/reticon/core.py",
             "app/assets/nav_icons/*"]


def blocked(names: list[str]) -> list[str]:
    return sorted({n for n in names for pat in FORBIDDEN if fnmatch.fnmatch(n, pat)})


def main() -> int:
    out = subprocess.run(["git", "diff", "--cached", "--name-only"], capture_output=True, text=True, check=True)
    names = [line.strip() for line in out.stdout.splitlines() if line.strip()]
    bad = blocked(names)
    if bad:
        print("BLOCKED: forbidden files are staged: " + ", ".join(bad))
        return 1
    print(f"ok: {len(names)} staged file(s), none forbidden")
    return 0


if __name__ == "__main__":
    sys.exit(main())
