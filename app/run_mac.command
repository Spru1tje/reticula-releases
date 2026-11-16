#!/usr/bin/env bash
# Reticula - first run creates .venv and installs; later runs just start the app
set -e
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
# macOS: if python3 is missing install it from python.org or "brew install python"
if [ ! -d .venv ]; then
  "$PY" -m venv .venv
  .venv/bin/python -m pip install --upgrade pip
  .venv/bin/python -m pip install -r requirements.txt
fi
exec .venv/bin/python tools/launch.py "$@"
