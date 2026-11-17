@echo off
REM Reticula - first run creates .venv and installs; later runs just start the app.
REM Python 3.11-3.13 is used (the scientific packages have no ready-made wheels for newer versions yet).
REM Set RETICULA_PY=3.12 (for example) to force one version.
cd /d "%~dp0"
if not exist .venv\Scripts\python.exe (
  set "RPY="
  if defined RETICULA_PY (
    py -%RETICULA_PY% -c "import sys" >nul 2>nul && set "RPY=py -%RETICULA_PY%"
  )
  if not defined RPY for %%V in (3.13 3.12 3.11) do if not defined RPY (
    py -%%V -c "import sys" >nul 2>nul && set "RPY=py -%%V"
  )
  if not defined RPY (
    python -c "import sys; sys.exit(0 if (3,11) <= sys.version_info[:2] <= (3,13) else 1)" >nul 2>nul && set "RPY=python"
  )
  if not defined RPY (
    echo Python 3.11, 3.12 or 3.13 not found - install it from https://www.python.org and tick "Add to PATH"
    pause
    exit /b 1
  )
  %RPY% -m venv .venv
  if errorlevel 1 (echo Could not create the virtual environment & pause & exit /b 1)
  .venv\Scripts\python.exe -m pip install --upgrade pip
  .venv\Scripts\python.exe -m pip install -r requirements.txt
  if errorlevel 1 (echo Installing the packages failed - see the messages above & rmdir /s /q .venv & pause & exit /b 1)
)
.venv\Scripts\python.exe tools\launch.py
pause
