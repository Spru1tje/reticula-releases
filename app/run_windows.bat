@echo off
REM Reticula - first run creates .venv and installs; later runs just start the app
cd /d "%~dp0"
if not exist .venv\Scripts\python.exe (
  py -3 -m venv .venv 2>nul || python -m venv .venv
  if errorlevel 1 (echo Python 3.10+ not found - install it from https://www.python.org and tick "Add to PATH" & pause & exit /b 1)
  .venv\Scripts\python.exe -m pip install --upgrade pip
  .venv\Scripts\python.exe -m pip install -r requirements.txt
)
.venv\Scripts\python.exe tools\launch.py
pause
