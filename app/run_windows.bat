@echo off
REM Reticula - first run creates .venv and installs; later runs just start the app.
REM Python 3.11-3.13 is used (the scientific packages have no ready-made wheels for newer versions yet).
REM Set RETICULA_PY=3.12 (for example) to force one version.
cd /d "%~dp0"
if exist .venv\Scripts\python.exe goto run
set "RPY="
if defined RETICULA_PY py -%RETICULA_PY% -c "import sys" >nul 2>nul && set "RPY=py -%RETICULA_PY%"
if defined RPY goto make
py -3.13 -c "import sys" >nul 2>nul && set "RPY=py -3.13"
if defined RPY goto make
py -3.12 -c "import sys" >nul 2>nul && set "RPY=py -3.12"
if defined RPY goto make
py -3.11 -c "import sys" >nul 2>nul && set "RPY=py -3.11"
if defined RPY goto make
python -c "import sys; sys.exit(0 if (3,11) <= sys.version_info[:2] <= (3,13) else 1)" >nul 2>nul && set "RPY=python"
if defined RPY goto make
echo Python 3.11, 3.12 or 3.13 not found - install it from https://www.python.org and tick "Add to PATH"
pause
exit /b 1
:make
%RPY% -m venv .venv
if errorlevel 1 goto fail_venv
.venv\Scripts\python.exe -m pip install --upgrade pip
.venv\Scripts\python.exe -m pip install -r requirements.txt
if errorlevel 1 goto fail_pip
:run
.venv\Scripts\python.exe tools\launch.py
pause
exit /b 0
:fail_venv
echo Could not create the virtual environment
pause
exit /b 1
:fail_pip
echo Installing the packages failed - see the messages above
rmdir /s /q .venv
pause
exit /b 1
