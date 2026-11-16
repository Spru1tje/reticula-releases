# PyInstaller spec - Reticula for Windows (one-folder, windowed). Run from the repository root:
#   pyinstaller --noconfirm --clean pkg/mfs_win.spec        -> dist/Reticula/Reticula.exe
import sys
from pathlib import Path
ROOT = Path(SPECPATH).parent
sys.path.insert(0, str(Path(SPECPATH)))
import spec_common as SC

APP = ROOT / "app"
datas = [(str(APP), "app"), (str(ROOT / "icon" / "reticula.ico"), ".")]
tess = ROOT / "tesseract"                                  # Tesseract OCR for scale-bar labels (copied there by the build workflow)
if tess.is_dir():
    datas.append((str(tess), "tesseract"))
d2, binaries, hidden = SC.collect(APP)
datas += d2
# the launcher looks for the icon next to itself
a = Analysis([str(ROOT / "mfs_launcher.py")], pathex=[str(APP)], datas=datas, binaries=binaries, hiddenimports=hidden,
             excludes=SC.EXCLUDES, noarchive=False)
a.binaries = SC.keep(a.binaries)
a.datas = SC.keep(a.datas)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="Reticula", icon=str(ROOT / "icon" / "reticula.ico"), console=False, upx=False)
coll = COLLECT(exe, a.binaries, a.datas, name="Reticula", upx=False)
