# PyInstaller spec - Reticula for Linux (one-folder build). Run from the repository root:
#   pyinstaller --noconfirm --clean pkg/mfs_linux.spec      -> dist/Reticula/Reticula
# Tesseract is not bundled on Linux: install it with the system package manager (e.g. apt install tesseract-ocr).
import re
import sys
from pathlib import Path
ROOT = Path(SPECPATH).parent
sys.path.insert(0, str(Path(SPECPATH)))
import spec_common as SC

APP = ROOT / "app"
VERSION = re.search(r'__version__\s*=\s*"([^"]+)"', (APP / "mfs" / "__init__.py").read_text(encoding="utf-8")).group(1)

datas = SC.app_datas(APP)
d2, binaries, hidden = SC.collect(APP)
datas += d2

a = Analysis([str(ROOT / "mfs_launcher.py")], pathex=[str(APP)], datas=datas, binaries=binaries, hiddenimports=hidden,
             excludes=SC.EXCLUDES, noarchive=False)
a.binaries = SC.keep(a.binaries)
a.datas = SC.keep(a.datas)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="Reticula", console=False, upx=False)
coll = COLLECT(exe, a.binaries, a.datas, name="Reticula", upx=False)
