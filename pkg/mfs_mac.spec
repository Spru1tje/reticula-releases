# PyInstaller spec - Reticula for macOS (one-folder .app bundle). Run from the repository root:
#   pyinstaller --noconfirm --clean pkg/mfs_mac.spec        -> dist/Reticula.app
import re
import sys
from pathlib import Path
ROOT = Path(SPECPATH).parent
sys.path.insert(0, str(Path(SPECPATH)))
import spec_common as SC

APP = ROOT / "app"
VERSION = re.search(r'__version__\s*=\s*"([^"]+)"', (APP / "mfs" / "__init__.py").read_text(encoding="utf-8")).group(1)
ICON = ROOT / "reticula.icns"                              # made by pkg/make_icns.sh
TESS = ROOT / "tesseract"                                  # made by pkg/bundle_tesseract.sh

datas = SC.app_datas(APP)
if ICON.is_file():
    datas.append((str(ICON), "."))
if TESS.is_dir():
    datas.append((str(TESS), "tesseract"))
d2, binaries, hidden = SC.collect(APP)
datas += d2

a = Analysis([str(ROOT / "mfs_launcher.py")], pathex=[str(APP)], datas=datas, binaries=binaries, hiddenimports=hidden,
             excludes=SC.EXCLUDES, noarchive=False)
a.binaries = SC.keep(a.binaries)
a.datas = SC.keep(a.datas)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="Reticula", icon=str(ICON) if ICON.is_file() else None,
          console=False, upx=False, disable_windowed_traceback=False)
coll = COLLECT(exe, a.binaries, a.datas, name="Reticula", upx=False)
app = BUNDLE(
    coll, name="Reticula.app", icon=str(ICON) if ICON.is_file() else None, bundle_identifier="be.ugent.reticula",
    info_plist={
        "CFBundleName": "Reticula", "CFBundleDisplayName": "Reticula",
        "CFBundleShortVersionString": VERSION, "CFBundleVersion": VERSION,
        "CFBundleGetInfoString": f"Reticula {VERSION}", "NSHumanReadableCopyright": "Built by Stijn Germonpré",
        "NSHighResolutionCapable": True, "LSMinimumSystemVersion": "13.0",
        "CFBundleDocumentTypes": [{"CFBundleTypeName": "Reticula project", "CFBundleTypeRole": "Editor", "CFBundleTypeExtensions": ["reticula"],
                                   "LSHandlerRank": "Owner", "CFBundleTypeDescription": "Reticula project file"}],
        "UTExportedTypeDeclarations": [{"UTTypeIdentifier": "be.ugent.reticula.project", "UTTypeDescription": "Reticula project",
                                        "UTTypeConformsTo": ["public.data", "public.archive"],
                                        "UTTypeTagSpecification": {"public.filename-extension": ["reticula"],
                                                                   "public.mime-type": "application/x-reticula-project"}}],
    },
)
