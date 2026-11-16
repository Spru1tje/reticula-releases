# PyInstaller spec - Reticula 2.11.0 macOS (one-folder .app bundle)
# Run from the repo root:  pyinstaller --noconfirm --clean pkg/mfs_mac.spec
from PyInstaller.utils.hooks import collect_data_files, copy_metadata, collect_submodules
from pathlib import Path

ROOT = Path(SPECPATH).parent  # repo root (this file lives in pkg/)
APP = ROOT / "app"
ICON = ROOT / "reticula.icns"
TESS = ROOT / "tesseract"

datas = [(str(APP), "app")]
if ICON.is_file():
    datas.append((str(ICON), "."))
if TESS.is_dir():
    datas.append((str(TESS), "tesseract"))
datas += collect_data_files("streamlit_image_coordinates")
datas += collect_data_files("streamlit")
datas += collect_data_files("pymatgen", includes=["**/*.json", "**/*.yaml", "**/*.json.gz", "**/*.csv"])
datas += collect_data_files("matplotlib")
for pkg in ["streamlit", "altair", "pandas", "numpy", "pyarrow", "protobuf", "click", "tornado", "starlette", "uvicorn",
            "packaging", "pillow", "matplotlib", "scipy", "pymatgen", "pymatgen-core", "monty", "spglib", "tifffile",
            "ncempy", "nmrglue", "openpyxl", "xlrd", "jsonschema", "narwhals", "pydeck", "watchdog", "requests",
            "websockets", "anyio", "h11", "jinja2", "uncertainties", "ruamel.yaml", "tqdm", "brukeropusreader",
            "pytesseract", "streamlit-image-coordinates"]:
    try:
        datas += copy_metadata(pkg)
    except Exception as e:
        print("no metadata for", pkg, e)

hidden = ["h5py", "hdf5plugin", "streamlit.web.cli", "streamlit.runtime.scriptrunner.magic_funcs",
          "tifffile", "nmrglue", "brukeropusreader", "ncempy.io", "xlrd", "openpyxl",
          "pymatgen.core", "pymatgen.analysis.diffraction.xrd", "pymatgen.io.cif",
          "scipy.optimize", "scipy.signal", "scipy.sparse.linalg", "scipy.spatial", "scipy.interpolate", "scipy.integrate",
          "matplotlib.backends.backend_pdf", "matplotlib.backends.backend_svg", "matplotlib.backends.backend_agg",
          "PIL.Image", "PIL.TiffImagePlugin", "pytesseract", "streamlit_image_coordinates",
          "mfs.semcal", "mfs.dsc", "mfs.psd", "mfs.kinetics", "mfs.xrd", "mfs.samplestore", "mfs.si_export",
          "xml.sax.saxutils", "zipfile", "mfs.spectra", "mfs.xps", "mfs.tga", "mpl_toolkits.axes_grid1.inset_locator",
          "mfs.particles", "mfs.edx", "mfs.ratios", "mfs.templates", "mfs.batchproc", "mfs.si_bundle", "mfs.sampledb",
          "sqlite3", "_sqlite3", "mfs.cif", "mfs.cifplot", "mfs.cifsi", "mfs.cifpxrd", "mfs.poretools", "mfs.hirshfeld",
          "mfs.pagestate", "mfs.uvvis", "mfs.uvvis_page", "mfs.rigaku", "mfs.particle_picker", "mfs.shapeprior",
          "scipy.signal._signaltools", "contourpy", "scipy.signal", "scipy.optimize", "streamlit.components.v1",
          "scipy.spatial", "scipy.spatial._qhull", "scipy.spatial._ckdtree", "scipy.stats", "scipy.ndimage",
          "mfs.nodes", "mfs.chem", "mfs.defects", "mfs.project", "mfs.globalstyle", "mfs.update", "mfs.samplelist",
          "mfs.style", "urllib.request", "http.client", "ssl", "tempfile", "tkinter", "_tkinter"]
hidden += collect_submodules("streamlit")
hidden += [m for m in collect_submodules("PIL") if "test" not in m]
hidden += collect_submodules("uvicorn") + collect_submodules("starlette")

a = Analysis(
    [str(ROOT / "mfs_launcher.py")],
    pathex=[str(APP)],
    datas=datas,
    hiddenimports=hidden,
    excludes=["pytest", "IPython", "jupyter", "notebook", "PyQt5", "PyQt6", "PySide2", "PySide6", "wx",
              "sympy", "plotly", "lxml", "pyarrow.flight", "pyarrow.parquet", "pyarrow.dataset", "pyarrow.substrait",
              "pyarrow.tests", "scipy.tests", "numpy.tests", "pandas.tests", "matplotlib.tests", "networkx.drawing",
              "sphinx", "pymatgen.ext", "pymatgen.io.vasp.outputs"],
    noarchive=False,
)

import re as _re
_DROP = _re.compile(
    r"(pyarrow[\\/](include|tests|.*(flight|s3|substrait|dataset|acero|parquet|orc|gcsfs|azurefs|hdfs|cuda)[^\\/]*\.(so|dylib|pyd))"
    r"|pyarrow[\\/].*\.(lib|pyx|pxd|pxi|h|hpp|cc)$|[\\/]tests?[\\/]|mpl-data[\\/]sample_data"
    r"|\.(lib|pyi|c|h)$|[\\/]lxml[\\/])",
    _re.I,
)

def _keep(toc):
    out = [e for e in toc if not _DROP.search(e[0].replace("\\", "/"))]
    print("size filter dropped", len(toc) - len(out), "entries")
    return out

a.binaries = _keep(a.binaries)
a.datas = _keep(a.datas)

pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, [],
    exclude_binaries=True,
    name="Reticula",
    icon=str(ICON) if ICON.is_file() else None,
    console=False,
    upx=False,
    disable_windowed_traceback=False,
)
coll = COLLECT(exe, a.binaries, a.datas, name="Reticula", upx=False)
app = BUNDLE(
    coll,
    name="Reticula.app",
    icon=str(ICON) if ICON.is_file() else None,
    bundle_identifier="be.ugent.reticula",
    info_plist={
        "CFBundleName": "Reticula",
        "CFBundleDisplayName": "Reticula",
        "CFBundleShortVersionString": "2.11.0",
        "CFBundleVersion": "2.11.0",
        "CFBundleGetInfoString": "Reticula 2.11.0 — Built by Stijn Germonpré",
        "NSHumanReadableCopyright": "Built by Stijn Germonpré",
        "NSHighResolutionCapable": True,
        "LSMinimumSystemVersion": "13.0",
        "CFBundleDocumentTypes": [
            {
                "CFBundleTypeName": "Reticula project",
                "CFBundleTypeRole": "Editor",
                "CFBundleTypeExtensions": ["reticula"],
                "LSHandlerRank": "Owner",
                "CFBundleTypeDescription": "Reticula project file",
            }
        ],
        "UTExportedTypeDeclarations": [
            {
                "UTTypeIdentifier": "be.ugent.reticula.project",
                "UTTypeDescription": "Reticula project",
                "UTTypeConformsTo": ["public.data", "public.archive"],
                "UTTypeTagSpecification": {
                    "public.filename-extension": ["reticula"],
                    "public.mime-type": "application/x-reticula-project",
                },
            }
        ],
    },
)
