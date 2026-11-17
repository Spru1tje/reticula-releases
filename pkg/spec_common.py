"""Shared pieces of the PyInstaller specs (Windows and macOS). Imported by pkg/mfs_win.spec and pkg/mfs_mac.spec."""
from __future__ import annotations
import ast
import importlib.util
import re
import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_data_files, collect_dynamic_libs, collect_submodules, copy_metadata

STDLIB = set(sys.stdlib_module_names)

# Python 3.12.0 miscompiles module-level comprehensions once PyInstaller rewrites the code objects: the frozen app then fails with
# "NameError: name 'obj' is not defined" in scipy.stats (15 of 40 self-tests fail). Fixed in Python 3.12.1.
if sys.version_info[:3] == (3, 12, 0):
    raise SystemExit("Do not build Reticula with Python 3.12.0 (the frozen scipy.stats is broken). Use Python 3.12.1 or newer, or 3.13.")


def app_imports(app_dir: Path) -> list[str]:
    """Every absolute module the app's own code imports (also inside functions). PyInstaller only analyses the launcher, and
    the app is shipped as data files, so anything the app imports lazily would otherwise be missing from the frozen build."""
    found: set[str] = set()
    for p in app_dir.rglob("*.py"):
        parts = set(p.relative_to(app_dir).parts)
        if parts & {"tests", "tools", "__pycache__"}:
            continue
        try:
            tree = ast.parse(p.read_text(encoding="utf-8"))
        except Exception:
            continue
        for n in ast.walk(tree):
            if isinstance(n, ast.Import):
                found.update(a.name for a in n.names)
            elif isinstance(n, ast.ImportFrom) and n.level == 0 and n.module:
                found.add(n.module)
                found.update(f"{n.module}.{a.name}" for a in n.names)
    out = []
    for m in sorted(found):
        top = m.split(".")[0]
        if top in ("mfs", "app", "tools"):
            continue
        try:
            if importlib.util.find_spec(m) is not None:
                out.append(m)
        except Exception:
            pass
    return out


METADATA = ["streamlit", "altair", "pandas", "numpy", "pyarrow", "protobuf", "click", "tornado", "starlette", "uvicorn", "packaging",
            "pillow", "matplotlib", "scipy", "pymatgen", "pymatgen-core", "monty", "spglib", "tifffile", "ncempy", "nmrglue", "openpyxl",
            "xlrd", "jsonschema", "narwhals", "pydeck", "watchdog", "requests", "websockets", "anyio", "h11", "jinja2", "uncertainties",
            "ruamel.yaml", "tqdm", "brukeropusreader", "pytesseract", "streamlit-image-coordinates", "streamlit-ketcher", "rdkit"]

ALWAYS_HIDDEN = ["h5py", "hdf5plugin", "streamlit.web.cli", "streamlit.runtime.scriptrunner.magic_funcs", "pymatgen.core",
                 "pymatgen.analysis.diffraction.xrd", "pymatgen.io.cif", "scipy.optimize", "scipy.signal", "scipy.sparse.linalg",
                 "scipy.spatial", "scipy.interpolate", "scipy.integrate", "scipy.spatial._qhull", "scipy.spatial._ckdtree", "scipy.stats",
                 "scipy.ndimage", "scipy.signal._signaltools", "matplotlib.backends.backend_pdf", "matplotlib.backends.backend_svg",
                 "matplotlib.backends.backend_agg", "mpl_toolkits.axes_grid1.inset_locator", "PIL.Image", "PIL.TiffImagePlugin", "contourpy",
                 "sqlite3", "_sqlite3", "ssl", "tempfile", "tkinter", "_tkinter", "streamlit.components.v1"]

EXCLUDES = ["pytest", "IPython", "jupyter", "notebook", "PyQt5", "PyQt6", "PySide2", "PySide6", "wx", "sympy", "plotly", "lxml",
            "pyarrow.flight", "pyarrow.parquet", "pyarrow.dataset", "pyarrow.substrait", "pyarrow.tests", "scipy.tests", "numpy.tests",
            "pandas.tests", "matplotlib.tests", "networkx.drawing", "sphinx", "pymatgen.ext", "pymatgen.io.vasp.outputs"]

DROP = re.compile(r"(pyarrow[\\/](include|tests|.*(flight|s3|substrait|dataset|acero|parquet|orc|gcsfs|azurefs|hdfs|cuda)[^\\/]*\.(dll|pyd|so|dylib))"
                  r"|pyarrow[\\/].*\.(lib|pyx|pxd|pxi|h|hpp|cc)$|[\\/]tests?[\\/]|mpl-data[\\/]sample_data|\.(lib|pyi|c|h)$|[\\/]lxml[\\/])", re.I)


def keep(toc):
    out = [e for e in toc if not DROP.search(e[0].replace("/", "\\")) and not DROP.search(e[0])]
    print("size filter dropped", len(toc) - len(out), "entries")
    return out


APP_PARTS = ("app.py", "mfs", "references", "example_data", "docs", "assets", ".streamlit", "README.md", "CHANGELOG.md")


def app_datas(app_dir: Path):
    """The program files that go into the build as data: only APP_PARTS, without caches. (Bundling the whole app/ folder would
    also take tests, tools and - on a developer machine - .venv, exports and private test data.)"""
    out = []
    for part in APP_PARTS:
        src = app_dir / part
        if src.is_file():
            out.append((str(src), "app"))
        elif src.is_dir():
            for f in sorted(src.rglob("*")):
                if f.is_file() and "__pycache__" not in f.parts and f.suffix not in (".pyc", ".pyo"):
                    out.append((str(f), str(Path("app") / f.relative_to(app_dir).parent)))
        else:
            raise FileNotFoundError(f"app part missing: {src}")
    return out


def collect(app_dir: Path):
    """datas, binaries, hiddenimports shared by both platforms."""
    datas, binaries = [], []
    for pkg in ("streamlit_image_coordinates", "streamlit", "matplotlib"):
        datas += collect_data_files(pkg)
    datas += collect_data_files("pymatgen", includes=["**/*.json", "**/*.yaml", "**/*.json.gz", "**/*.csv"])
    hidden = list(ALWAYS_HIDDEN)
    for pkg in ("streamlit_ketcher", "rdkit"):          # v2.13 NMR page: structure editor front-end + chemistry toolkit
        try:
            d, b, h = collect_all(pkg)
            datas += d; binaries += b; hidden += h
        except Exception as e:
            print("collect_all failed for", pkg, e)
    for pkg in METADATA:
        try:
            datas += copy_metadata(pkg)
        except Exception as e:
            print("no metadata for", pkg, e)
    hidden += app_imports(app_dir)
    hidden += collect_submodules("streamlit")
    hidden += [m for m in collect_submodules("PIL") if "test" not in m]
    hidden += collect_submodules("uvicorn") + collect_submodules("starlette")
    return datas, binaries, sorted(set(hidden))
