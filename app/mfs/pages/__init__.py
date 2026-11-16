"""Page modules of Reticula (split out of the former 3,800-line app.py in v2.13).

Each module holds the code of one or a few technique pages, moved verbatim. The code still uses the app's helper
functions (show, dl, files_input, newfig …) and the per-run sidebar values (WIDTH, ASPECT, FONT …) as plain globals;
`bind(globals())` is called by app.py once per run, after the sidebar is built, and writes exactly those names into
each module (the module's INJECTED tuple). A name that does not exist raises at once - not on a rarely used page."""
from __future__ import annotations

import importlib

MODULES = ("spectra", "nmr", "sorption", "thermal", "defects", "kinetics", "elemental", "microscopy", "bioassays", "structure", "workflow")


def modules():
    return [importlib.import_module(f"{__name__}.{m}") for m in MODULES]


def bind(ns: dict) -> None:
    mods = modules()
    for m in mods:
        for name in m.INJECTED:
            if name in ns:
                val = ns[name]
            else:                                   # defined in another page module (e.g. node_picker used by the DSC page)
                owner = next((o for o in mods if o is not m and name in vars(o) and name not in o.INJECTED), None)
                if owner is None:
                    raise NameError(f"mfs.pages.{m.__name__.rsplit('.', 1)[-1]} needs '{name}' from app.py, which does not define it")
                val = vars(owner)[name]
            setattr(m, name, val)
