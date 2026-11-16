"""v2.13: the page code moved out of app.py into mfs/pages/*.py. These checks need no browser:
every name a page module takes from the app is listed in INJECTED, and app.py really provides it."""
import ast
import builtins
import os
import symtable
import sys
from unittest.mock import MagicMock

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
_STUBBED = []
_BEFORE = set(sys.modules)
try:
    import streamlit  # noqa: F401
except ImportError:                                # the structural checks below do not need a real Streamlit
    for _m in ("streamlit", "streamlit.components", "streamlit.components.v1"):
        if _m not in sys.modules:
            sys.modules[_m] = MagicMock(); _STUBBED.append(_m)
from mfs import pages as PG  # noqa: E402
_MODS = PG.modules()                               # import every page module now ...
if _STUBBED:                                       # ... then leave no trace of the stub: other test files must see "no streamlit" as before
    for _m in list(sys.modules):
        if _m in _STUBBED or (_m not in _BEFORE and (_m == "mfs" or _m.startswith("mfs."))):
            sys.modules.pop(_m, None)
            _parent, _, _child = _m.rpartition(".")           # `from mfs import x` reads the attribute of the package object
            if _parent in sys.modules and hasattr(sys.modules[_parent], _child):
                delattr(sys.modules[_parent], _child)


def _app_names():
    tree = ast.parse(open(os.path.join(ROOT, "app.py"), encoding="utf-8").read())
    names = set()
    for n in tree.body:
        if isinstance(n, (ast.FunctionDef, ast.ClassDef)):
            names.add(n.name)
        elif isinstance(n, (ast.Assign, ast.AnnAssign)):
            for t in (n.targets if isinstance(n, ast.Assign) else [n.target]):
                for e in (t.elts if isinstance(t, ast.Tuple) else [t]):
                    if isinstance(e, ast.Name):
                        names.add(e.id)
        elif isinstance(n, (ast.Import, ast.ImportFrom)):
            names.update((a.asname or a.name).split(".")[0] for a in n.names)
    return names


def _undefined_globals(path):
    """Names that functions/classes in the module read as globals but the module itself never defines or imports."""
    src = open(path, encoding="utf-8").read()
    top = symtable.symtable(src, path, "exec")
    defined = {s.get_name() for s in top.get_symbols() if s.is_assigned() or s.is_imported() or s.is_namespace()}
    out = set()

    def walk(t):
        for s in t.get_symbols():
            if t.get_type() != "module" and s.is_referenced() and s.is_global() and not s.is_assigned():
                out.add(s.get_name())
        for c in t.get_children():
            walk(c)
    walk(top)
    return out - defined - set(dir(builtins))


@pytest.mark.parametrize("name", PG.MODULES)
def test_injected_list_is_complete_and_available(name):
    mod = next(m for m in _MODS if m.__name__.endswith("." + name))
    path = os.path.join(ROOT, "mfs", "pages", name + ".py")
    missing = _undefined_globals(path) - set(mod.INJECTED)
    assert not missing, f"{name}: names used but neither defined, imported nor in INJECTED: {sorted(missing)}"
    app = _app_names(); others = [o for o in _MODS if o is not mod]
    for n in mod.INJECTED:
        assert n in app or any(n in vars(o) and n not in o.INJECTED for o in others), f"{name}: '{n}' is provided by nobody"


def test_bind_fails_loudly_on_a_missing_name():
    ns = {n: object() for n in _app_names()}
    _orig = PG.modules
    PG.modules = lambda: _MODS                     # bind() must not re-import (and so must not need Streamlit) here
    try:
        PG.bind(ns)                                 # complete namespace works
        ns.pop("show")
        with pytest.raises(NameError):
            PG.bind(ns)
    finally:
        PG.modules = _orig


def test_every_page_of_the_dispatch_table_is_importable():
    src = open(os.path.join(ROOT, "app.py"), encoding="utf-8").read()
    import re
    fns = set(re.findall(r'"[^"]+": (page_\w+)', src))
    assert len(fns) >= 20
    have = set()
    for m in _MODS:
        have |= {k for k in vars(m) if k.startswith("page_")}
    have |= {n for n in _app_names() if n.startswith("page_")}
    assert fns <= have, fns - have
