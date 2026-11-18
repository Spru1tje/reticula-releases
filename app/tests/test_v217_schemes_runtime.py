"""v2.17 (Schemes & Icons, phase 1): the vendored icon library loads and renders from inside Reticula."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mfs.reticon.runtime import Library

HERE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "mfs")
RUNTIME = os.path.join(HERE, "icons_runtime", "specs.json")
CANVAS = os.path.join(HERE, "components", "scheme_canvas")


def test_library_loads_searches_and_renders_an_icon():
    lib = Library.load(RUNTIME)
    hits = lib.search("flask")
    assert hits, "search for 'flask' should find the flask icons"
    svg = lib.svg("round_bottom_flask", shape="base", style="illustration")
    assert svg.lstrip().startswith("<") and "<svg" in svg


def test_canvas_component_and_variants_are_static_files():
    assert os.path.exists(os.path.join(CANVAS, "index.html"))
    assert os.path.exists(os.path.join(CANVAS, "catalog.json"))
    assert os.path.isdir(os.path.join(CANVAS, "variants"))
