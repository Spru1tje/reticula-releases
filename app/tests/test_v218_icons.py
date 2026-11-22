"""v2.18: the Reticon library shipped with Reticula holds the 500 new Reticula icons, each with five variants (base + 5 = 6 shapes)."""
import json
import os

HERE = os.path.dirname(__file__)
RT = os.path.join(HERE, "..", "mfs", "icons_runtime")

# one name per family, so a broken family shows up here
SAMPLE = ["r22_vessel_00_empty_plain", "r22_chart_bars_00", "r22_rotor_04_14_tube", "r22_cuvette_10_50_screw", "xrd_pattern_1", "epr_spectrum_1", "tapered_tube_gel", "florence_flask_crystals", "glovebox_screen",
          "tube_furnace_knob", "nanocube_defect", "mesoporous_sphere_hollow", "process_straight_100", "process_tee_140"]


def _load(name):
    with open(os.path.join(RT, name), encoding="utf-8") as f:
        return json.load(f)


def test_library_has_the_new_icons():
    specs = _load("specs.json")["icons"]
    assert len(specs) >= 5000          # v2.20: about 5000 icons in total (built-in + v2.18 + v2.20 families)
    for n in SAMPLE:
        assert n in specs, n


def test_every_sampled_new_icon_has_five_variants():
    variants = _load("specs_variants.json")
    for n in SAMPLE:
        assert len(variants.get(n, {})) == 5, (n, sorted(variants.get(n, {})))


def test_new_icons_are_vector_only():
    specs = _load("specs.json")["icons"]
    allowed = {"path", "rect", "circle", "ellipse", "line", "polygon", "polyline"}
    for n in SAMPLE:
        ops = specs[n]["shapes"]["base"]["ops"]
        tags = {o["tag"] for o in ops}
        assert tags <= allowed, (n, tags - allowed)     # no images, no masks, no text
