"""v2.19.0 Reticon: the ten new material-science icons (five shape variants each) have known element counts and render to PNG."""
from pathlib import Path

import resvg_py

from mfs.reticon import materials_2190 as M
from mfs.reticon.runtime import Library


def _count(icon, tag, role=None):
    return sum(1 for o in icon.ops if o.tag == tag and (role is None or o.role == role))


def test_simple_cubic_base_has_eight_corner_atoms():
    ic = M.simple_cubic("base")
    assert _count(ic, "circle", "a") == 8


def test_face_centred_base_adds_five_face_atoms_to_the_eight_corners():
    ic = M.face_centred("base")
    assert _count(ic, "circle", "a") == 8
    assert _count(ic, "circle", "s") == 5


def test_mesoporous_base_has_sixteen_pores_cut_from_the_body():
    ic = M.mesoporous("base")
    assert _count(ic, "circle", "x") == 16


def test_cylindrical_pore_base_has_two_channels():
    ic = M.cylindrical_pore("base")
    assert _count(ic, "circle", "x") == 2


def test_nanowire_base_has_three_wires():
    ic = M.nanowire("base")
    assert _count(ic, "rect") == 3


def test_every_new_icon_has_five_variants_and_renders_to_png(tmp_path):
    lib = Library.load(Path(M.__file__).resolve().parents[2] / "mfs" / "icons_runtime" / "specs.json").load_variants()
    for name in M.BUILDERS:
        assert set(lib.icons[name]["shapes"]) == set(M.VARIANTS), name
        for v in M.VARIANTS:
            svg = lib.svg(name, shape=v)
            png = bytes(resvg_py.svg_to_bytes(svg_string=svg, width=96, height=96))
            assert png[:8] == b"\x89PNG\r\n\x1a\n", (name, v)
            assert len(png) > 200, (name, v)


def _new_names():
    from pathlib import Path as _P
    terms = [l.strip() for l in open(_P(__file__).resolve().parents[2] / "docs" / "icons" / "terms-new.txt", encoding="utf-8") if l.strip()]
    return list(M.BUILDERS) + terms[:290]


def test_three_hundred_new_icons_with_five_variants_each():
    names = _new_names()
    assert len(names) == 300 and len(set(names)) == 300
    lib = Library.load(Path(M.__file__).resolve().parents[2] / "mfs" / "icons_runtime" / "specs.json").load_variants()
    for n in names:
        assert n in lib.icons, n
        assert lib.icons[n]["category"] == "materials", n
        assert set(lib.icons[n]["shapes"]) == set(M.VARIANTS), n


def test_all_new_variants_render_to_png():
    lib = Library.load(Path(M.__file__).resolve().parents[2] / "mfs" / "icons_runtime" / "specs.json").load_variants()
    for n in _new_names():
        for v in M.VARIANTS:
            png = bytes(resvg_py.svg_to_bytes(svg_string=lib.svg(n, shape=v), width=48, height=48))
            assert png[:8] == b"\x89PNG\r\n\x1a\n", (n, v)


def test_family_rules_for_the_three_corrected_terms():
    from mfs.reticon import families_2190 as F
    assert F.family_for("hysteresis_loop") == "curve"
    assert F.family_for("downconversion") == "curve"
    assert F.family_for("cryogel") == "foam"
