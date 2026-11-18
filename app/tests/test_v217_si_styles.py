"""v2.17: SI export in literature style (journal wording templates)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mfs import summary as SUM


def test_every_style_changes_heading_and_software_sentence():
    texts = {k: SUM.methods_text(["TGA", "PXRD"], "2.17.0", k) for k in SUM.SI_STYLES}
    assert len(set(texts.values())) == len(SUM.SI_STYLES)
    for k, t in texts.items():
        assert t.startswith("### " + SUM.SI_STYLES[k]["heading"])
        assert "Figuryte 2.17.0" in t and "**TGA.**" in t and "**PXRD.**" in t
        assert "  " not in t


def test_default_style_and_unknown_style_fall_back_to_acs():
    assert SUM.methods_text(["TGA"], "2.17.0") == SUM.methods_text(["TGA"], "2.17.0", "ACS")
    assert SUM.methods_text(["TGA"], "2.17.0", "nope") == SUM.methods_text(["TGA"], "2.17.0", "ACS")


def test_no_version_gives_no_double_space_and_no_pages_gives_empty():
    t = SUM.methods_text(["TGA"], "", "RSC")
    assert "Figuryte (open" in t and "  " not in t
    assert SUM.methods_text(["Unknown page"], "2.17.0", "ACS") == ""
