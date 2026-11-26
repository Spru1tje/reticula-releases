"""v2.22.7: internal consistency of the repo reference MIL-53_Cr.cif (COD 1502688). No external data.
Checks bond lengths that the CIF itself implies. The C–O carboxylate check is marked xfail on purpose:
the repo CIF has no C–O band of carboxylate length (see the finding in the report), so the test documents the gap
and will flag (strict xfail) when the CIF is corrected."""
import os, sys, itertools, collections
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest

CIF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "references", "library", "MIL-53_Cr.cif")


def _pair_distances(cutoff=2.3):
    gemmi = pytest.importorskip("gemmi")
    small = gemmi.make_small_structure_from_block(gemmi.cif.read_file(CIF)[0])
    c = small.cell
    sites = small.get_all_unit_cell_sites()
    fr = np.array([[s.fract.x, s.fract.y, s.fract.z] for s in sites])
    el = [s.element.name for s in sites]
    out = collections.defaultdict(list)
    for i in range(len(sites)):
        for j in range(len(sites)):
            for t in itertools.product([-1, 0, 1], repeat=3):
                if i == j and t == (0, 0, 0):
                    continue
                d = c.orthogonalize(gemmi.Fractional(*(fr[j] + np.array(t) - fr[i]))).length()
                if d <= cutoff:
                    out[tuple(sorted([el[i], el[j]]))].append(d)
    return {k: np.array(v) for k, v in out.items()}


def test_cr_o_distances_are_in_the_octahedral_range():
    d = _pair_distances()
    cro = d[("Cr", "O")]
    assert cro.size > 0
    assert 1.8 <= cro.min() and cro.max() <= 2.05


def test_c_c_distances_are_aromatic():
    d = _pair_distances()
    cc = d[("C", "C")]
    assert cc.size > 0
    assert 1.30 <= cc.min() and cc.max() <= 1.55


@pytest.mark.xfail(strict=True, reason="repo MIL-53_Cr.cif: no carboxylate C–O band (1.22–1.32 Å); C–O lengths are "
                                       "1.174 Å or >= 1.45 Å. The reference needs correcting before O–C validation.")
def test_c_o_carboxylate_band_exists():
    d = _pair_distances()
    co = d[("C", "O")]
    assert np.any((co >= 1.22) & (co <= 1.32))
