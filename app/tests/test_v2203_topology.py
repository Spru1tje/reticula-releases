"""v2.20.3: own topology library (pcu, fcu, dia) built from the textbook geometry; checked on UiO-66 (fcu)."""
import os
import numpy as np
import pytest

from mfs import topology as T


def test_library_counts():
    assert T.describe("pcu") == {"nodes": 1, "edges_per_cell": 3, "connectivity": 6.0}
    assert T.describe("fcu") == {"nodes": 4, "edges_per_cell": 24, "connectivity": 12.0}
    assert T.describe("dia") == {"nodes": 8, "edges_per_cell": 16, "connectivity": 4.0}


CIF = os.path.join(os.path.dirname(__file__), "..", "references", "UiO-66_Zr_COD4512074.cif")


@pytest.mark.skipif(not os.path.exists(CIF), reason="reference CIF missing")
def test_uio66_is_fcu():
    pymatgen = pytest.importorskip("pymatgen.core")
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    s = pymatgen.Structure.from_file(CIF)
    L = s.lattice.matrix
    def pd(a, b):
        d = a - b; d -= np.round(d); return np.linalg.norm(d @ L, axis=-1)
    def centroids(pts, cut):
        D = np.array([[pd(pts[i], pts[j]) for j in range(len(pts))] for i in range(len(pts))])
        lab = fcluster(linkage(squareform(D, checks=False), "single"), t=cut, criterion="distance")
        out = []
        for k in np.unique(lab):
            q = pts[lab == k]; ref = q[0]
            out.append((np.array([ref + ((p - ref) - np.round(p - ref)) for p in q]).mean(axis=0)) % 1)
        return np.array(out)
    zr_c = centroids(np.array([x.frac_coords for x in s if x.species_string.startswith("Zr")]), 4.2)
    lk_c = centroids(np.array([x.frac_coords for x in s if x.species_string.startswith("C")]), 1.7)
    nodes, edges = T.fcu()
    assert len(zr_c) == 4 and len(lk_c) == 24
    # origin: the cluster centre that fits the fcc nodes best (the cell origin is a choice)
    best = min(((np.sqrt(np.mean([min(pd(c - sh, n) ** 2 for n in nodes) for c in zr_c])), sh) for sh in zr_c), key=lambda t: t[0])
    sh = best[1]
    assert best[0] < 0.05, best[0]
    mids = np.array([((nodes[i] + nodes[j] + sf) / 2 + sh) % 1 for i, j, sf in edges])
    d = [min(pd(c, m) for m in mids) for c in lk_c]
    assert float(np.max(d)) < 0.05, float(np.max(d))
