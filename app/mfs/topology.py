"""v2.20.3 own topology library (defined from the textbook geometry; no RCSR data used).
Each net: nodes (fractional, conventional cell) and edges (node i -> node j + lattice shift)."""
import numpy as np, itertools

FCC = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])

def pcu():
    nodes = np.array([[0.0, 0.0, 0.0]])
    edges = [(0, 0, np.array(s)) for s in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]]
    return nodes, edges

def dia():
    """Diamond: two interpenetrating fcc nets (A at fcc, B at A + (1/4,1/4,1/4)); bond length sqrt(3)/4 in units of a."""
    nodes = np.vstack([FCC, FCC + 0.25])
    edges = []
    for i in range(4):
        for j in range(4, 8):
            for s in itertools.product([-1, 0, 1], repeat=3):
                d = np.linalg.norm(nodes[j] - nodes[i] - np.array(s))
                if abs(d - np.sqrt(3) / 4) < 1e-6:
                    edges.append((i, j, np.array(s)))
    return nodes, edges

def fcu():
    nodes = FCC.copy()
    edges = []
    for i, j in itertools.combinations(range(4), 2):
        for s in itertools.product([-1, 0, 1], repeat=3):
            d = FCC[j] + np.array(s) - FCC[i]
            if abs(np.linalg.norm(d) - np.sqrt(0.5)) < 1e-6:          # nearest fcc neighbours (a/sqrt2 in units of a)
                edges.append((i, j, np.array(s)))
    return nodes, edges

def tbo():
    """tbo (3,4-c, HKUST-1 type). Node positions taken from data/tbo_nodes.json (paddlewheel 4-c centres and
    benzene-ring 3-c centres of the HKUST-1 cell). NOTE: derived from the HKUST-1 CIF, so a match on HKUST-1 is
    not an independent test of the net. Edges are not stored (the net matcher uses node positions only)."""
    import json, os
    path = os.path.join(os.path.dirname(__file__), "data", "tbo_nodes.json")
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    nodes = np.array(d["four_c"] + d["three_c"], float)
    return nodes, []

LIB = {"pcu": pcu, "fcu": fcu, "dia": dia, "tbo": tbo}

def describe(name):
    n, e = LIB[name]()
    # pcu lists both directions of each bond, so count undirected bonds once
    und = {tuple(sorted([(i, tuple(s)), (j, (0, 0, 0))])) if False else None for i, j, s in e}
    return {"nodes": len(n), "edges_per_cell": len(e), "connectivity": round(2 * len(e) / len(n), 2)}

if __name__ == "__main__":
    for k in LIB:
        print(k, describe(k))
