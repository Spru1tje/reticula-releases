"""Snapshot every defect method's numbers + formula for all bundled Zr/Hf inputs (+ variants), to prove the generalisation keeps them identical."""
import sys, json, glob, os, copy, math
import pathlib; ROOT = str(pathlib.Path(__file__).resolve().parents[2]); sys.path.insert(0, ROOT)
from mfs import defects as D, report as REP
def snap(inp):
    res = D.run_all(inp); out = {}
    for m in res["methods"]:
        out[m.key] = {k: (None if v is None else (round(v, 9) if isinstance(v, float) else v)) for k, v in
                      dict(L=m.L, x=m.x, m=m.m, Cl=m.Cl, OS=m.OS, M=m.M_formula, f=m.formula, st=m.status, msg=m.message,
                           steps="|".join(m.steps)).items()}
        out[m.key]["extra"] = {k: (round(v, 9) if isinstance(v, float) else v) for k, v in m.extra.items() if not isinstance(v, (list, dict))}
    out["cons"] = {k: (round(v, 9) if isinstance(v, float) else v) for k, v in res["consensus"].items()}
    return out
def cases():
  cases = {}
  for p in sorted(glob.glob(os.path.join(ROOT, "example_data", "Defects", "defect_input_*.json"))):
    d = json.load(open(p, encoding="utf-8")); inp = REP.input_from_dict(d)
    if inp is None:
        raise SystemExit("no dict_to_input")
    for metal in ("Zr", "Hf", "mixed"):
        i2 = copy.deepcopy(inp); i2.metal = metal
        if metal == "mixed": i2.f_hf = 0.3
        cases[os.path.basename(p) + "|" + metal] = snap(i2)
        for fw in ("UiO-66", "NU-1000", "MOF-808"):
            i3 = copy.deepcopy(i2); i3.framework = fw; i3.CN = None
            cases[os.path.basename(p) + "|" + metal + "|" + fw] = snap(i3)
  return json.loads(json.dumps(cases, sort_keys=True, default=str))
if __name__ == "__main__":
    json.dump(cases(), open(sys.argv[1], "w", encoding="utf-8"), indent=0, sort_keys=True, default=str)
