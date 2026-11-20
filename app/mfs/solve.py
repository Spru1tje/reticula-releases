r"""v2.18.4: crystal-structure solution from a laboratory powder pattern of an unknown compound - the decision tree of the field, automated.

The route (each step reports what it decided and why; any step can be given by the user instead):
  0. DATA CHECK        resolution, number of peaks, range. A silicon (or other) reference pattern gives the zero shift and the instrument profile.
  1. UNIT CELL         peak picking -> indexing (several candidate cells) -> a Le Bail fit of every candidate on the real pattern; the cell with the lowest
                       Rwp wins, a smaller cell wins over a supercell that fits no better.
  2. SPACE GROUP       the reflections that a space group forbids must be absent: isolated forbidden reflections are tested against the Le Bail intensities
                       (systematic absences); the groups that agree are ranked, centrosymmetric first.
  3. STRUCTURE SOLUTION charge flipping with the space group imposed, many random starts in parallel; the density peaks give the atoms.
  4. ATOMS             peaks -> symmetry orbits -> elements (from the composition you give) by the R factor of the amplitudes; spurious peaks are removed.
  5. REFINEMENT        Rietveld refinement of the model with soft bond restraints (the existing refinement module) - the structure you get is a starting
                       model for a careful refinement, not a finished publication structure.
Outputs: a CIF, the Le Bail / Rietveld pattern files, a report of every decision. Honest limits are listed in the result (`warnings`)."""
from __future__ import annotations

import math
import os
import time

import numpy as np

from mfs import chargeflip as CF
from mfs import rietveld as RV
from mfs import spacegroups as SG
from mfs import xtal

_ASSIGNABLE_MIN_Z = 3                                   # hydrogen and helium are not placed


class Steps(list):
    def add(self, key, title, status, text, data=None):
        self.append({"key": key, "title": title, "status": status, "text": text, "data": data or {}})


# ------------------------------------------------------------------------------------------------------------------------------ 0. data
def data_check(x, y, wavelength):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    dmin = wavelength / (2 * math.sin(math.radians(x.max() / 2)))
    step = float(np.median(np.diff(x)))
    from mfs import phaseid as PID
    pk = PID.exp_peaks(x, y, 0.015)[0]
    notes = []
    if dmin > 1.6:
        notes.append(f"The pattern stops at d = {dmin:.2f} Å (2θ = {x.max():.0f}°). For a structure solution data to d ≈ 1.2 Å (Cu Kα: 2θ ≈ 80°) or better are needed; the result will only show the heaviest atoms.")
    if len(pk) < 12:
        notes.append(f"Only {len(pk)} peaks were found. Too few for a reliable cell and structure: measure longer / to a higher angle.")
    if step > 0.03:
        notes.append(f"The step ({step:.3f}°) is coarse; 0.01–0.02° is better for overlapping peaks.")
    return {"dmin": dmin, "step": step, "n_peaks": int(len(pk)), "peaks": np.asarray(pk, float), "notes": notes, "tt_max": float(x.max()), "tt_min": float(x.min())}


# ------------------------------------------------------------------------------------------------------------------------------ 1. cell
def lebail_fit(x, y, cell, radiation="Cu Kα1", bg_order=6, instrument=None, offset="displacement", glass_bg=False, progress=None, fix_cell=False):
    ph = RV.Phase.from_cell(tuple(cell), "unknown phase")
    R = RV.Refinement(np.asarray(x, float), np.asarray(y, float), [ph], radiation, bg_order, "lebail", 0.0, instrument=instrument, offset=offset,
                      glass_bg=glass_bg, progress=progress)
    if fix_cell:
        R.fix_cell = True
    res = R.run()
    return R, res


def physical_cores():
    """Physical cores (hyper-threads add little to numpy / scipy work): half of the logical count on hyper-threaded CPUs, at least 1."""
    n = os.cpu_count() or 2
    return max(1, n // 2) if n >= 4 else n


def _lebail_job(a):
    """One candidate cell: Le Bail fit; module level so that a process pool can run it. A bad candidate returns (None, error) and does not stop the others."""
    x, y, cell, radiation, instrument, offset = a
    try:
        return lebail_fit(x, y, cell, radiation, 6, instrument, offset)
    except Exception as e:
        return None, {"error": str(e)}


def pick_cell(x, y, wavelength, peaks, radiation="Cu Kα1", systems=("cubic", "tetragonal", "hexagonal", "orthorhombic"), n_candidates=4, instrument=None,
              threads=None, offset="displacement", progress=None, steps=None, cell_given=None):
    """Indexing + Le Bail of the candidates; returns {cell, system, R, res, candidates:[...]} for the winner."""
    from mfs import indexing as IX
    prog = progress or (lambda f, t: None)
    if cell_given is not None:
        cands = [{"system": RV.detect_system(tuple(cell_given)), "cell": tuple(cell_given), "M20": None, "volume": RV.cell_volume(tuple(cell_given)), "n_unindexed": 0, "quality": "given",
                  "centring": "?"}]
    else:
        prog(0.0, "indexing")
        cands = IX.index(peaks, wavelength, systems=systems, max_solutions=max(n_candidates + 2, 6), progress=lambda f, t: prog(0.6 * f, "indexing – " + t))
        cands = cands[: max(n_candidates, 1)]
    # large cells: each Le Bail fit costs far more (many reflections); only the best-indexed candidates are fitted
    if cell_given is None and cands and RV.cell_volume(tuple(cands[0]["cell"])) > 2000.0:
        cands = cands[:2]
    if not cands:
        if steps is not None:
            steps.add("cell", "Unit cell", "failed", "No unit cell could be found from the peak positions. Possible reasons: impurity peaks, a cell of low symmetry (monoclinic: tick the option; triclinic is not searched), "
                                                       "a zero shift of the pattern (use a silicon reference), or too few peaks.")
        return None
    # measured on the development PC (MOF-size cells): sequential fits were faster than a process pool (worker start-up re-imports the modules)
    args = [(x, y, c["cell"], radiation, instrument, offset) for c in cands]
    out = [_lebail_job(a) for a in args]
    rows = []
    for c, (R, res) in zip(cands, out):
        if R is None:
            continue
        cell = res["phases"][0]["cell"]
        rows.append({**c, "cell_fit": tuple(float(v) for v in cell), "Rwp": float(res["stats"]["Rwp"]), "volume_fit": float(res["phases"][0]["volume"]), "R": R, "res": res})
    if not rows:
        return None
    best_rwp = min(r["Rwp"] for r in rows)
    ok = [r for r in rows if r["Rwp"] <= best_rwp * 1.10 + 0.005]
    win = min(ok, key=lambda r: r["volume_fit"])                          # a supercell that fits no better than the small cell is not wanted
    prog(0.7, "cell chosen")
    if steps is not None:
        tbl = [{"#": k + 1, "system": r["system"], "centring": r.get("centring", ""), "a (Å)": round(r["cell_fit"][0], 4), "b (Å)": round(r["cell_fit"][1], 4), "c (Å)": round(r["cell_fit"][2], 4),
                "β (°)": round(r["cell_fit"][4], 3), "volume (Å³)": round(r["volume_fit"], 1), "M20": r.get("M20"), "Rwp Le Bail (%)": round(100 * r["Rwp"], 2), "chosen": r is win}
               for k, r in enumerate(rows)]
        steps.add("cell", "Unit cell", "ok", f"{len(rows)} candidate cell(s) from indexing were fitted with Le Bail on the real pattern (no structure needed). "
                  f"Chosen: {win['system']}, V = {win['volume_fit']:.1f} Å³, Rwp = {100 * win['Rwp']:.1f} % – the lowest Rwp, and the smallest volume among cells that fit about equally well.", {"table": tbl})
    return {"cell": win["cell_fit"], "system": win["system"], "R": win["R"], "res": win["res"], "candidates": rows, "centring": win.get("centring"), "Rwp": win["Rwp"]}


# ------------------------------------------------------------------------------------------------------------------------------ 2. space group
def reflection_groups(R, res, wavelength):
    """The Le Bail reflection groups of the fitted phase: members (hkl), area, isolation, detection threshold, LP-corrected squared amplitude."""
    ph = R.phases[0]
    v = R._values()
    cell = R.cell(0, v)
    Gs = np.linalg.inv(RV._cell_metric(cell))
    inv_d2 = np.einsum("ki,ij,kj->k", ph.hkl, Gs, ph.hkl)
    d = 1.0 / np.sqrt(inv_d2)
    th = np.arcsin(np.clip(wavelength / (2 * d), 0, 0.9995))
    tt = 2 * np.degrees(th) + v["zero"] + v["disp"] * np.cos(th)
    tan = np.tan(th)
    H, _ = RV._tch(v["U"] * tan ** 2 + v["V"] * tan + v["W"], v["X"] / np.cos(th) + v["Y"] * tan)
    s = float(res["phases"][0]["scale"])
    c_lb = np.asarray(R.I_free[0], float) * s                         # Le Bail area (scaled counts x degree)
    resid = R.y - (res["calc"] * R.scale_y)
    # Matched-filter area of every group: the Le Bail area of an absent or weak reflection is biased upwards (noise below the background is clipped away);
    # adding the projection of the final residual on the group's own profile removes the bias: a = c + sum(P r) / sum(P^2)  (P = profile of unit area).
    I = c_lb.copy()
    sig_a = np.full(len(c_lb), np.inf)
    sig_pt = 1.4826 * float(np.median(np.abs(resid - np.median(resid)))) + 1e-12
    rad = list(R.rad)
    tt_k = np.empty((len(rad), len(c_lb)))
    for ri, (lam, wl) in enumerate(rad):
        sn = np.clip(lam / (2.0 * d), 0, 0.9995)
        th_ = np.arcsin(sn)
        tt_k[ri] = 2 * np.degrees(th_) + v["zero"] + v["disp"] * np.cos(th_)
    wsum = sum(w for _, w in rad)
    num = np.zeros(len(c_lb))
    den = np.zeros(len(c_lb))
    for k in range(len(c_lb)):
        lo_ = np.searchsorted(R.x, tt_k[:, k].min() - 3 * H[k])
        hi_ = np.searchsorted(R.x, tt_k[:, k].max() + 3 * H[k])
        if hi_ - lo_ < 3:
            continue
        xs = R.x[lo_:hi_]
        eta_k = RV._tch(np.array([v["U"] * tan[k] ** 2 + v["V"] * tan[k] + v["W"]]), np.array([v["X"] / np.cos(th[k]) + v["Y"] * tan[k]]))[1][0]
        P = np.zeros(len(xs))
        for ri, (lam, wl) in enumerate(rad):
            P += wl / wsum * RV._pv(xs - tt_k[ri, k], H[k], eta_k)
        d2 = float(np.sum(P * P)) * 1.0
        if d2 <= 0:
            continue
        num[k] = float(np.sum(P * resid[lo_:hi_]))
        den[k] = d2
        local = resid[lo_:hi_]
        sg = 1.4826 * float(np.median(np.abs(local - np.median(local)))) if len(local) > 6 else sig_pt
        sig_a[k] = max(sg, 0.3 * sig_pt) / math.sqrt(d2)
    ok_ = den > 0
    I[ok_] = c_lb[ok_] + num[ok_] / den[ok_]
    order = np.argsort(m_ := ph.m_grp)
    bounds = np.searchsorted(m_[order], np.arange(len(ph.hkl) + 1))
    members = [ph.m_hkl[order[bounds[k]:bounds[k + 1]]] for k in range(len(ph.hkl))]
    srt = np.argsort(tt)
    gap = np.full(len(tt), 99.0)
    tts = tt[srt]
    dif = np.diff(tts)
    gap[srt[1:]] = np.minimum(gap[srt[1:]], dif)
    gap[srt[:-1]] = np.minimum(gap[srt[:-1]], dif)
    out = []
    floor = 0.003 * float(np.max(I[ok_])) if ok_.any() else 0.0                  # a bump below 0.3 % of the strongest line is profile imperfection, not a reflection
    x0, x1 = float(R.x.min()), float(R.x.max())
    for k in range(len(ph.hkl)):
        thr = max(4.0 * float(sig_a[k]), floor) if np.isfinite(sig_a[k]) else float("inf")  # 4 standard deviations of the matched-filter area
        lp = float(xtal.lorentz_polarisation(tt[k], R.mono))
        inside = (tt[k] - 2.5 * H[k] > x0) and (tt[k] + 2.5 * H[k] < x1)               # a line at the edge of the pattern cannot be tested
        out.append({"hkl": members[k], "I": float(I[k]), "thr": float(thr), "iso": bool(gap[k] > 1.2 * H[k] and inside), "d": float(d[k]), "tt": float(tt[k]),
                    "H": float(H[k]), "G": float(max(I[k], 0.0) / max(lp, 1e-9)), "n_members": int(len(members[k])), "sigma": float(sig_a[k])})
    return out


def space_group_stage(groups, system, steps=None, max_list=12):
    ranked = SG.rank(system, groups)
    cons = [r for r in ranked if r["consistent"]]
    if steps is not None:
        top = ranked[:max_list]
        tbl = [{"space group": f"{r['label']} (No. {r['number']})", "centro": "yes" if r["centro"] else "no", "forbidden reflections found absent": r["n_abs"], "forbidden but present": r["n_viol"],
                "agrees": "yes" if r["consistent"] else "no"} for r in top]
        n_iso = sum(1 for g in groups if g["iso"])
        steps.add("sg", "Space group", "ok" if cons else "warning",
                  f"{n_iso} isolated reflection groups were used. A space group is kept when every reflection it forbids is absent in the pattern. "
                  + (f"{len(cons)} group(s) agree; the most informative one (most confirmed absences, centrosymmetric first) is tried first." if cons else
                     "No group agrees fully (impurity or Le Bail artefacts?): the best-scoring groups are tried."), {"table": tbl})
    return ranked


# ------------------------------------------------------------------------------------------------------------------------------ 3. charge flipping
def build_problem(R, groups, ranked_item, dmin_cf, normalise=True):
    """Charge-flipping problem for one space group: reflections (full sphere) of the Le Bail groups down to dmin_cf with the amplitudes from the intensities."""
    ph = R.phases[0]
    cell = R.cell(0, R._values())
    keep = [k for k, g in enumerate(groups) if g["d"] >= dmin_cf and g["I"] > 0]
    hkl, grp = [], []
    for newk, k in enumerate(keep):
        h = groups[k]["hkl"]
        hkl.append(np.vstack([h, -h]))
        grp.append(np.full(2 * len(h), newk))
    hkl = np.rint(np.vstack(hkl)).astype(int)
    grp = np.concatenate(grp)
    Rm, tm = SG.operations(ranked_item["number"], ranked_item["perm"])
    forb = SG.forbidden(Rm, tm, hkl)
    Rp, tp = SG.point_operations(ranked_item["number"], ranked_item["perm"])
    s2 = np.array([1.0 / groups[k]["d"] ** 2 / 4.0 for k in keep])
    Gtot = np.array([groups[k]["G"] for k in keep])
    n_allowed = np.bincount(grp[~forb], minlength=len(keep)).astype(float)
    mean_f2 = Gtot / np.maximum(n_allowed, 1.0)
    if normalise:
        E2 = CF.normalise_amplitudes(s2, mean_f2, np.maximum(n_allowed, 1.0))
    else:
        E2 = mean_f2
    P = CF.Problem(cell, hkl, grp, E2 * np.maximum(n_allowed, 1.0), dmin_cf, sym=(Rp, tp), forbidden=forb)
    P.raw_G = Gtot
    P.keep = keep
    P.s2 = s2
    return P


# ------------------------------------------------------------------------------------------------------------------------------ 4. atoms
def orbits_from_peaks(pk, cell, Rm, tm, tol=0.45, min_rel=0.10):
    """Group density peaks into symmetry orbits. pk = [(x, y, z, height)] strongest first. Returns [{'x': (3,), 'h': height, 'mult': n, 'images': (n,3)}]."""
    Gm = RV._cell_metric(cell)
    pts = np.array([p[:3] for p in pk])
    hts = np.array([p[3] for p in pk])
    used = np.zeros(len(pts), bool)
    orbs = []
    hmax = hts.max() if len(hts) else 1.0
    for i in range(len(pts)):
        if used[i] or hts[i] < min_rel * hmax:
            continue
        x_i = _snap_to_special_position(pts[i], Rm, tm, Gm)
        imgs = np.mod(np.einsum("nij,j->ni", Rm, x_i) + tm, 1.0)
        # distinct images
        u = []
        for im in imgs:
            if not any(_pbc_dist(im, q, Gm) < 0.25 for q in u):
                u.append(im)
        u = np.array(u)
        members = []
        for j in range(len(pts)):
            if not used[j] and np.min([_pbc_dist(pts[j], q, Gm) for q in u]) < tol:
                members.append(j)
        for j in members:
            used[j] = True
        orbs.append({"x": x_i, "h": float(np.mean(hts[members])), "mult": int(len(u)), "images": u})
    return orbs


def _snap_to_special_position(x, Rm, tm, Gm, tol=0.6, rounds=3):
    """A density peak sits close to, not exactly on, a special position (mirror plane, axis, centre): the operations that map it onto itself within `tol` Å are
    averaged, which puts it exactly on the position. A general position stays as it is."""
    x = np.asarray(x, float).copy()
    for _ in range(rounds):
        imgs = np.einsum("nij,j->ni", Rm, x) + tm
        dlt = imgs - x
        dlt -= np.round(dlt)
        dist = np.sqrt(np.einsum("ni,ij,nj->n", dlt, Gm, dlt))
        sel = dist < tol
        if sel.sum() <= 1:
            break
        step = dlt[sel].mean(axis=0)
        x = x + step
        if np.linalg.norm(step) < 1e-6:
            break
    return np.mod(x, 1.0)


def _pbc_dist(a, b, Gm):
    d = np.asarray(a) - np.asarray(b)
    d -= np.round(d)
    return float(np.sqrt(d @ Gm @ d))


class AtomModel:
    """Orbits with elements; fast R factor of the group amplitudes."""

    def __init__(self, orbs, cell, hkl_full, grp, G_raw, s2, elements, n_expected=None):
        self.orbs = orbs
        self.cell = cell
        self.elements = list(elements)
        self.grp = grp
        self.G_raw = np.asarray(G_raw, float)
        self.s2g = np.asarray(s2, float)
        Gs = np.linalg.inv(RV._cell_metric(cell))
        self.hkl = np.asarray(hkl_full, float)
        inv = np.einsum("ki,ij,kj->k", self.hkl, Gs, self.hkl)
        self.s2r = inv / 4.0
        self.fe = {el: xtal.f_xray(el, self.s2r) for el in self.elements}
        # orbit phase sums S_o(h) = sum over images exp(2 pi i h.x)
        self.S = np.stack([np.exp(2j * np.pi * (self.hkl @ o["images"].T)).sum(axis=1) for o in orbs], 1) if orbs else np.zeros((len(self.hkl), 0), complex)
        self.n_expected = n_expected

    def rfactor(self, assign, active, return_scale=False):
        idx = [i for i in range(len(self.orbs)) if active[i]]
        if not idx:
            return 1.0
        F = np.zeros(len(self.hkl), complex)
        for i in idx:
            F += self.fe[assign[i]] * self.S[:, i]
        f2 = np.abs(F) ** 2
        Gc = np.bincount(self.grp, weights=f2, minlength=len(self.G_raw))
        ok = (self.G_raw > 0) & (Gc > 0)
        if ok.sum() < 3:
            return 1.0
        o, c = np.sqrt(self.G_raw[ok]), np.sqrt(Gc[ok])
        s2 = self.s2g[ok]
        best, bB, bK = 9.0, 0.0, 1.0
        for B in (0.0, 1.0, 2.0, 3.5, 5.0, 8.0):
            cb = c * np.exp(-B * s2)
            K = float(np.sum(o * cb) / max(np.sum(cb * cb), 1e-30))
            r = float(np.sum(np.abs(o - K * cb)) / np.sum(o))
            if r < best:
                best, bB, bK = r, B, K
        self._B, self._K = bB, bK
        return (best, bB, bK) if return_scale else best


def _score(model, assign, active, counts, mult, w_missing=0.35):
    r = model.rfactor(assign, active)
    got = {}
    for i in range(len(assign)):
        if active[i]:
            got[assign[i]] = got.get(assign[i], 0.0) + mult[i]
    tot = max(sum(counts.values()), 1e-9)
    miss = sum(abs(counts[e] - got.get(e, 0.0)) for e in counts) / tot
    return r + w_missing * miss, r, miss


def assign_one(model: AtomModel, counts, order, mult, sweeps=3):
    """Elements to orbits for one assumed number of formula units: the heaviest element takes the highest peaks until its count is reached, then the next element, ...;
    then single-orbit changes (exchange of two elements, swap of an active and an inactive orbit, deactivation) are accepted when the score (R factor of the amplitudes +
    0.35 x the fraction of atoms that does not agree with the composition) improves."""
    els = sorted(counts, key=lambda e: -xtal.atomic_number(e))
    n = len(model.orbs)
    assign = [els[-1]] * n
    active = np.zeros(n, bool)
    free = list(order)
    for e in els:
        need = counts[e]
        got = 0.0
        for i in list(free):
            if got >= need * 0.999:
                break
            if got + mult[i] <= need * 1.02 + 0.01 or got == 0.0 and mult[i] <= need * 1.5:
                assign[i] = e
                active[i] = True
                got += mult[i]
                free.remove(i)
    best, r, miss = _score(model, assign, active, counts, mult)
    for _ in range(sweeps):
        improved = False
        idx_active = [i for i in order if active[i]]
        # change the element of an active orbit
        for i in idx_active:
            for e in els:
                if e == assign[i]:
                    continue
                old = assign[i]
                assign[i] = e
                sc, _, _ = _score(model, assign, active, counts, mult)
                if sc < best - 1e-4:
                    best, improved = sc, True
                else:
                    assign[i] = old
        # exchange an active orbit for an inactive one of the same multiplicity
        for i in idx_active:
            for j in order:
                if active[j] or mult[j] != mult[i]:
                    continue
                active[i], active[j] = False, True
                assign[j] = assign[i]
                sc, _, _ = _score(model, assign, active, counts, mult)
                if sc < best - 1e-4:
                    best, improved = sc, True
                    break
                active[i], active[j] = True, False
        # activate / deactivate
        for i in order:
            if active[i]:
                active[i] = False
                sc, _, _ = _score(model, assign, active, counts, mult)
                if sc < best - 1e-4:
                    best, improved = sc, True
                else:
                    active[i] = True
            else:
                active[i] = True
                for e in els:
                    assign[i] = e
                    sc, _, _ = _score(model, assign, active, counts, mult)
                    if sc < best - 1e-4:
                        best, improved = sc, True
                        break
                else:
                    active[i] = False
        if not improved:
            break
    sc, r, miss = _score(model, assign, active, counts, mult)
    return {"assign": assign, "active": active, "score": sc, "R": r, "missing": miss}


def assign_elements(model: AtomModel, composition, z_candidates):
    """Assign elements to the orbits for every candidate number of formula units Z in the cell; the best score wins. composition: {element: amount per formula unit} (H ignored)."""
    comp = {e: float(v) for e, v in composition.items() if xtal.atomic_number(e) >= _ASSIGNABLE_MIN_Z and v > 0}
    if not comp or not model.orbs:
        return None
    mult = np.array([o["mult"] for o in model.orbs], float)
    order = list(np.argsort(-np.array([o["h"] for o in model.orbs])))
    best = None
    for Zf in z_candidates:
        counts = {e: Zf * v for e, v in comp.items()}
        a = assign_one(model, counts, order, mult)
        a["Z"] = Zf
        if best is None or a["score"] < best["score"] - 1e-9:
            best = a
    return best


def _orbit_images(x, Rm, tm, tol=0.25, Gm=None):
    imgs = np.mod(np.einsum("nij,j->ni", Rm, x) + tm, 1.0)
    u = []
    for im in imgs:
        if not any(_pbc_dist(im, q, Gm) < tol for q in u):
            u.append(im)
    return np.array(u)


def refine_positions(model: AtomModel, assign, active, Rm, tm, max_shift_A=0.6, max_nfev=25):
    """Least-squares refinement of the orbit coordinates against the group amplitudes (R factor residual). The symmetrised structure factor of an orbit is
    (m / n_ops) x the sum over ALL operations, so an orbit near a special position stays symmetric. The result is re-snapped to special positions."""
    from scipy.optimize import least_squares
    idx = [i for i in range(len(model.orbs)) if active[i]]
    if not idx:
        return
    Gm = RV._cell_metric(model.cell)
    n_ops = len(Rm)
    x0 = np.array([model.orbs[i]["x"] for i in idx], float)
    mult = np.array([model.orbs[i]["mult"] for i in idx], float)
    sg = np.sqrt(np.maximum(model.G_raw, 0.0))
    okg = sg > 0
    amax = max_shift_A / np.array(model.cell[:3])
    fe = [model.fe[assign[i]] for i in idx]

    def struct_f2(X, Bv=None):
        F = np.zeros(len(model.hkl), complex)
        for k in range(len(idx)):
            imgs = np.einsum("nij,j->ni", Rm, X[k]) + tm
            S = np.exp(2j * np.pi * (model.hkl @ imgs.T)).sum(axis=1) * (mult[k] / n_ops)
            F += fe[k] * S
        return F

    def resid(p):
        X = p.reshape(-1, 3)
        F = struct_f2(X)
        Gc = np.bincount(model.grp, weights=np.abs(F) ** 2, minlength=len(model.G_raw))
        c = np.sqrt(np.maximum(Gc, 0.0))
        best = None
        for B in (0.0, 1.0, 2.0, 3.5, 5.0):
            cb = c * np.exp(-B * model.s2g)
            K = float(np.sum(sg[okg] * cb[okg]) / max(np.sum(cb[okg] ** 2), 1e-30))
            r = (sg[okg] - K * cb[okg])
            ss = float(np.sum(r * r))
            if best is None or ss < best[0]:
                best = (ss, r)
        return np.concatenate([best[1] / max(float(np.sum(sg[okg])), 1e-30) * 50.0, 0.02 * (p - x0.ravel()) / amax_flat])
    amax_flat = np.tile(amax, len(idx))
    lo = (x0 - amax).ravel() if False else (x0 - amax[None, :]).ravel()
    hi = (x0 + amax[None, :]).ravel()
    try:
        sol = least_squares(resid, x0.ravel(), bounds=(lo, hi), max_nfev=max_nfev, diff_step=1e-3)
    except Exception:
        return
    X = sol.x.reshape(-1, 3)
    for k, i in enumerate(idx):
        xs = _snap_to_special_position(X[k], Rm, tm, Gm, tol=0.35)
        imgs = _orbit_images(xs, Rm, tm, 0.25, Gm)
        if abs(len(imgs) - model.orbs[i]["mult"]) > 0:
            xs = np.mod(X[k], 1.0)
            imgs = _orbit_images(xs, Rm, tm, 0.25, Gm)
            if abs(len(imgs) - model.orbs[i]["mult"]) > 0:
                continue                                       # the multiplicity would change: keep the old position
        model.orbs[i]["x"] = xs
        model.orbs[i]["images"] = imgs
        model.S[:, i] = np.exp(2j * np.pi * (model.hkl @ imgs.T)).sum(axis=1)


def _extend_model(model: AtomModel, new_orbs):
    model.orbs = list(model.orbs) + list(new_orbs)
    cols = [np.exp(2j * np.pi * (model.hkl @ o["images"].T)).sum(axis=1) for o in new_orbs]
    model.S = np.hstack([model.S, np.stack(cols, 1)]) if cols else model.S


def difference_peaks(P, model: AtomModel, assign, active, Rm, tm, min_rel=0.22, min_dist=0.9):
    """Difference Fourier map (|Fo| - |Fc|, phases of the model): positive peaks are atoms that the model still misses. |Fo| of every reflection is the group amplitude
    shared in the proportion of the model's own |Fc|^2 (the usual Le Bail partition). Returns new symmetry orbits (not on an existing atom)."""
    r, B, K = model.rfactor(assign, active, return_scale=True)
    F = np.zeros(len(model.hkl), complex)
    for i in range(len(model.orbs)):
        if active[i]:
            F += model.fe[assign[i]] * model.S[:, i]
    Fc = F * np.exp(-B * model.s2r)
    a2 = np.abs(Fc) ** 2
    grp = model.grp
    tot = np.bincount(grp, weights=a2, minlength=len(model.G_raw))
    members = np.bincount(grp, minlength=len(model.G_raw)).astype(float)
    share = np.where(tot[grp] > 1e-12, a2 / np.maximum(tot[grp], 1e-12), 1.0 / np.maximum(members[grp], 1.0))
    Fo = np.sqrt(np.maximum(model.G_raw[grp] * share, 0.0)) / max(K, 1e-12)
    delta = (Fo - np.abs(Fc)) * np.where(P.forbidden, 0.0, 1.0)
    coef = delta * np.exp(1j * np.angle(Fc + 1e-30))
    grid = np.zeros(P.n, complex)
    grid.reshape(-1)[P.flat] = coef
    rho = np.fft.ifftn(grid).real
    pk = CF.peaks(rho, model.cell, n_max=40, min_rel=min_rel, min_dist=0.7)
    if not pk:
        return []
    orbs = orbits_from_peaks(pk, model.cell, Rm, tm, tol=0.5, min_rel=min_rel)
    Gm = RV._cell_metric(model.cell)
    have = [im for o, on in zip(model.orbs, active) if on for im in o["images"]]
    out = []
    for o in orbs:
        if have and min(_pbc_dist(o["x"], q, Gm) for q in have) < min_dist:
            continue
        out.append(o)
    return out


def complete_model(P, model: AtomModel, res, counts, Rm, tm, max_cycles=6, log=None):
    """Structure completion: while the composition is not reached or the R factor is high, add the strongest difference-map peak as the element that lacks most
    (accepted only when the score improves). Works on a copy of the assignment."""
    assign = list(res["assign"])
    active = list(res["active"])
    mult = np.array([o["mult"] for o in model.orbs], float)
    els = sorted(counts, key=lambda e: -xtal.atomic_number(e))
    score, r, miss = _score(model, assign, np.array(active), counts, mult)
    for cyc in range(max_cycles):
        got = {}
        for i, o in enumerate(model.orbs):
            if active[i]:
                got[assign[i]] = got.get(assign[i], 0.0) + o["mult"]
        lacking = [e for e in els if got.get(e, 0.0) < counts[e] - 0.5]
        if not lacking and r < 0.08:
            break
        cands = difference_peaks(P, model, assign, np.array(active), Rm, tm)
        if not cands:
            break
        n0 = len(model.orbs)
        _extend_model(model, cands[:6])
        assign = assign + [els[-1]] * (len(model.orbs) - n0)
        active = active + [False] * (len(model.orbs) - n0)
        mult = np.array([o["mult"] for o in model.orbs], float)
        best = (score, None, None)
        for j in range(n0, len(model.orbs)):
            for e in (lacking or els):
                if got.get(e, 0.0) + model.orbs[j]["mult"] > counts[e] + 0.5:
                    continue
                assign[j] = e
                active[j] = True
                sc, _, _ = _score(model, assign, np.array(active), counts, mult)
                active[j] = False
                if sc < best[0] - 1e-4:
                    best = (sc, j, e)
        if best[1] is None:
            model.orbs = model.orbs[:n0]
            model.S = model.S[:, :n0]
            assign, active = assign[:n0], active[:n0]
            break
        sc, j, e = best
        assign[j], active[j] = e, True
        # drop the other trial orbits
        keep = list(range(n0)) + [j]
        model.orbs = [model.orbs[k] for k in keep]
        model.S = model.S[:, keep]
        assign = [assign[k] for k in keep]
        active = [active[k] for k in keep]
        mult = np.array([o["mult"] for o in model.orbs], float)
        score, r, miss = _score(model, assign, np.array(active), counts, mult)
        if log is not None:
            log.append(f"completion {cyc + 1}: +{e} at {np.round(model.orbs[-1]['x'], 3).tolist()} (x{model.orbs[-1]['mult']}) -> R(F) {100 * r:.1f} %")
    return {"assign": assign, "active": np.array(active), "score": score, "R": r, "missing": miss, "Z": res.get("Z")}


def z_candidates_for(cell, composition, z_formula=None, vmin=7.0, vmax=45.0, limit=6):
    """Plausible numbers of formula units per cell: the volume per non-hydrogen atom between vmin and vmax (Å^3; dense oxides ~ 8-15, porous frameworks up to ~40)."""
    if z_formula:
        return [int(z_formula)]
    n = sum(v for e, v in composition.items() if xtal.atomic_number(e) >= _ASSIGNABLE_MIN_Z)
    V = RV.cell_volume(tuple(cell))
    zs = [z for z in range(1, 97) if vmin <= V / (z * max(n, 1e-9)) <= vmax]
    zs.sort(key=lambda z: abs(math.log(V / (z * n) / 14.0)))
    return sorted(zs[:limit]) or [1]

# ------------------------------------------------------------------------------------------------------------------------------ model -> structure
def structure_from_orbits(orbs, assign, active, cell):
    from pymatgen.core import Lattice, Structure
    species, coords = [], []
    for o, a, on in zip(orbs, assign, active):
        if not on:
            continue
        for im in o["images"]:
            species.append(a)
            coords.append(im)
    return Structure(Lattice.from_parameters(*cell), species, coords, coords_are_cartesian=False)


def cif_text(structure, symprec=0.05, comment=""):
    import warnings
    warnings.filterwarnings("ignore")
    from pymatgen.io.cif import CifWriter
    txt = str(CifWriter(structure, symprec=symprec))
    if comment:
        txt = "# " + comment.replace("\n", "\n# ") + "\n" + txt
    return txt


def parse_formula(text):
    """'Zr6O4(OH)4(C8H4O4)6' or 'Zr O C' -> {element: amount}. Hydrogen is kept here and ignored later."""
    import re
    import warnings
    warnings.filterwarnings("ignore")
    from pymatgen.core import Composition
    t = (text or "").strip()
    if not t:
        return {}
    if re.fullmatch(r"[A-Z][a-z]?([\s,;]+[A-Z][a-z]?)*", t):                  # only element symbols: one of each
        return {e: 1.0 for e in re.split(r"[\s,;]+", t)}
    return {str(k): float(v) for k, v in Composition(t).get_el_amt_dict().items()}


def report_text(out, name="pattern", version=""):
    """Plain-text report of every decision (for the ZIP and the SI)."""
    L = [f"Reticula {version} - structure solution from a powder pattern: {name}", "=" * 72, ""]
    for s in out["steps"]:
        L.append(f"[{s['status'].upper()}] {s['title']}")
        L.append("  " + s["text"])
        tb = s.get("data", {}).get("table")
        if tb:
            keys = list(tb[0].keys())
            L.append("  " + " | ".join(keys))
            for r in tb[:30]:
                L.append("  " + " | ".join(str(r.get(k, "")) for k in keys))
        L.append("")
    if out.get("ok"):
        c = out["cell"]
        L += [f"Cell: a = {c[0]:.4f}, b = {c[1]:.4f}, c = {c[2]:.4f} A, alpha = {c[3]:.2f}, beta = {c[4]:.2f}, gamma = {c[5]:.2f} deg", f"Space group tried: {out['space_group']['label']} (No. {out['space_group']['number']})",
              f"R(F) of the Le Bail amplitudes: {100 * out['R_F']:.1f} %", ""]
    L += ["Warnings and limits:"] + ["  - " + w for w in out.get("warnings", [])]
    return "\n".join(L)


# ------------------------------------------------------------------------------------------------------------------------------ pipeline
def solve(x, y, radiation="Cu Kα1", composition=None, z_formula=None, cell_given=None, space_group_given=None, instrument=None, zero_shift=0.0,
          systems=("cubic", "tetragonal", "hexagonal", "orthorhombic"), dmin_cf=None, n_trials=12, n_iter=300, n_sg_try=2, n_sg_max=6, good_R=0.10, threads=None, progress=None, refine=True, seed=0,
          offset="displacement"):
    """The whole decision tree. composition: {'Zr': 6, 'O': 32, 'C': 48} (per formula unit; hydrogen is ignored) or None (then the atoms are only located and named by height).
    Returns a dict: steps (the decision log), cell, space_group, structure (pymatgen), cif (text), R_F, refinement (Rietveld result or None), warnings."""
    t0 = time.time()
    prog = progress or (lambda f, t: None)
    steps = Steps()
    warns = []
    x = np.asarray(x, float) - float(zero_shift)
    y = np.asarray(y, float)
    lam = RV.RADIATION[radiation][0][0]
    # 0 -----------------------------------------------------------------------------------------------------------------
    prog(0.01, "checking the data")
    dc = data_check(x, y, lam)
    steps.add("data", "Data check", "ok" if not dc["notes"] else "warning",
              f"{dc['n_peaks']} peaks, 2θ {dc['tt_min']:.1f}–{dc['tt_max']:.1f}° (d down to {dc['dmin']:.2f} Å), step {dc['step']:.3f}°."
              + (f" Zero shift {zero_shift:+.3f}° from the reference sample was applied." if zero_shift else " No reference sample: a zero shift or sample displacement is refined in the Le Bail fit.")
              + ((" " + " ".join(dc["notes"])) if dc["notes"] else ""))
    warns += dc["notes"]
    # 1 -----------------------------------------------------------------------------------------------------------------
    pc = pick_cell(x, y, lam, dc["peaks"], radiation, systems, 4, instrument, threads, offset, lambda f, t: prog(0.02 + 0.30 * f, t), steps, cell_given)
    if pc is None:
        return {"steps": steps, "ok": False, "warnings": warns + ["no unit cell"], "seconds": time.time() - t0}
    R, res = pc["R"], pc["res"]
    cell = pc["cell"]
    system = pc["system"]
    # 2 -----------------------------------------------------------------------------------------------------------------
    prog(0.34, "space group from the systematic absences")
    groups = reflection_groups(R, res, lam)
    if space_group_given:
        ranked = [r for r in SG.rank(system, groups) if r["number"] == int(space_group_given)] or SG.rank(system, groups)[:1]
        steps.add("sg", "Space group", "ok", f"Given: No. {space_group_given}.")
    else:
        ranked = space_group_stage(groups, system, steps)
    cons = [r for r in ranked if r["consistent"]]
    trylist = (cons or ranked)[: max(n_sg_max, 1)]
    # 3 + 4 -------------------------------------------------------------------------------------------------------------
    dmin = dmin_cf or max(dc["dmin"], 1.4)
    comp = {k: v for k, v in (composition or {}).items() if xtal.atomic_number(k) >= _ASSIGNABLE_MIN_Z}
    best = None
    tried = []
    log_lines = []
    for si, sgi in enumerate(trylist):
        if best is not None and best["R_F"] < good_R and si >= max(n_sg_try, 1):
            break                                                       # a good model was found in the first groups: the lower-ranked groups are not needed
        f0 = 0.38 + 0.4 * si / len(trylist)
        f1 = 0.38 + 0.4 * (si + 1) / len(trylist)
        prog(f0, f"charge flipping in {sgi['label']}")
        try:
            P = build_problem(R, groups, sgi, dmin)
        except Exception as e:
            tried.append({"space group": sgi["label"], "error": str(e)})
            continue
        results = CF.run_trials(P, n_trials=n_trials, n_iter=n_iter, threads=threads, seed0=seed + 1000 * si, progress=lambda f, t: prog(f0 + (f1 - f0) * f, f"{sgi['label']}: {t}"))
        Rm_, tm_ = SG.operations(sgi["number"], sgi["perm"])
        zc = z_candidates_for(cell, comp, z_formula) if comp else []
        for ti, r in enumerate(results[: max(3, n_trials // 3)]):
            pk = CF.peaks(r["rho"], cell, min_rel=0.05, n_max=300)
            if len(pk) < 1:
                continue
            orbs = orbits_from_peaks(pk, cell, Rm_, tm_, min_rel=0.05)
            if not orbs or not comp:
                continue
            import copy
            for Zf in zc:                                                 # every plausible number of formula units is carried through the completion: the final score decides
                orbs_z = copy.deepcopy(orbs)
                m = AtomModel(orbs_z, cell, P.hkl, P.group, P.raw_G, P.s2, list(comp), None)
                a = assign_elements(m, comp, [Zf])
                if a is None:
                    continue
                counts = {e: Zf * v for e, v in comp.items()}
                lines = []
                refine_positions(m, a["assign"], a["active"], Rm_, tm_)
                a = complete_model(P, m, a, counts, Rm_, tm_, log=lines)
                refine_positions(m, a["assign"], a["active"], Rm_, tm_)
                mult_ = np.array([o["mult"] for o in m.orbs], float)
                a["score"], a["R"], a["missing"] = _score(m, a["assign"], np.array(a["active"]), counts, mult_)
                a["Z"] = Zf
                tried.append({"space group": sgi["label"], "start": ti, "R charge flipping (%)": round(100 * r["R_final"], 1), "formula units Z": Zf,
                              "atoms placed": int(sum(o["mult"] for o, on in zip(m.orbs, a["active"]) if on)), "R(F) (%)": round(100 * a["R"], 1),
                              "composition mismatch (%)": round(100 * float(a["missing"]), 1)})
                if best is None or a["score"] < best["score"] - 1e-9:
                    best = {"R_F": a["R"], "score": a["score"], "Z": Zf, "sg": sgi, "orbs": list(m.orbs), "assign": list(a["assign"]), "active": np.array(a["active"]), "trial": r, "P": P}
                    log_lines = lines
    if best is None:
        steps.add("cf", "Structure solution (charge flipping)", "failed",
                  ("No composition was given: the atoms cannot be named. Give the elements of the compound (for example Zr, O, C) so that the peaks can be assigned." if not comp else
                   "Charge flipping did not give an atomic model that explains the intensities."), {"table": tried})
        return {"steps": steps, "ok": False, "cell": cell, "system": system, "warnings": warns, "groups": groups, "ranked": ranked, "R": R, "lebail": res, "seconds": time.time() - t0}
    steps.add("cf", "Structure solution (charge flipping)", "ok",
              f"{n_trials} random starts per space group, run in parallel; the density peaks of the best starts were turned into atoms by the R factor of the intensities. "
              f"Best: {best['sg']['label']}, R(F) = {100 * best['R_F']:.1f} % over the Le Bail amplitudes (about 20–35 % is typical for a correct heavy-atom model from powder data).", {"table": tried})
    struct = structure_from_orbits(best["orbs"], best["assign"], best["active"], cell)
    comment = ("Reticula 2.18.4 structure solution from powder data (charge flipping). Atom positions are a STARTING MODEL: light atoms may be missing, hydrogen is not placed; "
               "refine it with restraints before use. Space group " + best["sg"]["label"])
    txt = cif_text(struct, 0.05, comment)
    steps.add("atoms", "Atoms", "ok", f"{len(struct)} atoms in the cell: " + ", ".join(f"{e} {n}" for e, n in sorted(struct.composition.get_el_amt_dict().items(), key=lambda kv: -kv[1])) +
              f". {int(best['active'].sum())} symmetry-independent positions.", {})
    out = {"steps": steps, "ok": True, "cell": cell, "system": system, "space_group": best["sg"], "structure": struct, "cif": txt, "R_F": best["R_F"], "warnings": warns,
           "groups": groups, "ranked": ranked, "R": R, "lebail": res, "orbs": best["orbs"], "assign": best["assign"], "active": best["active"], "seconds": time.time() - t0}
    # 5 -----------------------------------------------------------------------------------------------------------------
    if refine:
        prog(0.84, "Rietveld refinement of the model")
        try:
            ph = RV.Phase.from_cif(txt.encode("utf-8"), "solved structure")
            res_r = RV.quick_refine(x, y, [ph], radiation, "rietveld", None, 6, None, 0.0, None, offset, None, True, 1.0, instrument, False, None)
            res_r["_phase_objs"] = [ph]
            out["refinement"] = res_r
            steps.add("rietveld", "Rietveld refinement", "ok", f"Rwp = {100 * res_r['stats']['Rwp']:.1f} %, GoF = {res_r['stats']['chi2'] ** 0.5:.2f}. Cell, scale, background, peak shape, atom positions (soft bond restraints) "
                      "and one displacement parameter per element were refined on the real pattern. Rwp below about 15 % with a smooth difference curve means the model describes the pattern; "
                      "a large Rwp means missing atoms (light atoms, guests in the pores) or a wrong space group.", {})
        except Exception as e:
            steps.add("rietveld", "Rietveld refinement", "failed", f"Not possible: {e}")
    warns.append("A structure from powder data is a model: check the bond lengths, the coordination of the metal and the pore content (guest molecules and solvent are not found). "
                 "Light atoms and hydrogen are often not located.")
    out["warnings"] = warns
    prog(1.0, "done")
    return out
