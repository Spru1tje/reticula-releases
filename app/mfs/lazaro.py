r"""v2.18.4: the multifaceted TGA method of Abánades Lázaro (Eur. J. Inorg. Chem. 2020, 4284-4294, DOI 10.1002/ejic.202000656) for the exact composition of
defective, multivariate and composite MOFs.

The idea (written per metal atom, then multiplied by the metals per node):

    R_exp = M_w(structure) / M_w(residue)            (the experimental structure : residue mass ratio, from the TGA curve)

The structure is built from a fixed building unit (SBU), x linkers, modulators in a known molar ratio to the linker (from 1H NMR, HPLC, ICP ...), optionally
fixed counter-ions (Cl- from EDX/ICP) and ONE balancing species (OH-, O2- or Cl-) that closes the charge balance:

    charge:   q_avail = x * sum_i r_i z_i  + sum fixed n z  + y z_bal
    mass:     R * (M_res + x * sum_i r_i f_i M_i) = M_sbu + x * sum_i r_i M_i + sum fixed n M + y M_bal          (f_i: fraction of component i left in the residue)

Both equations are linear in x, so the composition follows exactly (no assumption about which decomposition step belongs to which species):

    x = [ R M_res - M_sbu - sum fixed n M - (q_avail - sum fixed n z) M_bal/z_bal ]  /  [ sum r_i M_i - (sum r_i z_i) M_bal/z_bal - R sum r_i f_i M_i ]

Two stages are used: the desolvated structure (about 200 degrees C: SBU still hydroxylated, modulators present) and the dehydroxylated structure (DH MOF, just
before the linker burns off: modulators gone, SBU dehydroxylated, charge closed by O2-). Composites: an inert second material (oxide / metal) is removed with its
share of the residue; a material that burns away before the linker is not in the DH mass and its weight percent follows from the start ratio; a material that
burns away together with the linker is removed with its known weight fraction.
"""
from __future__ import annotations

from dataclasses import dataclass

from .chem import M_CL, M_O, M_OH


@dataclass
class Part:
    name: str
    M: float          # molar mass of the bound (deprotonated) form
    z: float          # anionic charge it takes from the node (0 neutral, negative = cation)
    r: float = 1.0    # molar amount relative to the main linker
    f: float = 0.0    # fraction of its mass that stays in the residue (partial decomposition)


def solve_per_metal(R, M_res, M_sbu, q_avail, linkers, mods=(), fixed=(), balance=None):
    """x = main linker per metal (and every other amount) from the experimental ratio R.
    linkers, mods: lists of Part (r relative to the main linker, which is linkers[0] with r = 1); fixed: [(Part, n_per_metal)];
    balance: Part(M, z) that closes the charge (None = no closure: the charge left is reported)."""
    parts = list(linkers) + list(mods)
    A = sum(p.r * p.M for p in parts)
    C = sum(p.r * p.z for p in parts)
    P = sum(p.r * p.f * p.M for p in parts)
    Mf = sum(n * p.M for p, n in fixed)
    Qf = sum(n * p.z for p, n in fixed)
    if balance is not None and balance.z:
        k = balance.M / balance.z
        num = R * M_res - M_sbu - Mf - (q_avail - Qf) * k
        den = A - C * k - R * P
    else:
        num = R * M_res - M_sbu - Mf
        den = A - R * P
    if abs(den) < 1e-12:
        raise ValueError("the denominator vanishes (the components have no mass left to count)")
    x = num / den
    y = (q_avail - Qf - x * C) / balance.z if (balance is not None and balance.z) else 0.0
    left = q_avail - Qf - x * C - (y * balance.z if balance is not None and balance.z else 0.0)
    amounts = {p.name: x * p.r for p in parts}
    mass = M_sbu + Mf + x * A + (y * balance.M if balance is not None else 0.0)
    R_calc = mass / (M_res + x * P)
    return {"x": x, "amounts": amounts, "y": y, "charge_left": left, "mass": mass, "R_calc": R_calc, "A": A, "C": C, "P": P}


def forward_R(x, M_res, M_sbu, q_avail, linkers, mods=(), fixed=(), balance=None):
    """The experimental ratio that a given composition would give (used for tests and for the 'expected R' of the ideal structure)."""
    parts = list(linkers) + list(mods)
    A = sum(p.r * p.M for p in parts); C = sum(p.r * p.z for p in parts); P = sum(p.r * p.f * p.M for p in parts)
    Mf = sum(n * p.M for p, n in fixed); Qf = sum(n * p.z for p, n in fixed)
    y = (q_avail - Qf - x * C) / balance.z if (balance is not None and balance.z) else 0.0
    return (M_sbu + Mf + x * A + (y * balance.M if balance is not None else 0.0)) / (M_res + x * P)


def composite_weight_fraction(M_struct_per_metal, M_res, R_start):
    """Mass fraction of a second material that burns away completely BEFORE the linker (it is in the start mass, not in the DH mass):
    the MOF part is M_struct, its residue M_res, so  R_start = M_struct / (S_MOF * M_res)  with S_MOF = mass fraction of the MOF  =>  w = 1 - M_struct / (R_start M_res)."""
    return 1.0 - M_struct_per_metal / (R_start * M_res)


# ------------------------------------------------------------------------------------------------ the two methods on a DefectInput
BALANCE = {"OH": ("OH-", M_OH, 1.0), "O": ("O2-", M_O, 2.0), "Cl": ("Cl-", M_CL, 1.0)}


def _balance(key):
    if key in BALANCE:
        n, M, z = BALANCE[key]
        return Part(n, M, z)
    return None


def _collect(inp, stage):
    """(linkers, mods, message). Linkers: the main linker first (r = 1). Ratios come from the NMR integrals (a_i = A_i / H_i).
    Start stage: every species with an integral. DH stage: linkers + the species flagged 'at the plateau'."""
    ratios = inp.nmr_ratios()
    ml = inp.main_linker()
    a0 = ratios.get(ml.name)
    linkers = [Part(ml.name, ml.M, ml.z, 1.0, ml.residue_frac)]
    mods = []
    others = [s for s in inp.species if s.name != ml.name and s.name in ratios]
    if a0 is None:
        if others:
            return None, None, "needs the NMR integral of the main linker (it sets the scale of the other species)"
        return linkers, mods, ""
    for s in others:
        if stage == "dh" and s.role != "linker" and not s.at_plateau:
            continue                                   # modulators that burned off before the linker are not in the DH structure
        p = Part(s.name, s.M, s.z, ratios[s.name] / a0, s.residue_frac)
        (linkers if s.role == "linker" else mods).append(p)
    return linkers, mods, ""


def _stage_R(inp, W, stage):
    """R_exp of the MOF part at a stage, after the composite correction. W = mass at the stage (same scale as W_end)."""
    Wend = inp.tga_W_end
    kind = getattr(inp, "composite_kind", "none")
    frac = inp.mof_residue_fraction if getattr(inp, "mof_residue_fraction", None) else inp.framework_oxide_fraction()
    note = ""
    if kind == "before_dh" and stage == "start":
        return None, "a second material that burns away before the linker is still in the mass at this stage: use (i2) for the MOF and the composite weight percent below"
    if kind in ("with_dh", "after_dh") and getattr(inp, "composite_wt_pct", None):
        w = inp.composite_wt_pct / 100.0
        return W * (1 - w) / Wend, f"composite removed with its known weight fraction {100 * w:.1f} %"
    if frac < 1:
        return (W - Wend * (1 - frac)) / (Wend * frac), f"inert second material removed with its residue share ({100 * (1 - frac):.1f} % of the residue)"
    return W / Wend, note


def _tex_num(v, n=3):
    return f"{v:.{n}f}"


def method_i(inp):
    """Returns [i1, i2] MethodResult objects (status 'n/a' with a message when the data are missing)."""
    from .defects import REFS, MethodResult, open_sites, build_formula  # noqa: F401  (build_formula kept for symmetry)
    out = []
    nM = inp.nM(); Mres = inp.MMO2(); Q = inp.Q(); q_avail = Q / nM
    nd = inp.node_def()
    ml = inp.main_linker()
    ref = [("Abánades Lázaro, I. A multifaceted thermogravimetric analysis method for the exact determination of the composition of multifunctional "
            "metal-organic framework materials. Eur. J. Inorg. Chem. 2020, 4284-4294.", "10.1002/ejic.202000656")]
    sym = inp.M_symbol().replace("$_{", "").replace("}$", "").replace("{", "").replace("}", "")

    # ------------------------------------------------------------------ i1: the desolvated structure (about 200 degrees C)
    r1 = MethodResult("i1", "(i1) Lázaro 2020: experimental ratio structure : residue at ~200 °C + NMR + charge balance")
    r1.refs = ref
    r1.equations = [
        r"R_\mathrm{exp} = \frac{M_w(\mathrm{structure})}{M_w(\mathrm{residue})} = \frac{M_\mathrm{SBU} + x\sum_i r_iM_i + \sum n_jM_j + y\,M_\mathrm{bal}}{M_\mathrm{res} + x\sum_i r_if_iM_i}",
        r"q_\mathrm{avail} = x\sum_i r_iz_i + \sum n_jz_j + y\,z_\mathrm{bal}\qquad(\mathrm{charge\ balance,\ per\ metal})",
        r"x = \frac{R_\mathrm{exp}M_\mathrm{res} - M_\mathrm{SBU} - \sum n_jM_j - (q_\mathrm{avail}-\sum n_jz_j)\,M_\mathrm{bal}/z_\mathrm{bal}}{\sum_i r_iM_i - \left(\sum_i r_iz_i\right)M_\mathrm{bal}/z_\mathrm{bal} - R_\mathrm{exp}\sum_i r_if_iM_i}",
        r"L = n_M\,x,\qquad m = n_M\,x\,r_\mathrm{cap}\qquad(n_M:\ \mathrm{metals\ per\ node})",
    ]
    r1.symbols = {"R_exp": "mass of the desolvated structure divided by the residue mass (W_200 / W_end, composite-corrected)",
                  "M_SBU": "building unit per metal, e.g. ZrO4/6(OH)4/6 for UiO-66", "M_res": "residue per metal (ZrO2, HfO2, ZnO, TiO2 …)",
                  "r_i": "molar amount of component i relative to the main linker (NMR)", "z_i": "charge of component i", "f_i": "fraction of component i left in the residue (0 = burns away)",
                  "n_j": "fixed counter-ion per metal (Cl- from EDX/ICP)", "y": "balancing species per metal", "q_avail": "charge available per metal (Q / n_M)"}
    r1.assumptions = ["The structure at the chosen temperature is desolvated (pores empty) and the node is not yet dehydroxylated.",
                      "Every charge that the linkers, modulators and measured counter-ions do not balance is balanced by ONE species (default OH-; set in the inputs). "
                      "Unlike method (b) no coordinated water is added: the structure is the desolvated one.",
                      "The composition is exact for the chosen model: it does not depend on which decomposition step you assign to which species.",
                      "The nature of the balancing species must be checked with another technique (ICP, EDX, FT-IR): the charge left is reported."]
    r1.pitfalls = ["Residual solvent or water at the chosen temperature is counted as structure (R too high). Hold the temperature or choose the plateau of the TGA curve.",
                   "A modulator that is missing from the NMR table is counted as linker or as the balancing species.",
                   "A second material that burns away before the linker is in the mass at this stage: use (i2) and the composite weight percent."]
    # ------------------------------------------------------------------ i2: the dehydroxylated structure (DH MOF)
    r2 = MethodResult("i2", "(i2) Lázaro 2020: experimental ratio DH structure : residue (just before the linker burns) + O2- compensation")
    r2.refs = ref
    r2.equations = [
        r"R_\mathrm{exp}^\mathrm{DH} = \frac{M_w(\mathrm{DH\ structure})}{M_w(\mathrm{residue})} = \frac{M_w[\mathrm{ZrO(L)_xO_{1-x}}]}{M_w[\mathrm{ZrO_2}]}\quad(\mathrm{per\ metal,\ UiO\text{-}66})",
        r"x = \frac{R^\mathrm{DH}_\mathrm{exp}M_\mathrm{res} - M_\mathrm{SBU}^\mathrm{DH} - \dots}{\sum_i r_iM_i - \left(\sum_i r_iz_i\right)M_\mathrm{O}/2 - R^\mathrm{DH}_\mathrm{exp}\sum_i r_if_iM_i}\qquad\left(\mathrm{UiO\text{-}66:}\ x = \frac{R^\mathrm{DH}M_{ZrO_2} - M_{ZrO} - M_O}{M_L - M_O}\right)",
    ]
    r2.symbols = dict(r1.symbols)
    r2.symbols["M_SBU^DH"] = "dehydroxylated building unit per metal, e.g. ZrO for UiO-66"
    r2.assumptions = ["At the plateau the modulators have burned away (species flagged 'at plateau' stay), the node is dehydroxylated (M6O4(OH)4 → M6O6) and every charge a missing linker leaves is "
                      "closed by an oxo (O2-): ZrO(L)xO(1-x). Option 'no compensation': ZrO(L)x.",
                      "The oxygen that the linker's carboxylate gives to the residue is counted exactly (the reason why this method has no error where the older plateau methods have up to 21 %).",
                      "A chloride that survives to the plateau is not included here: use (i1) if Cl matters."]
    r2.pitfalls = ["Overlapping decomposition steps (a modulator that burns together with the linker) are no problem here: list that species as 'at plateau' with its NMR ratio.",
                   "The plateau must be the mass just before the linker starts to burn; read it at the DTG minimum between the steps."]
    kind = getattr(inp, "composite_kind", "none")
    Wend = inp.tga_W_end
    if not Wend:
        r1.message = r2.message = "needs the TGA residue mass (W_end)"
        return [r1, r2]
    bal1 = _balance(getattr(inp, "lazaro_balance", "OH") or "OH")
    L0 = inp.L0()

    def _node_result(r, stage, R, Msbu, bal, linkers, mods, fixed):
        sol = solve_per_metal(R, Mres, Msbu, q_avail, linkers, mods, fixed, bal)
        x = sol["x"]
        L = nM * x
        am = {k: nM * v for k, v in sol["amounts"].items() if k != ml.name}
        sp = {s.name: s for s in inp.species}
        m = sum(v for k, v in am.items() if sp.get(k) is not None and sp[k].role == "cap")
        clv = nM * sum(n for p, n in fixed if p.name.startswith("Cl"))
        if bal is not None and bal.name.startswith("Cl"):
            clv += nM * sol["y"]
        r.L, r.x, r.m, r.Cl = L, L0 - L, m, clv
        oth = {k: v for k, v in am.items() if sp.get(k) is not None and sp[k].z > 0}
        r.OS = open_sites(inp, L, m, clv, sum(sp[k].z * v for k, v in oth.items() if k in sp) if oth else m)
        r.status = "ok" if (abs(sol["charge_left"]) < 0.02 or bal is None) and x > 0 and sol["y"] > -0.02 else "warn"
        parts_txt = f"({ml.name.split(' (')[0]}){L:.2f}" + "".join(f"({k.split(' (')[0]}){v:.2f}" for k, v in am.items() if v > 1e-4)
        if bal is not None and sol["y"] * nM > 1e-3:
            parts_txt += f"({bal.name.rstrip('-').rstrip('2')}){sol['y'] * nM:.2f}"
        r.formula = nd["core"].format(M=sym) + parts_txt + (" (DH)" if stage == "dh" else "")
        r.M_formula = nM * sol["mass"]
        r.extra = {"R_exp": R, "R_calc": sol["R_calc"], "charge_left_per_node": nM * sol["charge_left"], "y_balance_per_node": nM * sol["y"],
                   "M_structure_per_metal": sol["mass"], "x_per_metal": x}
        if sol["y"] < -0.02 and bal is not None:
            r.message = f"negative balancing species ({nM * sol['y']:.2f} per node): the composition is not consistent with this model (check the NMR table and the stage)"
        elif abs(sol["charge_left"]) > 0.02 and bal is None:
            r.message = f"charge left per node: {nM * sol['charge_left']:+.2f} (not balanced: add the balancing species or measure it)"
        return sol

    # --- i1
    if not inp.tga_W_200:
        r1.message = "needs the TGA mass at ~200 °C (desolvated structure)"
    else:
        R1, note1 = _stage_R(inp, inp.tga_W_200, "start")
        lk, md, msg = _collect(inp, "start")
        if R1 is None:
            r1.message = note1
        elif lk is None:
            r1.message = msg
        else:
            fixed = []
            cl_node = None
            if inp.cl_per_node is not None:
                cl_node = inp.cl_per_node
            elif inp.cl_wt_pct is not None:
                cl_node = 0.0
                for _ in range(4):                                  # Cl wt% needs the structure mass: iterate
                    fixed_try = [(Part("Cl-", M_CL, 1.0), cl_node / nM)]
                    s_ = solve_per_metal(R1, Mres, inp.M_clust() / nM, q_avail, lk, md, fixed_try, bal1)
                    cl_node = inp.cl_wt_pct / 100 * nM * s_["mass"] / M_CL
            if cl_node:
                fixed = [(Part("Cl-", M_CL, 1.0), cl_node / nM)]
            sol = _node_result(r1, "start", R1, inp.M_clust() / nM, bal1, lk, md, fixed)
            r1.steps = [
                f"R_\\mathrm{{exp}} = \\frac{{W_{{200}}}}{{W_\\mathrm{{end}}}} = \\frac{{{inp.tga_W_200:.2f}}}{{{Wend:.2f}}}"
                + (f"\\ \\to\\ {R1:.4f}\\ \\ (\\mathrm{{{note1.split(' (')[0].replace('%', chr(92) + '%')}}})" if note1 else f" = {R1:.4f}"),
                f"M_\\mathrm{{SBU}} = M_\\mathrm{{node}}/{nM} = {inp.M_clust() / nM:.2f},\\quad M_\\mathrm{{res}} = {Mres:.2f},\\quad q_\\mathrm{{avail}} = {Q:g}/{nM} = {q_avail:.3f}",
                f"\\sum r_iM_i = {sol['A']:.2f},\\quad \\sum r_iz_i = {sol['C']:.3f}" + (f",\\quad \\mathrm{{Cl}} = {cl_node:.2f}\\ \\mathrm{{per\\ node}}" if cl_node else ""),
                f"x = {sol['x']:.4f}\\ \\mathrm{{per\\ metal}}\\ \\Rightarrow\\ L = {nM}\\times{sol['x']:.4f} = {r1.L:.3f}\\ \\mathrm{{per\\ node}},\\quad x_\\mathrm{{missing}} = {L0:g} - {r1.L:.3f} = {r1.x:.3f}",
                f"\\mathrm{{check:}}\\ R_\\mathrm{{calc}} = {sol['R_calc']:.4f}\\ (\\mathrm{{must\\ equal\\ }}R_\\mathrm{{exp}}),\\quad \\mathrm{{charge\\ left}} = {nM * sol['charge_left']:+.3f}\\ \\mathrm{{per\\ node}}"]

    # --- i2
    if not inp.tga_W_plateau:
        r2.message = "needs the TGA mass at the dehydroxylated plateau"
    else:
        R2, note2 = _stage_R(inp, inp.tga_W_plateau, "dh")
        lk, md, msg = _collect(inp, "dh")
        if lk is None:
            r2.message = msg
        else:
            if nd.get("plateau_dehydroxylated", True):
                Msbu_dh = (nM * inp.MM() + (nd.get("core_O", 0) + nd.get("core_OH", 0) / 2) * M_O) / nM
            else:
                Msbu_dh = inp.M_core() / nM
            bal2 = None if getattr(inp, "lazaro_dh_compensation", "O") == "none" else _balance("O")
            sol = _node_result(r2, "dh", R2, Msbu_dh, bal2, lk, md, [])
            r2.steps = [
                f"R_\\mathrm{{exp}}^\\mathrm{{DH}} = \\frac{{W_\\mathrm{{plat}}}}{{W_\\mathrm{{end}}}} = \\frac{{{inp.tga_W_plateau:.2f}}}{{{Wend:.2f}}}"
                + (f"\\ \\to\\ {R2:.4f}" if note2 else f" = {R2:.4f}"),
                f"M_\\mathrm{{SBU}}^\\mathrm{{DH}} = {Msbu_dh:.2f},\\quad M_\\mathrm{{res}} = {Mres:.2f},\\quad q_\\mathrm{{avail}} = {q_avail:.3f}",
                f"x = {sol['x']:.4f}\\ \\mathrm{{per\\ metal}}\\ \\Rightarrow\\ L = {nM}\\times{sol['x']:.4f} = {r2.L:.3f}\\ \\mathrm{{per\\ node}},\\quad x_\\mathrm{{missing}} = {L0:g} - {r2.L:.3f} = {r2.x:.3f}",
                f"\\mathrm{{check:}}\\ R_\\mathrm{{calc}} = {sol['R_calc']:.4f},\\quad \\mathrm{{O^{{2-}}\\ compensation}} = {nM * sol['y']:.3f}\\ \\mathrm{{per\\ node}}"]
            # a second material that burned away before the linker: its weight percent from the start ratio
            if kind == "before_dh" and inp.tga_W_200:
                lk1, md1, m1 = _collect(inp, "start")
                if lk1 is not None:
                    A1 = sum(p.r * p.M for p in lk1 + md1); C1 = sum(p.r * p.z for p in lk1 + md1)
                    y1 = (q_avail - sol["x"] * C1) / bal1.z if bal1 is not None else 0.0
                    Mstruct = inp.M_clust() / nM + sol["x"] * A1 + (y1 * bal1.M if bal1 is not None else 0.0)
                    w = composite_weight_fraction(Mstruct, Mres, inp.tga_W_200 / Wend)
                    r2.extra["composite_wt_pct"] = 100 * w
                    r2.steps.append(f"w_\\mathrm{{composite}} = 1 - \\frac{{M_\\mathrm{{MOF}}}}{{R_\\mathrm{{exp}}M_\\mathrm{{res}}}} = 1 - \\frac{{{Mstruct:.2f}}}{{{inp.tga_W_200 / Wend:.4f}\\times{Mres:.2f}}} = {w:.4f}\\ \\Rightarrow\\ {100 * w:.1f}\\ \\mathrm{{wt\\%}}")
                    r2.message = (r2.message + " · " if r2.message else "") + f"composite that burns before the linker: {100 * w:.1f} wt % of the desolvated material"
    out += [r1, r2]
    return out
