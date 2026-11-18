"""
Hand-made shape variants (real alternative drawings) for the most-used icons.
Registered into variants.OVERRIDES; each function gets the base Icon and returns a new Icon.
"""
from __future__ import annotations

import copy

from .core import Icon, Op, polar
from .library2 import rot_rect
from . import variants as V

_REG: dict[str, list] = {}


def over(name: str, key: str, label: str):
    def deco(fn):
        def build(base: Icon) -> Icon:
            ic = Icon(base.name, f"{base.label} ({label})", base.category, base.tags)
            fn(ic, base)
            return ic
        _REG.setdefault(name, []).append((key, build))
        return fn
    return deco


def keep(ic: Icon, base: Icon, drop=()):
    """Copy the base drawing into ic, skipping op indexes in `drop`."""
    ic.ops += [copy.deepcopy(o) for i, o in enumerate(base.ops) if i not in drop]
    ic.anchors = copy.deepcopy(base.anchors)


def moved(base: Icon, tf: str, drop=()):
    return [V.with_transform(o, tf) for i, o in enumerate(base.ops) if i not in drop]


# ----------------------------------------------------------------- glassware
@over("round_bottom_flask", "two-neck", "two-neck")
def _(ic, b):
    keep(ic, b)
    ic.spath("M40.6 25.6 L46.6 14.6", 3)
    ic.spath("M45.4 28.6 L51.4 17.6", 3)
    ic.line(44.8, 13.6, 53.2, 18.6, 3)
    ic.anchor("side_neck", 49, 15, "up")


@over("round_bottom_flask", "three-neck", "three-neck")
def _(ic, b):
    keep(ic, b)
    for s in (1, -1):
        x = lambda v: 32 + s * (v - 32)
        ic.spath(f"M{x(40.6)} 25.6 L{x(46.6)} 14.6", 3)
        ic.spath(f"M{x(45.4)} 28.6 L{x(51.4)} 17.6", 3)
        ic.line(x(44.8), 13.6, x(53.2), 18.6, 3)
    ic.anchor("left_neck", 15, 15, "up")
    ic.anchor("right_neck", 49, 15, "up")


@over("erlenmeyer_flask", "stoppered", "stoppered")
def _(ic, b):
    keep(ic, b)
    ic.rect(25, 1, 14, 9, 2.5)


@over("erlenmeyer_flask", "stirring", "with stir bar")
def _(ic, b):
    keep(ic, b, drop=(3, 4))
    ic.cut().spath("M22 44 Q32 49 42 44", 1.6)
    ic.acc().rect(26.5, 49.5, 11, 3, 1.5)


@over("beaker", "stirring", "with stir bar")
def _(ic, b):
    keep(ic, b)
    ic.cut().spath("M23.5 34 Q32 41 40.5 34", 1.6)
    ic.acc().rect(26, 49, 12, 3, 1.5)


@over("beaker", "covered", "covered (watch glass)")
def _(ic, b):
    keep(ic, b)
    ic.spath("M7 8.5 Q30 1 53 8.5", 2.6, glass=True)


@over("test_tube", "stoppered", "stoppered")
def _(ic, b):
    keep(ic, b, drop=(1,))
    ic.rect(21.5, 2, 21, 9, 2.5)


@over("test_tube", "boiling", "boiling")
def _(ic, b):
    keep(ic, b)
    for x, y, r in ((31, 44, 1.4), (33.4, 39, 1.1), (30.6, 34.5, 1.2), (33, 31.5, 0.9)):
        ic.acc().circle(x, y, r)
    for x in (16, 48):
        ic.acc().spath(f"M{x} 34 q-2.5 -3 0 -6 q2.5 -3 0 -6", 2.2)


@over("vial", "open", "open, cap beside")
def _(ic, b):
    ic.ops += moved(b, "translate(-8 0)", drop=(0, 1, 2, 3, 4))
    rot_rect(ic, 50, 47, 18, 12, 3, 90)
    for y in (42, 47, 52):
        ic.cut().line(45.5, y, 54.5, y, 1.8)


@over("vial", "septum", "septum + needle")
def _(ic, b):
    keep(ic, b)
    ic.cut().rect(27, 7, 10, 3, 1)
    ic.acc().line(32, 0.5, 32, 30, 1.6)
    ic.rect(28.5, 0, 7, 3, 1)


@over("volumetric_flask", "open", "unstoppered")
def _(ic, b):
    keep(ic, b, drop=(0,))
    ic.line(26.5, 6, 37.5, 6, 2.6)


@over("volumetric_flask", "filled", "filled to the mark")
def _(ic, b):
    keep(ic, b, drop=(3,))
    ic.sec().circle(32, 43.3, 10.6)
    ic.sec().rect(30.6, 15, 2.8, 22, 0.8)


@over("measuring_cylinder", "full", "full")
def _(ic, b):
    keep(ic, b, drop=(1,))
    ic.ops.insert(1, Op("s", "rect", {"x": 26.5, "y": 10, "width": 11, "height": 41.5, "rx": 1}))


@over("measuring_cylinder", "stoppered", "stoppered")
def _(ic, b):
    keep(ic, b)
    ic.rect(25, 0.5, 14, 6.5, 2)


@over("separatory_funnel", "draining", "stopcock open, draining")
def _(ic, b):
    keep(ic, b)
    ic.ops[-2].attrs["transform"] = "rotate(90 32 52)"
    ic.sec().circle(32, 62.6, 1.3)


@over("separatory_funnel", "single-phase", "one phase")
def _(ic, b):
    keep(ic, b, drop=(2, 3))
    ic.ops.insert(2, Op("s", "path", {"d": "M19.6 24 H44.4 C45 30 43 35.5 38.6 39.6 C36.5 41.5 34.5 42.8 32 44 "
                                            "C25.5 40 21.2 35.5 20 30.6 C19.4 28 19.3 26 19.6 24 Z"}))


@over("falcon_tube", "pellet", "pellet after centrifugation")
def _(ic, b):
    keep(ic, b)
    ic.acc().path("M29.6 52 H34.4 L32.8 55 H31.2 Z")


@over("falcon_tube", "open", "open, cap beside")
def _(ic, b):
    ic.ops += moved(b, "translate(-8 0)", drop=(0,))
    rot_rect(ic, 50, 50, 22, 12, 3, 90)


@over("microtube", "closed", "closed")
def _(ic, b):
    keep(ic, b, drop=(2,))
    ic.rect(20, 4, 24, 5, 2)
    ic.rect(41, 6, 6, 4, 1.5)


@over("microtube", "pellet", "pellet")
def _(ic, b):
    keep(ic, b)
    ic.acc().path("M30 51 H34 L32.9 54.6 Q32 56 31.1 54.6 Z")


@over("petri_dish", "side", "side view with lid")
def _(ic, b):
    ic.srect(6, 36, 52, 14, 2, 2.6, glass=True)
    ic.sec().rect(9, 41, 46, 7, 1)
    for x in (16, 27, 38, 47):
        ic.acc().ellipse(x, 40.5, 3, 1.4)
    ic.srect(4, 30, 56, 9, 2, 2.6, glass=True)


@over("petri_dish", "crystals", "crystallisation dish")
def _(ic, b):
    keep(ic, b, drop=(2, 3, 4, 5, 6))
    for (x, y, s) in ((24, 26, 4), (38, 25, 3), (36, 39, 4.5), (24, 40, 3), (44, 33, 2.5)):
        ic.acc().poly([(x, y - s), (x + s * .87, y - s / 2), (x + s * .87, y + s / 2), (x, y + s),
                       (x - s * .87, y + s / 2), (x - s * .87, y - s / 2)])


@over("cuvette", "capped", "capped")
def _(ic, b):
    keep(ic, b)
    ic.rect(21.5, 1, 21, 7, 2)


@over("cuvette", "in-beam", "in the light beam")
def _(ic, b):
    ic.acc().line(2, 38, 22, 38, 3.4)
    ic.ops += [Op("x", "rect", {"x": 21, "y": 4, "width": 22, "height": 56, "rx": 3})]
    keep(ic, b)
    ic.acc().line(42, 38, 61, 38, 2)


@over("reagent_bottle", "clear", "clear glass")
def _(ic, b):
    ic.rect(25, 3, 14, 9, 2)
    ic.spath("M26 12 H38 V17 Q50 20 50 31 V54 Q50 59 45 59 H19 Q14 59 14 54 V31 Q14 20 26 17 Z", 3, glass=True)
    ic.sec().rect(17.5, 36, 29, 19.5, 3)
    ic.cut().rect(20, 24, 24, 9, 1.5)


@over("reagent_bottle", "dropper", "dropper bottle")
def _(ic, b):
    ic.acc().ellipse(32, 8, 5, 6)
    ic.rect(27, 13, 10, 6, 1.5)
    ic.sec().path("M27 19 H37 V23 Q46 25 46 34 V54 Q46 59 41 59 H23 Q18 59 18 54 V34 Q18 25 27 23 Z")
    ic.cut().rect(23, 35, 18, 12, 2)


# ----------------------------------------------------------------- synthesis equipment
@over("autoclave", "open", "opened, liner out")
def _(ic, b):
    ic.ops += moved(b, "translate(-9 4) scale(0.86)", drop=(0, 1, 2, 3))
    rot_rect(ic, 20, 8, 40, 9, 2.5, -14)
    ic.rect(42, 22, 18, 34, 3.5)
    ic.rect(40.5, 18, 21, 5, 1.5)
    ic.sec().rect(45, 31, 12, 22, 2)


@over("autoclave", "liner", "Teflon liner only")
def _(ic, b):
    ic.rect(19, 15, 26, 44, 5)
    ic.rect(16.5, 8, 31, 8, 2.5)
    ic.cut().line(19, 19, 45, 19, 1.4)
    ic.cut().rect(23, 25, 18, 30, 3)
    ic.sec().rect(25.5, 34, 13, 18.5, 2)


@over("oven", "heating", "heating")
def _(ic, b):
    keep(ic, b)
    for x in (20, 28, 36):
        ic.acc().spath(f"M{x} 26 q-2.5 -2.5 0 -5 q2.5 -2.5 0 -5", 2)


@over("oven", "vacuum", "vacuum oven")
def _(ic, b):
    keep(ic, b, drop=(1, 2, 3, 4))
    ic.ops.insert(1, Op("x", "circle", {"cx": 26.5, "cy": 30.5, "r": 14}))
    ic.ops.insert(2, Op("s", "rect", {"x": 15, "y": 34, "width": 23, "height": 3.5, "rx": 1}))
    ic.acc().circle(50.5, 45, 4.5)
    ic.cut().line(50.5, 45, 53, 42.5, 1.4)


@over("centrifuge", "open", "lid open")
def _(ic, b):
    ic.poly([(12, 25), (17, 3), (47, 3), (52, 25)])
    ic.cut().poly([(17, 22), (21, 7), (43, 7), (47, 22)])
    ic.rect(7, 27, 50, 29, 5)
    ic.ellipse(32, 27, 22, 6)
    ic.cut().ellipse(32, 27, 18, 4.4)
    ic.sec().ellipse(32, 27, 11, 2.7)
    ic.cut().rect(14, 37, 16, 9, 2)
    ic.cut().circle(44, 42, 5)


@over("centrifuge", "micro", "microcentrifuge")
def _(ic, b):
    ic.rect(10, 34, 44, 22, 8)
    ic.path("M12 35 Q32 12 52 35 Z")
    ic.cut().circle(32, 28, 6.5)
    ic.sec().circle(32, 28, 3.2)
    ic.cut().rect(18, 42, 14, 7, 1.5)
    ic.cut().circle(42, 45.5, 3.4)


@over("glovebox", "in-use", "hands in the gloves")
def _(ic, b):
    keep(ic, b)
    for s in (1, -1):
        x = lambda v: 32 + s * (v - 32)
        ic.sec().path(f"M{x(18.5)} 38 C{x(18)} 32 {x(19)} 27 {x(21)} 23 L{x(25)} 23.5 C{x(25.5)} 28 {x(24.5)} 33 {x(24)} 38 Z")


@over("glovebox", "antechamber", "with antechamber")
def _(ic, b):
    ic.rect(2, 10, 48, 39, 4)
    ic.cut().rect(7, 15, 38, 15, 2.5)
    ic.cut().circle(17, 40, 5.5)
    ic.cut().circle(36, 40, 5.5)
    ic.sec().circle(17, 40, 3)
    ic.sec().circle(36, 40, 3)
    ic.rect(50, 22, 13, 18, 5)
    ic.cut().circle(56.5, 31, 4.5)
    ic.circle(56.5, 31, 2)
    ic.rect(6, 49, 5, 9, 1.5)
    ic.rect(41, 49, 5, 9, 1.5)


@over("hotplate_stirrer", "flask", "with round-bottom flask")
def _(ic, b):
    ic.spath("M29 2 V10.4 A12 12 0 1 0 35 10.4 V2", 3, glass=True)
    ic.sec().path("M23.1 23 A9.2 9.2 0 0 0 40.9 23 Z")
    ic.acc().rect(28.5, 27.5, 7, 2.2, 1.1)
    keep(ic, b, drop=(0, 1, 2))


@over("hotplate_stirrer", "heating", "heating")
def _(ic, b):
    keep(ic, b)
    ic.ops[3] = Op("a", ic.ops[3].tag, ic.ops[3].attrs)
    for x in (11, 53):
        ic.acc().spath(f"M{x} 30 q-2.5 -3 0 -6 q2.5 -3 0 -6", 2.2)


@over("analytical_balance", "draft-shield", "with draft shield")
def _(ic, b):
    keep(ic, b)
    ic.srect(11, 6, 42, 37, 2, 2.2, glass=True)
    ic.line(30, 6, 30, 43, 1.4)


@over("analytical_balance", "powder", "weighing powder")
def _(ic, b):
    keep(ic, b)
    ic.acc().path("M26 27 Q32 19 38 27 Z")


@over("sonicator", "probe", "probe sonicator")
def _(ic, b):
    ic.rect(25, 2, 14, 11, 2.5)
    ic.path("M27.5 13 H36.5 L34.5 25 H29.5 Z")
    ic.spath("M16 24 V52 Q16 58 22 58 H42 Q48 58 48 52 V24", 3, glass=True)
    ic.sec().rect(19.5, 34, 25, 20.5, 2)
    ic.rect(30.5, 25, 3, 22, 1)
    for x, y in ((28, 48), (36, 49), (32, 51.5), (26.5, 52), (37.5, 52.5)):
        ic.acc().circle(x, y, 1.1)


@over("sonicator", "flask", "bath with flask")
def _(ic, b):
    keep(ic, b, drop=(4, 5))
    ic.cut().circle(32, 32, 13)
    ic.cut().rect(27, 4, 10, 22, 2)
    ic.spath("M29.5 6 V22.3 A10 10 0 1 0 34.5 22.3 V6", 3, glass=True)
    ic.sec().path("M25 33 A7.5 7.5 0 0 0 39 33 Z")


@over("tube_furnace", "split", "split (opened)")
def _(ic, b):
    ic.rect(11, 6, 42, 14, 4)
    ic.rect(11, 37, 42, 14, 4)
    for x in (19, 25, 31, 37, 43):
        ic.acc().line(x, 16, x, 18.5, 1.6)
        ic.acc().line(x, 38.5, x, 41, 1.6)
    ic.sec().rect(2.5, 26, 59, 6, 3)
    ic.rect(15, 51, 5, 8, 1.5)
    ic.rect(44, 51, 5, 8, 1.5)
    ic.line(13, 20, 13, 37, 2)


@over("tube_furnace", "gas-flow", "with gas flow")
def _(ic, b):
    keep(ic, b, drop=(10, 11, 12))
    for x0 in (1, 50):
        ic.acc().line(x0, 20, x0 + 10, 20, 2.4)
        ic.acc().poly([(x0 + 10, 17), (x0 + 13.5, 20), (x0 + 10, 23)])


@over("micropipette", "dispensing", "dispensing a drop")
def _(ic, b):
    ic.ops += moved(b, "translate(0 -2.5)")
    ic.sec().path("M32 58.5 C34 61 34.2 63.5 32 63.5 C29.8 63.5 30 61 32 58.5 Z")


@over("micropipette", "eject", "ejecting the tip")
def _(ic, b):
    keep(ic, b, drop=(5,))
    ic.sec().path("M29.2 48.5 H34.8 L32.8 61 H31.2 Z")
    ic.ops[-1].attrs["transform"] = "translate(5 1) rotate(22 32 52)"
    ic.acc().spath("M42 46 l4 2 M42 51 l4 0", 1.8)


@over("gas_cylinder", "regulator", "with regulator")
def _(ic, b):
    keep(ic, b)
    ic.rect(35, 5, 12, 6, 2)
    ic.acc().circle(46, 3.6, 3.4)
    ic.acc().circle(52, 9, 3.4)
    ic.line(47, 10, 60, 10, 2)


@over("gas_cylinder", "lecture", "small lecture bottle")
def _(ic, b):
    ic.rect(25, 18, 14, 4, 2)
    ic.rect(29, 22, 6, 6)
    ic.rect(23, 28, 18, 32, 7)
    ic.cut().rect(23, 38, 18, 3)
    ic.acc().circle(32, 24.5, 1.4)


# ----------------------------------------------------------------- materials & proteins
@over("nanoparticle", "bare", "bare (no ligands)")
def _(ic, b):
    keep(ic, b, drop=tuple(range(24)))


@over("nanoparticle", "aggregate", "aggregated")
def _(ic, b):
    for cx, cy, r in ((22, 24, 12), (40, 22, 10), (30, 41, 12), (47, 40, 9)):
        ic.sec().circle(cx, cy, r)
        ic.cut().circle(cx, cy, r * .78)
        ic.circle(cx, cy, r * .66)
        ic.cut().circle(cx - r * .3, cy - r * .3, r * .16)


@over("mof_framework", "guests", "with guest molecules")
def _(ic, b):
    keep(ic, b)
    for x, y, r in ((30, 33, 3.4), (35.5, 36, 2.4), (27, 38.5, 2.2)):
        ic.acc().circle(x, y, r)


@over("mof_framework", "interpenetrated", "interpenetrated")
def _(ic, b):
    keep(ic, b)
    f = [(22, 33), (44, 33), (44, 55), (22, 55)]
    bk = [(x + 10, y - 10) for x, y in f]
    for i in range(4):
        for p, q in ((f[i], f[(i + 1) % 4]), (bk[i], bk[(i + 1) % 4]), (f[i], bk[i])):
            ic.acc().line(*p, *q, 2)
    for p in f + bk:
        ic.acc().circle(*p, 3)


@over("octahedral_sbu", "linkers", "with linkers")
def _(ic, b):
    pts = [(32, 8), (32, 56), (10, 34), (26, 40), (54, 30), (38, 24)]
    for (x, y) in pts:
        dx, dy = x - 32, y - 32
        n = max(1e-6, (dx * dx + dy * dy) ** .5)
        ic.sec().line(x, y, x + dx / n * 8, y + dy / n * 8, 2.6)
    ic.ops += moved(b, "translate(4 4) scale(0.875)")


@over("octahedral_sbu", "missing-linker", "missing linker, open site")
def _(ic, b):
    pts = [(32, 8), (32, 56), (10, 34), (26, 40), (38, 24)]
    for (x, y) in pts:
        dx, dy = x - 32, y - 32
        n = max(1e-6, (dx * dx + dy * dy) ** .5)
        ic.sec().line(x, y, x + dx / n * 8, y + dy / n * 8, 2.6)
    ic.ops += moved(b, "translate(4 4) scale(0.875)")
    ic.acc().scircle(51.25, 30.25, 5.5, 1.8)


@over("protein_ribbon", "substrate", "with bound substrate")
def _(ic, b):
    keep(ic, b, drop=(1,))
    ic.acc().poly([polar(44, 46, 6, 30 + 60 * k) for k in range(6)])
    ic.cut().circle(44, 46, 2.2)


@over("protein_ribbon", "sheet", "β-sheet")
def _(ic, b):
    for y, d in ((14, 1), (30, -1), (46, 1)):
        x0, x1 = (8, 50) if d == 1 else (56, 14)
        ic.sec().poly([(x0, y - 3.5), (x1, y - 3.5), (x1, y - 7), (x1 + d * 7, y), (x1, y + 7), (x1, y + 3.5), (x0, y + 3.5)])
    ic.spath("M56 17.5 C62 20 62 27 56 26.5", 2.2)
    ic.spath("M8 33.5 C2 36 2 43 8 42.5", 2.2)


# ----------------------------------------------------------------- instruments & lines
@over("well_plate_96", "lid", "with lid")
def _(ic, b):
    keep(ic, b)
    ic.srect(2.5, 10.5, 59, 43, 5, 2.2, glass=True)


@over("well_plate_96", "384", "384-well")
def _(ic, b):
    ic.rect(4, 13, 56, 38, 4)
    for r in range(8):
        for c in range(14):
            ic.cut().circle(9.5 + c * 3.85, 18 + r * 3.9, 1.25)


@over("schlenk_line", "cold-trap", "with cold trap")
def _(ic, b):
    keep(ic, b)
    ic.line(56, 27, 56, 33, 2.6)
    ic.rect(50, 33, 12, 25, 4)
    ic.cut().rect(53, 36, 6, 19, 2)
    ic.sec().rect(54.5, 42, 3, 12, 1)


@over("schlenk_line", "two-flasks", "two flasks")
def _(ic, b):
    keep(ic, b)
    ic.spath("M11 42 V46 A7 7 0 1 0 17 46 V42", 2.6, glass=True)
    ic.sec().path("M9.2 51 A5 5 0 0 0 18.8 51 Z")


@over("electron_microscope", "tem", "TEM (tall column)")
def _(ic, b):
    ic.rect(25, 1, 14, 7, 2.5)
    ic.rect(27, 8, 10, 40)
    for y in (12, 20, 28, 36):
        ic.rect(21, y, 22, 5, 2)
    ic.cut().line(32, 4, 32, 46, 2)
    ic.rect(10, 46, 44, 13, 3)
    ic.cut().ellipse(32, 52.5, 12, 3.6)
    ic.acc().ellipse(32, 52.5, 8, 2)


@over("electron_microscope", "beam", "beam on the sample")
def _(ic, b):
    ic.acc().poly([(26, 2), (38, 2), (33, 34), (31, 34)])
    ic.rect(8, 36, 48, 8, 2)
    ic.sec().rect(8, 44, 48, 14, 2)
    for a in (-60, -30, 210, 240):
        x2, y2 = polar(32, 35, 20, a)
        ic.sec().spath(f"M32 34 L{x2:.1f} {y2:.1f}", 1.8)
    ic.acc().circle(32, 35, 2.4)


def register():
    for name, pairs in _REG.items():
        assert len(pairs) == 2, (name, len(pairs))
        V.OVERRIDES[name] = (pairs[0], pairs[1])
    return len(_REG)
