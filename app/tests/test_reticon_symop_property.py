"""Property test for the CIF symmetry parser (_symop): fixed random seed, so the run is repeatable.

1. Random garbage strings: the parser either returns three numbers or raises ValueError. No other exception.
2. Random valid operations: the result equals the exact value, computed with fractions.Fraction.
Both copies of generators.py are checked.
"""
import importlib
import random
from fractions import Fraction

import pytest

MODULES = ["mfs.reticon.generators", "mfs.icons_runtime.reticon.generators"]
SEED = 20261127
POINT = (Fraction(1, 10), Fraction(2, 10), Fraction(3, 10))
ALPHABET = "xyz0123456789/+-.()*_ ,'\"abcdefg;:=!<>"


def _garbage(rng: random.Random) -> str:
    return "".join(rng.choice(ALPHABET) for _ in range(rng.randint(0, 14)))


def _term(rng: random.Random):
    """A random valid term and its exact value at POINT."""
    kind = rng.choice(["var", "num", "frac", "varfrac"])
    if kind == "var":
        name = rng.choice("xyz")
        return name, POINT["xyz".index(name)]
    if kind == "num":
        n = rng.randint(0, 9)
        return str(n), Fraction(n)
    if kind == "frac":
        a, b = rng.randint(0, 9), rng.randint(1, 9)
        return f"{a}/{b}", Fraction(a, b)
    name = rng.choice("xyz")
    b = rng.randint(1, 9)
    return f"{name}/{b}", POINT["xyz".index(name)] / b


def _component(rng: random.Random):
    text, value = "", Fraction(0)
    sign = rng.choice(["", "-"])
    t, v = _term(rng)
    text, value = sign + t, (-v if sign else v)
    for _ in range(rng.randint(0, 2)):
        op = rng.choice(["+", "-"])
        t, v = _term(rng)
        text += op + t
        value = value + v if op == "+" else value - v
    return text, value


@pytest.mark.parametrize("modname", MODULES)
def test_garbage_gives_numbers_or_value_error_only(modname):
    G = importlib.import_module(modname)
    rng = random.Random(SEED)
    for _ in range(2000):
        s = _garbage(rng)
        try:
            f = G._symop(s + ",y,z")
        except ValueError:
            continue
        out = f(0.1, 0.2, 0.3)
        assert len(out) == 3 and all(isinstance(v, float) for v in out), s


@pytest.mark.parametrize("modname", MODULES)
def test_valid_operations_match_exact_values(modname):
    G = importlib.import_module(modname)
    rng = random.Random(SEED + 1)
    for _ in range(500):
        parts = [_component(rng) for _ in range(3)]
        op = ",".join(t for t, _ in parts)
        got = G._symop(op)(*(float(v) for v in POINT))
        for g, (_, want) in zip(got, parts):
            assert abs(g - float(want)) < 1e-9, (op, got)
