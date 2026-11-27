"""CIF symmetry operations are parsed, never evaluated: good strings give the right coordinates,
hostile or malformed strings raise ValueError. Both copies of generators.py are checked.
"""
import importlib
import pytest

MODULES = ["mfs.reticon.generators", "mfs.icons_runtime.reticon.generators"]
POINT = (0.1, 0.2, 0.3)


@pytest.mark.parametrize("modname", MODULES)
@pytest.mark.parametrize("op,expected", [
    ("x,y,z", (0.1, 0.2, 0.3)),
    ("-x+1/2,y+1/2,-z", (0.4, 0.7, -0.3)),
    ("1/2+x,1/2-y,z+1/2", (0.6, 0.3, 0.8)),
    ("x/2,y,z", (0.05, 0.2, 0.3)),
])
def test_good_operations_give_the_right_coordinates(modname, op, expected):
    G = importlib.import_module(modname)
    got = G._symop(op)(*POINT)
    assert all(abs(a - b) < 1e-9 for a, b in zip(got, expected)), (op, got, expected)


HOSTILE = [
    "__import__('os').system('echo x'),y,z",
    "().__class__,y,z",
    "open('/etc/passwd'),y,z",
    "x**2,y,z",
    "x*2,y,z",
    "2x,y,z",
    "x.real,y,z",
    "1e5,y,z",
    "lambda: 1,y,z",
    "x;y,y,z",
    "x,y",            # only two components
    "x,y,z,x",        # four components
    "x,y,z/0",        # division by zero
    ",y,z",           # empty component
]


@pytest.mark.parametrize("modname", MODULES)
@pytest.mark.parametrize("op", HOSTILE)
def test_hostile_or_malformed_strings_raise(modname, op):
    G = importlib.import_module(modname)
    with pytest.raises(ValueError):
        G._symop(op)


@pytest.mark.parametrize("modname", MODULES)
def test_no_eval_in_generators(modname):
    src = importlib.import_module(modname).__file__
    text = open(src, encoding="utf-8").read()
    assert "eval(" not in text
