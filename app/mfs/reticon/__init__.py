__version__ = "4.2.0"
from .core import Icon, Style, Palette, render_svg
from .library import REGISTRY
from . import library2  # noqa: F401  (registers extended set)
from . import library3  # noqa: F401  (original redraws)
from . import library4  # noqa: F401  (v4 content)
from . import library5  # noqa: F401  (v5 batch 1)
from . import library6  # noqa: F401  (v5 batch 2)
from . import library7  # noqa: F401  (v5 batch 3)
from . import library8  # noqa: F401  (v5 batch 4: generated sets)
from . import overrides as _ov
_ov.register()

__all__ = ["Icon", "Style", "Palette", "render_svg", "REGISTRY", "get_icon", "all_icons"]


from . import anchors as _anchors


def get_icon(name: str) -> Icon:
    return _anchors.apply(REGISTRY[name]())


def all_icons() -> list[Icon]:
    return [_anchors.apply(b()) for b in REGISTRY.values()]
