"""A sans-serif TrueType font for the example-image scripts on any OS (Linux: Liberation/DejaVu, Windows: Arial, macOS: Helvetica)."""
from PIL import ImageFont

_REGULAR = ("LiberationSans-Regular.ttf", "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", "arial.ttf", "Arial.ttf",
            "/System/Library/Fonts/Helvetica.ttc", "/System/Library/Fonts/Supplemental/Arial.ttf", "DejaVuSans.ttf")
_BOLD = ("LiberationSans-Bold.ttf", "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", "arialbd.ttf", "Arial Bold.ttf",
         "/System/Library/Fonts/Supplemental/Arial Bold.ttf", "DejaVuSans-Bold.ttf")


def sans(size, bold=False):
    for name in (_BOLD if bold else ()) + _REGULAR:
        try:
            return ImageFont.truetype(name, size)
        except Exception:
            continue
    try:
        return ImageFont.load_default(size)
    except TypeError:                      # Pillow < 10.1
        return ImageFont.load_default()
