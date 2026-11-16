"""Synthetic SEM images with burned-in info strips in common vendor layouts (for the calibration tests)."""
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont

FONT = next((p for p in ["/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
                         "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", r"C:\Windows\Fonts\arial.ttf"] if os.path.exists(p)), None)


def _texture(w, h, seed=0):
    rng = np.random.default_rng(seed)
    a = rng.normal(110, 25, (h // 8 + 1, w // 8 + 1))
    a = np.array(Image.fromarray(np.clip(a, 0, 255).astype(np.uint8)).resize((w, h), Image.BICUBIC), float)
    for _ in range(60):                       # bright crystals
        x, y, r = rng.integers(0, w), rng.integers(0, h), rng.integers(8, 40)
        yy, xx = np.ogrid[:h, :w]
        a[(abs(xx - x) + abs(yy - y)) < r] = rng.integers(170, 235)
    return np.clip(a + rng.normal(0, 8, a.shape), 0, 255).astype(np.uint8)


def make(style, w=1024, h=884, seed=0):
    """Returns (uint8 image, truth dict(nm_per_px, bar_px, label, mag))."""
    strip = 80 if style != "overlay" else 0
    img = Image.new("L", (w, h + strip), 0 if style not in ("white_strip",) else 255)
    img.paste(Image.fromarray(_texture(w, h, seed)), (0, 0))
    d = ImageDraw.Draw(img)
    f = ImageFont.truetype(FONT, 22)
    y = h + 20
    if style == "zeiss":            # label above bar, left; Mag = 10.00 K X
        L = 300; d.text((40, y - 6), "1 µm", fill=255, font=f); d.rectangle([40, y + 26, 40 + L - 1, y + 33], fill=255)
        d.text((420, y), "EHT = 5.00 kV", fill=255, font=f); d.text((640, y), "Mag = 10.00 K X", fill=255, font=f)
        return np.array(img), dict(nm_per_px=1000 / L, bar_px=L, label=(1.0, "µm"), mag=10000.0)
    if style == "fei_box":          # framed bar with end ticks, '500 nm' below, 'x50 000' mag field
        L = 220; x0 = 700
        d.text((30, y), "HV 5.00 kV", fill=255, font=f); d.text((200, y), "mag 50 000 x", fill=255, font=f)
        d.rectangle([x0, y + 4, x0 + L - 1, y + 10], fill=255)
        d.rectangle([x0, y - 2, x0 + 2, y + 16], fill=255); d.rectangle([x0 + L - 3, y - 2, x0 + L - 1, y + 16], fill=255)
        d.text((x0 + L // 2 - 35, y + 22), "500 nm", fill=255, font=f)
        return np.array(img), dict(nm_per_px=500 / L, bar_px=L, label=(500.0, "nm"), mag=50000.0)
    if style == "white_strip":      # black bar on a white strip (Hitachi-like), '10.0kx'
        L = 410; x0 = 500
        d.text((20, y), "5.0kV 8.1mm x10.0k SE", fill=0, font=f)
        d.rectangle([x0, y + 8, x0 + L - 1, y + 13], fill=0); d.text((x0 + L + 12, y), "5.00um", fill=0, font=f)
        return np.array(img), dict(nm_per_px=5000 / L, bar_px=L, label=(5.0, "µm"), mag=10000.0)
    if style == "overlay":          # bar drawn on the image itself, no strip
        L = 180; x0 = w - L - 40; yb = h - 50
        d.rectangle([x0 - 10, yb - 36, x0 + L + 10, yb + 14], fill=0)
        d.rectangle([x0, yb, x0 + L - 1, yb + 6], fill=255); d.text((x0 + L // 2 - 30, yb - 32), "200 nm", fill=255, font=f)
        return np.array(img), dict(nm_per_px=200 / L, bar_px=L, label=(200.0, "nm"), mag=None)
    raise ValueError(style)


STYLES = ["zeiss", "fei_box", "white_strip", "overlay"]
