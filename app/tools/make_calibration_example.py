"""Build example_data/SEM_calibration/JEOL-style_noMetadata_x25000.tif: an SEM image WITHOUT pixel-size metadata
(plain TIFF), JEOL-style black info strip with 'x25,000' and a 1 µm bar. Truth: 1000/252 = 3.968 nm/px (bar = 252 px)."""
import os, sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from mfs import readers as R

src = R.read_sem(open(os.path.join(ROOT, "example_data/SEM/Sulfo-MOF-801_SEM.tif"), "rb"), "x.tif")
img = src["image"][: src["databar_row"]]
W, H, STRIP = 1280, 960, 64
body = np.array(Image.fromarray(img).resize((W, H - STRIP), Image.BICUBIC))
canvas = Image.new("L", (W, H), 0)
canvas.paste(Image.fromarray(body), (0, 0))
d = ImageDraw.Draw(canvas)
font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", 24)
y = H - STRIP + 18
d.text((14, y), "SEI", fill=255, font=font)
d.text((90, y), "5.0kV", fill=255, font=font)
d.text((210, y), "x25,000", fill=255, font=font)
BAR = 252                                  # 1 µm at 3.96875 nm/px
x0 = 400
d.rectangle([x0, y + 4, x0 + BAR - 1, y + 9], fill=255)       # bar, label to its right like JEOL
d.text((x0 + BAR + 16, y - 2), "1µm", fill=255, font=font)
d.text((780, y), "WD 8.0mm", fill=255, font=font)
d.text((1010, y), "MOF  SEM", fill=255, font=font)
out = os.path.join(ROOT, "example_data/SEM_calibration/JEOL-style_noMetadata_x25000.tif")
canvas.save(out, compression="tiff_deflate")
print(out, "truth nm/px", 1000 / BAR)
