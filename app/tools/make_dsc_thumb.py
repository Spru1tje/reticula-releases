"""Gallery thumbnail for the DSC page (assets/gallery/16_DSC.jpg) from the bundled TRIOS example."""
import io, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
from mfs import dsc, style
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d = dsc.read_dsc(os.path.join(R, "example_data", "DSC", "MOF-PEG_TA-TRIOS_EXAMPLE.txt"))
cs = [c for c in d["curves"]][:2]
m = d["meta"].get("mass_mg")
items, peaks, tgs = [], {}, {}
for i, c in enumerate(cs):
    c["id"] = i
    items.append((("heating" if c["kind"] == "heating" else "cooling"), c, dsc.to_specific(c, m)))
    peaks[i] = dsc.evaluate_all(c, "linear", m)
    if c["kind"] == "heating":
        w = dsc.auto_tg_window(c, peaks[i]); g = dsc.glass_transition(c, *w, mass_mg=m)
        if g: tgs[i] = [g]
style.apply_style(8)
fig, ax = plt.subplots(figsize=(3.2, 2.4))
dsc.draw_dsc(ax, items, "T", True, 0.6, peaks=peaks, tgs=tgs, annotate=True)
ax.legend(frameon=False, loc="upper right", fontsize=6)
fig.tight_layout(); b = io.BytesIO(); fig.savefig(b, dpi=150, format="png"); b.seek(0)
im = Image.open(b).convert("RGB"); im.thumbnail((480, 360))
bg = Image.new("RGB", (480, 360), "white"); bg.paste(im, ((480 - im.width) // 2, (360 - im.height) // 2))
bg.save(os.path.join(R, "assets", "gallery", "16_DSC.jpg"), quality=90); print("ok")
