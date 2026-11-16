"""Multi-image SEM panels: every annotation (scale-bar label, magnification, panel letter) lies inside its own axes."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mfs import plots as P, style

E = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "example_data", "SEM")


def _ims():
    from PIL import Image
    out = []
    for f in sorted(os.listdir(E))[:4]:
        a = np.asarray(Image.open(os.path.join(E, f)).convert("L"))
        out.append({"image": a, "pixel_nm": 2.5, "crop_row": int(a.shape[0] * 0.9), "mag_text": "×25 000", "title": "MOF"})
    return out


def _outside(fig):
    fig.canvas.draw(); r = fig.canvas.get_renderer(); bad = []
    for ax in fig.axes:
        ab = ax.get_window_extent(r)
        for t in ax.texts:
            if not t.get_text():
                continue
            bb = t.get_window_extent(r)
            if not (bb.x0 >= ab.x0 and bb.x1 <= ab.x1 and bb.y0 >= ab.y0 and bb.y1 <= ab.y1):
                bad.append(t.get_text())
    return bad


def test_labels_inside_axes_all_layouts():
    ims = _ims()
    for ncols, width, base in [(4, 3.25, 9), (2, 3.25, 10), (4, 7.0, 12), (1, 3.25, 8), (3, 7.0, 7)]:
        style.apply_style(); plt.rcParams["font.size"] = base
        fig = P.fig_sem_panel(ims[: max(ncols, 2)], ncols, width)
        assert not _outside(fig), (ncols, width, _outside(fig))
        fss = {round(t.get_fontsize(), 2) for ax in fig.axes for t in ax.texts if t.get_text()}
        assert len(fss) == 1                                   # one uniform font across panels
        plt.close(fig)


def test_font_override_and_relative_thickness():
    ims = _ims()
    fig = P.fig_sem_panel(ims[:2], 2, 7.0, font_size=8)
    assert all(abs(t.get_fontsize() - 8) < 1e-6 for ax in fig.axes for t in ax.texts if t.get_text())
    assert P.panel_font_size(0.8, 0.8) == 4.5 and P.panel_font_size(5, 5) == 9.0
    plt.close(fig)
