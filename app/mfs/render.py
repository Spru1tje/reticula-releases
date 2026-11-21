"""Renderer helper: text measurements need a renderer, but not every matplotlib canvas has get_renderer()
(e.g. a plain FigureCanvasBase). Fall back to an Agg renderer of the same size and resolution."""
from __future__ import annotations


def renderer_for(fig):
    canvas = fig.canvas
    if hasattr(canvas, "get_renderer"):
        return canvas.get_renderer()
    from matplotlib.backends.backend_agg import RendererAgg
    return RendererAgg(fig.bbox.width, fig.bbox.height, fig.dpi)
