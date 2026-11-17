"""v2.14: progress bar with an estimate of the time left for the particle sizer (segmentation of one or many images).

Streamlit blocks the script while a segmentation runs, so the bar is moved by a small worker thread. The estimate is
seconds per megapixel of the image for each algorithm, learned from the runs of this session (exponential average; a
first guess before the first run). Results that come straight from the cache return at once and show no bar."""
from __future__ import annotations

import threading
import time

import streamlit as st

DEFAULT_S_PER_MP = {"adaptive": 6.0, "watershed": 3.0, "auto": 1.0}
SHOW_AFTER_S = 0.4                      # a cached result is back before this: nothing is drawn


def _ctx():
    try:
        from streamlit.runtime.scriptrunner import add_script_run_ctx, get_script_run_ctx
        return add_script_run_ctx, get_script_run_ctx()
    except Exception:
        return None, None


class SegProgress:
    """One bar for a whole batch: 'Image 2 of 5 · 14 s elapsed · about 31 s left'."""

    def __init__(self, n_images: int, label: str = "Segmenting particles"):
        self.n, self.label = max(1, int(n_images)), label
        self.box = st.empty()
        self.t0 = time.perf_counter()
        self.done = 0
        self.rates = st.session_state.setdefault("_seg_rate", {})
        self._pending = 0.0             # estimated seconds of the images still to do (set per call)

    @staticmethod
    def megapixels(img) -> float:
        a = getattr(img, "shape", None)
        return (float(a[0]) * float(a[1]) / 1e6) if a is not None and len(a) >= 2 else 1.0

    def estimate(self, algo: str, mp: float) -> float:
        return self.rates.get(algo, DEFAULT_S_PER_MP.get(algo, 2.0)) * max(mp, 0.05)

    def call(self, img, algo: str, fn, *args, **kw):
        """Run fn(img, *args) with the bar; returns its result. algo ∈ adaptive | watershed | auto."""
        mp = self.megapixels(img)
        est = self.estimate(algo, mp)
        left_after = (self.n - self.done - 1) * est
        stop = threading.Event()
        t_start = time.perf_counter()
        add_ctx, ctx = _ctx()

        def _loop():
            while not stop.wait(0.25):
                el = time.perf_counter() - t_start
                if el < SHOW_AFTER_S:
                    continue
                e = max(est, el * 1.15)                          # never claims to be done early
                frac = (self.done + min(0.97, el / e)) / self.n
                left = max(0.0, e - el) + left_after
                try:
                    self.box.progress(min(frac, 0.99), text=f"{self.label} · image {self.done + 1} of {self.n} · {el:.0f} s elapsed · "
                                                            f"about {left:.0f} s left (estimate from earlier runs)")
                except Exception:
                    return
        th = threading.Thread(target=_loop, daemon=True)
        if add_ctx and ctx is not None:
            add_ctx(th, ctx)
        th.start()
        try:
            return fn(img, *args, **kw)
        finally:
            stop.set(); th.join(timeout=1.0)
            dt = time.perf_counter() - t_start
            if dt > SHOW_AFTER_S:                                 # a real run (not a cache hit): learn the speed
                r = dt / max(mp, 0.05)
                self.rates[algo] = r if algo not in self.rates else 0.6 * self.rates[algo] + 0.4 * r
            self.done += 1

    def finish(self):
        try:
            self.box.empty()
        except Exception:
            pass
