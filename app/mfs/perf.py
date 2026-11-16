"""v2.11 speed helpers: a content fingerprint of a matplotlib figure (stable across reruns), a bounded cache of preview
PNGs keyed on it, a cheap preview interpolation for big images (exports keep the full-quality one), trace downsampling
for interactive display and a file-content hash for st.cache_data keys."""
from __future__ import annotations
import hashlib, io, pickle, threading
from collections import OrderedDict
import numpy as np

_LOCK = threading.Lock()
_CACHE: "OrderedDict[str, bytes]" = OrderedDict()
MAX_ITEMS, MAX_BYTES = 96, 160 * 2 ** 20
STATS = {"hits": 0, "misses": 0, "render_s": 0.0}


class _Sink:
    def __init__(self):
        self.h = hashlib.blake2b(digest_size=16)

    def write(self, b):
        self.h.update(b)


def _arr_digest(a):
    a = np.ascontiguousarray(a)
    if a.dtype == object:
        return hashlib.blake2b(repr(a.tolist()).encode(), digest_size=16).hexdigest()
    try:
        return hashlib.blake2b(memoryview(a).cast("B"), digest_size=16).hexdigest()
    except (TypeError, ValueError):
        return hashlib.blake2b(a.tobytes(), digest_size=16).hexdigest()


class _P(pickle.Pickler):
    def reducer_override(self, obj):
        from matplotlib.transforms import TransformNode
        from matplotlib.figure import Figure
        if isinstance(obj, np.ndarray):                      # arrays by content digest (no copy into the pickle)
            return (np.ndarray, (obj.shape, obj.dtype.str, _arr_digest(obj)))
        if isinstance(obj, TransformNode):                   # '_parents' holds id()s -> unstable
            st = dict(obj.__getstate__()); st.pop("_parents", None)
            return (type(obj), (), st)
        if isinstance(obj, Figure):                          # the pyplot figure number changes every run
            st = dict(obj.__getstate__()); st.pop("_number", None)
            return (Figure, (), st)
        return NotImplemented


def fingerprint(fig, *extra):
    """Hash of everything that is drawn (artists, data, styles) + the rcParams + `extra`; None if not picklable."""
    import matplotlib as mpl
    s = _Sink()
    try:
        _P(s, protocol=4).dump((fig, extra))
        s.write(repr(sorted((k, repr(v)) for k, v in mpl.rcParams.items())).encode())
    except Exception:
        return None
    return s.h.hexdigest()


def _cache_put(k, v):
    with _LOCK:
        _CACHE[k] = v; _CACHE.move_to_end(k)
        while len(_CACHE) > MAX_ITEMS or sum(len(x) for x in _CACHE.values()) > MAX_BYTES:
            _CACHE.popitem(last=False)


def fast_images(fig, interp="antialiased"):
    """Context: big images drawn with a cheap interpolation (lanczos on a 1–4 Mpx SEM image costs ~0.3 s per panel at
    200 dpi); restored afterwards so exports keep the original."""
    from contextlib import contextmanager

    @contextmanager
    def cm():
        saved = []
        ims = list(getattr(fig, "images", []))
        for ax in fig.axes:                                  # not fig.findobj(): that builds every tick (~50 ms/axes)
            ims.extend(getattr(ax, "images", []))
        for im in ims:
            a = im.get_array()
            if a is not None and getattr(a, "size", 0) > 250_000 and im.get_interpolation() not in ("nearest", "none", interp):
                saved.append((im, im.get_interpolation())); im.set_interpolation(interp)
        try:
            yield
        finally:
            for im, i in saved:
                im.set_interpolation(i)
    return cm()


def fast_lines(fig, max_points=8000):
    """Context: lines with more than `max_points` points are drawn min/max-decimated (peaks kept) for the on-screen
    preview only; exports and CSVs use the full data."""
    from contextlib import contextmanager

    @contextmanager
    def cm():
        saved = []
        for ax in fig.axes:
            for ln in ax.get_lines():
                try:
                    x, y = ln.get_xdata(orig=True), ln.get_ydata(orig=True)
                    if np.ndim(x) == 1 and len(x) > max_points and len(x) == len(y) and ln.get_marker() in ("None", None, "", " "):
                        xs, ys = downsample(np.asarray(x, float), np.asarray(y, float), max_points)
                        saved.append((ln, x, y)); ln.set_data(xs, ys)
                except Exception:
                    continue
        try:
            yield
        finally:
            for ln, x, y in saved:
                ln.set_data(x, y)
    return cm()


def preview_png(fig, dpi=200, lock=None):
    """PNG bytes for the on-screen preview; identical figures (same fingerprint) are not rendered again."""
    import time
    k = fingerprint(fig, dpi)
    if k is not None:
        with _LOCK:
            b = _CACHE.get(k)
            if b is not None:
                _CACHE.move_to_end(k); STATS["hits"] += 1
                return b
    t = time.perf_counter()
    buf = io.BytesIO()
    from contextlib import nullcontext
    with (lock or nullcontext()), fast_images(fig), fast_lines(fig):
        fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    b = buf.getvalue()
    STATS["misses"] += 1; STATS["render_s"] += time.perf_counter() - t
    if k is not None:
        _cache_put(k, b)
    return b


def downsample(x, y, max_points=4000):
    """Min/max decimation for display (keeps peaks, vectorised): <= ~max_points points. Exports use the full data."""
    x = np.asarray(x); y = np.asarray(y)
    n = len(x)
    if n <= max_points or max_points < 4:
        return x, y
    nb = max_points // 2
    k = int(np.ceil(n / nb))
    m = (n // k) * k
    yy = np.where(np.isfinite(y[:m]), y[:m], np.nan).reshape(-1, k)
    base = np.arange(yy.shape[0]) * k
    with np.errstate(all="ignore"):
        lo = np.nan_to_num(yy, nan=np.inf).argmin(1); hi = np.nan_to_num(yy, nan=-np.inf).argmax(1)
    idx = np.sort(np.concatenate([base + lo, base + hi, np.arange(m, n)]))
    idx = np.unique(idx)
    return x[idx], y[idx]


def file_hash(f) -> str:
    """Content hash of an uploaded/bundled file object or bytes (st.cache_data key)."""
    if isinstance(f, (bytes, bytearray, memoryview)):
        b = bytes(f)
    else:
        try:
            b = f.getvalue()
        except AttributeError:
            pos = f.tell(); f.seek(0); b = f.read(); f.seek(pos)
    return hashlib.blake2b(b, digest_size=16).hexdigest()


# ----------------------------------------------------------------------------- Streamlit: O(1) widget state lookups
class _LazyFiltered(dict):
    """Stand-in for SessionState.filtered_state. Streamlit's widgets only ask `key in state` / `state[key]`, but building
    the full dict costs O(all keys) per widget (~1–2 ms each with a few hundred keys). Membership and lookup are answered
    directly; anything else (iteration, len, items, ...) materialises the real dict."""
    __slots__ = ("_safe",)

    def __init__(self, safe):
        dict.__init__(self); self._safe = safe

    def _full(self):
        with self._safe._lock:
            return _ORIG_FILTERED.fget(self._safe)

    def __contains__(self, key):
        from streamlit.runtime.state.session_state import _is_internal_key
        from streamlit.runtime.state.common import is_element_id
        if not isinstance(key, str) or _is_internal_key(key) or is_element_id(key):
            return key in self._full()
        with self._safe._lock:
            try:
                self._safe._state[key]
                return True
            except KeyError:
                return False

    def __getitem__(self, key):
        if key in self:
            with self._safe._lock:
                return self._safe._state[key]
        raise KeyError(key)

    def get(self, key, default=None):
        return self[key] if key in self else default

    def __iter__(self):
        return iter(self._full())

    def __len__(self):
        return len(self._full())

    def keys(self):
        return self._full().keys()

    def items(self):
        return self._full().items()

    def values(self):
        return self._full().values()

    def __repr__(self):
        return repr(self._full())

    def __eq__(self, other):
        return self._full() == other

    def copy(self):
        return dict(self._full())


_ORIG_FILTERED = None


def install_fast_widget_state():
    """Patch SafeSessionState.filtered_state (guarded: any import problem leaves Streamlit untouched)."""
    global _ORIG_FILTERED
    try:
        from streamlit.runtime.state.safe_session_state import SafeSessionState
        from streamlit.runtime.state.session_state import _is_internal_key   # noqa: F401  (must exist)
        from streamlit.runtime.state.common import is_element_id             # noqa: F401
        cur = SafeSessionState.__dict__.get("filtered_state")
        if not isinstance(cur, property) or getattr(cur.fget, "_mfs_fast", False):
            return
        _ORIG_FILTERED = cur
        def fget(self):
            return _LazyFiltered(self)
        fget._mfs_fast = True
        SafeSessionState.filtered_state = property(fget)
    except Exception:
        pass
