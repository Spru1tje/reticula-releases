"""v2.14: parallel particle segmentation (one image per worker process) with an own result cache and real progress.

* Several images are segmented at the same time on all CPU cores (a 16-image batch that took ~50 s in series takes about as long
  as the slowest single image on a many-core PC).
* The pool is created once and kept (starting a worker costs 1-3 s: it imports numpy / scipy).
* Any failure of the pool (a worker killed, no permission to start processes, frozen app without freeze_support) falls back to
  the serial path - never to an error: the result is identical, only slower.
* Results are cached by a hash of the image and the parameters (LRU, 64 entries) - the same role st.cache_data had before.
* In the frozen Windows/macOS app multiprocessing needs multiprocessing.freeze_support() in mfs_launcher.py (done) and the
  workers need the app folder on sys.path (passed by the initialiser)."""
from __future__ import annotations

import hashlib
import os
import sys
import threading
import time
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed

MAX_CACHE = 64
_LOCK = threading.Lock()
_POOL = None
_POOL_WORKERS = 0
_CACHE: "OrderedDict[str, tuple]" = OrderedDict()
_DISABLED = os.environ.get("MFS_NO_PARALLEL") == "1"
LAST = {"mode": "", "workers": 0, "seconds": 0.0, "n": 0, "error": ""}      # what the last batch did (for the UI caption / tests)


def n_workers():
    try:
        c = os.cpu_count() or 2
    except Exception:
        c = 2
    return max(1, min(8, c - 1))


def _init_worker(paths):
    for p in paths:
        if p and p not in sys.path:
            sys.path.insert(0, p)
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):      # one thread per worker: no oversubscription
        os.environ.setdefault(k, "1")


def _work(kind, img, args):
    from mfs import segjobs
    return segjobs.run(kind, img, args)


def _get_pool():
    global _POOL, _POOL_WORKERS
    with _LOCK:
        if _POOL is None:
            _POOL_WORKERS = n_workers()
            _POOL = ProcessPoolExecutor(max_workers=_POOL_WORKERS, initializer=_init_worker, initargs=([p for p in sys.path if p],))
        return _POOL


def shutdown():
    global _POOL
    with _LOCK:
        if _POOL is not None:
            try:
                _POOL.shutdown(wait=False, cancel_futures=True)
            except Exception:
                pass
            _POOL = None


def key(kind, img, args):
    h = hashlib.md5()
    h.update(kind.encode())
    try:
        h.update(img.tobytes())
        h.update(repr(img.shape).encode())
    except Exception:
        h.update(repr(img).encode())
    h.update(repr(args).encode())
    return h.hexdigest()


def _put(k, v):
    with _LOCK:
        _CACHE[k] = v
        _CACHE.move_to_end(k)
        while len(_CACHE) > MAX_CACHE:
            _CACHE.popitem(last=False)


def _get(k):
    with _LOCK:
        v = _CACHE.get(k)
        if v is not None:
            _CACHE.move_to_end(k)
        return v


def cached(kind, img, args):
    return _get(key(kind, img, args)) is not None


def clear_cache():
    with _LOCK:
        _CACHE.clear()


def get(kind, img, args):
    """One result: from the cache or computed here (serial)."""
    from mfs import segjobs
    k = key(kind, img, args)
    v = _get(k)
    if v is None:
        v = segjobs.run(kind, img, args)
        _put(k, v)
    return v


def prefetch(jobs, progress=None, parallel=True):
    """Make sure every job (kind, img, args) is in the cache; the missing ones are computed in parallel.
    progress(done, total, seconds) is called after every finished image (done counts the cached ones too). Returns LAST."""
    t0 = time.perf_counter()
    todo = []
    for kind, img, args in jobs:
        k = key(kind, img, args)
        if _get(k) is None:
            todo.append((k, kind, img, args))
    total, done = len(jobs), len(jobs) - len(todo)
    LAST.update({"mode": "cache" if not todo else "", "workers": 0, "n": len(todo), "error": "", "seconds": 0.0})
    if progress and todo:
        progress(done, total, 0.0)
    if not todo:
        return LAST
    ran_parallel = False
    if parallel and not _DISABLED and len(todo) > 1:
        try:
            pool = _get_pool()
            futs = {pool.submit(_work, kind, img, args): k for k, kind, img, args in todo}
            for fu in as_completed(futs):
                _put(futs[fu], fu.result())
                done += 1
                if progress:
                    progress(done, total, time.perf_counter() - t0)
            ran_parallel = True
            LAST.update({"mode": "parallel", "workers": _POOL_WORKERS})
        except Exception as e:                           # broken pool / cannot start processes: finish serially
            LAST["error"] = f"{type(e).__name__}: {e}"
            shutdown()
            todo = [(k, kind, img, args) for k, kind, img, args in todo if _get(k) is None]
            done = total - len(todo)
    if not ran_parallel:
        from mfs import segjobs
        for k, kind, img, args in todo:
            if _get(k) is None:
                _put(k, segjobs.run(kind, img, args))
            done += 1
            if progress:
                progress(done, total, time.perf_counter() - t0)
        LAST["mode"] = LAST["mode"] or "serial"
    LAST["seconds"] = time.perf_counter() - t0
    return LAST
