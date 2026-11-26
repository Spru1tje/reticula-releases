"""Synthetic test of touching and separated rods in the rod measurement (mfs.rods.segment_rods).

Two 40 x 8 nm rods end to end, 1 nm per pixel, flat background.
  * gap 6 px (separated): two objects; each length measured 42.4 nm (pinned within 1 nm), i.e. +6 % over 40 nm.
  * gap 0 px (touching): ONE object of 90.9 nm (pinned within 1 nm). The two rods together are 80 nm in the
    synthetic image, so the chain is measured as one object and its length is overestimated. This is a KNOWN
    LIMITATION of the current method: a chain of touching rods is not split and is not reported as two rods.
No real micrograph is used; this file is synthetic only.
"""
import numpy as np
from scipy import ndimage as ndi

from mfs import rods as RO


def _image(rods, size=200):
    yy, xx = np.mgrid[0:size, 0:size].astype(float)
    im = np.full((size, size), 200.0)
    for cx, cy, L, W, th in rods:
        t = np.deg2rad(th)
        dx, dy = xx - cx, yy - cy
        u = dx * np.cos(t) + dy * np.sin(t)
        v = -dx * np.sin(t) + dy * np.cos(t)
        im[(u / (L / 2)) ** 2 + (v / (W / 2)) ** 2 <= 1] = 60.0
    return ndi.gaussian_filter(im, 0.7)


def _pair(gap):
    x1 = 100 - (40 + gap) / 2
    x2 = 100 + (40 + gap) / 2
    return _image([(x1, 100, 40, 8, 0), (x2, 100, 40, 8, 0)])


def test_separated_rods_are_two_objects_with_pinned_length():
    _, objs = RO.segment_rods(_pair(6), 1.0, min_area_px=25)
    assert len(objs) == 2
    assert all(abs(o["length_nm"] - 42.4) < 1.0 for o in objs)


def test_touching_rods_are_one_object_known_limitation():
    _, objs = RO.segment_rods(_pair(0), 1.0, min_area_px=25)
    assert len(objs) == 1                      # chain is not split (known limitation)
    assert abs(objs[0]["length_nm"] - 90.9) < 1.0   # pinned measured length of the merged chain
