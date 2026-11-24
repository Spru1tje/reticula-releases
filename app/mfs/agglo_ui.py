"""v2.21: Streamlit tagger for agglomerates. Click an object on the micrograph overlay, tag it 'agglomerate' or 'rod',
and the learned model (mfs.agglo_learn) is refitted and applied to the measurements. Shows the cross-validated quality
and which features drive the decision, so the user can judge the model before trusting its counts."""
from __future__ import annotations

import copy
import math

import numpy as np
import streamlit as st

from mfs import agglo_learn as AL
from mfs import rods as R


def _overlay_image(im, lab, objs):
    """The overlay of mfs.rods (red rods, green agglomerates, blue film edge) as a PIL image."""
    from PIL import Image
    from scipy import ndimage as ndi
    edge = (lab > 0) & (lab != ndi.grey_erosion(lab, size=3))
    rgb = np.stack([np.asarray(im)] * 3, -1).astype(np.uint8)
    film = {o["id"] for o in objs if o["film_edge"]}
    agg = {o["id"] for o in objs if o.get("agglomerate")}
    rgb[edge] = [255, 0, 0]
    rgb[edge & np.isin(lab, list(agg))] = [0, 200, 0]
    rgb[edge & np.isin(lab, list(film))] = [0, 120, 255]
    return Image.fromarray(rgb)


@st.cache_data(show_spinner=False, max_entries=16)
def segment_cached(im, nm_per_px, min_area_px=25):
    """Segmentation and measurement of one micrograph, cached: a click in the tagger reruns the script but does not
    segment again (the segmentation is the slow part, 0.3 s to several seconds)."""
    lab, objs = R.segment_rods(im, nm_per_px, min_area_px=min_area_px)
    return lab, objs


def _nearest(objs, x, y, max_px=8.0):
    best, bd = None, max_px
    for o in objs:
        d = math.hypot(o["cx"] - x, o["cy"] - y)
        if d < bd:
            best, bd = o, d
    return best


def render_tagger(im, lab, objs, nm_per_px, image_name: str, store_path=None):
    """Render the tagger for one micrograph. Returns the summary of the model application (learned or rule)."""
    from streamlit_image_coordinates import streamlit_image_coordinates as sic

    tags = AL.load_tags(store_path)
    model = AL.fit(tags)
    objs = copy.deepcopy(objs)                    # the cached objects stay untouched; the flags are set on this copy
    summ = AL.apply(objs, model, nm_per_px, fallback=R.flag_agglomerates)

    st.markdown("**Agglomerate tag** — click an object on the picture, then choose what it is.")
    c1, c2 = st.columns([3, 2])
    with c1:
        click = sic(_overlay_image(im, lab, objs), key=f"agglo_click_{image_name}")
    with c2:
        cnt = AL.counts(tags)
        st.caption(f"Tags stored: {cnt['agglomerate']} agglomerate · {cnt['rod']} rod/particle")
        cv = AL.cross_validate(tags)
        if model is None:
            need = max(0, AL.MIN_PER_CLASS - min(cnt["agglomerate"], cnt["rod"]))
            st.info(f"Not enough tags yet: the rule is used. Tag at least {AL.MIN_PER_CLASS} of each kind "
                    f"({need} more of the rarer kind).")
        else:
            st.success(f"Learned model active ({model['n']} tags). "
                       f"Cross-validated accuracy {cv['accuracy']:.0%}, balanced {cv['balanced_accuracy']:.0%} "
                       f"(agglomerate recall {cv['sensitivity_agglomerate']:.0%}, rod recall {cv['specificity_rod']:.0%}).")
            st.caption("Coefficients (standardised; positive = pushes towards agglomerate): "
                       + ", ".join(f"{k} {v:+.2f}" for k, v in model["coef"].items()))
        st.caption(f"Objects now counted as agglomerate: {summ['n_agglomerate']} of {len(objs)} ({summ['method']}).")

        if click and "x" in click:
            o = _nearest(objs, click["x"], click["y"])
            if o is None:
                st.warning("No object at that spot. Click on the centre of a particle.")
            else:
                st.write(f"Object at ({o['cx']:.0f}, {o['cy']:.0f}) px: length {o['length_nm']:.1f} nm, "
                         f"width {o['width_nm']:.1f} nm, solidity {o['solidity']:.2f}, "
                         f"p(agglomerate) {o.get('p_agglomerate', float('nan')):.2f}")
                b1, b2, b3 = st.columns(3)
                if b1.button("Agglomerate", key=f"agg_{image_name}_{o['id']}"):
                    tags = AL.add_tag(tags, image_name, o, nm_per_px, _median_area(objs, nm_per_px), 1)
                    AL.save_tags(tags, store_path); st.rerun()
                if b2.button("Rod / particle", key=f"rod_{image_name}_{o['id']}"):
                    tags = AL.add_tag(tags, image_name, o, nm_per_px, _median_area(objs, nm_per_px), 0)
                    AL.save_tags(tags, store_path); st.rerun()
                if b3.button("Remove tag", key=f"rm_{image_name}_{o['id']}"):
                    tags = AL.remove_tag(tags, image_name, o)
                    AL.save_tags(tags, store_path); st.rerun()
    return summ


def _median_area(objs, nm_per_px):
    if not objs:
        return 1.0
    return float(np.median([o["area_px"] * nm_per_px ** 2 for o in objs]))
