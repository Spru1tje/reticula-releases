"""
Semantic connection points (in 64-grid coordinates) for icons that are connected in schemes.
Automatic anchors (top / bottom / left / right / center of the drawing box) exist for every icon;
these named ones are added on top (same id = override). Directions: up | down | left | right | out.
"""
ANCHORS = {
    "round_bottom_flask": [("neck", 32, 5, "up")],
    "erlenmeyer_flask": [("neck", 32, 5, "up")],
    "volumetric_flask": [("neck", 32, 2, "up")],
    "test_tube": [("mouth", 32, 6, "up")],
    "beaker": [("mouth", 32, 7, "up"), ("spout", 11, 8, "left")],
    "separatory_funnel": [("top", 32, 1, "up"), ("outlet", 32, 61, "down")],
    "buchner_filtration": [("funnel", 32, 6, "up"), ("vacuum", 49, 47, "right")],
    "liebig_condenser": [("top", 32, 2, "up"), ("bottom", 32, 62, "down"),
                         ("water_in", 14, 48, "left"), ("water_out", 50, 14, "right")],
    "reflux_setup": [("top", 32, 4, "up"), ("water_in", 22, 27, "left"), ("water_out", 42, 10, "right")],
    "schlenk_line": [("vacuum", 4, 13.5, "left"), ("gas", 4, 23.5, "left"), ("tap1", 14, 42, "down"),
                     ("tap2", 28, 42, "down"), ("tap3", 42, 42, "down")],
    "schlenk_flask": [("neck", 32, 5, "up"), ("side_arm", 52, 11, "up")],
    "tube_furnace": [("inlet", 2.5, 30.5, "left"), ("outlet", 61.5, 30.5, "right")],
    "flow_reactor": [("inlet", 2, 22, "left"), ("outlet", 62, 44, "right")],
    "syringe_pump": [("outlet", 60, 27, "right")],
    "vacuum_pump": [("inlet", 17, 14, "up"), ("exhaust", 58, 33, "right")],
    "gas_cylinder": [("valve", 32, 3, "up")],
    "lc_column": [("inlet", 6, 6, "left"), ("outlet", 58, 58, "right")],
    "esi_source": [("capillary", 2, 32, "left"), ("ms_inlet", 52, 32, "right")],
    "lc_ms": [("lc_out", 26, 44, "right"), ("ms_in", 32, 44, "left")],
    "rotary_evaporator": [("vacuum", 54, 4, "up")],
    "centrifugal_filter": [("top", 32, 3, "up")],
    "micropipette": [("tip", 32, 61, "down")],
    "pasteur_pipette": [("tip", 32, 60, "down")],
    "uv_vis": [("beam_in", 3, 32, "left"), ("beam_out", 61, 32, "right")],
    "cuvette": [("beam_in", 23.5, 38, "left"), ("beam_out", 40.5, 38, "right")],
    "laser": [("beam", 60, 32, "right")],
    "potentiostat": [("cell", 48, 36, "up")],
    "sorption_analyzer": [("port1", 14, 46, "down"), ("port2", 24, 46, "down"), ("port3", 34, 46, "down")],
    "arrow_forward": [("start", 6, 32, "left"), ("end", 58, 32, "right")],
    "arrow_double": [("start", 6, 32, "left"), ("end", 58, 32, "right")],
    "arrow_curved": [("start", 10, 46, "down"), ("end", 54, 36, "down")],
    "arrow_branch": [("start", 4, 32, "left"), ("end1", 56, 16, "right"), ("end2", 56, 48, "right")],
}


def apply(icon):
    for aid, x, y, d in ANCHORS.get(icon.name, ()):
        icon.anchors = [a for a in icon.anchors if a["id"] != aid] + [{"id": aid, "x": x, "y": y, "dir": d}]
    return icon
