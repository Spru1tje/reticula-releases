"""v2.19.4 CIE chart: the sRGB gamut test and the displayed colour of a chromaticity, with known answers."""
from mfs import plextras as PX


def test_white_point_is_inside_and_neutral():
    rgb, inside = PX.xy_to_srgb(0.3127, 0.3290)          # D65
    assert inside
    assert all(abs(c - rgb[0]) < 0.01 for c in rgb)      # neutral grey / white


def test_sRGB_primaries_are_on_the_gamut_boundary_and_show_pure_colours():
    r, inside_r = PX.xy_to_srgb(0.64, 0.33)
    g, inside_g = PX.xy_to_srgb(0.30, 0.60)
    assert inside_r and inside_g
    assert r[0] > 0.99 and r[1] < 0.02 and r[2] < 0.02    # pure red
    assert g[1] > 0.99 and g[0] < 0.02 and g[2] < 0.02    # pure green


def test_saturated_spectral_point_outside_sRGB_is_flagged():
    _, inside = PX.xy_to_srgb(0.70, 0.30)                  # beyond the red primary
    assert not inside


def test_draw_cie_runs_and_marks_points(tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    PX.draw_cie(ax, [(0.3127, 0.3290), (0.70, 0.30)], labels=["D65", "outside"])
    labels = [line.get_label() for line in ax.get_lines()]
    assert "sRGB gamut" in labels
    assert "D65" in labels and "outside" in labels
    plt.close(fig)
