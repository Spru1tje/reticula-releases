"""v2.17: EJIC 2020 rule (Lázaro et al., doi:10.1002/ejic.202000656): R_exp uses the mass at the START of the decomposition profile."""
import numpy as np

from mfs import analysis as A


def _curve():
    """Synthetic MOF TGA: flat 100 % up to ~400 °C, then one linker-combustion step (Gaussian DTG peak at 450 °C), residue 55 %."""
    T = np.linspace(25, 800, 1550)
    rate = 0.9 * np.exp(-0.5 * ((T - 450) / 18.0) ** 2)            # % per °C
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (rate[1:] + rate[:-1]) * np.diff(T))])
    return T, 100.0 - cum


def test_onset_is_where_the_mass_starts_to_fall():
    T, w = _curve()
    a = A.tga_analysis(T, w)
    assert a["T_onset"] is not None
    # the foot of the peak at 450 °C (sigma 18 °C): the mass starts to fall around 400-420 °C
    assert 390 <= a["T_onset"] <= 430
    assert 90.0 < a["W_onset"] < 99.0      # the tangent crosses the plateau where the mass has just started to fall


def test_onset_differs_from_the_flattest_point_rule_is_an_explicit_choice():
    T, w = _curve()
    a = A.tga_analysis(T, w)
    # both rules are reported; the EJIC rule is the onset, the flattest stretch is kept as the default elsewhere
    assert a["T_plateau"] is not None and a["T_onset"] is not None


def test_no_fall_means_no_onset():
    T = np.linspace(25, 800, 1550)
    w = np.full_like(T, 100.0)
    a = A.tga_analysis(T, w)
    assert a["T_onset"] is None
