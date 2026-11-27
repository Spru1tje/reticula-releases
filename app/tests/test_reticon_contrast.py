"""Contrast of the Reticon palette tokens on a white background (WCAG relative luminance).

Ink (outlines, text-like marks) and accent (the one element an icon is about) must reach 3:1 on paper.
Muted is a decorative token (axes, supports, shadow) and is deliberately NOT held to 3:1; it is
documented here so the exception is visible. Dark theme is not defined yet, so it is not checked.
"""
from pathlib import Path

PAPER = "#FFFFFF"
TOKENS = {"ink": "#1F2D3D", "accent": "#2E6F95"}
DECORATIVE = {"muted": "#8FA3B5"}
MIN_RATIO = 3.0


def _lum(hexcode: str) -> float:
    channels = [int(hexcode[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    lin = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in channels]
    return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]


def contrast(a: str, b: str) -> float:
    hi, lo = sorted([_lum(a), _lum(b)], reverse=True)
    return (hi + 0.05) / (lo + 0.05)


def test_ink_and_accent_reach_three_to_one_on_paper():
    for name, colour in TOKENS.items():
        ratio = contrast(colour, PAPER)
        assert ratio >= MIN_RATIO, f"{name} {colour} on paper is {ratio:.2f}:1 (need {MIN_RATIO}:1)"


def test_contrast_helper_matches_known_values():
    assert round(contrast("#000000", "#FFFFFF"), 1) == 21.0
    assert round(contrast("#FFFFFF", "#FFFFFF"), 1) == 1.0


def test_detects_low_contrast():
    assert contrast("#DDDDDD", PAPER) < MIN_RATIO
