"""The dark palette in the sidebar is applied only on a dark (or unknown) Streamlit theme; on a light theme it is skipped."""
import logging
import types

import pytest

from mfs import sidebar_nav as SN


def _set_theme(monkeypatch, theme_type):
    ctx = types.SimpleNamespace(theme=types.SimpleNamespace(type=theme_type))
    monkeypatch.setattr(SN.st, "context", ctx)


def test_light_theme_skips_dark_palette(monkeypatch):
    _set_theme(monkeypatch, "light")
    assert SN._streamlit_theme_is_dark() is False


def test_dark_theme_keeps_dark_palette(monkeypatch):
    _set_theme(monkeypatch, "dark")
    assert SN._streamlit_theme_is_dark() is True


def test_unknown_theme_keeps_dark_palette(monkeypatch):
    _set_theme(monkeypatch, None)
    assert SN._streamlit_theme_is_dark() is True


def test_missing_context_keeps_dark_palette_and_logs(monkeypatch, caplog):
    class _Broken:
        @property
        def theme(self):
            raise RuntimeError("no script context")

    monkeypatch.setattr(SN.st, "context", _Broken())
    with caplog.at_level(logging.DEBUG, logger="mfs.sidebar_nav"):
        assert SN._streamlit_theme_is_dark() is True
    assert any("theme not available" in r.getMessage() for r in caplog.records)
