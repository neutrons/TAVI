"""Contain the entry point for the application."""

import neutrons_standard

neutrons_standard.init("tavi")

try:
    from ._version import __version__  # noqa: F401
except ImportError:
    __version__ = "unknown"
