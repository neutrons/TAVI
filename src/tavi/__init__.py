"""Contain the entry point for the application."""

import neutrons_standard

neutrons_standard.init("tavi")

from .library import init as init_library  # noqa: E402

init_library()

try:
    from ._version import __version__  # noqa: F401
except ImportError:
    __version__ = "unknown"
