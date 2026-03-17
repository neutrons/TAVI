"""Contain the entry point for the application."""

import neutrons_standard

neutrons_standard.init("tavi")

from .library import init as init_library  # noqa: E402

try:
    init_library()
except ModuleNotFoundError as e:
    # Unable to init submodule. Still valid for ci, proceeding.
    pass

try:
    from ._version import __version__  # noqa: F401
except ImportError:
    __version__ = "unknown"
