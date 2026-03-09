"""Contain the entry point for the application."""

import neutrons_standard

neutrons_standard.init("tavi")


from tavi.meta.logging import init_logging  # noqa: E402

init_logging()

try:
    from ._version import __version__  # noqa: F401
except ImportError:
    __version__ = "unknown"
