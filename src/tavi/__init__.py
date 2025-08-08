"""
Contains the entry point for the application
"""
try:
    from ._version import __version__  # noqa: F401
except ImportError:
    __version__ = "unknown"

def tavi():  # noqa: N802
    """Start Class"""

    from .tavimain import tavi  # noqa: E501, N813

    return tavi()