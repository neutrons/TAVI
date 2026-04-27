"""Contain the entry point for the application."""

# NOTE: WARNING, These MUST be imported before neutrons_standard due to an interaction between qtpy
#       and ruamel.yaml that causes some parsing warnings.
import qtpy.QtCore  # isort: skip
import qtpy.QtWidgets  # isort: skip

import neutrons_standard  # isort: skip


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
