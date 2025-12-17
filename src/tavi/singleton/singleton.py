"""Singleton wrapper."""
from functools import wraps
from typing import List

def Singleton(orig_cls) -> None:
    """Singleton wrapper function."""
    orig_new = orig_cls.__new__
    orig_init = orig_cls.__init__
    instance = None
    initialized = False

    @wraps(orig_cls.__init__)
    def __init__(self, *args, **kwargs):
        """Initialization."""
        nonlocal initialized
        if initialized:
            return
        initialized = True
        orig_init(self, *args, **kwargs)

    @wraps(orig_cls.__new__)
    def __new__(cls, *args, **kwargs):  # noqa: ARG001
        """Handle instance."""
        nonlocal instance
        if instance is None:
            # this needs to work with object.__new__, which only has only the `cls` arg
            instance = orig_new(cls)  # , *args, **kwargs)
        return instance
