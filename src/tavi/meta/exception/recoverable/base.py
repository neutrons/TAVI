"""Base class for recoverable exceptions."""

from tavi.meta.exception.tavi_exception import TaviError


class RecoverableError(TaviError):
    """Type that all recoverable exceptions must inherit."""

    pass
