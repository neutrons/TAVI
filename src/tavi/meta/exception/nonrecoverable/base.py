"""Base class for non-recoverable exceptions."""

from tavi.meta.exception.tavi_exception import TaviError


class NonRecoverableError(TaviError):
    """Type that all non-recoverable exceptions must inherit."""

    pass
