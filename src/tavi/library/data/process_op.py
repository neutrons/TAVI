"""Base class and shared helpers for the operations that derive ProcessedScan objects."""

import abc
from typing import Sequence

from tavi.library.data.scan import (
    UUID,
    ProcessedScan,
)
from tavi.library.data.tavi_data import TaviData

#: Name the rebinning operations give the uncertainty column they compute, matching the
#: ``error_name`` a ``PlotSeries`` is built with.
ERROR_COLUMN = "error"


def _norm_suffix(norm: str, multiplication: float) -> str:
    """
    Spell a normalization out for a friendly name, so two normalizations of one scan stay apart.

    Args:
        norm: Name of the column the y column was put on a per-unit footing against.
        multiplication: Scale applied to the y column.

    Returns:
        The suffix, naming the norm column and, unless it is 1, the scale.

    """
    return f"_norm{norm}" if multiplication == 1 else f"_norm{norm}x{multiplication:g}"


class ProcessOps(metaclass=abc.ABCMeta):
    """
    Base class for process operations.

    An operation resolves its origin uuids against a ``TaviData`` pool and
    returns a new ``ProcessedScan``. Producing is not storing: registering the
    result is the model layer's job.

    Subclasses take whatever else they need in their own ``__init__``, check
    their inputs in ``validate`` and build the result in ``exec``. Calling
    ``validate`` is left to ``exec``, not enforced here.
    """

    def __init__(self, tavi_data: TaviData, uuids: Sequence[UUID]) -> None:
        """
        Store what the operation runs on.

        Args:
            tavi_data: The pool ``uuids`` are looked up in.
            uuids: Origin scan uuids, in the order the operation takes them.

        """
        self.tavi_data = tavi_data
        self.uuids = uuids

    @property
    def size(self) -> int:
        """Number of uuids, can be used for GUI generation."""
        return len(self.uuids)

    @abc.abstractmethod
    def validate(self) -> None:
        """Raise if the operation's inputs cannot produce a scan."""

    @abc.abstractmethod
    def exec(self) -> ProcessedScan:
        """Run the operation and return the derived scan, registered nowhere."""
