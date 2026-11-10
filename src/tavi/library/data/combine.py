"""One entry point for combining loaded scans, appending or rebinning them as the parameters describe."""

from typing import Optional, Sequence

from tavi.library.data.algorithms.append_ops import AppendOp
from tavi.library.data.algorithms.bin_by_tol import ToleranceRebinOps
from tavi.library.data.algorithms.equal_step_bin import EqualStepRebinOps
from tavi.library.data.enum.rebin_mode import RebinMode
from tavi.library.data.process_op import ProcessOps
from tavi.library.data.scan import (
    UUID,
    ProcessedScan,
    RawScan,
)
from tavi.library.data.tavi_data import TaviData


class CombineOps(ProcessOps):
    """
    Combine one or more raw scans into a single ``ProcessedScan``, registered in ``tavi_data``.

    Which combination is run is read off the parameters rather than named: ``tol`` asks
    for a tolerance rebin, ``step`` for an equal-step one, and neither for a plain
    append. The combination itself is delegated to ``ToleranceRebinOps``,
    ``EqualStepRebinOps`` or ``AppendOp`` - this class resolves the mode, refuses
    parameters that mode has no use for, and stores what comes back.

    Unlike every other ``ProcessOps``, ``exec`` registers its result in
    ``tavi_data.processed_scans``. Combining is how a derived scan comes to exist in a
    project, and leaving that write to each caller is what has kept the store empty.
    """

    def __init__(
        self,
        tavi_data: TaviData,
        uuids: Sequence[UUID],
        columns: Sequence[str],
        norm: str = "",
        tol: Optional[float] = None,
        step: Optional[float] = None,
        start: Optional[float] = None,
        finish: Optional[float] = None,
        multiplication: float = 1.0,
    ) -> None:
        """
        Store the origins, the columns to carry over, and whichever rebin the caller is asking for.

        ``tol`` and ``step`` select the mode and are mutually exclusive; ``start`` and
        ``finish`` are the equal-step grid's outer edges, and ``norm`` and
        ``multiplication`` the weighting both rebins apply. Supplying none of them
        appends the origins and nothing more.
        """
        super().__init__(tavi_data, uuids)
        self.columns = columns
        self.norm = norm
        self.tol = tol
        self.step = step
        self.start = start
        self.finish = finish
        self.multiplication = multiplication

    @property
    def mode(self) -> RebinMode:
        """The combination the parameters describe, for a caller deciding whether to run this op."""
        if self.tol is not None and self.step is not None:
            raise ValueError(
                f"Combine takes either a tolerance or a step, not both, got tol={self.tol!r}, step={self.step!r}."
            )
        if self.tol is not None:
            return RebinMode.TOLERANCE
        if self.step is not None:
            return RebinMode.EQUAL_STEP
        return RebinMode.NONE

    def validate(self) -> None:
        """Check the parameters belong to the mode they resolve to, then everything that mode's op requires."""
        mode = self.mode

        # the rebins read min()/max() off the appended x column, and an append of nothing
        # yields empty columns rather than failing - neither is a combination.
        if not self.uuids:
            raise ValueError("Combine needs at least one scan to combine.")

        if mode is RebinMode.NONE:
            # silently dropping these would hand back an unnormalized scan to a caller who
            # asked for a normalized one, with nothing in the result saying so.
            if self.norm:
                raise ValueError(
                    f"Combine normalization column {self.norm!r} needs a rebin to weigh by it; "
                    "pass a tolerance or a step."
                )
            if self.multiplication != 1:
                raise ValueError(
                    f"Combine multiplication {self.multiplication!r} needs a rebin to scale; "
                    "pass a tolerance or a step."
                )

        if mode is not RebinMode.EQUAL_STEP:
            # edges are positions on a grid, and only a step lays one down.
            if self.start is not None or self.finish is not None:
                raise ValueError("Combine start and finish only applies to an equal-step rebin; pass a step.")

        # the rebins refuse an already-reduced origin themselves, but the append does not;
        # checked here so that one façade means raw scans whichever mode it resolves to.
        for uuid in self.uuids:
            scan = self.tavi_data.fetch_by_uuid(uuid)
            if not isinstance(scan, RawScan):
                raise ValueError(
                    f"Combine can only operate on raw data, but {scan.tavimeta.friendly_name!r} "
                    f"({uuid.value}) is already processed."
                )

        self._resolve_delegate().validate()

    def exec(self) -> ProcessedScan:
        """
        Run the combination the parameters describe and register the result in ``tavi_data.processed_scans``.

        The result is whatever the delegated op built - its columns, its ``error``
        column where it computes one, its ``prov.contributing_scans`` and its friendly
        name - and is returned as well as stored, so a caller need not look it back up.
        """
        self.validate()
        processed_scan = self._resolve_delegate().exec()
        self.register(processed_scan)
        return processed_scan

    def register(self, processed_scan: ProcessedScan) -> None:
        """
        Write a combined scan into ``tavi_data.processed_scans`` under its own uuid.

        Kept apart from ``exec`` so that the store write has one named place: a caller
        holding a scan this op produced earlier - one rebuilt from a project file, say -
        can register it without re-running the combination.
        """
        self.tavi_data.processed_scans[processed_scan.uuid] = processed_scan

    def _resolve_delegate(self) -> ProcessOps:
        """
        Build the op the mode calls for, reading the parameters as they stand.

        ``mode`` alone decides which one, so there is a single place the parameters are
        read as a choice of operation, and asking it here is what refuses a tolerance
        and a step together however this is reached.

        Built afresh on every call rather than cached, so that reassigning a parameter
        between ``validate`` and ``exec`` takes effect as it does on the ops this
        delegates to - a cache keyed on nothing would leave ``mode`` reporting one
        combination while ``exec`` ran another. Construction only stores values.
        """
        mode = self.mode
        if mode is RebinMode.TOLERANCE:
            return ToleranceRebinOps(
                self.tavi_data,
                self.uuids,
                self.columns,
                norm=self.norm,
                tol=self.tol,
                multiplication=self.multiplication,
            )
        if mode is RebinMode.EQUAL_STEP:
            return EqualStepRebinOps(
                self.tavi_data,
                self.uuids,
                self.columns,
                norm=self.norm,
                step=self.step,
                start=self.start,
                finish=self.finish,
                multiplication=self.multiplication,
            )
        return AppendOp(self.tavi_data, self.uuids, self.columns)
