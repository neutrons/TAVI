"""Model for combining raw scans into one derived scan."""

from typing import Optional

from tavi.backend.model.interface.combine_model_interface import CombineModelInterface
from tavi.library.data.combine import CombineOps
from tavi.library.data.combine_request import CombineRequest
from tavi.library.data.enum.rebin_mode import RebinMode
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.scan import ProcessedScan
from tavi.library.data.tavi_data import TaviData
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.exception_event import ReportErrorEvent
from tavi.meta.event.type.model_event import AddProcessedScanEvent
from tavi.meta.exception.nonrecoverable.base import NonRecoverableError


class CombineModel(CombineModelInterface):
    """
    Turns a user's combination request into a ``CombineOps`` run, and announces what it produced.

    The combining itself belongs to the library: ``CombineOps`` resolves which of append,
    tolerance rebin and equal-step rebin the parameters describe, refuses the parameters that
    mode has no use for, and registers its result in ``TaviData.processed_scans``. This model's
    own work is the translation either side of that - the typed text of a request into numbers
    and a column list, and a failure anywhere in the run into one message the user sees.

    It holds the project's whole ``TaviData`` rather than a store handle, because a combination
    resolves its origins across every store and writes the result back into one.
    """

    def __init__(self, tavi_data: TaviData) -> None:
        """Initialize with the project's TaviData - the scans a combination draws from and registers into."""
        super().__init__()

        self._tavi_data = tavi_data
        self._event_broker = EventBroker()

    def combine(self, request: CombineRequest) -> ModelResponse:
        """Run the combination ``request`` describes, register the result, and announce it."""
        try:
            ops = self._build_ops(request)
            processed_scan = ops.exec()
        except (ValueError, KeyError) as error:
            # Every refusal the ops raise is about what the user asked for - an unusable number, a
            # column no origin carries, a rebin parameter the chosen mode has no use for - so each
            # is reported to them rather than raised into a worker thread nothing is watching.
            return self._report_error(_message_of(error))

        self._publish_add_processed_scan(processed_scan)
        return ModelResponse(code=ResponseCode.OK)

    def _build_ops(self, request: CombineRequest) -> CombineOps:
        """Translate a request into the op it describes, raising ValueError on anything unusable."""
        tolerance = step = start = finish = None
        # The mode decides which numbers are read at all, so a stale value left in a field the
        # user is not using can neither fail the parse nor reach CombineOps as a second mode.
        if request.rebin_mode is RebinMode.TOLERANCE:
            tolerance = _number("Tolerance", request.tolerance)
        elif request.rebin_mode is RebinMode.EQUAL_STEP:
            step = _number("Step", request.step)
            start = _optional_number("Start", request.start)
            finish = _optional_number("Stop", request.finish)

        return CombineOps(
            self._tavi_data,
            request.source_scan_uuids,
            self._columns(request),
            norm=request.normalization.strip(),
            tol=tolerance,
            step=step,
            start=start,
            finish=finish,
            multiplication=_optional_number("Multiplication", request.multiplication) or 1.0,
        )

    @staticmethod
    def _columns(request: CombineRequest) -> list[str]:
        """The columns to carry into the result: x first, y second, then the normalization column."""
        x_axis = request.x_axis.strip()
        y_axis = request.y_axis.strip()
        if not x_axis or not y_axis:
            raise ValueError("Combine needs both an x-axis and a y-axis column.")

        columns = [x_axis, y_axis]
        # The append drops every column it was not asked for, so a norm column left out here would
        # be gone by the time a rebin needs to weigh by it.
        normalization = request.normalization.strip()
        if normalization and normalization not in columns:
            columns.append(normalization)
        return columns

    def _publish_add_processed_scan(self, scan: ProcessedScan) -> None:
        self._event_broker.publish(
            AddProcessedScanEvent(
                uuid=scan.uuid,
                friendly_name=scan.tavimeta.friendly_name,
                friendly_path=scan.tavimeta.friendly_path,
            )
        )

    def _report_error(self, message: str) -> ModelResponse:
        """Surface a combination failure to the user instead of failing silently."""
        self._event_broker.publish(ReportErrorEvent(error=NonRecoverableError(message, "")))
        return ModelResponse(code=ResponseCode.ERROR, message=message)


def _message_of(error: Exception) -> str:
    """The message an exception carries, without the quotes ``str()`` puts around a KeyError's."""
    return str(error.args[0]) if error.args else str(error)


def _number(label: str, text: str) -> float:
    """Parse a field the chosen mode requires, naming the field if it is not a number."""
    value = _optional_number(label, text)
    if value is None:
        raise ValueError(f"{label} is required for this rebin.")
    return value


def _optional_number(label: str, text: str) -> Optional[float]:
    """Parse a field that may be left blank, naming the field if what is there is not a number."""
    text = text.strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        raise ValueError(f"{label} '{text}' is not numeric.") from None
