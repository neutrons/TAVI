"""Model for combining raw scans into one derived scan."""

import math
from typing import Optional

import numpy as np

from tavi.backend.model.interface.combine_model_interface import CombineModelInterface
from tavi.library.data.algorithms.histogram import POSITION_TOL, Histogram, measured_positions
from tavi.library.data.combine import CombineOps
from tavi.library.data.combine_request import CombineRequest
from tavi.library.data.enum.rebin_mode import RebinMode
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.scan import UUID, HistogramData, ProcessedScan
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

    def _build_ops(self, request: CombineRequest) -> CombineOps | Histogram:
        """Translate a request into the op it describes, raising ValueError on anything unusable."""
        if request.rebin_mode is RebinMode.HISTOGRAM:
            return self._build_histogram(request)

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

    def _build_histogram(self, request: CombineRequest) -> Histogram:
        """Translate a request into a ``Histogram`` over ``[x_axis, y_axis, z_axis, normalization]``, the columns as named."""
        x_axis = request.x_axis.strip()
        z_axis = request.z_axis.strip()
        y_axis = request.y_axis.strip()
        normalization = request.normalization.strip()
        if not x_axis or not z_axis or not y_axis:
            raise ValueError("Histogram needs an x-axis, a y-axis and a z-axis column.")
        # every bin is counts over summed norm, so there is no histogram without one.
        if not normalization:
            raise ValueError("Histogram needs a normalization column.")

        # Histogram lays grids by bin count, so the panel's Start/Stop/Step become explicit
        # edges here - which also keeps a typed edge exactly where the user put it.
        x_edges = self._histogram_edges(
            request.source_scan_uuids, x_axis, "X", request.x_start, request.x_stop, request.x_step
        )
        y_edges = self._histogram_edges(
            request.source_scan_uuids, y_axis, "Y", request.y_start, request.y_stop, request.y_step
        )

        return Histogram(
            self._tavi_data,
            request.source_scan_uuids,
            [x_axis, y_axis, z_axis, normalization],
            bins=[x_edges, y_edges],
            multiplication=_optional_number("Multiplication", request.multiplication) or 1.0,
        )

    def _histogram_edges(
        self, uuids: list[UUID], column: str, label: str, start_text: str, stop_text: str, step_text: str
    ) -> np.ndarray:
        """
        Lay one map axis's bin edges from Start/Stop/Step, as the 1D equal-step rebin lays its grid.

        A blank Start is read off the data half a step below its smallest value, centring the first
        bin on it. A typed Stop is an edge to cover: the last edge is the first one at or past it.
        A blank Stop instead ends the grid at the bin holding the largest value.
        """
        start = _optional_number(f"{label} Start", start_text)
        stop = _optional_number(f"{label} Stop", stop_text)
        step = _optional_number(f"{label} Step", step_text)
        if step is None:
            return self._position_edges(uuids, column, label, start, stop)
        if not math.isfinite(step) or step <= 0:
            raise ValueError(f"{label} Step must be a positive, finite number, got {step_text.strip()!r}.")

        largest = None
        if start is None or stop is None:
            values = self._column_values(uuids, column)
            if start is None:
                start = float(values.min()) - step / 2
            if stop is None:
                largest = float(values.max())
                stop = largest + step / 2
        if stop <= start:
            raise ValueError(f"{label} Stop must be greater than {label} Start, got start={start!r}, stop={stop!r}.")

        if largest is not None:
            # Covering the padded max + step/2 as if it were an edge adds a bin whenever the data's
            # span is not a whole number of steps - readback jitter (l at 1.6499 and 1.8001, step
            # 0.15) makes that nearly always - and that bin is empty. Stop at the one holding max.
            count = math.floor((largest - start) / step + 1e-9) + 1
        else:
            # float noise can put a whole number of steps just past an integer (0.15 / 0.05 is
            # 3.0000000000000004), which a bare ceil would turn into an extra, empty bin.
            count = math.ceil((stop - start) / step - 1e-9)
        return start + step * np.arange(max(count, 1) + 1)

    def _position_edges(
        self, uuids: list[UUID], column: str, label: str, start: Optional[float], stop: Optional[float]
    ) -> np.ndarray:
        """
        Lay one bin per measured position along one map axis, for a Step left blank.

        Values within ``POSITION_TOL`` are one position, so readback jitter does not split a
        setpoint. Each inner edge sits midway between neighbouring positions - so unevenly spaced
        setpoints still get a bin each - and each outer edge mirrors its neighbour, centring the
        outer bins. A typed Start or Stop limits which values count and replaces that outer edge.
        """
        if start is not None and stop is not None and stop <= start:
            raise ValueError(f"{label} Stop must be greater than {label} Start, got start={start!r}, stop={stop!r}.")
        values = self._column_values(uuids, column)
        if start is not None:
            values = values[values >= start]
        if stop is not None:
            values = values[values <= stop]
        if values.size == 0:
            raise ValueError(f"{label} Start/Stop hold none of the {column!r} values.")

        positions = measured_positions(values, POSITION_TOL)
        inner = [(high + next_low) / 2 for (_, high), (next_low, _) in zip(positions, positions[1:])]
        (first_low, first_high), (last_low, last_high) = positions[0], positions[-1]
        if inner:
            low_edge = first_low - (inner[0] - first_high)
            high_edge = last_high + (last_low - inner[-1])
        else:
            # one position leaves no spacing to size its bin by; numpy's own choice for a single value.
            centre = (first_low + first_high) / 2
            low_edge, high_edge = centre - 0.5, centre + 0.5
        if start is not None:
            low_edge = start
        if stop is not None:
            high_edge = stop
        return np.array([low_edge, *inner, high_edge])

    def _column_values(self, uuids: list[UUID], column: str) -> np.ndarray:
        """Every finite value of ``column`` across the origins, for reading a blank grid edge off the data."""
        if not uuids:
            raise ValueError("Histogram needs at least one scan to bin.")
        arrays = []
        for uuid in uuids:
            scan = self._tavi_data.fetch_by_uuid(uuid)
            if column not in scan.data.data:
                raise KeyError(f"Column {column!r} is not in {scan.tavimeta.friendly_name!r}.")
            arrays.append(np.asarray(scan.data.data[column], dtype=float))
        values = np.concatenate(arrays)
        values = values[np.isfinite(values)]
        if values.size == 0:
            raise ValueError(f"Histogram column {column!r} has no finite values to lay a grid on.")
        return values

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

    def _publish_add_processed_scan(self, scan: ProcessedScan | HistogramData) -> None:
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
