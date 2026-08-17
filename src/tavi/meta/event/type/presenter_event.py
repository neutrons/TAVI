"""Events that Presenters emit."""

from typing import Optional

from tavi.library.data.plot import Plot
from tavi.library.data.scan import UUID, RawScan, Scan
from tavi.meta.event.event_interface import Event


class FocusEvent(Event):
    """Event to focus on specific items by UUID."""

    ids: list[UUID]


class DownstreamReadyEvent(Event):
    """Notify upstream that consumers are ready at startup."""

    pass


class RawScanFocusEvent(Event):
    """Event to plot a list of raw scans."""

    scans: list[RawScan]


class SavePlotEvent(Event):
    """Request to save a plot into the project's data store."""

    plot: Plot


class PlotFocusEvent(Event):
    """
    Event to render a list of plots.

    ``scans`` carries the Scan objects every series across ``plots`` points at (by
    ``source_scan_uuid``) — ``RawScan`` today, but any ``Scan`` (e.g. a future
    ``ProcessedScan``) may be referenced. Deep-copied by ``EventBroker.publish`` like every
    other event field. Presenters/views resolve each series against this snapshot — they
    never hold a live handle into a model's scan storage.
    """

    plots: list[Plot]
    scans: dict[UUID, Scan] = {}


class FocusActivePlotEvent(Event):
    """
    Request to switch which currently-focused plot is "active", by uuid.

    Unlike ``FocusEvent``, this does not ask a model to re-resolve or re-broadcast the
    whole focused batch (no ``PlotFocusEvent``, no view re-render) — only the single plot
    named by ``uuid`` is looked up, and the result is announced via ``ActivePlotChangedEvent``.
    """

    uuid: UUID


class ActivePlotChangedEvent(Event):
    """
    Event announcing which single plot is currently "active" (selected via the plotter's plot dropdown).

    Published by the plotter presenter both on a fresh ``PlotFocusEvent`` (defaulting to the first
    plot) and whenever the dropdown selection changes. ``scans`` is the same snapshot carried by the
    triggering ``PlotFocusEvent``, so consumers (e.g. the data widget) can resolve the active plot's
    series without reaching into any model's live storage.
    """

    plot: Optional[Plot] = None
    scans: dict[UUID, Scan] = {}
