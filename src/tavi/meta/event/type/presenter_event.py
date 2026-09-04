"""Events that Presenters emit."""

from typing import Optional

from tavi.library.data.plot import Plot, PlotSeries
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
    Request to make one already-focused series active, by its source scan's uuid.

    Handled by both ``TaviProjectModel`` and ``PlotModel`` — each searches its own focused plots
    (``TaviData.plots`` vs. an unsaved preview in ``PlotModel._last_plots``) for a series whose
    ``source_scan_uuid`` matches, and no-ops otherwise. A series' source scan is what identifies
    it - not its containing ``Plot``'s uuid - so one entry can be picked out of an otherwise-fused,
    multi-series saved plot exactly as it would be among several single-series preview plots. The
    publisher does not need to know which model currently owns the matching series.
    """

    uuid: UUID


class ActivePlotChangedEvent(Event):
    """
    Event announcing the scan (and series) backing whichever single series is currently "active".

    Selected via the plotter's "Current Plot" dropdown - one entry per series, not per Plot, so a
    fused multi-series plot still offers each of its series individually. Carries the ``Scan``
    itself rather than the ``Plot`` and a snapshot to resolve it against, since a ``Plot`` may be
    an unsaved preview with nowhere persistent to live; consumers that only display data (e.g.
    the data widget) care about the scan, not the plot's save state. ``None`` means no series is
    currently active.

    ``series`` is carried alongside the scan so the plotter can resync its own axis/preset fields
    to whichever series just became active - a plain default-axis scan lookup wouldn't reflect a
    per-series edit (e.g. "Apply All" off).
    """

    scan: Optional[Scan] = None
    series: Optional[PlotSeries] = None
