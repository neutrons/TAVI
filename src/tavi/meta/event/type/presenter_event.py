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
    Request to make one already-focused plot active, by uuid.

    Handled by both ``TaviProjectModel`` and ``PlotModel`` — each checks whether ``uuid`` is one
    of its own (``TaviData.plots`` vs. an unsaved preview in ``PlotModel._last_plots``) and
    no-ops otherwise. Saved-plot uuids (fresh ``uuid4()`` on save) and preview-plot uuids (the
    same uuid as the ``RawScan`` they preview) never collide, so exactly one model ever acts on
    a given uuid. The publisher does not need to know which kind of plot is currently focused.
    """

    uuid: UUID


class ActivePlotChangedEvent(Event):
    """
    Event announcing the scan backing whichever single plot is currently "active".

    Selected via the plotter's plot dropdown. Carries the ``Scan`` itself — the active plot's first contributing scan (see
    ``first_contributing_scan``) — rather than the ``Plot`` and a snapshot to resolve it against.
    A ``Plot`` may be an unsaved preview with nowhere persistent to live; consumers that only
    display data (e.g. the data widget) care about the scan, not the plot's save state, and a
    plot's first contributing scan always exists in the current setup. ``None`` means no plot is
    currently active.
    """

    scan: Optional[Scan] = None
