"""Events that Presenters emit."""

from typing import Optional

from tavi.library.data.fit_entry import FitCurve, FitEntry, FitResultSummary
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
    """
    Event to plot a list of raw scans.

    ``also_plots`` carries any saved ``Plot``s focused in the same tree multiselect (see
    ``TaviProjectModel._handle_focus_event``) - ``PlotModel`` folds them into the same preview
    batch and publishes one merged ``PlotFocusEvent``, so a scan+plot multiselect overlays
    both instead of one clobbering the other.
    """

    scans: list[RawScan]
    also_plots: list[Plot] = []


class FitFocusEvent(Event):
    """
    Announce that a list of fits is now focused. Mirrors RawScanFocusEvent/PlotFocusEvent.

    UI-facing only (``PlotterPresenter``/``FittingPresenter``): they clear stale state and mark
    these uuids as pending a curve. The backend recompute trigger is the separate
    ``FitRecomputeEvent``, published right after this one - keeping them distinct (rather than
    having ``FitModel`` also subscribe to this one) guarantees the UI has already marked its
    pending state by the time ``FitModel``'s resulting ``FitComputedEvent`` arrives, regardless of
    subscriber registration order.

    ``exclusive`` is False when this batch's original selection also included raw scans or
    plots (see ``TaviProjectModel._handle_focus_event``) - subscribers should overlay onto
    that still-focused scan/plot canvas rather than clearing it.
    """

    fits: list[FitEntry]
    exclusive: bool = True
    scans: dict[UUID, Scan] = {}
    """The Scan each focused fit's own series points at, same contract as ``PlotFocusEvent.scans``.
    Every fit carries the full ``PlotSeries`` it was fit against - scan, x/y columns and
    normalization alike - so this is what lets re-focusing a fit redraw the data underneath its
    curve, and repopulate the Data File tab, instead of leaving them blank."""


class FitRecomputeEvent(Event):
    """Ask FitModel to recompute each of these fits against its source series' current data."""

    fits: list[FitEntry]


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


class PeakParamsSuggestedEvent(Event):
    """
    Announce a heuristic initial guess for one peak's amplitude/center/FWHM, from real data.

    ``source_scan_uuid`` lets the presenter drop a stale reply that arrives after the user has
    since switched to a different active series - the same staleness concern ``FitComputedEvent``
    guards against.
    """

    source_scan_uuid: UUID
    amplitude: float
    center: float
    fwhm: float


class BackgroundParamsSuggestedEvent(Event):
    """
    Announce a heuristic initial guess for the background's slope/intercept, from real data.

    Mirrors ``PeakParamsSuggestedEvent``, ``source_scan_uuid`` staleness guard included. Kept a
    separate event rather than extra fields on that one, so a peak guess and a background guess
    can never partly overwrite each other's fields in the panel.
    """

    source_scan_uuid: UUID
    slope: float
    intercept: float


class FitComponentsVisibilityChangedEvent(Event):
    """
    Announce that the fitting panel's "Plot Separately" checkbox was toggled.

    Purely a display preference, so it carries no fit identity: every drawn fit's components
    show or hide together. ``FitCurve`` always carries its components, whether or not they are
    currently shown, so ``PlotterPresenter`` answers this by toggling artists already on the
    canvas rather than asking ``FitModel`` to recompute anything.
    """

    visible: bool


class FitComputedEvent(Event):
    """
    Announce a fit result, whether freshly performed or recomputed after being reselected.

    Two independent subscribers act on this: ``PlotterPresenter``/``Plot1DView`` render ``curve``
    alongside the data it was fit against and ``FittingPresenter`` reflects ``result`` in the
    fitting panel; ``TaviProjectModel`` persists ``fit`` (the spec only - never ``curve``) into
    ``TaviData.fits`` and announces it in the project tree via ``FitAppendEvent`` the first time
    it sees this uuid. Unlike ``SavePlotEvent`` (one consumer), this one intentionally has several.
    """

    fit: FitEntry
    curve: FitCurve
    result: FitResultSummary
