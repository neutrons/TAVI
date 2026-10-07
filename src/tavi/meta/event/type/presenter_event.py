"""Events that Presenters emit."""

from typing import Optional

from tavi.library.data.fit_entry import FitEntry, FitMember, FitOutcome
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


class ClearFocusEvent(Event):
    """
    Announce that whatever was focused no longer is - the first half of every new selection.

    Published by ``TaviProjectModel`` first thing on every ``FocusEvent``, ahead of the focus events
    it routes to, so every subscriber has dropped its state for the old selection before any
    ``ActivePlotChangedEvent``/``SyncFitEvent`` from the new one arrives. It says nothing about what
    a subscriber should do with that; ``FittingPresenter``, for one, resets its panel and forgets
    which fit covers each scan. ``PlotFocusEvent`` can't serve as the cue - ``PlotModel``
    republishes it for every field edit, which is a refresh of the same focus, not a new one.
    """

    pass


class FitFocusEvent(Event):
    """
    Announce that a list of fits is now focused. Mirrors RawScanFocusEvent/PlotFocusEvent.

    UI-facing only (``PlotterPresenter``/``FittingPresenter``): they clear stale state and mark
    these uuids as pending a curve. The backend recompute trigger is the separate
    ``FitRecomputeEvent``, published right after this one - keeping them distinct (rather than
    having ``FitModel`` also subscribe to this one) guarantees the UI has already marked its
    pending state by the time ``FitModel``'s resulting ``SyncFitEvent`` arrives, regardless of
    subscriber registration order.

    ``exclusive`` is False when this batch's original selection also included raw scans or
    plots (see ``TaviProjectModel._handle_focus_event``) - subscribers should overlay onto
    that still-focused scan/plot canvas rather than clearing it.
    """

    fits: list[FitEntry]
    exclusive: bool = True
    scans: dict[UUID, Scan] = {}
    """The Scan every member of every focused fit points at, same contract as ``PlotFocusEvent.scans``.
    Every member carries the full ``PlotSeries`` it was fit against - scan, x/y columns and
    normalization alike - so this is what lets re-focusing a fit redraw the data underneath its
    curves, and repopulate the Data File tab, instead of leaving them blank."""


class FitRecomputeEvent(Event):
    """Ask FitModel to recompute every member of each of these fits against its source series' current data."""

    fits: list[FitEntry]


class SavePlotEvent(Event):
    """Request to save a plot into the project's data store."""

    plot: Plot


class SaveFitEvent(Event):
    """
    Request to save a freshly performed fit's members into the project's data store.

    Published by ``FitModel`` after Perform Fit, just before the matching ``SyncFitEvent`` - never for
    a recompute, which changes nothing worth saving. ``TaviProjectModel`` merges ``members`` into
    ``TaviData.fits[fit_uuid]`` by scan, so the members this leaves out stay as they were.
    ``FitWindowPresenter`` also takes it as the cue that the following sync is a fresh fit.
    """

    fit_uuid: UUID
    members: list[FitMember]
    """In fit order. Every member is a different series."""


class UndoFitMemberEvent(Event):
    """
    Ask ``TaviProjectModel`` to roll one member of a fit back to the state before its last save.

    Published by ``FittingPresenter`` for the active series only - undo is for dialing one member
    in, so it is offered only while Perform Fit would refit that member alone.
    """

    fit_uuid: UUID
    source_scan_uuid: UUID


class RedoFitMemberEvent(Event):
    """Ask ``TaviProjectModel`` to reapply the member state the last ``UndoFitMemberEvent`` rolled back."""

    fit_uuid: UUID
    source_scan_uuid: UUID


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
    since switched to a different active series - the same staleness concern ``SyncFitEvent``
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


class ShowScanTitleChangedEvent(Event):
    """
    Announce that the plotter's "Show Title" checkbox was toggled.

    A display preference, so it names no scan: every series is labelled the same way.
    ``PlotModel`` owns the answer - ``PlotSeries.scan_name`` is set when a preview plot is
    built - so it relabels what's currently focused and republishes ``PlotFocusEvent`` rather
    than the view relabelling artists it doesn't own the names of.
    """

    show_title: bool


class SyncFitEvent(Event):
    """
    Sync the UI to freshly computed members of a fit - the "sync" step of focus -> calculate -> sync.

    Published by ``FitModel`` both after Perform Fit and after a saved fit is recomputed for display
    (``FitRecomputeEvent``). Display only: ``PlotterPresenter``/``Plot1DView`` render each outcome's
    ``curve`` alongside the data it was fit against, ``FittingPresenter`` reflects the active series'
    member in the fitting panel, and ``FitWindowPresenter`` fills the windows of a fit it was just told
    was saved (``SaveFitEvent``). Persisting the fit is ``SaveFitEvent``'s job, not this one's.

    ``outcomes`` may cover only some of the fit's members - refitting one series of a sequential fit
    syncs just that one - which is why this names the fit by uuid rather than carrying a whole
    ``FitEntry``.
    """

    fit_uuid: UUID
    outcomes: list[FitOutcome]
    """In fit order. Every outcome is a different series."""


class ApplyAllChangedEvent(Event):
    """
    Announce the plotter's "Apply All" checkbox state.

    ``FittingPresenter`` reads it as "fit every focused series, each seeded from the one before" (a
    sequential fit) when checked, versus "fit only the series picked in the Current Plot dropdown"
    when not - the same scope the checkbox already gives the plotter's own field edits.
    """

    apply_all: bool


class SyncFitSpecEvent(Event):
    """
    Sync the fitting panel to one saved member of a fit - its spec and last result - without refitting it.

    Published by ``FitModel.sync_fit_spec`` when a series is picked in the "Current Plot" dropdown
    and the panel already knows a fit covering it, so switching between plots is cheap.
    """

    fit_uuid: UUID
    member: FitMember
