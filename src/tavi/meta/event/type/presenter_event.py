"""Events that Presenters emit."""

from tavi.library.data.fit_entry import FitEntry, FitMember, FitOutcome
from tavi.library.data.plot import Plot, PlotSeries
from tavi.library.data.scan import UUID, RawScan, Scan
from tavi.meta.event.event_interface import Event


class FocusEvent(Event):
    """Event to focus on specific items by UUID."""

    ids: list[UUID]


class StartApplicationEvent(Event):
    """Announce that the application's views are built, so models can push their initial state to them."""

    pass


class FocusRawScanEvent(Event):
    """
    Add raw scans to the focus.

    ``PlotModel`` answers with a ``FocusPlotEvent`` of one preview plot per scan. Like every
    ``Focus*`` event this only adds: a new selection is cleared first by ``ClearFocusEvent``, so a
    scan+plot multiselect publishes this and ``FocusPlotEvent`` side by side without either one
    clobbering the other.
    """

    scans: list[RawScan]


class ClearFocusEvent(Event):
    """
    Announce that whatever was focused no longer is - the first half of every new selection.

    Published by ``TaviProjectModel`` first thing on every ``FocusEvent``, ahead of the ``Focus*``
    events it routes to, so every subscriber has dropped its state for the old selection before
    any of the new one arrives. The stage is part of the focus, so it is cleared with it - no
    ``ClearStageEvent`` follows. It says nothing about what a subscriber should do with that;
    ``FittingPresenter``, for one, resets its panel.
    """

    pass


class FocusFitEvent(Event):
    """
    Add fits to the focus. Like every ``Focus*`` event this only adds - see ``ClearFocusEvent``.

    ``PlotModel`` adds the series each member was fit against, if not already focused, by
    publishing a ``FocusPlotEvent`` - so a fit's data is drawn under its curve. ``PlotterPresenter``
    and ``FittingPresenter`` mark these uuids as pending a curve. The recompute trigger is the
    separate ``RecomputeFitEvent``, published right after this one - keeping them distinct (rather
    than having ``FitModel`` also subscribe to this one) guarantees the UI has already marked its
    pending state by the time ``FitModel``'s resulting ``SyncFitEvent`` arrives, regardless of
    subscriber registration order.
    """

    fits: list[FitEntry]
    scans: dict[UUID, Scan] = {}
    """The Scan every member of every focused fit points at, same contract as ``FocusPlotEvent.scans``."""


class RecomputeFitEvent(Event):
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


class FocusPlotEvent(Event):
    """
    Add plots to the focus. Like every ``Focus*`` event this only adds - see ``ClearFocusEvent``.

    ``scans`` carries the Scan objects every series across ``plots`` points at (by
    ``source_scan_uuid``) — ``RawScan`` today, but any ``Scan`` (e.g. a future
    ``ProcessedScan``) may be referenced. Deep-copied by ``EventBroker.publish`` like every
    other event field. Presenters/views resolve each series against this snapshot — they
    never hold a live handle into a model's scan storage.
    """

    plots: list[Plot]
    scans: dict[UUID, Scan] = {}


class SyncPlotEvent(Event):
    """
    Sync the UI to new content for the plots already focused - same focus, updated series.

    Published by ``PlotModel`` after a field edit or a scan removal, carrying every focused plot
    as it now stands. Kept apart from ``FocusPlotEvent`` because the reactions differ: a refresh
    replaces what's drawn but must not reset anything a new selection would.
    """

    plots: list[Plot]
    scans: dict[UUID, Scan] = {}


class SyncPeakParamsEvent(Event):
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


class SyncBackgroundParamsEvent(Event):
    """
    Announce a heuristic initial guess for the background's slope/intercept, from real data.

    Mirrors ``SyncPeakParamsEvent``, ``source_scan_uuid`` staleness guard included. Kept a
    separate event rather than extra fields on that one, so a peak guess and a background guess
    can never partly overwrite each other's fields in the panel.
    """

    source_scan_uuid: UUID
    slope: float
    intercept: float


class SetFitComponentsVisibleEvent(Event):
    """
    Announce that the fitting panel's "Plot Separately" checkbox was toggled.

    Purely a display preference, so it carries no fit identity: every drawn fit's components
    show or hide together. ``FitCurve`` always carries its components, whether or not they are
    currently shown, so ``PlotterPresenter`` answers this by toggling artists already on the
    canvas rather than asking ``FitModel`` to recompute anything.
    """

    visible: bool


class SyncFitEvent(Event):
    """
    Sync the UI to freshly computed members of a fit - the "sync" step of focus -> calculate -> sync.

    Published by ``FitModel`` both after Perform Fit and after a saved fit is recomputed for display
    (``RecomputeFitEvent``). Display only: ``PlotterPresenter``/``Plot1DView`` render each outcome's
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


class ClearStageEvent(Event):
    """
    Announce that nothing is staged any more. Mirrors ``ClearFocusEvent`` for the stage.

    The stage is the part of the focus that edits and fits apply to. Published by
    ``PlotterPresenter`` before it restages - when "Apply All" is toggled, or when another series
    is picked in the "Current Plot" dropdown while it's off.
    """

    pass


class StageSeriesEvent(Event):
    """
    Add focused series to the stage, by source scan uuid. Like ``Focus*`` events this only adds.

    Staging is the user's choice of what edits and fits apply to: every focused series with
    "Apply All" checked, only the one picked in the "Current Plot" dropdown without it. Published
    by ``PlotterPresenter``; ``PlotModel`` tracks the stage and answers with ``SyncStageEvent``.
    """

    source_scan_uuids: list[UUID]


class SyncStageEvent(Event):
    """
    Sync the UI to the staged series as they now stand, and the scans backing them.

    Published by ``PlotModel``, which tracks the stage: after ``StageSeriesEvent``, when a field edit
    changes a staged series, and when a removal drops one. The first series leads - it is the one
    the plotter's fields, the data tab and the fitting panel show. Empty means nothing is staged.
    """

    series: list[PlotSeries]
    scans: dict[UUID, Scan] = {}


class SyncFitSpecEvent(Event):
    """
    Sync the fitting panel to one saved member of a fit - its spec and last result - without refitting it.

    Published by ``FitModel.sync_fit_spec`` when a series is picked in the "Current Plot" dropdown
    and the panel already knows a fit covering it, so switching between plots is cheap.
    """

    fit_uuid: UUID
    member: FitMember
