"""New Raw Scan Event module."""

from tavi.library.data.scan import UUID
from tavi.meta.event.event_interface import Event


class RawScanAppendEvent(Event):
    """Indicates a new RawScan has been added to the Project."""

    uuid: UUID
    friendly_name: str
    friendly_path: str


class PlotAppendEvent(Event):
    """Indicates a new Plot has been added to the Project."""

    uuid: UUID
    friendly_name: str
    friendly_path: str


class FitAppendEvent(Event):
    """Indicates a new Fit has been added to the Project. Mirrors PlotAppendEvent."""

    uuid: UUID
    friendly_name: str
    friendly_path: str


class RawScanRemoveEvent(Event):
    """Indicates a RawScan has been removed from the Project and is no longer in TaviData."""

    uuid: UUID


class PlotRemoveEvent(Event):
    """Indicates a Plot has been removed from the Project and is no longer in TaviData."""

    uuid: UUID


class FitRemoveEvent(Event):
    """Indicates a Fit has been removed from the Project and is no longer in TaviData."""

    uuid: UUID


class SyncFitHistoryEvent(Event):
    """
    Announce whether one member of a fit can now be undone or redone.

    Published by ``TaviProjectModel`` whenever a member's history moves - a save pushes a step, an
    undo/redo moves one between stacks - so ``FittingPresenter`` can enable its Undo/Redo buttons
    without asking the model, keeping only these flags per member.
    """

    fit_uuid: UUID
    source_scan_uuid: UUID
    can_undo: bool
    can_redo: bool


class RestoreFitMemberEvent(Event):
    """
    Announce that one member of a fit was rolled back or forward to an earlier saved state.

    Published by ``TaviProjectModel`` just before the ``FitRecomputeEvent`` that redraws the member,
    the same way ``FitFocusEvent`` precedes it: ``FitWindowPresenter`` marks that member's open
    window as pending here, so the resulting ``SyncFitEvent`` refreshes it regardless of subscriber
    registration order.
    """

    fit_uuid: UUID
    source_scan_uuid: UUID


class RawScanLoadingEvent(Event):
    """loading raw data event."""

    raw_scan_uuid: list[str]


class SyncRecentProjects(Event):
    """Update list of recent projects."""

    recent_projects: list[str]
