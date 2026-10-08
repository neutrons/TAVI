"""New Raw Scan Event module."""

from tavi.library.data.scan import UUID
from tavi.meta.event.event_interface import Event


class AddRawScanEvent(Event):
    """Indicates a new RawScan has been added to the Project."""

    uuid: UUID
    friendly_name: str
    friendly_path: str


class AddPlotEvent(Event):
    """Indicates a new Plot has been added to the Project."""

    uuid: UUID
    friendly_name: str
    friendly_path: str


class AddFitEvent(Event):
    """Indicates a new Fit has been added to the Project. Mirrors AddPlotEvent."""

    uuid: UUID
    friendly_name: str
    friendly_path: str


class RemoveRawScanEvent(Event):
    """Indicates a RawScan has been removed from the Project and is no longer in TaviData."""

    uuid: UUID


class RemovePlotEvent(Event):
    """Indicates a Plot has been removed from the Project and is no longer in TaviData."""

    uuid: UUID


class RemoveFitEvent(Event):
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

    Published by ``TaviProjectModel`` just before the ``RecomputeFitEvent`` that redraws the member,
    the same way ``FocusFitEvent`` precedes it: ``FitWindowPresenter`` marks that member's open
    window as pending here, so the resulting ``SyncFitEvent`` refreshes it regardless of subscriber
    registration order.
    """

    fit_uuid: UUID
    source_scan_uuid: UUID


class SyncRecentProjectsEvent(Event):
    """Update list of recent projects."""

    recent_projects: list[str]
