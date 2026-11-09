"""New Raw Scan Event module."""

from tavi.library.data.filter import ScanFilter
from tavi.library.data.scan import UUID
from tavi.meta.event.event_interface import Event


class AddRawScanEvent(Event):
    """Indicates a RawScan is listed in the Project: newly loaded, or passing a newly applied filter."""

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


class AddProcessedScanEvent(Event):
    """
    Indicates a new ProcessedScan - a combination of raw scans - has been added to the Project.

    Published by ``CombineModel`` once the combination is registered in ``TaviData.processed_scans``.
    Mirrors ``AddPlotEvent``: the uuid and the name to list it under, never the scan itself.
    """

    uuid: UUID
    friendly_name: str
    friendly_path: str


class AddFilterEvent(Event):
    """
    Indicates ``FilterModel`` has saved a new filter to the Project. Filters are not listed in the project tree.

    Carries the whole filter, not just its name, because the filter panel shows its criteria too.
    """

    filter: ScanFilter


class SyncFilterEvent(Event):
    """
    Announce how a saved filter stands after an edit to it was saved into the Project.

    Published by ``FilterModel`` when ``save_filter`` replaces a filter already saved; a new filter is
    announced with ``AddFilterEvent`` instead. The filter panel uses it to drop its unsaved-changes marker.
    """

    filter: ScanFilter


class ClearRawScansEvent(Event):
    """
    Announce that no raw scan is listed any more - the first half of relisting them, paired with ``AddRawScanEvent``.

    Published by ``TaviProjectModel`` when the applied filters change, after ``ClearFocusEvent`` and before
    one ``AddRawScanEvent`` per raw scan that passes them. The scans stay in ``TaviData``: nothing has left
    the project, which is what ``RemoveRawScanEvent`` would mean.
    """

    pass


class RemoveRawScanEvent(Event):
    """Indicates a RawScan has been removed from the Project and is no longer in TaviData."""

    uuid: UUID


class RemovePlotEvent(Event):
    """Indicates a Plot has been removed from the Project and is no longer in TaviData."""

    uuid: UUID


class RemoveFitEvent(Event):
    """Indicates a Fit has been removed from the Project and is no longer in TaviData."""

    uuid: UUID


class RemoveFilterEvent(Event):
    """Indicates ``FilterModel`` has removed a filter from the Project and it is no longer in TaviData."""

    uuid: UUID


class SyncItemNameEvent(Event):
    """
    Announce that an item already listed in the Project is now known by a different name.

    Published by ``TaviProjectModel.rename_item`` once the new name is stored. One event covers
    plots and fits alike: a renamed row is a renamed row, and the tree relabels it the same way
    whichever it is. The tree waits for this rather than relabelling the row itself, so it can
    never show a name the project does not hold.
    """

    uuid: UUID
    friendly_name: str


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


class ApplyFilterEvent(Event):
    """
    Add a filter to those narrowing which scans the project lists.

    Published by ``FilterModel``. Like ``Focus*`` events this only adds: every applied filter narrows
    the list further, so replacing the applied filter is a ``ClearFilterEvent`` followed by this one.
    ``TaviProjectModel`` answers by relisting the raw scans that pass - see ``ClearRawScansEvent``.
    """

    filter: ScanFilter


class ClearFilterEvent(Event):
    """
    Announce that no filter narrows the project's scan list any more. Mirrors ``ClearFocusEvent`` for filters.

    Published by ``FilterModel`` when the applied filter is cleared, removed from the project, or about to
    be replaced. ``TaviProjectModel`` answers by relisting every raw scan - see ``ClearRawScansEvent``.
    """

    pass


class SyncFilterValidationEvent(Event):
    """
    Announce why a filter cannot be applied - or, with ``errors`` empty, that it can.

    Published by ``FilterModel`` for ``validate_filter``, and for ``apply_filter`` when the selected filter
    is invalid and so was not applied.
    """

    filter_uuid: UUID
    errors: list[str]
