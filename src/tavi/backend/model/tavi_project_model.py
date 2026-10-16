"""Tavi Project."""

import logging
from collections.abc import Callable
from typing import Optional

from neutrons_standard.config import Resource
from neutrons_standard.decorators.singleton import Singleton
from ruamel.yaml import YAML

from tavi.backend.model.fit_history import FitHistory
from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface
from tavi.backend.model.plot_resolver import scans_for_plots
from tavi.library.data.fit_entry import FitEntry, FitMember
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.plot import Plot
from tavi.library.data.scan import UUID, RawScan
from tavi.library.data.tavi_data import TaviData
from tavi.library.storage.controller.raw_scan_load_controller import RawScanLoadController
from tavi.library.storage.interface.filestore_interface import Filestore
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import (
    AddFitEvent,
    AddPlotEvent,
    AddRawScanEvent,
    RemoveFitEvent,
    RemovePlotEvent,
    RemoveRawScanEvent,
    RestoreFitMemberEvent,
    SyncFitHistoryEvent,
    SyncRecentProjectsEvent,
)
from tavi.meta.event.type.presenter_event import (
    ClearFocusEvent,
    FocusEvent,
    FocusFitEvent,
    FocusPlotEvent,
    FocusRawScanEvent,
    RecomputeFitEvent,
    SaveFitEvent,
    SavePlotEvent,
    StartApplicationEvent,
)

logger = logging.getLogger(__name__)


@Singleton
class TaviProjectModel(TaviProjectInterface):
    """Tavi project class."""

    def __init__(self, filestore: Filestore) -> None:
        """Init tavi data."""
        self.filestore = filestore
        self.tavi_data: TaviData = TaviData(raw_scans={}, plots={}, fits={})
        self._event_broker: EventBroker = EventBroker()
        self.raw_scan_load_controller: RawScanLoadController = RawScanLoadController()
        # Session-only: earlier states of each fit member, for undo/redo. Kept beside TaviData
        # rather than in it, so dialing a fit in never grows the saved project.
        self._fit_history = FitHistory()

        self._event_broker.register(StartApplicationEvent, self.sync_on_ready)
        self._event_broker.register(FocusEvent, self._handle_focus_event)
        self._event_broker.register(SavePlotEvent, self._handle_save_plot_event)
        self._event_broker.register(SaveFitEvent, self._handle_save_fit_event)

    def get_plots_handle(self) -> dict:
        """Return reference to the plots dict."""
        return self.tavi_data.plots

    def get_raw_scans_handle(self) -> dict:
        """Return reference to the raw_scans dict."""
        return self.tavi_data.raw_scans

    def get_fits_handle(self) -> dict:
        """Return reference to the fits dict."""
        return self.tavi_data.fits

    def load_raw_scan_from_folder(self, folder: str) -> ModelResponse:
        """Load a folder containing raw scans, skipping any whose file content is already in the project."""
        raw_scans: list[RawScan] = self.raw_scan_load_controller.load_folder(folder)
        events = []
        for scan in raw_scans:
            # A raw scan's uuid is the md5 of its file text: an identical uuid is an unchanged
            # file, so the stored copy is kept (preserving edits to its writable tavimeta). A
            # file grown since the last load hashes differently and comes in as a new scan.
            if scan.uuid in self.tavi_data.raw_scans:
                logger.info(f"Skipping already-loaded scan {scan.tavimeta.friendly_name!r} ({scan.uuid.value}).")
                continue
            self.tavi_data.raw_scans[scan.uuid] = scan
            events.append(
                AddRawScanEvent(
                    uuid=scan.uuid, friendly_name=scan.tavimeta.friendly_name, friendly_path=scan.tavimeta.friendly_path
                )
            )

        for event in events:
            self._event_broker.publish(event)

        return ModelResponse(code=ResponseCode.OK)

    def remove_items(self, uuids: list[UUID]) -> ModelResponse:
        """Drop the named items, and whatever they orphan, from the project and announce each removal."""
        purged = self.tavi_data.purge(uuids)
        self._fit_history.forget(fit_uuids=set(purged.fits), scan_uuids=set(purged.raw_scans))

        for uuid in purged.raw_scans:
            self._event_broker.publish(RemoveRawScanEvent(uuid=uuid))
        for uuid in purged.plots:
            self._event_broker.publish(RemovePlotEvent(uuid=uuid))
        for uuid in purged.fits:
            self._event_broker.publish(RemoveFitEvent(uuid=uuid))

        return ModelResponse(code=ResponseCode.OK)

    def sync_on_ready(self, _: StartApplicationEvent) -> None:
        """Sync with downstream when its ready."""
        self.emit_sync_recent_projects()

    def emit_sync_recent_projects(self) -> None:
        """Notify consumers of latest recent projects."""
        recent_projects = self._get_recent_projects()
        e = SyncRecentProjectsEvent(recent_projects=recent_projects)
        self._event_broker.publish(e)

    def _get_recent_projects(self) -> list[str]:
        # TODO: Demo purposes only. Remove this line when settings.yaml is actually used.
        self.filestore.write_user_data_file("settings.yaml", Resource.read("default_settings.yml"))
        raw_settings_yml = self.filestore.read_user_data_file("settings.yaml")
        yaml = YAML()
        settings = yaml.load(raw_settings_yml)
        settings_dict = dict(settings)
        return settings_dict["TAVI"]["recent"]["projects"]

    def _handle_save_plot_event(self, e: SavePlotEvent) -> None:
        """Record a presenter-submitted plot in ``tavi_data`` and announce it."""
        self.tavi_data.plots[e.plot.uuid] = e.plot
        run_names = "_".join(series.run_name for series in e.plot.series)
        friendly_name = f"{run_names}_Plot"
        self._event_broker.publish(AddPlotEvent(uuid=e.plot.uuid, friendly_name=friendly_name, friendly_path=""))

    def _handle_save_fit_event(self, e: SaveFitEvent) -> None:
        """
        Merge a fit's members into ``tavi_data`` and announce it, mirroring ``_handle_save_plot_event``.

        Each member replaces the stored member for the same scan, so refitting one series of a
        sequential fit leaves its other members as they were. The state each one replaces becomes
        an undo step for that member. Only a new fit is announced to the project tree - refitting an
        existing one would otherwise add its uuid to the tree twice.
        """
        existing = self.tavi_data.fits.get(e.fit_uuid)
        if existing is not None:
            for member in e.members:
                previous = existing.member_for(member.source_scan_uuid)
                # A failed fit has nothing on screen worth going back to.
                if previous is not None and previous.result is not None:
                    self._fit_history.record(e.fit_uuid, previous)
            self.tavi_data.fits[e.fit_uuid] = existing.with_members(e.members)
        else:
            fit = FitEntry(uuid=e.fit_uuid, name=self._fit_name(e.members), members=e.members)
            self.tavi_data.fits[fit.uuid] = fit
            self._event_broker.publish(AddFitEvent(uuid=fit.uuid, friendly_name=fit.name, friendly_path=""))
        for member in e.members:
            self._publish_fit_history(e.fit_uuid, member.source_scan_uuid)

    def undo_fit_member(self, fit_uuid: UUID, source_scan_uuid: UUID) -> ModelResponse:
        """Roll one member of a fit back to the state before its last save."""
        self._step_fit_history(fit_uuid, source_scan_uuid, self._fit_history.undo)
        return ModelResponse(code=ResponseCode.OK)

    def redo_fit_member(self, fit_uuid: UUID, source_scan_uuid: UUID) -> ModelResponse:
        """Reapply the member state the last undo rolled back."""
        self._step_fit_history(fit_uuid, source_scan_uuid, self._fit_history.redo)
        return ModelResponse(code=ResponseCode.OK)

    def _step_fit_history(
        self, fit_uuid: UUID, source_scan_uuid: UUID, step: Callable[[UUID, FitMember], Optional[FitMember]]
    ) -> None:
        """
        Swap one member for the state ``step`` returns, then redraw just that member.

        The restored member is recomputed from its spec rather than replayed, the same as selecting
        a saved fit in the tree - TAVI never caches fitted curves.
        """
        fit = self.tavi_data.fits.get(fit_uuid)
        current = fit.member_for(source_scan_uuid) if fit is not None else None
        if current is None:
            return
        restored = step(fit_uuid, current)
        if restored is None:
            return
        fit = fit.with_members([restored])
        self.tavi_data.fits[fit_uuid] = fit
        self._publish_fit_history(fit_uuid, source_scan_uuid)
        self._event_broker.publish(RestoreFitMemberEvent(fit_uuid=fit_uuid, source_scan_uuid=source_scan_uuid))
        self._event_broker.publish(RecomputeFitEvent(fits=[fit.model_copy(update={"members": [restored]})]))

    def _publish_fit_history(self, fit_uuid: UUID, source_scan_uuid: UUID) -> None:
        key = (fit_uuid, source_scan_uuid)
        self._event_broker.publish(
            SyncFitHistoryEvent(
                fit_uuid=fit_uuid,
                source_scan_uuid=source_scan_uuid,
                can_undo=self._fit_history.can_undo(key),
                can_redo=self._fit_history.can_redo(key),
            )
        )

    def _fit_name(self, members: list[FitMember]) -> str:
        """Name a new fit after the run it covers, or the first and last runs of a sequential one."""
        first, last = members[0].series.run_name, members[-1].series.run_name
        if len(members) == 1:
            return f"{first}_Fit"
        return f"{first}-{last}_Fit"

    def _handle_focus_event(self, e: FocusEvent) -> None:
        """
        Clear the old focus, then route a ``FocusEvent``'s ids to one additive ``Focus*`` event per item type.

        Only ``TaviData`` knows what type each tree uuid is, which is why the tree publishes one
        generic ``FocusEvent`` and this model narrows it.
        """
        # First, so every subscriber has let go of the old selection before this one's chain reaches it.
        self._event_broker.publish(ClearFocusEvent())
        ids = e.ids
        raw_scans: list[RawScan] = []
        plots: list[Plot] = []
        fits: list[FitEntry] = []
        for uuid in ids:
            # FocusEvent ids come from the project tree, which only ever lists uuids
            # TaviData actually owns — fetch_by_uuid raising here means tree/TaviData are
            # out of sync, a bug worth surfacing loudly rather than silently dropping the item.
            inst = self.tavi_data.fetch_by_uuid(uuid)
            if isinstance(inst, RawScan):
                raw_scans.append(inst)
            if isinstance(inst, Plot):
                plots.append(inst)
            if isinstance(inst, FitEntry):
                fits.append(inst)

        # A focused Plot's own attached fits (stamped on it by PlotModel.save_focused_plots when
        # it was saved) are re-focused too, even though the user only selected the Plot itself -
        # this is how re-selecting a saved plot brings its fit curves back, not just its data.
        attached_fits = {fit.uuid: fit for fit in fits}
        for plot in plots:
            for fit_uuid in plot.fits:
                if fit_uuid not in attached_fits and fit_uuid in self.tavi_data.fits:
                    attached_fits[fit_uuid] = self.tavi_data.fits[fit_uuid]
        fits = list(attached_fits.values())

        if raw_scans:
            self._event_broker.publish(FocusRawScanEvent(scans=raw_scans))
        if plots:
            scans = scans_for_plots(plots, self.tavi_data.raw_scans)
            self._event_broker.publish(FocusPlotEvent(plots=plots, scans=scans))
        if fits:
            # FocusFitEvent (UI state) is published fully - every subscriber done - before
            # RecomputeFitEvent (backend trigger), so PlotterPresenter/FittingPresenter have
            # already marked these uuids pending by the time FitModel's resulting
            # SyncFitEvent arrives - see FocusFitEvent's docstring.
            # Carry each fit's own source scan along, so the presenter can redraw the data the
            # fit was made against without reaching into this model's storage.
            fit_scans = {
                member.source_scan_uuid: self.tavi_data.raw_scans[member.source_scan_uuid]
                for fit in fits
                for member in fit.members
                if member.source_scan_uuid in self.tavi_data.raw_scans
            }
            self._event_broker.publish(FocusFitEvent(fits=fits, scans=fit_scans))
            self._event_broker.publish(RecomputeFitEvent(fits=fits))
