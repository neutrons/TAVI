"""Tavi Project."""

from neutrons_standard.config import Resource
from neutrons_standard.decorators.singleton import Singleton
from ruamel.yaml import YAML

from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface
from tavi.backend.model.plot_resolver import find_series_by_source, scans_for_plots
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.plot import Plot
from tavi.library.data.scan import UUID, RawScan
from tavi.library.data.tavi_data import TaviData
from tavi.library.storage.controller.raw_scan_load_controller import RawScanLoadController
from tavi.library.storage.interface.filestore_interface import Filestore
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import (
    PlotAppendEvent,
    PlotRemoveEvent,
    RawScanAppendEvent,
    RawScanRemoveEvent,
    SyncRecentProjects,
)
from tavi.meta.event.type.presenter_event import (
    ActivePlotChangedEvent,
    DownstreamReadyEvent,
    FocusActivePlotEvent,
    FocusEvent,
    PlotFocusEvent,
    RawScanFocusEvent,
    SavePlotEvent,
)


@Singleton
class TaviProjectModel(TaviProjectInterface):
    """Tavi project class."""

    def __init__(self, filestore: Filestore) -> None:
        """Init tavi data."""
        self.filestore = filestore
        self.tavi_data: TaviData = TaviData(raw_scans={}, plots={})
        self._event_broker: EventBroker = EventBroker()
        self.raw_scan_load_controller: RawScanLoadController = RawScanLoadController()

        self._event_broker.register(DownstreamReadyEvent, self.sync_on_ready)
        self._event_broker.register(FocusEvent, self._handle_focus_event)
        self._event_broker.register(FocusActivePlotEvent, self._handle_active_plot_focus_event)
        self._event_broker.register(SavePlotEvent, self._handle_save_plot_event)

    def get_plots_handle(self) -> dict:
        """Return reference to the plots dict."""
        return self.tavi_data.plots

    def get_raw_scans_handle(self) -> dict:
        """Return reference to the raw_scans dict."""
        return self.tavi_data.raw_scans

    def load_raw_scan_from_folder(self, folder: str) -> ModelResponse:
        """Load a folder containing raw scans."""
        raw_scans: list[RawScan] = self.raw_scan_load_controller.load_folder(folder)
        events = []
        for scan in raw_scans:
            self.tavi_data.raw_scans[scan.uuid] = scan
            events.append(
                RawScanAppendEvent(
                    uuid=scan.uuid, friendly_name=scan.tavimeta.friendly_name, friendly_path=scan.tavimeta.friendly_path
                )
            )

        for event in events:
            self._event_broker.publish(event)

        return ModelResponse(code=ResponseCode.OK)

    def remove_items(self, uuids: list[UUID]) -> ModelResponse:
        """
        Drop raw scans and plots from the project and announce each removal.

        A saved plot holds no data of its own, only a ``source_scan_uuid`` per series, so a
        series whose scan is being removed can no longer be resolved and is dropped with it.
        The plot itself survives as long as it has a series left — removing one run should not
        destroy the rest of a fused, multi-series plot. Unknown uuids are ignored: the tree may
        ask twice for the same item (e.g. a folder and a scan inside it both selected).
        """
        requested = list(dict.fromkeys(uuids))
        removed_scans = [uuid for uuid in requested if uuid in self.tavi_data.raw_scans]
        removed_plots = [uuid for uuid in requested if uuid in self.tavi_data.plots]

        orphaned = set(removed_scans)
        pruned_plots = {}
        for plot_uuid, plot in self.tavi_data.plots.items():
            if plot_uuid in removed_plots:
                continue
            surviving = [series for series in plot.series if series.source_scan_uuid not in orphaned]
            if len(surviving) == len(plot.series):
                continue
            if surviving:
                pruned_plots[plot_uuid] = plot.model_copy(update={"series": surviving})
            else:
                removed_plots.append(plot_uuid)

        # Mutate in place: PlotModel holds this same dict by reference (get_raw_scans_handle).
        self.tavi_data.plots.update(pruned_plots)
        for uuid in removed_scans:
            del self.tavi_data.raw_scans[uuid]
        for uuid in removed_plots:
            del self.tavi_data.plots[uuid]

        for uuid in removed_scans:
            self._event_broker.publish(RawScanRemoveEvent(uuid=uuid))
        for uuid in removed_plots:
            self._event_broker.publish(PlotRemoveEvent(uuid=uuid))

        return ModelResponse(code=ResponseCode.OK)

    def sync_on_ready(self, _: DownstreamReadyEvent) -> None:
        """Sync with downstream when its ready."""
        self.emit_sync_recent_projects()

    def emit_sync_recent_projects(self) -> None:
        """Notify consumers of latest recent projects."""
        recent_projects = self._get_recent_projects()
        e = SyncRecentProjects(recent_projects=recent_projects)
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
        run_names = "_".join(series.scan_name for series in e.plot.series)
        friendly_name = f"{run_names}_Plot"
        self._event_broker.publish(PlotAppendEvent(uuid=e.plot.uuid, friendly_name=friendly_name, friendly_path=""))

    def _handle_active_plot_focus_event(self, e: FocusActivePlotEvent) -> None:
        """
        Resolve one series, by its source scan's uuid, across every currently-saved plot.

        This is how a single series is picked out of an otherwise-fused, multi-series saved
        plot for "Current Plot" browsing/editing. ``uuid`` may belong to a series living in an
        unsaved preview plot instead (``PlotModel``'s to handle) — a miss here just means this
        uuid isn't currently one of ours, not a bug.
        """
        match = find_series_by_source(list(self.tavi_data.plots.values()), e.uuid)
        if match is None:
            return
        _, series = match
        scan = self.tavi_data.raw_scans[series.source_scan_uuid]
        self._event_broker.publish(ActivePlotChangedEvent(scan=scan, series=series))

    def _handle_focus_event(self, e: FocusEvent) -> None:
        """Route a ``FocusEvent`` to type-specific downstream events."""
        ids = e.ids
        raw_scans: list[RawScan] = []
        plots: list[Plot] = []
        for uuid in ids:
            # FocusEvent ids come from the project tree, which only ever lists uuids
            # TaviData actually owns — fetch_by_uuid raising here means tree/TaviData are
            # out of sync, a bug worth surfacing loudly rather than silently dropping the item.
            inst = self.tavi_data.fetch_by_uuid(uuid)
            if isinstance(inst, RawScan):
                raw_scans.append(inst)
            if isinstance(inst, Plot):
                plots.append(inst)

        if raw_scans:
            self._event_broker.publish(RawScanFocusEvent(scans=raw_scans))
        if plots:
            scans = scans_for_plots(plots, self.tavi_data.raw_scans)
            self._event_broker.publish(PlotFocusEvent(plots=plots, scans=scans))
