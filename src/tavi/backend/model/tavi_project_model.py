"""Tavi Project."""

from neutrons_standard.config import Resource
from neutrons_standard.decorators.singleton import Singleton
from ruamel.yaml import YAML

from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface
from tavi.backend.model.plot_resolver import scans_for_plots
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.plot import Plot
from tavi.backend.model.plot_resolver import scans_for_plots
from tavi.library.data.scan import RawScan
from tavi.library.data.tavi_data import TaviData
from tavi.library.storage.controller.raw_scan_load_controller import RawScanLoadController
from tavi.library.storage.interface.filestore_interface import Filestore
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import RawScanAppendEvent, SyncRecentProjects
from tavi.meta.event.type.presenter_event import DownstreamReadyEvent, FocusEvent, PlotFocusEvent, RawScanFocusEvent


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

    def _handle_focus_event(self, e: FocusEvent) -> None:
        """Route a ``FocusEvent`` to type-specific downstream events."""
        ids = e.ids
        raw_scans: list[RawScan] = []
        plots: list[Plot] = []
        for uuid in ids:
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
