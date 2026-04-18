"""Tavi Project."""

from neutrons_standard.decorators.singleton import Singleton

from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.scan import RawScan
from tavi.library.data.tavi_data import TaviData
from tavi.library.storage.controller.raw_scan_load_controller import RawScanLoadController
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import RawScanAppendEvent


@Singleton
class TaviProjectModel(TaviProjectInterface):
    """Tavi project class."""

    def __init__(self) -> None:
        """Init tavi data."""
        self.tavi_data: TaviData = TaviData(raw_scans={})
        self._event_broker: EventBroker = EventBroker()
        self.raw_scan_load_controller: RawScanLoadController = RawScanLoadController()

    def load_raw_scan_from_folder(self, folder: str) -> None:
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
