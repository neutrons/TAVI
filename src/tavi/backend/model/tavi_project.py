"""Tavi Project."""

from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface
from tavi.event_broker.event_broker import EventBroker
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.meta.decorators.singleton import Singleton


@Singleton
class TaviProject(TaviProjectInterface):
    def __init__(self) -> None:
        """Init tavi data."""
        self._event_broker = EventBroker()

    def send(self, event):
        self._event_broker.publish(event)

    def load_raw_scan_from_folder(self, folder: str) -> None:
        print("folder director received by model:", folder)
        # TO DO
        # Implement load raw scan from folder logic
        # raw_scan_loading_event = RawScanLoadingEvent(raw_scan_uuid = ...)
        # self.send(raw_scan_loading_event)
        return ModelResponse(code=ResponseCode.OK, message="TODO: implement loading backend")
