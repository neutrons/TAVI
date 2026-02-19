"""Tavi Project."""

from tavi.event_broker.event_broker import EventBroker
from tavi.event_broker.event_type import Event

from tavi.backend.model.interface.tavi_project_interface import TaviProjectInterface
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.meta.decorators.singleton import Singleton


@Singleton
class TaviProjectModel(TaviProjectInterface):
    """Tavi project class."""

    def __init__(self) -> None:
        """Init tavi data."""
        self._event_broker = EventBroker()

    def send(self, event: Event) -> None:
        """Send pre-register event to event broker."""
        self._event_broker.publish(event)

    def load_raw_scan_from_folder(self, folder: str) -> None:
        """Load a folder containing raw scans."""
        print("folder director received by model:", folder)
        # TO DO
        # Implement load raw scan from folder logic
        # raw_scan_loading_event = RawScanLoadingEvent(raw_scan_uuid = ...)
        # self.send(raw_scan_loading_event)
        return ModelResponse(code=ResponseCode.OK, message="TODO: implement loading backend")
