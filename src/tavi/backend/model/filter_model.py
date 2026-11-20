"""Filter model module."""

from typing import Optional

from tavi.backend.model.interface.filter_model_interface import FilterModelInterface
from tavi.library.data.filter import ScanFilter, sample_filter
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.scan import UUID, RawScan
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import (
    AddFilterEvent,
    ApplyFilterEvent,
    ClearFilterEvent,
    RemoveFilterEvent,
    SyncFilterEvent,
    SyncFilterKeySuggestionsEvent,
    SyncFilterValidationEvent,
)


class FilterModel(FilterModelInterface):
    """Owns the saved filters and which is selected and applied, and announces what is applied."""

    def __init__(self, filters: dict[UUID, ScanFilter], raw_scans: dict[UUID, RawScan]) -> None:
        """Initialize with a handle into TaviData's filter and raw scans storages, the Sample filter selected, and nothing applied."""
        super().__init__()
        # TaviData's own dict, edited in place: filters are saved with the project, but no other
        # model reads or writes them, so there is nothing to route through TaviProjectModel.
        self._filters = filters
        # TaviData's dict, treated as read-only from this model
        self._raw_scans = raw_scans
        self._sample = sample_filter()
        self._selected_uuid: UUID = self._sample.uuid
        self._applied_uuid: Optional[UUID] = None
        self._metadata_keys: dict[str, int] = dict()
        self._data_keys: dict[str, int] = dict()

        self._event_broker = EventBroker()

    def save_filter(self, scan_filter: ScanFilter) -> ModelResponse:
        """Save ``scan_filter`` into the project, replacing any saved filter with the same uuid."""
        if scan_filter.uuid == self._sample.uuid:
            return ModelResponse(code=ResponseCode.ERROR, message="The Sample filter cannot be changed.")
        is_new = scan_filter.uuid not in self._filters
        # Saving an edit to the applied filter does not reapply it: what is applied is what the
        # user last chose to apply, until they apply again.
        self._filters[scan_filter.uuid] = scan_filter
        if is_new:
            self._event_broker.publish(AddFilterEvent(filter=scan_filter))
        else:
            self._event_broker.publish(SyncFilterEvent(filter=scan_filter))
        return ModelResponse(code=ResponseCode.OK)

    def remove_filter(self, filter_uuid: UUID) -> ModelResponse:
        """Drop a saved filter from the project, clearing it first if it is the one applied."""
        if filter_uuid == self._sample.uuid:
            return ModelResponse(code=ResponseCode.ERROR, message="The Sample filter cannot be removed.")
        if filter_uuid not in self._filters:
            return ModelResponse(code=ResponseCode.ERROR, message=f"No saved filter {filter_uuid.value}.")
        if filter_uuid == self._applied_uuid:
            self.clear_filter()
        if filter_uuid == self._selected_uuid:
            self._selected_uuid = self._sample.uuid
        del self._filters[filter_uuid]
        self._event_broker.publish(RemoveFilterEvent(uuid=filter_uuid))
        return ModelResponse(code=ResponseCode.OK)

    def select_filter(self, filter_uuid: UUID) -> ModelResponse:
        """Make a saved filter the one ``apply_filter`` applies."""
        if self._lookup(filter_uuid) is None:
            return ModelResponse(code=ResponseCode.ERROR, message=f"No saved filter {filter_uuid.value}.")
        self._selected_uuid = filter_uuid
        return ModelResponse(code=ResponseCode.OK)

    def validate_filter(self, scan_filter: ScanFilter) -> ModelResponse:
        """Check ``scan_filter`` can be applied and publish the outcome as a SyncFilterValidationEvent."""
        errors = scan_filter.validate_filter()
        self._event_broker.publish(SyncFilterValidationEvent(filter_uuid=scan_filter.uuid, errors=errors))
        return self._validation_response(errors)

    def apply_filter(self, scan_filter: ScanFilter) -> ModelResponse:
        """
        Save, select and apply ``scan_filter`` in place of any applied one, if it is valid.

        One call rather than save, select and apply in turn: each Proxy call runs on its own worker
        thread, so separate calls could apply before the panel's latest edits were saved. The Sample
        filter is applied as built in, whatever is passed for it.
        """
        if scan_filter.uuid == self._sample.uuid:
            scan_filter = self._sample
        else:
            self.save_filter(scan_filter)
        self._selected_uuid = scan_filter.uuid

        errors = scan_filter.validate_filter()
        if errors:
            self._event_broker.publish(SyncFilterValidationEvent(filter_uuid=scan_filter.uuid, errors=errors))
            return self._validation_response(errors)

        self.clear_filter()
        self._applied_uuid = scan_filter.uuid
        self._event_broker.publish(ApplyFilterEvent(filter=scan_filter))
        return ModelResponse(code=ResponseCode.OK)

    def clear_filter(self) -> ModelResponse:
        """Stop filtering the project's scans, publishing ClearFilterEvent."""
        if self._applied_uuid is None:
            return ModelResponse(code=ResponseCode.OK)
        self._applied_uuid = None
        self._event_broker.publish(ClearFilterEvent())
        return ModelResponse(code=ResponseCode.OK)

    def _scan_keys_event(self) -> SyncFilterKeySuggestionsEvent:
        return SyncFilterKeySuggestionsEvent(
            metadata_keys=sorted(list(self._metadata_keys.keys())), data_keys=sorted(list(self._data_keys.keys()))
        )

    def increment_scan_keys(self, uuid: UUID) -> ModelResponse:
        """Increment key references on scan addition for use in key suggestion."""
        if uuid not in self._raw_scans:
            return ModelResponse(code=ResponseCode.ERROR, message="Scan not found on attempted key update.")
        scan = self._raw_scans[uuid]

        for metadata_key in scan.metadata.data.keys():
            reference_count = self._metadata_keys.setdefault(metadata_key, 0) + 1
            self._metadata_keys[metadata_key] = reference_count

        for data_key in scan.data.data.keys():
            reference_count = self._data_keys.setdefault(data_key, 0) + 1
            self._data_keys[data_key] = reference_count

        self._event_broker.publish(self._scan_keys_event())

        return ModelResponse(code=ResponseCode.OK)

    def decrement_scan_keys(self, uuid: UUID, scan: RawScan) -> ModelResponse:
        """Decrement key references on scan removal for use in key suggestion."""
        if uuid in self._raw_scans:
            return ModelResponse(code=ResponseCode.ERROR, message="Scan found on attempted key removal.")

        for metadata_key in scan.metadata.data.keys():
            reference_count = self._metadata_keys.setdefault(metadata_key, 1) - 1
            if reference_count > 0:
                self._metadata_keys[metadata_key] = reference_count
            if reference_count < 1:
                self._metadata_keys.pop(metadata_key)

        for data_key in scan.data.data.keys():
            reference_count = self._data_keys.setdefault(data_key, 1) - 1
            if reference_count > 0:
                self._data_keys[data_key] = reference_count
            if reference_count < 1:
                self._data_keys.pop(data_key)

        self._event_broker.publish(self._scan_keys_event())

        return ModelResponse(code=ResponseCode.OK)

    def _lookup(self, filter_uuid: UUID) -> Optional[ScanFilter]:
        if filter_uuid == self._sample.uuid:
            return self._sample
        return self._filters.get(filter_uuid)

    def _validation_response(self, errors: list[str]) -> ModelResponse:
        if errors:
            return ModelResponse(code=ResponseCode.ERROR, message=" ".join(errors))
        return ModelResponse(code=ResponseCode.OK)
