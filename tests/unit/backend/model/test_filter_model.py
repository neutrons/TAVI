"""Tests for FilterModel."""

from tavi.backend.model.filter_model import FilterModel
from tavi.library.data.enum.filter_operator import FilterOperator
from tavi.library.data.enum.filter_scope import FilterScope
from tavi.library.data.filter import FilterCriterion, ScanFilter, sample_filter, sample_filter_uuid
from tavi.library.data.model_response import ResponseCode
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
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


def make_filter(uuid_val="filter-001", value="10", name="Filter") -> ScanFilter:
    criterion = FilterCriterion(scope=FilterScope.META, key="scan", operator=FilterOperator.GREATER, value=value)
    return ScanFilter(uuid=UUID(value=uuid_val), name=name, criteria=[criterion])


def make_raw_scan(uuid_val="scan-001", metadata=None):
    if metadata is None:
        metadata = {"scan": 1, "time": 0}
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data={"qh": [1.0, 2.0], "en": [3.0, 4.0], "time": [0.0, 1.0]}),
        metadata=ScanMetadata(data=metadata),
        tavimeta=TaviMetadata(
            default_axis=("qh", "en"),
            normalization=("monitor", 1.0),
            friendly_name="test_scan",
            friendly_path="/exp1",
        ),
        prov=Provenance(raw_file="scan.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def record_filter_events():
    received = []
    broker = EventBroker()
    event_types = (AddFilterEvent, ApplyFilterEvent, ClearFilterEvent, RemoveFilterEvent, SyncFilterKeySuggestionsEvent, SyncFilterValidationEvent)
    for event_type in event_types:
        broker.register(event_type, received.append)
    return received


def make_model_with(*filters, scans=None) -> FilterModel:
    """Build a FilterModel over a filter store already holding ``filters``, as TaviData would."""
    if scans is None:
        scans = dict()
    return FilterModel({scan_filter.uuid: scan_filter for scan_filter in filters}, raw_scans=scans)


# ---------------------------------------------------------------------------
# __init__
# ---------------------------------------------------------------------------


def test_starts_with_sample_selected_and_nothing_applied():
    model = make_model_with()

    assert model._selected_uuid == sample_filter_uuid()
    assert model._applied_uuid is None


def test_sample_is_never_put_in_the_filter_store():
    filters = {}
    FilterModel(filters, dict()).apply_filter(sample_filter())

    assert filters == {}


# ---------------------------------------------------------------------------
# save_filter / select_filter
# ---------------------------------------------------------------------------


def test_save_filter_writes_new_filter_into_store_and_announces_it():
    filters = {}
    received = record_filter_events()

    response = FilterModel(filters, dict()).save_filter(make_filter())

    assert response.code == ResponseCode.OK
    assert filters == {UUID(value="filter-001"): make_filter()}
    assert received == [AddFilterEvent(filter=make_filter())]


def test_resaving_filter_replaces_it_without_adding_it_again():
    model = make_model_with(make_filter(value="1"))
    received = record_filter_events()

    model.save_filter(make_filter(value="2"))

    assert model._filters[UUID(value="filter-001")].criteria[0].value == "2"
    assert received == []


def test_resaving_filter_announces_how_it_now_stands():
    model = make_model_with(make_filter(value="1"))
    received = []
    EventBroker().register(SyncFilterEvent, received.append)

    model.save_filter(make_filter(value="2"))

    assert received == [SyncFilterEvent(filter=make_filter(value="2"))]


def test_saving_new_filter_does_not_sync_it():
    received = []
    EventBroker().register(SyncFilterEvent, received.append)

    make_model_with().save_filter(make_filter())

    assert received == []


def test_resaving_applied_filter_does_not_reapply_it():
    model = make_model_with(make_filter(value="1"))
    model.apply_filter(make_filter(value="1"))
    received = record_filter_events()

    model.save_filter(make_filter(value="2"))

    assert received == []


def test_save_filter_refuses_to_change_sample():
    model = make_model_with()
    edited = sample_filter().model_copy(update={"name": "Mine"})

    assert model.save_filter(edited).code == ResponseCode.ERROR
    assert model._filters == {}


def test_select_filter_sees_filters_saved_after_construction():
    filters = {}
    model = FilterModel(filters, dict())
    filters[UUID(value="filter-001")] = make_filter()

    assert model.select_filter(UUID(value="filter-001")).code == ResponseCode.OK


def test_select_filter_tracks_selection():
    model = make_model_with(make_filter())

    assert model.select_filter(UUID(value="filter-001")).code == ResponseCode.OK
    assert model._selected_uuid == UUID(value="filter-001")


def test_select_unknown_filter_keeps_selection():
    model = make_model_with()

    assert model.select_filter(UUID(value="nope")).code == ResponseCode.ERROR
    assert model._selected_uuid == sample_filter_uuid()


# ---------------------------------------------------------------------------
# remove_filter
# ---------------------------------------------------------------------------


def test_remove_filter_drops_it_from_store_and_announces_it():
    model = make_model_with(make_filter())
    received = record_filter_events()

    assert model.remove_filter(UUID(value="filter-001")).code == ResponseCode.OK

    assert model._filters == {}
    assert received == [RemoveFilterEvent(uuid=UUID(value="filter-001"))]


def test_remove_filter_refuses_sample():
    assert make_model_with().remove_filter(sample_filter_uuid()).code == ResponseCode.ERROR


def test_remove_unknown_filter_is_an_error():
    received = record_filter_events()

    assert make_model_with().remove_filter(UUID(value="nope")).code == ResponseCode.ERROR
    assert received == []


def test_removing_selected_filter_falls_back_to_sample():
    model = make_model_with(make_filter())
    model.select_filter(UUID(value="filter-001"))

    model.remove_filter(UUID(value="filter-001"))

    assert model._selected_uuid == sample_filter_uuid()


def test_removing_applied_filter_clears_it_before_announcing_removal():
    model = make_model_with(make_filter())
    model.apply_filter(make_filter())
    received = record_filter_events()

    model.remove_filter(UUID(value="filter-001"))

    assert received == [ClearFilterEvent(), RemoveFilterEvent(uuid=UUID(value="filter-001"))]
    assert model._applied_uuid is None


def test_removing_another_filter_leaves_applied_and_selected_alone():
    model = make_model_with(make_filter(), make_filter("filter-002"))
    model.apply_filter(make_filter())

    model.remove_filter(UUID(value="filter-002"))

    assert model._applied_uuid == model._selected_uuid == UUID(value="filter-001")


# ---------------------------------------------------------------------------
# validate_filter
# ---------------------------------------------------------------------------


def test_validate_filter_publishes_no_errors_for_valid_filter():
    received = record_filter_events()

    response = make_model_with().validate_filter(make_filter())

    assert response.code == ResponseCode.OK
    assert received == [SyncFilterValidationEvent(filter_uuid=UUID(value="filter-001"), errors=[])]


def test_validate_filter_publishes_errors_for_invalid_filter():
    received = record_filter_events()

    response = make_model_with().validate_filter(make_filter(value="ten"))

    assert response.code == ResponseCode.ERROR
    assert len(received) == 1
    assert received[0].filter_uuid == UUID(value="filter-001")
    assert "is not a number" in received[0].errors[0]


def test_validate_filter_does_not_save_or_apply():
    model = make_model_with()
    received = record_filter_events()

    model.validate_filter(make_filter())

    assert model._filters == {}
    assert not any(isinstance(e, ApplyFilterEvent) for e in received)


# ---------------------------------------------------------------------------
# apply_filter / clear_filter
# ---------------------------------------------------------------------------


def test_apply_filter_applies_and_selects_it():
    model = make_model_with(make_filter())
    received = record_filter_events()

    assert model.apply_filter(make_filter()).code == ResponseCode.OK

    assert received == [ApplyFilterEvent(filter=make_filter())]
    assert model._applied_uuid == model._selected_uuid == UUID(value="filter-001")


def test_apply_filter_saves_the_panels_edits_first():
    model = make_model_with(make_filter(value="1"))
    received = record_filter_events()

    model.apply_filter(make_filter(value="2"))

    assert model._filters[UUID(value="filter-001")].criteria[0].value == "2"
    assert received == [ApplyFilterEvent(filter=make_filter(value="2"))]


def test_apply_filter_clears_previously_applied_filter_first():
    model = make_model_with(make_filter())
    model.apply_filter(sample_filter())
    received = record_filter_events()

    model.apply_filter(make_filter())

    assert received == [ClearFilterEvent(), ApplyFilterEvent(filter=make_filter())]


def test_apply_sample_uses_the_built_in_sample_whatever_is_passed():
    model = make_model_with()
    received = record_filter_events()
    edited = sample_filter().model_copy(update={"name": "Mine", "criteria": []})

    model.apply_filter(edited)

    assert received == [ApplyFilterEvent(filter=sample_filter())]
    assert model._filters == {}


def test_apply_invalid_filter_publishes_validation_and_applies_nothing():
    model = make_model_with()
    received = record_filter_events()

    assert model.apply_filter(make_filter(value="")).code == ResponseCode.ERROR

    assert [type(e) for e in received] == [AddFilterEvent, SyncFilterValidationEvent]
    assert model._applied_uuid is None


def test_apply_invalid_filter_still_saves_and_selects_it():
    model = make_model_with()

    model.apply_filter(make_filter(value=""))

    assert UUID(value="filter-001") in model._filters
    assert model._selected_uuid == UUID(value="filter-001")


def test_apply_invalid_filter_leaves_applied_filter_in_place():
    model = make_model_with()
    model.apply_filter(sample_filter())

    model.apply_filter(make_filter(value=""))

    assert model._applied_uuid == sample_filter_uuid()


def test_clear_filter_publishes_clear_and_forgets_applied():
    model = make_model_with()
    model.apply_filter(sample_filter())
    received = record_filter_events()

    assert model.clear_filter().code == ResponseCode.OK

    assert received == [ClearFilterEvent()]
    assert model._applied_uuid is None


def test_clear_filter_with_nothing_applied_publishes_nothing():
    received = record_filter_events()

    make_model_with().clear_filter()

    assert received == []


def test_clear_filter_keeps_selection():
    model = make_model_with(make_filter())
    model.apply_filter(make_filter())

    model.clear_filter()

    assert model._selected_uuid == UUID(value="filter-001")


# ---------------------------------------------------------------------------
# scan_keys
# ---------------------------------------------------------------------------

def test_increment_scan_keys_stores_and_publishes():
    filters = dict()
    scan_a = make_raw_scan("scan-a", metadata={"scan_a_key": "some_value"})
    scan_b = make_raw_scan("scan-b", metadata={"scan_b_key": "some_value"})
    scans = {scan_a.uuid: scan_a, scan_b.uuid: scan_b}
    received = record_filter_events()
    model = FilterModel(filters, scans)

    response = model.increment_scan_keys(scan_a.uuid)
    assert response.code == ResponseCode.OK
    assert len(received) == 1
    assert isinstance(received[0], SyncFilterKeySuggestionsEvent)
    response = model.increment_scan_keys(scan_b.uuid)
    assert len(received) == 2
    assert isinstance(received[1], SyncFilterKeySuggestionsEvent)
    event = received[1]
    assert event.metadata_keys
    assert event.data_keys

    for key in scan_a.metadata.data.keys():
        assert key in event.metadata_keys
    for key in scan_b.metadata.data.keys():
        assert key in event.metadata_keys
    for key in scan_a.data.data.keys():
        assert key in event.data_keys
    for key in scan_b.data.data.keys():
        assert key in event.data_keys


def test_increment_refuses_bad_lookup():
    response = FilterModel(dict(), dict()).increment_scan_keys(UUID("scan_a"))
    assert response.code == ResponseCode.ERROR


def test_decrement_scan_keys_removes_and_publishes():
    filters = dict()
    scan_a = make_raw_scan("scan-a", metadata={"scan_a_key": "some_value", "scan": 0})
    scan_b = make_raw_scan("scan-b", metadata={"scan_b_key": "some_value", "scan": 1})
    scans = {scan_a.uuid: scan_a, scan_b.uuid: scan_b}
    received = record_filter_events()
    model = FilterModel(filters, scans)

    model.increment_scan_keys(scan_a.uuid)
    model.increment_scan_keys(scan_b.uuid)
    scans.pop(scan_b.uuid)
    response = model.decrement_scan_keys(scan_b.uuid, scan_b)
    assert response.code == ResponseCode.OK
    assert len(received) == 3
    assert isinstance(received[2], SyncFilterKeySuggestionsEvent)
    event = received[2]
    assert event.metadata_keys
    assert event.data_keys

    for key in scan_a.metadata.data.keys():
        assert key in event.metadata_keys
    assert "scan_b_key" not in event.metadata_keys
    for key in scan_a.data.data.keys():
        assert key in event.data_keys


def test_decrement_refuses_removal_of_present_scan():
    filters = dict()
    scan_a = make_raw_scan("scan-a", metadata={"scan_a_key": "some_value", "scan": 0})
    scans = {scan_a.uuid: scan_a}
    response = FilterModel(filters, scans).decrement_scan_keys(scan_a.uuid, scan_a)
    assert response.code == ResponseCode.ERROR
