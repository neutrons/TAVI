"""Tests for FilterModel."""

from tavi.backend.model.filter_model import FilterModel
from tavi.library.data.enum.filter_operator import FilterOperator
from tavi.library.data.filter import FilterCriterion, ScanFilter, sample_filter, sample_filter_uuid
from tavi.library.data.model_response import ResponseCode
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import (
    AddFilterEvent,
    ApplyFilterEvent,
    ClearFilterEvent,
    RemoveFilterEvent,
    SyncFilterValidationEvent,
)


def make_filter(uuid_val="filter-001", value="10", name="Filter") -> ScanFilter:
    criterion = FilterCriterion(key="scan", operator=FilterOperator.GREATER, value=value)
    return ScanFilter(uuid=UUID(value=uuid_val), name=name, criteria=[criterion])


def record_filter_events():
    received = []
    broker = EventBroker()
    event_types = (AddFilterEvent, ApplyFilterEvent, ClearFilterEvent, RemoveFilterEvent, SyncFilterValidationEvent)
    for event_type in event_types:
        broker.register(event_type, received.append)
    return received


def make_model_with(*filters) -> FilterModel:
    """Build a FilterModel over a filter store already holding ``filters``, as TaviData would."""
    return FilterModel({scan_filter.uuid: scan_filter for scan_filter in filters})


# ---------------------------------------------------------------------------
# __init__
# ---------------------------------------------------------------------------


def test_starts_with_sample_selected_and_nothing_applied():
    model = make_model_with()

    assert model._selected_uuid == sample_filter_uuid()
    assert model._applied_uuid is None


def test_sample_is_never_put_in_the_filter_store():
    filters = {}
    FilterModel(filters).apply_filter(sample_filter())

    assert filters == {}


# ---------------------------------------------------------------------------
# save_filter / select_filter
# ---------------------------------------------------------------------------


def test_save_filter_writes_new_filter_into_store_and_announces_it():
    filters = {}
    received = record_filter_events()

    response = FilterModel(filters).save_filter(make_filter())

    assert response.code == ResponseCode.OK
    assert filters == {UUID(value="filter-001"): make_filter()}
    assert received == [AddFilterEvent(filter=make_filter())]


def test_resaving_filter_replaces_it_without_announcing_again():
    model = make_model_with(make_filter(value="1"))
    received = record_filter_events()

    model.save_filter(make_filter(value="2"))

    assert model._filters[UUID(value="filter-001")].criteria[0].value == "2"
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
    model = FilterModel(filters)
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
