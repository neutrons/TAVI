"""Tests for tavi.library.data.filter."""

import pytest

from tavi.library.data.enum.filter_operator import FilterOperator
from tavi.library.data.filter import FilterCriterion, ScanFilter, sample_filter, sample_filter_uuid
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata


def make_raw_scan(metadata=None, uuid_val="scan-001") -> RawScan:
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data={"qh": [1.0], "en": [2.0]}),
        metadata=ScanMetadata(data=metadata or {}),
        tavimeta=TaviMetadata(default_axis=("qh", "en"), friendly_name="test_scan", friendly_path="/exp1"),
        prov=Provenance(raw_file="scan.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_criterion(key="scan", operator=FilterOperator.GREATER_EQUAL, value="0") -> FilterCriterion:
    return FilterCriterion(key=key, operator=operator, value=value)


def make_filter(*criteria, name="Filter") -> ScanFilter:
    return ScanFilter(name=name, criteria=list(criteria))


# ---------------------------------------------------------------------------
# FilterCriterion.matches
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("operator", "value", "expected"),
    [
        (FilterOperator.EQUAL, "42", True),
        (FilterOperator.EQUAL, "42.0", True),
        (FilterOperator.NOT_EQUAL, "42", False),
        (FilterOperator.GREATER, "41", True),
        (FilterOperator.GREATER, "42", False),
        (FilterOperator.LESS, "43", True),
        (FilterOperator.GREATER_EQUAL, "42", True),
        (FilterOperator.LESS_EQUAL, "41", False),
    ],
)
def test_numeric_operators_compare_as_numbers(operator, value, expected):
    assert make_criterion("scan", operator, value).matches({"scan": "42"}) is expected


@pytest.mark.parametrize(
    ("operator", "value", "expected"),
    [
        (FilterOperator.IS, "mnf3", True),
        (FilterOperator.IS, "MnF", False),
        (FilterOperator.IS_NOT, "MnF3", False),
        (FilterOperator.IS_NOT, "NiO", True),
        (FilterOperator.HAS, "nf", True),
        (FilterOperator.HAS_NOT, "nf", False),
        (FilterOperator.HAS_NOT, "Ni", True),
    ],
)
def test_text_operators_ignore_case(operator, value, expected):
    assert make_criterion("samplename", operator, value).matches({"samplename": "MnF3"}) is expected


@pytest.mark.parametrize("operator", list(FilterOperator))
def test_missing_key_never_matches(operator):
    assert make_criterion("temp", operator, "1").matches({"scan": "1"}) is False


def test_numeric_operator_on_non_numeric_metadata_does_not_match():
    assert make_criterion("scan", FilterOperator.GREATER, "0").matches({"scan": "n/a"}) is False


# ---------------------------------------------------------------------------
# FilterCriterion.validate_criterion
# ---------------------------------------------------------------------------


def test_valid_criterion_has_no_errors():
    assert make_criterion().validate_criterion() == []


def test_blank_key_is_invalid():
    assert make_criterion(key=" ").validate_criterion() == ["Criterion has no key."]


def test_blank_value_is_invalid():
    assert make_criterion(value="").validate_criterion() == ["Criterion on 'scan' has no value."]


def test_non_numeric_value_for_numeric_operator_is_invalid():
    errors = make_criterion(operator=FilterOperator.LESS, value="abc").validate_criterion()

    assert len(errors) == 1
    assert "is not a number" in errors[0]


def test_non_numeric_value_for_text_operator_is_valid():
    assert make_criterion(operator=FilterOperator.HAS, value="abc").validate_criterion() == []


# ---------------------------------------------------------------------------
# ScanFilter
# ---------------------------------------------------------------------------


def test_filter_matches_only_when_every_criterion_does():
    scan = make_raw_scan({"scan": "5", "samplename": "MnF3"})
    passing = make_criterion("scan", FilterOperator.GREATER, "1")
    failing = make_criterion("samplename", FilterOperator.IS, "NiO")

    assert make_filter(passing).matches(scan) is True
    assert make_filter(passing, failing).matches(scan) is False


def test_filter_without_criteria_matches_everything():
    assert make_filter().matches(make_raw_scan()) is True


def test_filter_validation_collects_every_criterion_error():
    scan_filter = make_filter(make_criterion(key=""), make_criterion(value=""))

    assert len(scan_filter.validate_filter()) == 2


def test_blank_filter_name_is_invalid():
    assert make_filter(name=" ").validate_filter() == ["Filter has no name."]


def test_filters_get_distinct_uuids():
    assert make_filter().uuid != make_filter().uuid


def test_sample_filter_passes_numbered_scans():
    scan_filter = sample_filter()

    assert scan_filter.uuid == sample_filter_uuid()
    assert scan_filter.name == "Sample"
    assert scan_filter.validate_filter() == []
    assert scan_filter.matches(make_raw_scan({"scan": "12"})) is True
