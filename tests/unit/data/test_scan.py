import pytest
from pydantic import ValidationError

from tavi.library.data.scan import Scan, Data, MetaData


def test_scan_constructs_with_valid_types():
    scan = Scan(
        data=Data(),
        metadata=MetaData(),
        error_messages=(),
        provenance=(),
    )

    assert isinstance(scan.data, Data)
    assert isinstance(scan.metadata, MetaData)
    assert scan.error_messages == ()
    assert scan.provenance == ()


def test_scan_rejects_wrong_data_type():
    with pytest.raises(ValidationError):
        Scan(
            data="not-rawdata",  # wrong
            metadata=MetaData(),
            error_messages=(),
            provenance=(),
        )


def test_scan_rejects_wrong_metadata_type():
    with pytest.raises(ValidationError):
        Scan(
            data=Data(),
            metadata=123,  # wrong
            error_messages=(),
            provenance=(),
        )


def test_scan_rejects_non_tuple_error_messages():
    with pytest.raises(ValidationError):
        Scan(
            data=Data(),
            metadata=MetaData(),
            error_messages="this is a str",  # wrong: expects tuple
            provenance=(),
        )


def test_scan_rejects_non_tuple_others():
    with pytest.raises(ValidationError):
        Scan(
            data=Data(),
            metadata=MetaData(),
            error_messages=(),
            provenance={"not": "a tuple"},  # wrong: expects tuple
        )


def test_scan_coerces_list_to_tuple_if_allowed_by_pydantic():
    scan = Scan(
        data=Data(),
        metadata=MetaData(),
        error_messages=("a", "b"),
        provenance=("raw",),
    )
    assert scan.error_messages == ("a", "b")
    assert scan.provenance == ("raw",)