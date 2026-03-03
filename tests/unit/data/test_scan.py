import pytest
from pydantic import ValidationError

from tavi.library.data.scan import Scan, RawData, RawMetaData


def test_scan_constructs_with_valid_types():
    scan = Scan(
        data=RawData(),
        metadata=RawMetaData(),
        error_messages=(),
        others=(),
    )

    assert isinstance(scan.data, RawData)
    assert isinstance(scan.metadata, RawMetaData)
    assert scan.error_messages == ()
    assert scan.others == ()


def test_scan_rejects_wrong_data_type():
    with pytest.raises(ValidationError):
        Scan(
            data="not-rawdata",  # wrong
            metadata=RawMetaData(),
            error_messages=(),
            others=(),
        )


def test_scan_rejects_wrong_metadata_type():
    with pytest.raises(ValidationError):
        Scan(
            data=RawData(),
            metadata=123,  # wrong
            error_messages=(),
            others=(),
        )


def test_scan_rejects_non_tuple_error_messages():
    with pytest.raises(ValidationError):
        Scan(
            data=RawData(),
            metadata=RawMetaData(),
            error_messages="this is a str",  # wrong: expects tuple
            others=(),
        )


def test_scan_rejects_non_tuple_others():
    with pytest.raises(ValidationError):
        Scan(
            data=RawData(),
            metadata=RawMetaData(),
            error_messages=(),
            others={"not": "a tuple"},  # wrong: expects tuple
        )


def test_scan_coerces_list_to_tuple_if_allowed_by_pydantic():
    scan = Scan(
        data=RawData(),
        metadata=RawMetaData(),
        error_messages=["a", "b"],
        others=[1, 2],
    )
    assert scan.error_messages == ("a", "b")
    assert scan.others == (1, 2)