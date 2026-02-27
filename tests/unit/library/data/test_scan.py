import dataclasses

import pytest

from tavi.library.data.scan import RawData, RawMetaData, Scan


def test_rawdata_is_dataclass_and_instantiable():
    assert dataclasses.is_dataclass(RawData)
    rd = RawData()
    assert isinstance(rd, RawData)


def test_rawmetadata_is_dataclass_and_instantiable():
    assert dataclasses.is_dataclass(RawMetaData)
    md = RawMetaData()
    assert isinstance(md, RawMetaData)


def test_scan_is_dataclass_and_fields_roundtrip():
    assert dataclasses.is_dataclass(Scan)

    data = RawData()
    metadata = RawMetaData()
    error_message = ("warning: motor stalled", "note: low counts")
    others = ("aux1", {"key": "value"}, 123)

    scan = Scan(
        data=data,
        metadata=metadata,
        error_message=error_message,
        others=others,
    )

    assert scan.data is data
    assert scan.metadata is metadata
    assert scan.error_message == error_message
    assert scan.others == others


def test_scan_equality_semantics_from_dataclass():
    data1 = RawData()
    meta1 = RawMetaData()

    s1 = Scan(data=data1, metadata=meta1, error_message=("a",), others=("b",))
    s2 = Scan(data=data1, metadata=meta1, error_message=("a",), others=("b",))
    s3 = Scan(data=RawData(), metadata=meta1, error_message=("a",), others=("b",))

    assert s1 == s2
    assert s1 != s3


def test_scan_type_annotations_exist():
    ann = Scan.__annotations__
    assert ann["data"] is RawData
    assert ann["metadata"] is RawMetaData
    assert ann["error_message"] is tuple
    assert ann["others"] is tuple
