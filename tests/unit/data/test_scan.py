# test_scan.py
import dataclasses
import pytest

# Adjust this import to match your project layout, e.g.:
# from tavi.backend.model.scan import RawData, RawMetaData, Scan
from tavi.library.data.scan import RawData, RawMetaData, Scan


def test_classes_are_dataclasses() -> None:
    assert dataclasses.is_dataclass(RawData)
    assert dataclasses.is_dataclass(RawMetaData)
    assert dataclasses.is_dataclass(Scan)


def test_rawdata_and_rawmetadata_can_be_instantiated() -> None:
    RawData()
    RawMetaData()


def test_scan_requires_all_fields() -> None:
    with pytest.raises(TypeError):
        Scan()  # missing required args


def test_scan_stores_fields_correctly() -> None:
    data = RawData()
    meta = RawMetaData()
    err = ("warning: low counts",)
    others = ("note: user tag", 123)

    scan = Scan(data=data, metadata=meta, error_message=err, others=others)

    assert scan.data is data
    assert scan.metadata is meta
    assert scan.error_message == err
    assert scan.others == others


def test_scan_equality_and_repr() -> None:
    data1, meta1 = RawData(), RawMetaData()
    data2, meta2 = RawData(), RawMetaData()

    s1 = Scan(data=data1, metadata=meta1, error_message=(), others=())
    s1_same = Scan(data=data1, metadata=meta1, error_message=(), others=())
    s2 = Scan(data=data2, metadata=meta2, error_message=(), others=())

    assert s1 == s1_same

    r = repr(s1)
    assert "Scan" in r
    assert "data" in r
    assert "metadata" in r
    assert "error_message" in r
    assert "others" in r


def test_scan_field_types_are_as_expected() -> None:
    fields = {f.name: f.type for f in dataclasses.fields(Scan)}
    assert fields["data"] is RawData
    assert fields["metadata"] is RawMetaData
    assert fields["error_message"] is tuple
    assert fields["others"] is tuple