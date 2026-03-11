"""Scan object."""

from typing import Dict, Optional, Tuple

from pydantic import BaseModel, Field
from pydantic.dataclasses import dataclass


@dataclass
class Data:
    """
    Multicolumn data from triple-axis scans.

    Should contain {column_name: list[float,...]} implemented as attributes:values
    """

    pass


@dataclass
class MetaData:
    """
    Meta data associated with a triple-axis scan.

    Should contain {name: metadata} and loader: loaderENUM. Each should be created as attributes during loading event.
    """

    pass


@dataclass
class TaviMetaData:
    """
    Tavi specific meta data.

    Args:
        default_axis:Tuple[str, str]
        nomarlization: Tuple[str, int] as (column_name, weight)

    """

    default_axis: Tuple[str, str]
    nomarlization: Optional[str] = None


@dataclass
class Provenance:
    """
    History of the Scan class..

    Args:
        raw_file: path of the raw scan file.
        contributing_scans: dict[uuid, weight]

    """

    raw_file: str
    contributing_scans: Dict[str, int]


@dataclass
class Scan:
    """
    Tavi scan data class.

    Represents a single scan within a Tavi project, containing both raw measurement
    data and associated metadata.

    Attributes:
        uuid: a string representing a universally unique identifier.
        data (Data): Numerical measurement arrays collected during the scan
            (e.g., motor positions, detector counts, temperatures).
        metadata (MetaData): Descriptive information about the scan
            (e.g., experiment details, instrument settings, sample information).
        tavimeta: tavi specific MetaData. Always have both read/write permission.
        prov (Provenance): History of the scan object.

    """

    uuid: str
    data: Data
    metadata: MetaData
    tavimeta: TaviMetaData
    prov: Provenance


class RawScan(BaseModel, Scan):
    """
    Raw scan class.

    Inherit from Scan but set uuid, data, metadata, prov as read-only. tavimeta is still writable.
    """

    uuid: str = Field(frozen=True)
    data: Data = Field(frozen=True)
    metadata: MetaData = Field(frozen=True)
    tavimeta: TaviMetaData = Field(frozen=False)
    prov: Provenance = Field(frozen=True)


class ComboScan(BaseModel, Scan):
    """Combined scan class. Same as Scan, all fields read/write-able."""

    pass
