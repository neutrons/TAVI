"""Scan object."""

from typing import Any, Dict, Optional, Tuple

from pydantic import BaseModel, Field
from pydantic.dataclasses import dataclass


@dataclass
class UUID:
    """uuid class."""

    value: str

    def __hash__(self) -> int:
        """Generate hash."""
        return hash(self.value)


class ScanData(BaseModel):
    """
    Multicolumn data from triple-axis scans.

    Should contain {column_name: list[float,...]}.
    """

    data: Dict[str, list[float]] = Field(default_factory=dict)

    def __getattr__(self, key: str) -> list[float]:
        """Allow access as ScanData.h etc."""
        if key in self.data:
            return self.data[key]
        raise AttributeError(
            f"{self.__class__.__name__!s} has no attribute {key!r}. Valid columns are: {list(self.data.keys())}"
        )

    def __dir__(self) -> list[str]:
        """Allow tab suggestion."""
        return sorted(set(super().__dir__()) | set(self.data))


class ScanMetadata(BaseModel):
    """
    Meta data associated with a triple-axis scan.

    Should contain {name: metadata} and loader: loaderENUM. Each should be created as attributes during loading event.
    """

    data: Dict[str, Any] = Field(default_factory=dict)

    def __getattr__(self, key: str) -> Any:
        """Allow access as ScanData.h etc."""
        if key in self.data:
            return self.data[key]
        raise AttributeError(
            f"{self.__class__.__name__!s} has no attribute {key!r}. Valid columns are: {list(self.data.keys())}"
        )

    def __dir__(self) -> list[str]:
        """Allow tab suggestion."""
        return sorted(set(super().__dir__()) | set(self.data))


@dataclass
class TaviMetadata:
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
    contributing_scans: Dict[UUID, int]


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

    uuid: UUID
    data: ScanData
    metadata: ScanMetadata
    tavimeta: TaviMetadata
    prov: Provenance


class RawScan(BaseModel, Scan):
    """
    Raw scan class.

    Inherit from Scan but set uuid, data, metadata, prov as read-only. tavimeta is still writable.
    """

    uuid: UUID = Field(frozen=True)
    data: ScanData = Field(frozen=True)
    metadata: ScanMetadata = Field(frozen=True)
    tavimeta: TaviMetadata = Field(frozen=False)
    prov: Provenance = Field(frozen=True)


class ComboScan(BaseModel, Scan):
    """Combined scan class. Same as Scan, all fields read/write-able."""

    pass
