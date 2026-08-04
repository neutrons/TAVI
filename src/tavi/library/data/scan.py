"""Scan object."""

from typing import Any, Dict, Optional, Tuple
from uuid import uuid4

from pydantic import BaseModel, Field
from pydantic.dataclasses import dataclass
from pydantic.fields import FieldInfo


@dataclass
class UUID:
    """uuid class."""

    value: str

    def __hash__(self) -> int:
        """Generate hash."""
        return hash(self.value)


def UUIDFactory() -> FieldInfo:
    """Return a Pydantic Field that generates a fresh UUID on instantiation."""
    return Field(default_factory=lambda: UUID(value=str(uuid4())))


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

    ``categories`` maps a front end display name to the key in ``data`` holding that category's fields, e.g.
    {"ORNL Metadata": "ORNL Metadata"}. It decouples how the frontend groups/labels metadata from how the
    loader organizes ``data``, while still allowing flat ``.`` access to fields nested under a category.
    """

    data: Dict[str, Any] = Field(default_factory=dict)
    categories: Dict[str, str] = Field(default_factory=dict)

    def __getattr__(self, key: str) -> Any:
        """Allow access as ScanData.h etc, searching top-level data first, then each category's fields."""
        if key in self.data:
            return self.data[key]
        for data_key in self.categories.values():
            category = self.data.get(data_key)
            if isinstance(category, dict) and key in category:
                return category[key]
        raise AttributeError(
            f"{self.__class__.__name__!s} has no attribute {key!r}. Valid columns are: {list(self.data.keys())}"
        )

    def __dir__(self) -> list[str]:
        """Allow tab suggestion."""
        keys = set(self.data)
        for data_key in self.categories.values():
            category = self.data.get(data_key)
            if isinstance(category, dict):
                keys |= set(category)
        return sorted(set(super().__dir__()) | keys)

    def by_category(self) -> Dict[str, Any]:
        """Return data grouped by front end category display name, for display purposes."""
        return {display_name: self.data.get(data_key, {}) for display_name, data_key in self.categories.items()}


@dataclass
class TaviMetadata:
    """
    Tavi specific meta data.

    Args:
        default_axis:Tuple[str, str]
        normalization: Tuple[str, int] as (column_name, weight)

    """

    default_axis: Tuple[str, str]
    friendly_name: str
    friendly_path: str
    normalization: Optional[Tuple[str, float]] = None


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
