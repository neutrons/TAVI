"""Scan object."""

from typing import Any, Dict, List, Optional, Tuple
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


def new_uuid() -> UUID:
    """Mint a fresh, unique UUID - the one place TAVI generates them, so the scheme can be swapped here alone."""
    return UUID(value=str(uuid4()))


def UUIDFactory() -> FieldInfo:
    """Return a Pydantic Field that generates a fresh UUID on instantiation."""
    return Field(default_factory=new_uuid)


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

    ``data`` stays a flat {field_name: value} mapping. ``categories`` maps a front end display name to the list
    of ``data`` keys that belong under it, e.g. {"ORNL Metadata": ["scan", "proposal", ...]}. This decouples how
    the frontend groups/labels metadata for display from the (flat) storage layout, without disturbing the
    single-level ``.`` access to fields.
    """

    data: Dict[str, Any] = Field(default_factory=dict)
    categories: Dict[str, List[str]] = Field(default_factory=dict)

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

    def by_category(self) -> Dict[str, Any]:
        """Return data grouped by front end category display name, for display purposes."""
        return {
            display_name: {key: self.data[key] for key in keys if key in self.data}
            for display_name, keys in self.categories.items()
        }


@dataclass
class TaviMetadata:
    """
    Tavi specific meta data.

    Args:
        default_axis:Tuple[str, str]
        normalization: (column_name, weight) for a scan normalized by one factor, or
            (column_name, grid) for a ``HistogramData``, the grid holding each bin's summed norm.

    """

    default_axis: Tuple[str, str]
    friendly_name: str
    friendly_path: str
    normalization: Optional[Tuple[str, float] | Tuple[str, List[List[float]]]] = None


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


class ProcessedScan(BaseModel, Scan):
    """
    Processed scan class. Same as Scan, all fields read/write-able.

    Support combined scans, rebined scans. Same as Scan class, multi-column data.
    """

    pass


class HistogramData(BaseModel, Scan):
    """
    A scan binned onto a regular 2D intensity map, every grid shaped ``(n_x, n_y)``. All fields read/write-able.

    Carries the ``Scan`` fields like any derived scan, plus the map: bin edges and centres per
    axis and the per-bin ``intensity``, ``error`` and point ``counts``, with ``NaN`` intensity
    and error where no point fell. Each bin's summed norm is ``tavimeta.normalization``, as
    ``(norm column, grid)``. Grids are nested lists of plain floats, as
    ``ScanData`` stores its columns, so the scan deep-copies and serializes like any other.
    The x and y column names live in ``tavimeta.default_axis``, where every scan keeps its
    axes; ``intensity_name`` names the column the intensity was binned from.
    """

    intensity_name: str
    x_edges: list[float]
    y_edges: list[float]
    x_centers: list[float]
    y_centers: list[float]
    intensity: list[list[float]]
    error: list[list[float]]
    counts: list[list[int]]
