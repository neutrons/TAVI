"""Scan object."""

from dataclasses import dataclass


@dataclass
class RawData:
    """Declare dataclass for type hints."""

    pass


@dataclass
class RawMetaData:
    """Declare dataclass for type hints."""

    pass


@dataclass
class Scan:
    """
    Tavi scan data class.

    Represents a single scan within a Tavi project, containing both raw measurement
    data and associated metadata.

    Attributes:
        data (RawData): Numerical measurement arrays collected during the scan
            (e.g., motor positions, detector counts, temperatures).
        metadata (RawMetaData): Descriptive information about the scan
            (e.g., experiment details, instrument settings, sample information).
        error_message (tuple): Messages or warnings associated with the scan,
            such as instrument errors or data quality issues.
        others (tuple): Miscellaneous or auxiliary information related to the scan
            that does not fit into "data", "metadata"` or "error_message".

    """

    data: RawData
    metadata: RawMetaData
    error_message: tuple
    others: tuple
