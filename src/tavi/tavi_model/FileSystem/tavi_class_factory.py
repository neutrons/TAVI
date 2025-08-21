from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np

TYPE_MAP = {
    "str": Optional[str],
    "float": Optional[float],
    "int": Optional[int],
    "bool": Optional[bool],
    "list": Optional[list],
    "ndarray": Optional[np.ndarray],
}

DEFAULT_TYPE = Optional[str]


@dataclass
class Scan:
    """
    Represents a single scan within a Tavi project, containing both raw measurement
    data and associated metadata.

    Attributes:
        data (RawData): Numerical measurement arrays collected during the scan
            (e.g., motor positions, detector counts, temperatures).
        metadata (RawMetaData): Descriptive information about the scan
            (e.g., experiment details, instrument settings, sample information).
        column_names (tuple): Ordered names of data columns in `data`, typically
            matching measurement channels or parameters.
        error_message (tuple): Messages or warnings associated with the scan,
            such as instrument errors or data quality issues.
        others (tuple): Miscellaneous or auxiliary information related to the scan
            that does not fit into `data` or `metadata`.
    """

    data: Any
    metadata: Any
    column_names: tuple
    error_message: tuple
    others: tuple
    ubconf: Any


@dataclass
class TaviProject:
    """
    Represents a complete Tavi project containing one or more scans.

    Attributes:
        scans (dict[str, Scan]): A mapping of scan identifiers (e.g., scan numbers or
            unique labels) to their corresponding `Scan` objects. Each `Scan` holds
            both the raw measurement data and the associated metadata for that scan.
    """

    scans: dict[str, Scan] = field(default_factory=dict)
