from dataclasses import dataclass
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
    ubconf: Optional[Any]
    error_message: tuple
