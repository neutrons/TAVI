import json
import os
from dataclasses import dataclass, field, make_dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Type

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


def generate_dataclass(class_name: str, schema: dict, use_slots=True) -> Type[Any]:
    """
    Dynamically create a dataclass from a schema.

    Args:
        class_name: Name of the generated class.
        schema: Dictionary mapping field names to type strings.
        use_slots: Whether to use __slots__ for the generated class.

    Returns:
        A new dataclass type.
    """
    fields = [(name, TYPE_MAP.get(type_str, DEFAULT_TYPE), field(default=None)) for name, type_str in schema.items()]
    return make_dataclass(class_name, fields, slots=use_slots)


def load_schema(filename: str | Path) -> Dict[str, Any]:
    """
    Load a JSON schema file into a Python dictionary.

    Args:
        filename: Path to the JSON file.

    Returns:
        The parsed JSON object as a dictionary.
    """
    path = Path(filename)
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


dir_path = os.path.dirname(os.path.realpath(__file__))
schema_dir = os.path.join(dir_path, "tavi_data_schema.json")
schema = load_schema(schema_dir)

RawMetaData = generate_dataclass("RawMetaData", schema["RawMetaData"])
RawData = generate_dataclass("RawData", schema["RawData"])
UbConf = generate_dataclass("UbConf", schema["UbConf"])


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

    data: RawData
    metadata: RawMetaData
    column_names: tuple
    error_message: tuple
    others: tuple


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
    ubconf: dict[str, UbConf] = field(default_factory=dict)
