import json
import logging
import os
from dataclasses import dataclass, field, make_dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Type

import numpy as np

import tavi.tavi_model.FileSystem.spice_reader as spice_reader

logger = logging.getLogger("TAVI")
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


def load_folder(dir):
    """
    Load SPICE data files from a directory into a TaviProject.

    This function reads each file in the specified directory using
    `spice_reader.read_spice_datafile()`, converts the numeric arrays and
    metadata into `RawData` and `RawMetaData` dataclass instances, and stores
    them as `Scan` objects in a `TaviProject`.

    Args:
        directory (str | Path): Path to the folder containing SPICE data files.

    Returns:
        TaviProject: A project object containing all loaded scans.
    """

    # Initialize TaviProject class
    tavi_project = TaviProject()
    data_dir = os.path.join(dir, "Datafiles")
    # load files into TaviProject scans
    for filename in os.listdir(data_dir):
        numeric_data, col_names, meta_data, others, error_message = spice_reader.read_spice_datafile(
            os.path.join(data_dir, filename)
        )
        rawdata = RawData()
        rawmetadata = RawMetaData()

        for col_name in col_names:
            attr_name = "Pt" if col_name == "Pt." else col_name
            if hasattr(rawdata, attr_name):
                try:
                    setattr(rawdata, attr_name, numeric_data[:, col_names.index(col_name)])
                # if the numeric_data only has 1 entry, we'd still like to parse it as a 1D numpy array
                except IndexError:
                    setattr(rawdata, attr_name, np.array([numeric_data[col_names.index(col_name)]]))
            else:
                logger.warning("New data entry detected, consider updating RawData entry in tavi_data_schema.json")

        for key, value in meta_data.items():
            # replace "-" or " " with "_" to consolidate with python attribute's format
            key = key.replace("-", "_").replace(" ", "_")
            if hasattr(rawmetadata, key):
                setattr(rawmetadata, key, value)
            else:
                logger.warning(
                    "New meta data entry detected, consider updating RawMetaData entry in tavi_data_schema.json"
                )

        scan = Scan(
            data=rawdata,
            metadata=rawmetadata,
            column_names=col_names,
            error_message=error_message,
            others=others,
        )
        tavi_project.scans[filename] = scan

    tmp_exist = False
    ub_dir = os.path.join(dir, "UBConf")
    for ub_filename in os.listdir(ub_dir):
        ub_conf = UbConf()
        if ub_filename == "tmp":
            tmp_exist = True
            continue
        ub_data = spice_reader.read_spice_ubconf(os.path.join(ub_dir, ub_filename))
        for key, value in ub_data.items():
            if hasattr(ub_conf, key):
                setattr(ub_conf, key, value)
            else:
                logger.warning("New UbConf found, consider updating UbConf entry in tavi_data_schema.json")
        tavi_project.ubconf[ub_filename] = ub_conf

    if tmp_exist:
        for ub_filename in os.listdir(os.path.join(ub_dir, "tmp")):
            ub_conf = UbConf()
            ub_data = spice_reader.read_spice_ubconf(os.path.join(ub_dir, "tmp", ub_filename))
            for key, value in ub_data.items():
                if hasattr(ub_conf, key):
                    setattr(ub_conf, key, value)
                else:
                    logger.warning("New UbConf found, consider updating UbConf entry in tavi_data_schema.json")
        tavi_project.ubconf["tmp-" + ub_filename] = ub_conf

    return tavi_project


if __name__ == "__main__":
    current_directory = os.getcwd()
    filepath = os.path.join(current_directory, "test_data", "exp424")
    tavi_project = load_folder(filepath)
    ubpath = os.path.join(current_directory, "test_data", "exp424", "UBConf", "UB02Jul2024_14108PM.ini")
    ub = spice_reader.read_spice_ubconf(ubpath)
    # print(ub)
    # print(tavi_project.scans["CG4C_exp0424_scan0073.dat"].metadata)
    print(tavi_project.ubconf["tmp-UB02Jul2024_21029PM.ini"])
