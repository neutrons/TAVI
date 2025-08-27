import logging
import os
from dataclasses import field, make_dataclass
from typing import Any, Iterable, Optional

import numpy as np

import tavi.tavi_model.FileSystem.spice_reader as spice_reader
from tavi.tavi_model.FileSystem.tavi_class_factory import Scan

logger = logging.getLogger("TAVI")


class LoadORNL:
    """
    Loader for ORNL SPICE triple-axis spectrometer data.

    This class provides utilities to read SPICE-formatted data files produced
    at ORNL instruments (e.g., CG4C, HB1, HB3, HB1A). It supports loading
    either an explicit list of data files or all files within a specified
    directory. Metadata, numeric arrays, and UB matrix configurations are
    parsed into structured dataclasses (`RawData`, `RawMetaData`, `UbConf`)
    and organized into a `TaviProject`.

    args:
        data_folder (str | PathLike):
            Directory containing SPICE data files. Used if `data_files` is not provided.
        data_files (Iterable[str] | None):
            Optional explicit list of filenames to load. If `None`, all files
            in `data_folder` are processed.
        ub_dir (str | PathLike | None):
            Optional path to a UBConf directory. If not provided, the loader
            searches for a `UBConf` folder adjacent to `data_folder`.

    Methods:
        score() -> float:
            Return a heuristic score for how likely the directory corresponds
            to ORNL data. Currently returns a placeholder value (100).
        load() -> TaviProject:
            Public entry point to load data files into a `TaviProject`.
        _load_files() -> TaviProject:
            Internal method implementing the actual parsing and data
            organization logic.
    """

    def __init__(
        self,
        data_folder: Optional[os.PathLike | str] = None,
        data_files: Optional[Iterable[str]] = None,
        ub_dir: Optional[os.PathLike | str] = None,
    ) -> None:
        self.data_folder = data_folder
        self.data_files = data_files
        self.ub_dir = ub_dir

    def score(self) -> float:
        return 100

    def load(self):
        """
        Load SPICE data files into a TaviProject.

        This method serves as the main entry point for reading triple-axis
        spectrometer (SPICE) data. It gathers data files from either:

        - `self.data_files` if a list of files is explicitly provided, or
        - all files in `self.data_folder` if no list is given.

        Each file is parsed into numeric arrays and metadata, dynamically mapped
        onto `RawData`, `RawMetaData`, and (if applicable) `UbConf` dataclasses.
        These are assembled into `Scan` objects, which are collected into a
        `TaviProject`.

        Returns:
            TaviProject: A project object containing all loaded scans, indexed by filename.
        """
        return self._load_files()

    def _load_files(self):
        """
        Internal method to load SPICE data files into a TaviProject.

        This function reads each SPICE data file from either the provided list
        of `data_files` or from the `data_folder` directory. For each file, it:

        1. Uses `spice_reader.read_spice_datafile()` to extract:
           - numeric data arrays
           - column names
           - metadata
           - additional information
           - error messages

        2. Dynamically generates `RawData` and `RawMetaData` dataclasses with
           attributes inferred from the file contents. Attribute names are
           sanitized (e.g., replacing `-`, `.`, or whitespace with `_`) and
           guarded against invalid Python identifiers.

        3. Handles one-dimensional numeric data gracefully by converting scalars
           to 1D numpy arrays.

        4. Processes metadata entries. If a "ubconf" key is present, the
           corresponding UB matrix configuration file is searched for in the
           `ub_dir`, or in a `UBConf` subfolder adjacent to the data folder.
           If found, the UB configuration is read using
           `spice_reader.read_spice_ubconf()` and stored in a dynamically
           generated `UbConf` dataclass.

        5. Constructs a `Scan` object containing:
           - the populated `RawData` instance
           - the populated `RawMetaData` instance
           - the UB configuration (if available)
           - column names
           - error messages
           - other auxiliary information

        6. Stores the `Scan` object in a `TaviProject`, keyed by the filename.

        Returns:
            TaviProject: A project object containing all scans indexed by filename.
        """
        tavi_project = {}
        list_of_files = self.data_files if self.data_files else os.listdir(self.data_folder)
        for filename in list_of_files:
            rawdata = make_dataclass("RawData", [], slots=True)
            rawmetadata = make_dataclass("RawMetaData", [], slots=True)
            ub_conf = make_dataclass("UbConf", [], slots=True)

            numeric_data, col_names, meta_data, others, error_message = spice_reader.read_spice_datafile(
                os.path.join(self.data_folder, filename)
            )

            # ------------------load data------------------
            for col_name in col_names:
                # guard against invalid format
                if col_name[0].isdigit():
                    col_name = "_" + col_name
                attr_name = col_name.replace("-", "_").replace(" ", "_").replace(".", "")
                try:
                    # adding a checker to avoid re-creating class attributes, make the code faster
                    if hasattr(rawdata, attr_name):
                        setattr(rawdata, attr_name, numeric_data[:, col_names.index(col_name)])
                    else:
                        rawdata = make_dataclass(
                            "RawData", fields=[(attr_name, np.ndarray, field(default=None))], bases=(rawdata,)
                        )
                        setattr(rawdata, attr_name, numeric_data[:, col_names.index(col_name)])

                # if the numeric_data only has 1 entry, we'd still like to parse it as a 1D numpy array
                except IndexError:
                    rawdata = make_dataclass(
                        "RawData", fields=[(attr_name, np.ndarray, field(default=None))], bases=(rawdata,)
                    )
                    setattr(rawdata, attr_name, np.array([numeric_data[col_names.index(col_name)]]))

            rawdata = make_dataclass(
                            "RawData", fields=[("column_names", list, field(default=None))], bases=(rawdata,)
                        )
            setattr(rawdata, "column_names", col_names)

            # ------------------load meta data and ubconf------------------
            for key, value in meta_data.items():
                # replace "-" or " " with "_" to consolidate with python attribute's format
                if key[0].isdigit():
                    key = "_" + key
                key = key.replace("-", "_").replace(" ", "_").replace(".", "")
                # adding a checker to avoid re-creating class attributes, make the code faster
                if hasattr(rawmetadata, key):
                    setattr(rawmetadata, key, value)
                else:
                    rawmetadata = make_dataclass(
                        "RawMetaData", fields=[(key, str, field(default=None))], bases=(rawmetadata,)
                    )
                    setattr(rawmetadata, key, value)

                match key:
                    case "ubconf":
                        ub_filename = value
                        if ub_filename:
                            # users don't need to set specific ub_directory, the loader will try to search for ub
                            if not self.ub_dir:
                                prev_dir = os.path.join(self.data_folder, os.pardir)
                                if "UBConf" in os.listdir(prev_dir):
                                    self.ub_dir = os.path.join(prev_dir, "UBConf")
                            # look into "UBConf" folder first
                            if ub_filename in os.listdir(self.ub_dir):
                                ub_data = spice_reader.read_spice_ubconf(os.path.join(self.ub_dir, ub_filename))
                                for key, value in ub_data.items():
                                    # adding a checker to avoid re-creating class attributes, make the code faster
                                    if hasattr(ub_conf, key):
                                        setattr(ub_conf, key, value)
                                    else:
                                        ub_conf = make_dataclass(
                                            "UbConf", fields=[(key, Any, field(default=None))], bases=(ub_conf,)
                                        )
                                        setattr(ub_conf, key, value)
                            # look into tmp folder just in case
                            elif ub_filename in os.listdir(os.path.join(self.ub_dir, "tmp")):
                                ub_data = spice_reader.read_spice_ubconf(os.path.join(self.ub_dir, ub_filename))
                                for key, value in ub_data.items():
                                    # adding a checker to avoid re-creating class attributes, make the code faster
                                    if hasattr(ub_conf, key):
                                        setattr(ub_conf, key, value)
                                    else:
                                        ub_conf = make_dataclass(
                                            "UbConf", fields=[(key, Any, field(default=None))], bases=(ub_conf,)
                                        )
                                        setattr(ub_conf, key, value)
                            else:
                                logger.warning("Can't find %s, please double check UBMatrix data", ub_filename)
            
            rawmetadata = make_dataclass(
                        "RawMetaData", fields=[("others", str, field(default=None))], bases=(rawmetadata,)
                    )
            setattr(rawmetadata, "others", others)

            scan = Scan(
                data=rawdata,
                metadata=rawmetadata,
                error_message=error_message,
                ubconf=ub_conf,
            )
            tavi_project[filename] = scan
        return tavi_project

