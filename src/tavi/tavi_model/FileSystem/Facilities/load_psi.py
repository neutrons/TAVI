import logging
import os
from typing import Iterable, Optional

logger = logging.getLogger("TAVI")


class LoadPSI:
    """
    Loader for NIST triple-axis spectrometer data.

    This class provides utilities to read SPICE-formatted data files produced
    at NIST. It supports loading
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
        data_folder: Optional[os.PathLike | str],
        data_files: Optional[Iterable[str]] = None,
        ub_dir: Optional[os.PathLike | str] = None,
    ) -> None:
        self.data_folder = data_folder
        self.data_files = data_files
        self.ub_dir = ub_dir

    def score(self) -> float:
        return 0

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
        pass
