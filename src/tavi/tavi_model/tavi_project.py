import os
from dataclasses import field
from typing import Any, Iterable, Optional

import numpy as np

from tavi.tavi_model.FileSystem.load_manager import LoadManager
from tavi.tavi_model.FileSystem.tavi_class_factory import Scan
from tavi.tavi_model.filter import Filter, Logic, Operations


class TaviProject:
    """
    A container class for managing a complete TAVI project, including
    raw TAS (triple-axis spectrometer) experimental data, processed data,
    filtering operations, fitting, and plotting management.

    The `TaviProject` serves as the central data structure for organizing
    experimental scans and coordinating workflows such as loading, filtering,
    combining, fitting, and visualizing data. It acts as a high-level interface
    to link raw data, metadata, and analysis pipelines in TAVI.

    Attributes:
        scans (dict[str, Scan]):
            A mapping from scan identifiers (e.g., filenames or scan numbers)
            to their corresponding `Scan` objects, which hold raw data and metadata.

        combined_data (dict[str, np.ndarray]):
            Stores merged or aggregated data from multiple scans for analysis.

        filtered_data (dict[str, np.ndarray]):
            Stores scan data after user-applied filters (e.g., Q/E cuts, masks).

        view_selected_data (dict[str, Any]):
            Holds user-selected data for visualization purposes.

        process_selected_data (dict[str, Any]):
            Holds user-selected data for processing workflows.

        fit_manager (dict[str, Any]):
            A container for fitting models, results, and related configurations.

        plot_manager (dict[str, Any]):
            A container for managing plotting sessions, styles, and figures.

    Methods:
        load_scans(data_folder, data_files, ub_dir=None, facility=None):
            Load scans from disk using the `LoadManager`. Populates `self.scans`.

        load_tavi():
            (Placeholder) Load a serialized TAVI project file.

        save_tavi():
            (Placeholder) Save the current TAVI project to disk.

        filter_scans():
            (Placeholder) Apply filters to raw scan data and store results.

        combine_data():
            (Placeholder) Combine multiple scans into `combined_data`.

        fit_data():
            (Placeholder) Perform fitting routines on selected datasets.

        plot_data():
            (Placeholder) Generate plots of scan or processed data.
    """

    def __init__(
        self,
        scans: dict[str, Scan] = field(default_factory=dict),
        combined_data: dict[str, np.array] = field(default_factory=dict),
        process_selected_data: list[str] = [],  # view select, process select
        view_selected_data: list[str] = [],  # view select, process select
        fit_manager: dict[str, Any] = field(default_factory=dict),
        plot_manager: dict[str, Any] = field(default_factory=dict),
    ):
        self.scans = scans
        self.combined_data = combined_data
        self.process_selected_data = process_selected_data
        self.view_selected_data = view_selected_data
        self.fit_manager = fit_manager
        self.plot_manager = plot_manager

    # --------------------Load Manager-------------------------------------
    def load_scans(
        self,
        data_folder: Optional[os.PathLike | str] = None,
        data_files: Optional[os.PathLike | str | Iterable[os.PathLike | str]] = None,
        ub_dir: Optional[os.PathLike] = None,
        facility: Optional[str] = None,
    ) -> None:
        """
        Load raw TAS (triple-axis spectrometer) scan files into the project.

        This method uses the `LoadManager` to read one or more scan files
        from the specified directory and populate the `scans` attribute with
        `Scan` objects containing raw data, metadata, and UB matrix information.

        Args:
            data_folder (os.PathLike | str, optional):
                Path to the directory containing the scan files.
                If not provided, defaults to the current working directory.

            data_files (os.PathLike | str | Iterable[os.PathLike | str], optional):
                One or more scan file names to load from the given `data_folder`.
                Can be a single filename, a `PathLike` object, or an iterable of file paths.

            ub_dir (os.PathLike, optional):
                Path to the directory containing UB matrix configuration files.
                Used to associate orientation matrices with the scans if available.

            facility (str, optional):
                Identifier for the facility (e.g., "ORNL", "ILL", "ISIS") that
                produced the data. Allows facility-specific loading behavior.

        Returns:
            None
                The method updates the `scans` attribute in place.
        """
        self.scans = LoadManager(
            data_folder=data_folder, data_files=data_files, ub_dir=ub_dir, facility=facility
        ).load()

    # TO DO
    def load_tavi():
        pass

    # TO DO
    def save_tavi():
        pass

    # TO DO
    def select_scans(
        self,
        conditions: Optional[list[tuple[str, Operations, str | float]]] = None,
        and_or: Optional[Logic] = None,
        category: Optional[str] = None,
        tol=0.01,
    ) -> None:
        filtered_data = Filter(self.scans, conditions=conditions, and_or=and_or).filter_data()
        match category:
            case "view":
                self.view_selected_data.append(filtered_data)
            case "model":
                self.process_selected_data.append(filtered_data)
            case _:
                self.process_selected_data.append(filtered_data)

    # TO DO
    def combine_data():
        pass

    # TO DO
    def fit_data():
        pass

    # TO DO
    def plot_data():
        pass


if __name__ == "__main__":
    current_directory = os.getcwd()
    filepath = os.path.join(current_directory, "test_data", "exp424", "Datafiles")
    # files = ["CG4C_exp0424_scan0041.dat", "CG4C_exp0424_scan0042.dat"]
    TaviProj = TaviProject()

    TaviProj.load_scans(filepath)

    filename = "CG4C_exp0424_scan0042.dat"
    TaviProj.select_scans(
        conditions=(["time", Operations.CONTAINS, "AM"], ["scan", Operations.IS, "42"]), and_or=Logic.AND
    )
    print(TaviProj.process_selected_data)
    # print(type(TaviProj.scans[filename].metadata.scan))
#     print(TaviProj.scans[filename].ubconf)
#     print(TaviProj.scans[filename].data.column_names)
#     print(TaviProj.scans[filename].error_message)
# print(TaviProj.scans[filename].metadata.time)
