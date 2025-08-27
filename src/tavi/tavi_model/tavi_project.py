import os
from dataclasses import field
from typing import Any, Iterable, Optional

import numpy as np

from tavi.tavi_model.FileSystem.load_manager import LoadManager
from tavi.tavi_model.FileSystem.tavi_class_factory import Scan


class TaviProject:
    """
    Represents a complete Tavi project containing one or more scans.

    Attributes:
        scans (dict[str, Scan]): A mapping of scan identifiers (e.g., scan numbers or
            unique labels) to their corresponding `Scan` objects. Each `Scan` holds
            both the raw measurement data and the associated metadata for that scan.
    """

    def __init__(
        self,
        scans: dict[str, Scan] = field(default_factory=dict),
        combined_data: dict[str, np.array] = field(default_factory=dict),
        filtered_data: dict[str, np.array] = field(default_factory=dict),
        view_selected_data: dict[str, Any] = field(default_factory=dict),  # view select, process select
        process_selected_data: dict[str, Any] = field(default_factory=dict),  # view select, process select
        fit_manager: dict[str, Any] = field(default_factory=dict),
        plot_manager: dict[str, Any] = field(default_factory=dict),
    ):
        self.scans = scans
        self.combined_data = combined_data
        self.filtered_data = filtered_data
        self.view_selected_data = view_selected_data
        self.process_selected_data = process_selected_data
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
        call back for loading scans
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
    def filter_scans():
        pass

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
    files = ["CG4C_exp0424_scan0041.dat", "CG4C_exp0424_scan0042.dat"]
    TaviProj = TaviProject()

    TaviProj.load_scans(filepath, files)

    filename = "CG4C_exp0424_scan0042.dat"
    print(TaviProj.scans[filename].data.h)
    print(TaviProj.scans[filename].ubconf)
    print(TaviProj.scans[filename].data.column_names)
    print(TaviProj.scans[filename].error_message)
    print(TaviProj.scans[filename].metadata.others)
