import os
from typing import Iterable, Optional

from tavi.tavi_model.combine_data import CombineManager
from tavi.tavi_model.FileSystem.load_manager import LoadManager
from tavi.tavi_model.filter import Filter, Logic, Operations
from tavi.tavi_model.tavi_data import TaviData


class TaviProject:
    """
    Central container for managing a full TAVI project workflow.

    The `TaviProject` class provides a high-level interface for loading,
    filtering, processing, and analyzing triple-axis spectrometer (TAS) data.
    It acts as the core object for managing raw scans, metadata, selected data,
    and later processing steps such as combining, fitting, and plotting.

    Internally, `TaviProject` owns a `TaviData` object, which stores and
    organizes all scan-related information (raw and processed).

    Typical usage involves:
        1. Loading raw scan files into the project with `load_scans`.
        2. Selecting subsets of scans or points with `select_scans`.
        3. Performing transformations such as combining or fitting.
        4. Visualizing or exporting results.

    Attributes:
        tavi_data (TaviData):
            The main data container that stores scan lists, processed data,
            and user selections (for viewing or modeling).

    Notes:
        - The `load_scans` method uses a `LoadManager` backend for parsing
          facility-specific file formats (e.g., ORNL, ILL, ISIS).
        - Scan filtering is handled by the `Filter` class with logical
          conditions and operations defined in `tavi_model.filter`.
        - Future methods (`load_tavi`, `save_tavi`, `combine_data`,
          `fit_data`, `plot_data`) will extend project persistence,
          advanced analysis, and visualization.
    """

    def __init__(self):
        self.tavi_data = TaviData()

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

        Attributes:
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
        self.tavi_data.rawdataptr = LoadManager(
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
        filter_name: Optional[str] = None,
        conditions: Optional[list[tuple[str, Operations, str | float]]] = None,
        and_or: Optional[Logic] = None,
        category: Optional[str] = None,
        tol=0.01,
    ) -> None:
        filtered_data = Filter(self.tavi_data.scan_list, conditions=conditions, and_or=and_or, tol=tol).filter_data()
        match category:
            case "view":
                self.tavi_data.show_selected_data[filter_name] = filtered_data
            case "model":
                self.tavi_data.process_selected_data = filtered_data
            case _:
                self.tavi_data.show_selected_data[filter_name] = filtered_data

    # TO DO
    def combine_data(
        self,
        target_list: list[str],
        background_list: Optional[list[str]] = [],
        axis: Optional[tuple[str, str]] = None,
        tol: Optional[float] = 0.01,
    ):
        target = [self.tavi_data.rawdataptr[scan_name] for scan_name in target_list]
        background = [self.tavi_data.rawdataptr[scan_name] for scan_name in background_list]
        combined_data_1d = CombineManager(target=target, background=background).combine_1d(axis)
        return combined_data_1d

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

    # filename = "CG4C_exp0424_scan0042.dat"
    # TaviProj.select_scans(
    #     filter_name="scan_contains_42", conditions=([["scan", Operations.CONTAINS, "42"]]), and_or=Logic.OR
    # )

    # TaviProj.select_scans(filter_name="filter2", conditions=([["scan", Operations.CONTAINS, "4"]]), and_or=Logic.OR)
    # print(TaviProj.tavi_data.show_selected_data)
    #   print(type(TaviProj.scans[filename].metadata.scan))
    #   print(TaviProj.scans[filename].ubconf)
    #   print(TaviProj.scans[filename].data.Pt)
    #   print(TaviProj.scans[filename].error_message)
    #   print(TaviProj.scans[filename].metadata.time)

    # -----------------------combine data---------------------------
    target = ["CG4C_exp0424_scan0042.dat", "CG4C_exp0424_scan0042.dat"]
    test_return = TaviProj.combine_data(target_list=target, axis=("e", "detector"))
    print(test_return[0])
    print(test_return[1])
    print(test_return[2])
    print(test_return[3])
    print(test_return[4])
    print(test_return[5])
    print(len(test_return))
    # import matplotlib.pyplot as plt
    # plt.plot(test_return[3], test_return[4], '.')
    # plt.show()