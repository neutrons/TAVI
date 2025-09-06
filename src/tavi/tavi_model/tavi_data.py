import os
from dataclasses import field
from typing import Any, Iterable, Optional

import numpy as np

from tavi.tavi_model.FileSystem.load_manager import LoadManager
from tavi.tavi_model.FileSystem.tavi_class_factory import Scan
from tavi.tavi_model.filter import Filter, Logic, Operations


class TaviData:
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
        scan_list: dict[str, Scan] = field(default_factory=dict),
        combined_data: dict[str, np.array] = field(default_factory=dict),
        process_selected_data: list[str] = [],  # mouse selection
        show_selected_data: dict = {},  # display 
        fit_manager: dict[str, Any] = field(default_factory=dict),
        plot_manager: dict[str, Any] = field(default_factory=dict),
    ):
        self.scan_list = scan_list
        self.combined_data = combined_data
        self.process_selected_data = process_selected_data
        self.show_selected_data = show_selected_data
        self.fit_manager = fit_manager
        self.plot_manager = plot_manager