# -*- coding: utf-8 -*-
import os
from typing import Optional

import h5py

from tavi.data.scan import Scan
from tavi.data.scan_group import ScanGroup


class TAVI(object):
    """TAVI data file manager.

    TAVI_data contains four possible categories, including
    - data, a list of 1D scans, raws data.
    - process_data, including combined scan, 2D maps or dispersion plot
    - fits, contains fitting info including model, parameters and reduced chi_squared
    - plots, contains one or more scans and/or fits

    Attributes:
        file_path: path to a tavi file


    """

    def __init__(self) -> None:
        """Initialization"""
        self.file_path: Optional[str] = None
        self.data: dict = {}
        self.processed_data: dict = {}
        self.fits: dict = {}
        self.plots: dict = {}

    def new_file(self, file_path: Optional[str] = None) -> None:
        """Create a new tavi file
        Note:
            create at file_path if given, otherwise create a temporary
            file named tavi_temp.h5 at the current folder.
        """
        if file_path is None:
            current_path = os.getcwd()
            self.file_path = current_path + "/tavi_temp.h5"
        else:
            self.file_path = file_path

        try:
            with h5py.File(self.file_path, "w", track_order=True) as root:
                root.create_group("data", track_order=True)
                root.create_group("processed_data", track_order=True)
                root.create_group("fits", track_order=True)
                root.create_group("plots", track_order=True)
        except OSError:
            print(f"Cannot create tavi file at {self.file_path}")

    def get_nexus_data_from_disk(self, path_to_hdf5_folder):
        """Copy hdf5 data from path_to_hdf5_folder into tavi."""
        # validate path
        if path_to_hdf5_folder[-1] != "/":
            path_to_hdf5_folder += "/"

        scan_list = os.listdir(path_to_hdf5_folder)
        scan_list = [scan for scan in scan_list if scan.startswith("scan")]
        scan_list.sort()

        with h5py.File(self.file_path, "r+", track_order=True) as tavi_file:
            for scan in scan_list:
                scan_name = scan.split(".")[0]  # e.g. "scan0001"
                scan_path = path_to_hdf5_folder + scan
                with h5py.File(scan_path, "r") as scan_file:
                    dataset_name = scan_file[scan_name].attrs["dataset_name"]
                    if dataset_name not in tavi_file["/data"]:
                        tavi_file.create_group("/data/" + dataset_name)
                    scan_file.copy(
                        source=scan_file["/" + scan_name], dest=tavi_file["/data/" + dataset_name], expand_soft=True
                    )

        self.load_data()

    # TODO
    def get_spice_data_from_disk(self, path_to_spice_folder):
        """Load hdf5 data from path_to_hdf5.

        Args:
            path_to_spice_folder (str): path to spice folder
            OVERWRITE (bool): overwrite exsiting data if Ture, oterwise append new scans in data.
                Do not change processed_data, fit or plot
        """

        self.load_data()

    # TODO
    def get_data_from_oncat(self, user_credentials, ipts_info, OVERWRITE=True):
        """Load data from ONCat based on user_credentials and ipts_info.

        Args:
            user_credentials (str): username and password.
            scipts_infoans (str): ipts number and exp number.
            OVERWRITE (bool): overwrite exsiting data if Ture, oterwise append.
                Do not change processed_data, fit or plot
        """
        self.load_data()

    def load_data(self):
        """Load tavi data into memory"""

        with h5py.File(self.file_path, "r") as tavi_file:
            for dataset_name in (tavi_data := tavi_file["data"]):
                scans = {}
                for scan_name in (dataset := tavi_data[dataset_name]):
                    nexus_entry = dataset[scan_name]
                    _, scan_entry = Scan.from_nexus_entry(nexus_entry)
                    scans.update({scan_name: scan_entry})

                self.data.update({dataset_name: scans})

    def open_file(self, tavi_file_path):
        """Open existing tavi file"""
        # TODO validate path
        self.file_path = tavi_file_path
        self.load_data()

        # TODO load processed_data fits, and plots

    # TODO
    def save_file(self, path_to_hdf5: Optional[str] = None):
        """Save current data to a hdf5 on disk at path_to_hdf5

        Args:
            path_to_hdf5 (str): path to hdf5 data file, ends with '.h5'
        """
        if path_to_hdf5 is not None:
            self.file_path = path_to_hdf5

        # TODO save processed_data fits and plots

    def delete_data_entry(self, entry_id):
        """Delete a date entry from file based on entry_id.

        If an data entry if deleted, also delete the fits, scan_groups, plots
        that contains the conresponding entry.

        Args:
            entry_id (str): ID of entry to be deteled.
        """
        pass

    def select_entries(self, conditions):
        """Return a list of data entry IDs which satisfis the conditions.

        All data points in a scan have to satisfy all conditions.

        Args:
            conditions (dict): scan[key] equals value
                or in the range of (vaule[0], value[1])
                for example: {"IPTS":1234, "temperature":(2,5)}

        Returns:
            tuple: entry IDs

        """
        pass

    def get_selected(self):
        return "scan0001"

    def generate_scan_group(
        self,
        signals=None,
        backgrounds=None,
        signal_axes=(None, None, None),
        background_axes=(None, None, None),
    ):
        """Generate a scan group."""
        sg = ScanGroup(signals, backgrounds, signal_axes, background_axes)

        return sg
