# -*- coding: utf-8 -*-
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

    def __init__(self):
        """Initialization"""
        self.file_path = None
        self.data = {}
        self.processed_data = {}
        self.fits = {}
        self.plots = {}

    def new_tavi_file(self, file_path):
        """Create a new tavi file"""
        self.file_path = file_path

        try:
            with h5py.File(file_path, "w", track_order=True) as root:
                root.create_group("data", track_order=True)
                root.create_group("processed_data", track_order=True)
                root.create_group("fits", track_order=True)
                root.create_group("plots", track_order=True)
        except OSError:
            print(f"Cannot create tavi file at {file_path}")

    def load_nexus_data_from_disk(self, path_to_hdf5):
        """Load hdf5 data from path_to_hdf5.

        Args:
            path_to_hdf5 (str): path to hdf5 data file
            OVERWRITE (bool): overwrite exsiting data if Ture, oterwise append new scans in data.
                Do not change processed_data, fit or plot
        """
        # TODO validate path

        with h5py.File(self.file_path, "a") as tavi_file, h5py.File(path_to_hdf5, "r") as data_file:
            # IPTS1234_HB3_exp567
            data_id = data_file.attrs["file_name"].split("/")[-1]
            grp = tavi_file["data"].create_group(data_id)

            scans = {}
            for entry in data_file:
                data_file.copy(source=data_file[entry], dest=grp, expand_soft=True)

                if entry[0:4] == "scan":
                    s = Scan(data_file[entry])
                    scans.update({entry: s})

            self.data.update({data_id: scans})

    # TODO
    def load_spice_data_from_disk(self, path_to_spice_folder):
        """Load hdf5 data from path_to_hdf5.

        Args:
            path_to_spice_folder (str): path to spice folder
            OVERWRITE (bool): overwrite exsiting data if Ture, oterwise append new scans in data.
                Do not change processed_data, fit or plot
        """

        pass

    def load_data_from_oncat(self, user_credentials, ipts_info, OVERWRITE=True):
        """Load data from ONCat based on user_credentials and ipts_info.

        Args:
            user_credentials (str): username and password.
            scipts_infoans (str): ipts number and exp number.
            OVERWRITE (bool): overwrite exsiting data if Ture, oterwise append.
                Do not change processed_data, fit or plot
        """
        pass

    def open_tavi_file(self, tavi_file_path):
        """Open existing tavi file"""
        # TODO validate path
        self.file_path = tavi_file_path

        with h5py.File(tavi_file_path, "a") as tavi_file:
            # load datasets in data folder
            for data_id in tavi_file["data"].keys():
                dataset = tavi_file["data"][data_id]
                scans = {}
                for entry in dataset:
                    if entry[0:4] == "scan":
                        s = Scan(dataset[entry])
                        scans.update({entry: s})
                self.data.update({data_id: scans})

            # TODO
            # load processed_data fits, and plots

    def save_tavi_file(self, path_to_hdf5):
        """Save current data to a hdf5 on disk at path_to_hdf5

        Args:
            path_to_hdf5 (str): path to hdf5 data file, ends with '.h5'
        """
        pass

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
