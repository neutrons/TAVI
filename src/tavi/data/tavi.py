# -*- coding: utf-8 -*-
import os
from typing import Optional, Union

import h5py

from tavi.data.nexus_entry import NexusEntry
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

    def __init__(self, tavi_file_path: Optional[str] = None) -> None:
        """Initialization"""
        self.file_path: Optional[str] = None
        self.data: dict = {}
        self.processed_data: dict = {}
        self.fits: dict = {}
        self.plots: dict = {}

        if tavi_file_path is not None:
            self.open_file(tavi_file_path)

    @staticmethod
    def load_data(scan_dict: dict):
        """Load data under coresponding exp_id"""
        data_dict: dict = {}
        for scan_name, nxentry in scan_dict.items():
            exp_id = nxentry["attrs"]["dataset_name"]
            if not data_dict.get(exp_id):  # first entry
                data_dict.update({exp_id: {scan_name: nxentry}})
            else:  # exp_id existing
                data_dict[exp_id].update({scan_name: nxentry})
        return data_dict

    def load_nexus_data_from_disk(self, path_to_hdf5_folder):
        """Copy hdf5 data from path_to_hdf5_folder into tavi."""
        # validate path
        if path_to_hdf5_folder[-1] != "/":
            path_to_hdf5_folder += "/"

        scan_list = os.listdir(path_to_hdf5_folder)
        scan_list = [scan for scan in scan_list if scan.startswith("scan")]
        scan_list.sort()

        scan_dict = {}
        for scan_file in scan_list:
            scan = NexusEntry.from_nexus(path_to_hdf5_folder + scan_file)
            scan_dict.update(scan)
        self.data.update(TAVI.load_data(scan_dict))

    def load_spice_data_from_disk(self, path_to_spice_folder):
        """Load hdf5 data from path_to_hdf5.

        Args:
            path_to_spice_folder (str): path to spice folder
        """
        scan_dict = NexusEntry.from_spice(path_to_spice_folder)
        self.data.update(TAVI.load_data(scan_dict))

    # TODO
    def load_data_from_oncat(self, user_credentials, ipts_info):
        """Load data from ONCat based on user_credentials and ipts_info.

        Args:
            user_credentials (str): username and password.
            scipts_infoans (str): ipts number and exp number.
        """
        pass

    def get_exp_id(self) -> list[str]:
        "return a list of exp id from a tavi file"
        exp_id_list = []
        with h5py.File(self.file_path, "r") as tavi_file:
            for exp_id in tavi_file["data"].keys():
                exp_id_list.append(exp_id)
        return exp_id_list

    def open_file(self, tavi_file_path):
        """Open an existing tavi file"""

        if not os.path.exists(tavi_file_path):
            raise ValueError(f"TAVI file does not exist at path {tavi_file_path}")

        self.file_path = tavi_file_path

        # load data
        data_dict = {}
        exp_id_list = self.get_exp_id()
        for exp_id in exp_id_list:
            nxentries = NexusEntry.from_nexus(tavi_file_path, prefix=f"/data/{exp_id}")
            data_dict.update({exp_id: nxentries})
        self.data = data_dict
        # TODO load processed_data fits, and plots

    def save(self, file_path: Optional[str] = None):
        """Save current data to a hdf5 on disk at path_to_hdf5

        Args:
            path_to_hdf5 (str): path to hdf5 data file, ends with '.h5'
        """

        if file_path is None:
            current_path = os.getcwd()
            self.file_path = current_path + "/tavi_temp.h5"
        else:
            self.file_path = file_path

        try:
            with h5py.File(self.file_path, "w", track_order=True) as root:
                grp_data = root.create_group("data", track_order=True)
                for exp_id, scans in self.data.items():
                    grp_data.create_group(name=exp_id, track_order=True)
                    for entry_name, nxentry in scans.items():
                        name = f"data/{exp_id}/{entry_name}"
                        nxentry.to_nexus(file_path, name=name)
                # TODO save pdata, fits and plots
                grp_pdata = root.create_group("processed_data", track_order=True)
                for entry_name, nxentry in self.processed_data.items():
                    name = f"processed_data/{entry_name}"
                    nxentry.to_nexus(file_path, name=name)

                grp_fits = root.create_group("fits", track_order=True)
                grp_plots = root.create_group("plots", track_order=True)
        except OSError:
            print(f"Cannot create tavi file at {self.file_path}")

    def get_scan(self, scan_num: Union[tuple[str, int], int]) -> Scan:
        """Get the scan at location /data/exp_id/scanXXXX, return a Scan instance

        Arguments:
            scan_path (str): exp_id /scan_name. exp_id is in the format of
                IPTSXXXXX_INSTRU_expXXXX, it is needed when
                more than one experiment is loaded as data
        Return:
            Scan: an instance of Scan class
        """
        if isinstance(scan_num, tuple):
            exp_id, scan_num = scan_num
            scan_name = f"scan{scan_num:04}"
        elif isinstance(scan_num, int):
            exp_id = next(iter(self.data))
            scan_name = f"scan{scan_num:04}"
        else:
            raise ValueError(f"scan_num={scan_num} needs to be tuple or int")
        return Scan(scan_name, self.data[exp_id][scan_name])

    def group_scans(
        self,
        scan_nums: Union[list[int], list[tuple[str, int]]],
        name: str = "",
    ) -> ScanGroup:
        """Group scans for further processing

        Arguments:
            scan_nums:
            scan_group_name:

        Returns:
            ScanGroup
        """
        first_exp_id = next(iter(self.data))
        scan_num_list = []
        for scan_item in scan_nums:
            if isinstance(scan_item, int):
                scan_num_list.append((first_exp_id, f"scan{scan_item:04}"))
            elif isinstance(scan_item, tuple):
                exp_id, scan_num = scan_item
                scan_num_list.append((exp_id, f"scan{scan_num:04}"))
            else:
                raise ValueError(f"scan_num={scan_item} needs to be tuple or int")

        scans = []
        for exp_id, scan_name in scan_num_list:
            scans.append(Scan(scan_name, self.data[exp_id][scan_name]))
        return ScanGroup(scans, name=name)
