# -*- coding: utf-8 -*-
import h5py

from tavi.data.scan import Scan
from tavi.data.scan_group import ScanGroup


class TAVI(object):
    """TAVI data file manager.

    TAVI_data contains four possible categories, including
    - data, a list of 1D scans, raws data.
    - process_data, including combined scan, 2D maps or dispersion plot
    - fit, contains fitting info including model, parameters and reduced chi_squared
    - plot, contains one or more scans and/or fits

    Attributes:
        hdf5_path: save path to hdf5
        self.scans: list of Scan instances

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
        with h5py.File(file_path, "w") as root:
            root.create_group("data")
            root.create_group("processed_data")
            root.create_group("fits")
            root.create_group("plots")

    def load_nexus_data_from_disk(self, path_to_hdf5):
        """Load hdf5 data from path_to_hdf5.

        Args:
            path_to_hdf5 (str): path to hdf5 data file
            OVERWRITE (bool): overwrite exsiting data if Ture, oterwise append new scans in data.
                Do not change processed_data, fit or plot
        """
        # TODO check if file exsits

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

    def load_spice_data_from_disk(self, path_to_spice_folder, OVERWRITE=True):
        """Load hdf5 data from path_to_hdf5.

        Args:
            path_to_spice_folder (str): path to spice folder
            OVERWRITE (bool): overwrite exsiting data if Ture, oterwise append new scans in data.
                Do not change processed_data, fit or plot
        """

        # exp_info = [
        #     "experiment",
        #     "experiment_number",
        #     "proposal",
        #     "users",
        #     "local_contact",
        # ]

        # p = Path(path_to_spice_folder)
        # scans = sorted((p / "Datafiles").glob("*"))
        # instrument, exp = scans[0].parts[-1].split("_")[0:2]

        # # read in exp_info from the first scan and save as attibutes of the file
        # _, _, headers, _ = read_spice(scans[0])
        # ipts = headers["proposal"]
        # data_id = f"IPTS{ipts}_{instrument}_{exp}"  # e.g. "IPTS1234_HB3_exp567"
        # data_entries = []

        # for scan in scans:  # ignoring unused keys
        #     spice_data, col_headers, headers, unused = read_spice(scan)
        #     scan_id = ((scan.parts[-1].split("_"))[-1]).split(".")[0]  # e.g. "scan0001"

        #     meta_data = {"scan_id": scan_id}

        #     for k, v in headers.items():
        #         if k not in exp_info:  # ignore common keys in single scans
        #             if "," in v and k != "scan_title":  # vectors
        #                 meta_data.update({k: np.array([float(v0) for v0 in v.split(",")])})
        #             elif v.replace(".", "").isnumeric():  # numebrs only
        #                 if v.isdigit():  # int
        #                     meta_data.update({k: int(v)})
        #                 else:  # float
        #                     meta_data.update({k: float(v)})
        #             # separate COM/FWHM and its errorbar
        #             elif k == "Center of Mass":
        #                 com, e_com = v.split("+/-")
        #                 meta_data.update({"COM": float(com)})
        #                 meta_data.update({"COM_err": float(e_com)})
        #             elif k == "Full Width Half-Maximum":
        #                 fwhm, e_fwhm = v.split("+/-")
        #                 meta_data.update({"FWHM": float(fwhm)})
        #                 meta_data.update({"FWHM_err": float(e_fwhm)})
        #             else:  # other crap, keep as is
        #                 if k not in exp_info:
        #                     meta_data.update({k: v})

        #     data = {}

        #     if spice_data.ndim == 1:  # empty data or 1 point only
        #         if len(spice_data):  # 1 point only
        #             for idx, col_header in enumerate(col_headers):
        #                 data.update({col_header: spice_data[idx]})
        #         else:  # empty
        #             pass
        #     else:  # nomarl data
        #         for idx, col_header in enumerate(col_headers):
        #             data.update({col_header: spice_data[:, idx]})

        #     s = Scan()
        #     s.set_metadata(meta_data)
        #     s.set_data(data)

        #     data_entries.append(s)

        # self.data.update({data_id: data_entries})

    def load_data_from_oncat(self, user_credentials, ipts_info, OVERWRITE=True):
        """Load data from ONCat based on user_credentials and ipts_info.

        Args:
            user_credentials (str): username and password.
            scipts_infoans (str): ipts number and exp number.
            OVERWRITE (bool): overwrite exsiting data if Ture, oterwise append.
                Do not change processed_data, fit or plot
        """
        pass

    def open_tavi_file(self, file_path):
        """Open existing tavi file"""
        self.file_path = file_path
        with h5py.File(file_path, "a") as tavi_file:
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

    def save_to_file(self, path_to_hdf5):
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

        # """Load data from spice folder.

        # Args:
        #     path_to_spice_folder (str): spice folder, ends with '/'
        #     path_to_hdf5 (str): path to hdf5 data file, ends with '.h5'
        # """
        # exp_info = [
        #     "experiment",
        #     "experiment_number",
        #     "proposal",
        #     "users",
        #     "local_contact",
        # ]

        # p = Path(path_to_spice_folder)

        # with h5py.File(path_to_hdf5, "w") as f:

        #     scans = sorted((p / "Datafiles").glob("*"))
        #     instrument_str, exp_str = scans[0].parts[-1].split("_")[0:2]

        #     # read in exp_info from the first scan and save as attibutes of the file
        #     _, _, headers, _ = read_spice(scans[0])
        #     ipts = headers["proposal"]
        #     exp_id = "IPTS" + ipts + "_" + instrument_str

        #     grp_data = f.create_group("data_" + exp_id)
        #     grp_processed_data = f.create_group("processed_data")
        #     grp_fit = f.create_group("fits")
        #     grp_plot = f.create_group("plots")

        #     for k, v in headers.items():
        #         if k in exp_info:
        #             grp_data.attrs[k] = v

        #     # read scans into dataset1
        #     for scan in scans:  # ignoring unused keys
        #         spice_data, col_headers, headers, unused = read_spice(scan)

        #         scan_num = ((scan.parts[-1].split("_"))[-1]).split(".")[0]
        #         scan_id = exp_str + "_" + scan_num
        #         scan_entry = grp_data.create_group(scan_id)
        #         scan_entry.attrs["scan_id"] = scan_id

        #         for k, v in headers.items():
        #             if k not in exp_info:  # ignore common keys in single scans
        #                 if "," in v and k != "scan_title":  # vectors
        #                     scan_entry.attrs[k] = np.array([float(v0) for v0 in v.split(",")])
        #                 elif v.replace(".", "").isnumeric():  # numebrs only
        #                     if v.isdigit():  # int
        #                         scan_entry.attrs[k] = int(v)
        #                     else:  # float
        #                         scan_entry.attrs[k] = float(v)
        #                 # separate COM/FWHM and its errorbar
        #                 elif k == "Center of Mass":
        #                     com, e_com = v.split("+/-")
        #                     scan_entry.attrs["COM"] = float(com)
        #                     scan_entry.attrs["COM_err"] = float(e_com)
        #                 elif k == "Full Width Half-Maximum":
        #                     fwhm, e_fwhm = v.split("+/-")
        #                     scan_entry.attrs["FWHM"] = float(fwhm)
        #                     scan_entry.attrs["FWHM_err"] = float(e_fwhm)
        #                 else:  # other crap, keep as is
        #                     if k not in exp_info:
        #                         scan_entry.attrs[k] = v

        #         if spice_data.ndim == 1:  # empty data or 1 point only
        #             if len(spice_data):  # 1 point only
        #                 for idx, col_header in enumerate(col_headers):
        #                     scan_entry.create_dataset(col_header, data=spice_data[idx])
        #             else:  # empty
        #                 pass
        #         else:  # nomarl data
        #             for idx, col_header in enumerate(col_headers):
        #                 scan_entry.create_dataset(col_header, data=spice_data[:, idx])

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
