import numpy as np
import h5py
from pathlib import Path
from scans import Scan


class TAVI_Data(object):
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
        self.file_name = None
        self.exp_id = None
        self.data = []
        self.processed_data = []
        self.fits = []
        self.plots = []

    def load_data_from_disk(self, path_to_hdf5, OVERWRITE=True):
        """Load hdf5 data from path_to_hdf5.

        Args:
            path_to_hdf5 (str): path to hdf5 data file
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

    @staticmethod
    def convert_spice_to_hdf5(path_to_spice_folder, path_to_hdf5):
        """Load data from spice folder.

        Args:
            path_to_spice_folder (str): spice folder, ends with '/'
            path_to_hdf5 (str): path to hdf5 data file, ends with '.h5'
        """
        exp_info = [
            "experiment",
            "experiment_number",
            "proposal",
            "users",
            "local_contact",
        ]

        p = Path(path_to_spice_folder)

        with h5py.File(path_to_hdf5, "w") as f:

            scans = sorted((p / "Datafiles").glob("*"))
            instrument_str, exp_str = scans[0].parts[-1].split("_")[0:2]

            # read in exp_info from the first scan and save as attibutes of the file
            _, _, headers, _ = TAVI_Data._read_spice(scans[0])
            ipts = headers["proposal"]
            exp_id = "IPTS" + ipts + "_" + instrument_str

            grp_data = f.create_group("data_" + exp_id)
            grp_processed_data = f.create_group("processed_data")
            grp_fit = f.create_group("fits")
            grp_plot = f.create_group("plots")

            for k, v in headers.items():
                if k in exp_info:
                    grp_data.attrs[k] = v

            # read scans into dataset1
            for scan in scans:  # ignoring unused keys
                spice_data, col_headers, headers, unused = TAVI_Data._read_spice(scan)

                scan_num = ((scan.parts[-1].split("_"))[-1]).split(".")[0]
                scan_id = exp_str + "_" + scan_num
                scan_entry = grp_data.create_group(scan_id)
                scan_entry.attrs["scan_id"] = scan_id

                for k, v in headers.items():
                    if k not in exp_info:  # ignore common keys in single scans
                        if "," in v and k != "scan_title":  # vectors
                            scan_entry.attrs[k] = np.array([float(v0) for v0 in v.split(",")])
                        elif v.replace(".", "").isnumeric():  # numebrs only
                            if v.isdigit():  # int
                                scan_entry.attrs[k] = int(v)
                            else:  # float
                                scan_entry.attrs[k] = float(v)
                        # separate COM/FWHM and its errorbar
                        elif k == "Center of Mass":
                            com, e_com = v.split("+/-")
                            scan_entry.attrs["COM"] = float(com)
                            scan_entry.attrs["COM_err"] = float(e_com)
                        elif k == "Full Width Half-Maximum":
                            fwhm, e_fwhm = v.split("+/-")
                            scan_entry.attrs["FWHM"] = float(fwhm)
                            scan_entry.attrs["FWHM_err"] = float(e_fwhm)
                        else:  # other crap, keep as is
                            if k not in exp_info:
                                scan_entry.attrs[k] = v

                if spice_data.ndim == 1:  # empty data or 1 point only
                    if len(spice_data):  # 1 point only
                        for idx, col_header in enumerate(col_headers):
                            scan_entry.create_dataset(col_header, data=spice_data[idx])
                    else:  # empty
                        pass
                else:  # nomarl data
                    for idx, col_header in enumerate(col_headers):
                        scan_entry.create_dataset(col_header, data=spice_data[:, idx])

    @staticmethod
    def _read_spice(file_name):
        """Reads an ascii generated by spice, and returns a header structure and a data table

        Args:
            file_name (str): a string containing the filename

        Returns:
            spice_data (numpy.array): an array containing all columns/rows
            headers (dict): a dictionary containing information from the commented lines.
            col_headers (list):
            unused (dict): not yet used in hdf5
        """
        with open(file_name, encoding="utf-8") as f:
            all_content = f.readlines()
            header_list = [line.strip() for line in all_content if "#" in line]
            col_name_index = header_list.index("# col_headers =") + 1
            col_names = header_list[col_name_index].strip("#").split()
            header_list.pop(col_name_index)
            header_list.pop(col_name_index - 1)

        spice_data = np.genfromtxt(file_name, comments="#")

        col_headers = col_names
        headers = {}
        unused = []

        for line in header_list:
            line = line.strip("# ")
            if "=" in line:  # empty field
                if line[-1] == "=":
                    unused.append(line[:-1])
                else:
                    parts = line.split("=")
                    headers[parts[0].strip()] = parts[1].strip()
            elif "scan completed" in line:  # last line
                parts = line.split(" ")
                headers["end_time"] = parts[0] + " " + parts[1]
            else:  # empty field
                unused.append(line)

        return spice_data, col_headers, headers, unused


if __name__ == "__main__":
    spice_folder = "./tests/test_data_folder/exp758/"
    h5_file_name = "./tests/test_data_folder/tavi_exp758_2.h5"
    TAVI_Data.convert_spice_to_hdf5(spice_folder, h5_file_name)

    data = TAVI_Data()
    data.load_data_from_disk(h5_file_name)
