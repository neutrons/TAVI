import numpy as np
import h5py
from pathlib import Path
from datetime import datetime
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
    def convert_spice_to_nexus(path_to_spice_folder, path_to_hdf5):
        """Load data from spice folder.

        Args:
            path_to_spice_folder (str): spice folder, ends with '/'
            path_to_nexus (str): path to hdf5 data file, ends with '.h5'
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
            f.attrs["ipts"] = ipts
            f.attrs["instrument"] = instrument_str
            f.attrs["exp"] = headers["experiment_number"]

            f.attrs["name"] = headers["experiment"]
            f.attrs["users"] = headers["users"]
            f.attrs["local_contact"] = headers["local_contact"]

            # convert SPICE scans into nexus entries
            for scan in scans:  # ignoring unused keys
                spice_data, col_headers, headers, unused = TAVI_Data._read_spice(scan)

                # pre-processing
                scan_num = ((scan.parts[-1].split("_"))[-1]).split(".")[0]  # e.g. "scan0001"
                if "scan_title" in unused:
                    headers["scan_title"] = ""
                # TODO timezone
                start_date_time = "{} {}".format(headers["date"], headers["time"])
                start_time = datetime.strptime(start_date_time, "%m/%d/%Y %I:%M:%S %p")

                # /entry
                nxentry = f.create_group(scan_num)
                nxentry.attrs["NX_class"] = "NXentry"

                nxentry.create_dataset(name="title", data=headers["scan_title"], maxshape=None)
                nxentry["title"].attrs["type"] = "NX_CHAR"

                nxentry.create_dataset(name="start_time", data=start_time.isoformat(), maxshape=None)
                nxentry["start_time"].attrs["type"] = "NX_DATE_TIME"

                if "end_time" in headers:  # last scan never finished
                    end_date_time = headers["end_time"]
                    end_time = datetime.strptime(end_date_time, "%m/%d/%Y %I:%M:%S %p")
                    nxentry.create_dataset(name="end_time", data=end_time.isoformat(), maxshape=None)
                    nxentry["end_time"].attrs["type"] = "NX_DATE_TIME"

                nxentry.create_dataset(name="definition", data="NXtas", maxshape=None)
                nxentry["definition"].attrs["type"] = "NX_CHAR"

                # /entry/SPICElogs
                nxentry.create_group("SPICElogs")
                nxentry["SPICElogs"].attrs["NX_class"] = "NXcollection"
                # metadata to attibutes
                for k, v in headers.items():
                    if k not in exp_info:  # ignore common keys in single scans
                        if "," in v and k != "scan_title":  # vectors
                            nxentry["SPICElogs"].attrs[k] = np.array([float(v0) for v0 in v.split(",")])
                        elif v.replace(".", "").isnumeric():  # numebrs only
                            if v.isdigit():  # int
                                nxentry["SPICElogs"].attrs[k] = int(v)
                            else:  # float
                                nxentry["SPICElogs"].attrs[k] = float(v)
                        # separate COM/FWHM and its errorbar
                        elif k == "Center of Mass":
                            com, e_com = v.split("+/-")
                            nxentry["SPICElogs"].attrs["COM"] = float(com)
                            nxentry["SPICElogs"].attrs["COM_err"] = float(e_com)
                        elif k == "Full Width Half-Maximum":
                            fwhm, e_fwhm = v.split("+/-")
                            nxentry["SPICElogs"].attrs["FWHM"] = float(fwhm)
                            nxentry["SPICElogs"].attrs["FWHM_err"] = float(e_fwhm)
                        else:  # other crap, keep as is
                            if k not in exp_info:
                                nxentry["SPICElogs"].attrs[k] = v
                # motor position table to datasets
                if spice_data.ndim == 1:  # empty data or 1 point only
                    if len(spice_data):  # 1 point only
                        for idx, col_header in enumerate(col_headers):
                            nxentry["SPICElogs"].create_dataset(col_header, data=spice_data[idx])
                    else:  # empty
                        pass
                else:  # nomarl data
                    for idx, col_header in enumerate(col_headers):
                        nxentry["SPICElogs"].create_dataset(col_header, data=spice_data[:, idx])

                # /entry/user

                # /entry/instrument
                nxentry.create_group("instrument")
                nxentry["instrument"].attrs["NX_class"] = "NXinstrument"
                nxentry["instrument"].attrs["name"] = instrument_str
                nxentry["instrument"].attrs["sense"] = headers["sense"]

                # /entry/instrument/source
                nxentry["instrument"].create_group("source")
                nxentry["instrument/source"].attrs["NX_class"] = "NXsource"
                nxentry["instrument/source"].attrs["name"] = "HFIR"
                # Effective distance from sample Distance as seen by radiation from sample.
                # This number should be negative to signify that it is upstream of the sample.
                nxentry["instrument/source"].attrs["distance"] = -0.0
                nxentry["instrument/source"].attrs["type"] = "Reactor Neutron Source"
                nxentry["instrument/source"].attrs["probe"] = "neutron"
                # particle beam size in x and y
                nxentry["instrument/source"].attrs["sigma_x"] = 0.0
                nxentry["instrument/source"].attrs["sigma_y"] = 0.0

                # /entry/instrument/collimators
                nxentry["instrument"].create_group("collimators")
                nxentry["instrument/collimators"].attrs["types"] = "Soller"
                nxentry["instrument/collimators"].attrs["soller_angles"] = headers["collimation"]
                nxentry["instrument/collimators"].attrs["divergence_x"] = [0, 0, 0, 0]
                nxentry["instrument/collimators"].attrs["divergence_y"] = [0, 0, 0, 0]

                # /entry/instrument/slits
                nxslits = nxentry["instrument"].create_group("slits")
                match instrument_str:
                    case "HB1":
                        nxslits.attrs["slita_bt"] = spice_data[:, col_headers.index("slita_bt")]
                        nxslits.attrs["slita_lf"] = spice_data[:, col_headers.index("slita_lf")]
                        nxslits.attrs["slita_rt"] = spice_data[:, col_headers.index("slita_tt")]
                        nxslits.attrs["slita_tp"] = spice_data[:, col_headers.index("slita_tp")]
                        nxslits.attrs["slitb_bt"] = spice_data[:, col_headers.index("slitb_bt")]
                        nxslits.attrs["slitb_lf"] = spice_data[:, col_headers.index("slitb_lf")]
                        nxslits.attrs["slitb_rt"] = spice_data[:, col_headers.index("slitb_tt")]
                        nxslits.attrs["slitb_tp"] = spice_data[:, col_headers.index("slitb_tp")]
                    case "CG4C":
                        nxslits.attrs["bbb"] = spice_data[:, col_headers.index("bbb")]
                        nxslits.attrs["bbl"] = spice_data[:, col_headers.index("bbl")]
                        nxslits.attrs["bbr"] = spice_data[:, col_headers.index("bbr")]
                        nxslits.attrs["bbt"] = spice_data[:, col_headers.index("bbt")]
                        nxslits.attrs["bab"] = spice_data[:, col_headers.index("bab")]
                        nxslits.attrs["bal"] = spice_data[:, col_headers.index("bal")]
                        nxslits.attrs["bar"] = spice_data[:, col_headers.index("bar")]
                        nxslits.attrs["bat"] = spice_data[:, col_headers.index("bat")]

                # /entry/instrument/masks
                nxmask = nxentry["instrument"].create_group("mask")

                # /entry/instrument/monochromator
                nxmono = nxentry["instrument"].create_group("monochromator")
                nxmono.attrs["ei"] = spice_data[:, col_headers.index("ei")]
                # nxmono.attrs["usage"] = "Bragg"
                nxmono.attrs["type"] = headers["monochromator"]
                nxmono.attrs["focal_length"] = spice_data[:, col_headers.index("focal_length")]
                nxmono.attrs["m1"] = spice_data[:, col_headers.index("m1")]
                nxmono.attrs["m2"] = spice_data[:, col_headers.index("m2")]
                nxmono.attrs["marc"] = spice_data[:, col_headers.index("marc")]
                nxmono.attrs["mtrans"] = spice_data[:, col_headers.index("mtrans")]
                nxmono.attrs["mfocus"] = spice_data[:, col_headers.index("mfocus")]

                # /entry/instrument/analyzer
                nxano = nxentry["instrument"].create_group("analyzer")
                nxano.attrs["type"] = headers["analyzer"]
                nxano.attrs["ef"] = spice_data[:, col_headers.index("ef")]
                nxano.attrs["a1"] = spice_data[:, col_headers.index("a1")]
                nxano.attrs["a2"] = spice_data[:, col_headers.index("a2")]
                nxano.attrs["afocus"] = spice_data[:, col_headers.index("afocus")]

                if instrument_str == "CG4C":  # focused analyzer for CTAX
                    for i in range(8):  # qm1--qm8, xm1 -- xm8
                        nxano.attrs[f"qm{i+1}"] = spice_data[:, col_headers.index(f"qm{i+1}")]
                        nxano.attrs[f"xm{i+1}"] = spice_data[:, col_headers.index(f"xm{i+1}")]

                # /entry/instrument/detector
                nxdetector = nxentry["instrument"].create_group("detector")
                nxdetector.attrs["NX_class"] = "NXdetector"
                cts = spice_data[:, col_headers.index("detector")]
                nxdetector.attrs["data"] = np.array([int(ct) for ct in cts])

                # TODO HB1 polarized experiment

                # /entry/monitor
                nxmonitor = nxentry.create_group("monitor")
                nxmonitor.attrs["NX_class"] = "NXmonitor"
                if headers["preset_type"] == "normal":
                    nxmonitor.attrs["mode"] = headers["preset_channel"]
                    nxmonitor.attrs["preset"] = float(headers["preset_value"])
                    # all three recorded regardless of preset channel
                    nxmonitor.attrs["time"] = spice_data[:, col_headers.index("time")]
                    nxmonitor.attrs["monitor"] = spice_data[:, col_headers.index("monitor")]
                    nxmonitor.attrs["mcu"] = spice_data[:, col_headers.index("mcu")]

                elif headers["preset_type"] == "countfile":
                    # TODO HB1 polarized experiment
                    pass
                else:  # unknown
                    pass

                # /entry/sample
                nxsample = nxentry.create_group("sample")
                nxsample.attrs["name"] = headers["samplename"]
                nxsample.attrs["type"] = headers["sampletype"]
                nxsample.attrs["mosiac"] = headers["samplemosaic"]

                lcs = headers["latticeconstants"]
                nxsample.attrs["unit_cell"] = np.array([float(lc) for lc in lcs.split(",")])

                ub = headers["ubmatrix"]
                nxsample.attrs["ub_matrix"] = np.array([float(element) for element in ub.split(",")])
                pn = headers["plane_normal"]
                nxsample.attrs["plane_normal"] = np.array([float(element) for element in pn.split(",")])
                nxsample.attrs["ub_conf_file"] = headers["ubconf"]
                nxsample.attrs["ub_mode"] = int(headers["mode"])
                # TODO more UB info in the UBConf files

                nxsample.attrs["s1"] = spice_data[:, col_headers.index("s1")]
                nxsample.attrs["s2"] = spice_data[:, col_headers.index("s2")]
                nxsample.attrs["sgl"] = spice_data[:, col_headers.index("sgl")]
                nxsample.attrs["sgu"] = spice_data[:, col_headers.index("sgu")]
                nxsample.attrs["stl"] = spice_data[:, col_headers.index("stl")]
                nxsample.attrs["stu"] = spice_data[:, col_headers.index("stu")]

                nxsample.attrs["qh"] = spice_data[:, col_headers.index("h")]
                nxsample.attrs["qk"] = spice_data[:, col_headers.index("k")]
                nxsample.attrs["ql"] = spice_data[:, col_headers.index("l")]
                nxsample.attrs["en"] = spice_data[:, col_headers.index("e")]

                # temperture
                temperatue_str = [
                    "coldtip",
                    "tsample",
                    "temp_a",
                    "vti",
                    "sample",
                    "temp",
                    "dr_tsample",
                    "dr_temp",
                ]
                for t in temperatue_str:
                    if t in col_headers:
                        nxsample.attrs[t] = spice_data[:, col_headers.index(t)]

                # TODO field
                # TODO pressure

                # /entry/data
                nxdata = nxentry.create_group("data")
                nxdata.attrs["NX_class"] = "NXdata"

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
                    unused.append(line[:-2])  # remove  " ="
                else:
                    parts = line.split("=")
                    headers[parts[0].strip()] = parts[1].strip()
            elif "completed" in line or "stopped" in line:  # last line
                parts = line.split(" ")
                headers["end_time"] = parts[3] + " " + parts[0] + " " + parts[1]
            else:  # empty field
                unused.append(line)

        return spice_data, col_headers, headers, unused


if __name__ == "__main__":
    spice_folder = "./tests/test_data_folder/exp416/"
    # h5_file_name = "./tests/test_data_folder/tavi_exp758.h5"
    nexus_file_name = "./tests/test_data_folder/nexus_exp416.h5"
    TAVI_Data.convert_spice_to_nexus(spice_folder, nexus_file_name)

    # data = TAVI_Data()
    # data.load_data_from_disk(h5_file_name)
