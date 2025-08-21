import logging
import os
from dataclasses import field, make_dataclass
from typing import Iterable, Optional

import numpy as np

import tavi.tavi_model.FileSystem.spice_reader as spice_reader
from tavi.tavi_model.FileSystem.tavi_class_factory import Scan, TaviProject

logger = logging.getLogger("TAVI")


class LoadORNL:
    def __init__(
        self,
        data_folder: Optional[os.PathLike | str],
        data_file: Optional[os.PathLike | str | Iterable[os.PathLike | str]] = None,
    ) -> None:
        self.data_folder = data_folder
        self.data_files = data_file
        super().__init__()

    def score(self) -> float:
        # if it's below a thresh hold of 5 then we choose load_ornl
        # score = 0
        # try:
        #     filename = os.listdir(dir)[0]
        # except:
        #     logger.error("No file in directory, check directory contains triple-axis data files!")
        #     return -1
        # instrument_names = ["CG4C", "HB1", "HB3", "HB1A"]

        # # if instrument name contains HFIR instruments
        # for instrument_name in instrument_names:
        #     if instrument_name in filename:
        #         score += 10

        # with open(os.path.join(dir, filename), encoding="utf-8") as f:
        #     all_content = f.readlines()
        # headers = [line.strip() for line in all_content if "#" in line]

        # temporarily return a 100 score. Need more time in determining the scoring system
        return 100

    def load(self) -> TaviProject:
        """
        Load SPICE data files from either a directory or a list of data_files if specified into a TaviProject.

        This function reads each file in the specified directory using
        `spice_reader.read_spice_datafile()`, converts the numeric arrays and
        metadata into `RawData` and `RawMetaData` dataclass instances, and stores
        them as `Scan` objects in a `TaviProject`.

        Args:
            directory (str | Path): Path to the folder containing SPICE data files.

        Returns:
            TaviProject: A project object containing all loaded scans.
        """

        # Initialize TaviProject class
        if self.data_folder and not self.data_files:
            tavi_project = self._load_folder()
        # load files into TaviProject scans

        # tmp_exist = False
        # ub_dir = os.path.join(dir, "UBConf")
        # for ub_filename in os.listdir(ub_dir):
        #     ub_conf = UbConf()
        #     if ub_filename == "tmp":
        #         tmp_exist = True
        #         continue
        #     ub_data = spice_reader.read_spice_ubconf(os.path.join(ub_dir, ub_filename))
        #     for key, value in ub_data.items():
        #         if hasattr(ub_conf, key):
        #             setattr(ub_conf, key, value)
        #         else:
        #             logger.warning("New UbConf found, consider updating UbConf entry in tavi_data_schema.json")
        #     tavi_project.ubconf[ub_filename] = ub_conf

        # if tmp_exist:
        #     for ub_filename in os.listdir(os.path.join(ub_dir, "tmp")):
        #         ub_conf = UbConf()
        #         ub_data = spice_reader.read_spice_ubconf(os.path.join(ub_dir, "tmp", ub_filename))
        #         for key, value in ub_data.items():
        #             if hasattr(ub_conf, key):
        #                 setattr(ub_conf, key, value)
        #             else:
        #                 logger.warning("New UbConf found, consider updating UbConf entry in tavi_data_schema.json")
        #         tavi_project.ubconf["tmp-" + ub_filename] = ub_conf

        return tavi_project

    def _load_folder(self):
        tavi_project = TaviProject()

        for filename in os.listdir(self.data_folder):
            numeric_data, col_names, meta_data, others, error_message = spice_reader.read_spice_datafile(
                os.path.join(self.data_folder, filename)
            )
            rawdata = make_dataclass("RawData", [], slots=True)
            rawmetadata = make_dataclass("RawMetaData", [], slots=True)

            for col_name in col_names:
                # guard against invalid format
                if col_name[0].isdigit():
                    col_name = "_" + col_name
                attr_name = col_name.replace("-", "_").replace(" ", "_").replace(".", "")
                try:
                    rawdata = make_dataclass(
                        "RawData", fields=[(attr_name, np.ndarray, field(default=None))], bases=(rawdata,)
                    )
                    setattr(rawdata, attr_name, numeric_data[:, col_names.index(col_name)])

                # if the numeric_data only has 1 entry, we'd still like to parse it as a 1D numpy array
                except ValueError:
                    # rawdata = make_dataclass("rawdata", fields = [attr_name, np.ndarray, field(default = np.array([numeric_data[col_names.index(col_name)]]))], bases = (rawdata,))
                    rawdata = make_dataclass(
                        "RawData", fields=[attr_name, np.ndarray, field(default=None)], bases=(rawdata,)
                    )
                    setattr(rawdata, attr_name, np.array([numeric_data[col_names.index(col_name)]]))

            for key, value in meta_data.items():
                # replace "-" or " " with "_" to consolidate with python attribute's format
                if key[0].isdigit():
                    key = "_" + key
                key = key.replace("-", "_").replace(" ", "_").replace(".", "")
                if hasattr(rawmetadata, key):
                    setattr(rawmetadata, key, value)
                rawmetadata = make_dataclass(
                    "RawMetaData", fields=[(attr_name, str, field(default=None))], bases=(rawmetadata,)
                )
                setattr(rawmetadata, key, value)

            scan = Scan(
                data=rawdata,
                metadata=rawmetadata,
                column_names=col_names,
                error_message=error_message,
                others=others,
                ubconf=None,
            )
            tavi_project.scans[filename] = scan
            return tavi_project


if __name__ == "__main__":
    current_directory = os.getcwd()
    filepath = os.path.join(current_directory, "test_data", "exp424", "Datafiles")
    ornl = LoadORNL(filepath)
    tavi_project = ornl.load()
    print(tavi_project.scans["CG4C_exp0424_scan0073.dat"].data.q)
