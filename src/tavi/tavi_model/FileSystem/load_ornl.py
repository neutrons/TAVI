import logging
import os
from dataclasses import field, make_dataclass
from typing import Any, Iterable, Optional

import numpy as np

import tavi.tavi_model.FileSystem.spice_reader as spice_reader
from tavi.tavi_model.FileSystem.tavi_class_factory import Scan, TaviProject

logger = logging.getLogger("TAVI")


class LoadORNL:
    def __init__(
        self,
        data_folder: Optional[os.PathLike | str],
        data_file: Optional[os.PathLike | str | Iterable[os.PathLike | str]] = None,
        ub_dir: Optional[os.PathLike | str] = None,
    ) -> None:
        self.data_folder = data_folder
        self.data_files = data_file
        self.ub_dir = ub_dir

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
        # TO DO
        # load a list of file(s)
        return tavi_project

    def _load_folder(self):
        tavi_project = TaviProject()

        for filename in os.listdir(self.data_folder):
            rawdata = make_dataclass("RawData", [], slots=True)
            rawmetadata = make_dataclass("RawMetaData", [], slots=True)
            ub_conf = make_dataclass("UbConf", [], slots=True)

            numeric_data, col_names, meta_data, others, error_message = spice_reader.read_spice_datafile(
                os.path.join(self.data_folder, filename)
            )

            # load data
            for col_name in col_names:
                # guard against invalid format
                if col_name[0].isdigit():
                    col_name = "_" + col_name
                attr_name = col_name.replace("-", "_").replace(" ", "_").replace(".", "")
                try:
                    # adding a checker to avoid re-creating class attributes, make the code faster
                    if hasattr(rawdata, attr_name):
                        setattr(rawdata, attr_name, numeric_data[:, col_names.index(col_name)])
                    else:
                        rawdata = make_dataclass(
                            "RawData", fields=[(attr_name, np.ndarray, field(default=None))], bases=(rawdata,)
                        )
                        setattr(rawdata, attr_name, numeric_data[:, col_names.index(col_name)])

                # if the numeric_data only has 1 entry, we'd still like to parse it as a 1D numpy array
                except IndexError:
                    rawdata = make_dataclass(
                        "RawData", fields=[(attr_name, np.ndarray, field(default=None))], bases=(rawdata,)
                    )
                    setattr(rawdata, attr_name, np.array([numeric_data[col_names.index(col_name)]]))

            for key, value in meta_data.items():
                # replace "-" or " " with "_" to consolidate with python attribute's format
                if key[0].isdigit():
                    key = "_" + key
                key = key.replace("-", "_").replace(" ", "_").replace(".", "")
                # adding a checker to avoid re-creating class attributes, make the code faster
                if hasattr(rawmetadata, key):
                    setattr(rawmetadata, key, value)
                else:
                    rawmetadata = make_dataclass(
                        "RawMetaData", fields=[(key, str, field(default=None))], bases=(rawmetadata,)
                    )
                    setattr(rawmetadata, key, value)

                match key:
                    case "ubconf":
                        ub_filename = value
                        if ub_filename:
                            # users don't need to set specific ub_directory, the loader will try to search for ub
                            if not self.ub_dir:
                                prev_dir = os.path.join(self.data_folder, os.pardir)
                                if "UBConf" in os.listdir(prev_dir):
                                    self.ub_dir = os.path.join(prev_dir, "UBConf")
                            # look into "UBConf" folder first
                            if ub_filename in os.listdir(self.ub_dir):
                                ub_data = spice_reader.read_spice_ubconf(os.path.join(self.ub_dir, ub_filename))
                                for key, value in ub_data.items():
                                    # adding a checker to avoid re-creating class attributes, make the code faster
                                    if hasattr(ub_conf, key):
                                        setattr(ub_conf, key, value)
                                    else:
                                        ub_conf = make_dataclass(
                                            "UbConf", fields=[(key, Any, field(default=None))], bases=(ub_conf,)
                                        )
                                        setattr(ub_conf, key, value)
                            # look into tmp folder just in case
                            elif ub_filename in os.listdir(os.path.join(self.ub_dir, "tmp")):
                                ub_data = spice_reader.read_spice_ubconf(os.path.join(self.ub_dir, ub_filename))
                                for key, value in ub_data.items():
                                    # adding a checker to avoid re-creating class attributes, make the code faster
                                    if hasattr(ub_conf, key):
                                        setattr(ub_conf, key, value)
                                    else:
                                        ub_conf = make_dataclass(
                                            "UbConf", fields=[(key, Any, field(default=None))], bases=(ub_conf,)
                                        )
                                        setattr(ub_conf, key, value)
                            else:
                                logger.warning("Can't find %s, please double check UBMatrix data", ub_filename)
            scan = Scan(
                data=rawdata,
                metadata=rawmetadata,
                column_names=col_names,
                error_message=error_message,
                others=others,
                ubconf=ub_conf,
            )
            tavi_project.scans[filename] = scan
        return tavi_project


if __name__ == "__main__":
    current_directory = os.getcwd()
    filepath = os.path.join(current_directory, "test_data", "exp424", "Datafiles")
    ornl = LoadORNL(filepath)
    tavi_project = ornl.load()
    print(tavi_project.scans["CG4C_exp0424_scan0041.dat"].ubconf)
