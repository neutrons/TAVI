import os
from pathlib import Path
import logging
import numpy as np
import tavi.tavi_model.FileSystem.spice_reader as spice_reader
from tavi.tavi_model.FileSystem.tavi_class import TaviProject, RawData, RawMetaData, UbConf,Scan

logger = logging.getLogger("TAVI")

def load_folder(dir):
    """
    Load SPICE data files from a directory into a TaviProject.

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
    tavi_project = TaviProject()
    data_dir = os.path.join(dir, "Datafiles")
    # load files into TaviProject scans
    for filename in os.listdir(data_dir):
        numeric_data, col_names, meta_data, others, error_message = spice_reader.read_spice_datafile(
            os.path.join(data_dir, filename)
        )
        rawdata = RawData()
        rawmetadata = RawMetaData()

        for col_name in col_names:
            attr_name = "Pt" if col_name == "Pt." else col_name
            if hasattr(rawdata, attr_name):
                try:
                    setattr(rawdata, attr_name, numeric_data[:, col_names.index(col_name)])
                # if the numeric_data only has 1 entry, we'd still like to parse it as a 1D numpy array
                except IndexError:
                    setattr(rawdata, attr_name, np.array([numeric_data[col_names.index(col_name)]]))
            else:
                logger.warning("New data entry detected, consider updating RawData entry in tavi_data_schema.json")

        for key, value in meta_data.items():
            # replace "-" or " " with "_" to consolidate with python attribute's format
            key = key.replace("-", "_").replace(" ", "_")
            if hasattr(rawmetadata, key):
                setattr(rawmetadata, key, value)
            else:
                logger.warning(
                    "New meta data entry detected, consider updating RawMetaData entry in tavi_data_schema.json"
                )

        scan = Scan(
            data=rawdata,
            metadata=rawmetadata,
            column_names=col_names,
            error_message=error_message,
            others=others,
        )
        tavi_project.scans[filename] = scan

    tmp_exist = False
    ub_dir = os.path.join(dir, "UBConf")
    for ub_filename in os.listdir(ub_dir):
        ub_conf = UbConf()
        if ub_filename == "tmp":
            tmp_exist = True
            continue
        ub_data = spice_reader.read_spice_ubconf(os.path.join(ub_dir, ub_filename))
        for key, value in ub_data.items():
            if hasattr(ub_conf, key):
                setattr(ub_conf, key, value)
            else:
                logger.warning("New UbConf found, consider updating UbConf entry in tavi_data_schema.json")
        tavi_project.ubconf[ub_filename] = ub_conf

    if tmp_exist:
        for ub_filename in os.listdir(os.path.join(ub_dir, "tmp")):
            ub_conf = UbConf()
            ub_data = spice_reader.read_spice_ubconf(os.path.join(ub_dir, "tmp", ub_filename))
            for key, value in ub_data.items():
                if hasattr(ub_conf, key):
                    setattr(ub_conf, key, value)
                else:
                    logger.warning("New UbConf found, consider updating UbConf entry in tavi_data_schema.json")
        tavi_project.ubconf["tmp-" + ub_filename] = ub_conf

    return tavi_project


if __name__ == "__main__":
    current_directory = os.getcwd()
    filepath = os.path.join(current_directory, "test_data", "exp424")
    tavi_project = load_folder(filepath)
    ubpath = os.path.join(current_directory, "test_data", "exp424", "UBConf", "UB02Jul2024_14108PM.ini")
    ub = spice_reader.read_spice_ubconf(ubpath)
    # print(ub)
    # print(tavi_project.scans["CG4C_exp0424_scan0073.dat"].metadata)
    print(tavi_project.ubconf["tmp-UB02Jul2024_21029PM.ini"])
