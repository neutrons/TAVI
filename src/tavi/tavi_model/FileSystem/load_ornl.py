import logging
import os

import numpy as np

import tavi.tavi_model.FileSystem.spice_reader as spice_reader
from tavi.tavi_model.FileSystem.tavi_class import RawData, RawMetaData, Scan, TaviProject, UbConf

logger = logging.getLogger("TAVI")

def score_ornl(dir):
    
    # if it's below a thresh hold of 5 then we choose load_ornl
    score = 0
    try: 
        filename = os.listdir(dir)[0]
    except:
        logger.error("No file in directory, check directory contains triple-axis data files!")
        return -1
    instrument_names = ["CG4C", "HB1", "HB3", "HB1A"]
    for instrument_name in instrument_names:
        if instrument_name in filename:
            score += 10
    with open(os.path.join(dir, filename), encoding="utf-8") as f:
        all_content = f.readlines()
    headers = [line.strip() for line in all_content if "#" in line]
    print(headers[0])

def load_ornl(dir):
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
    filepath = os.path.join(current_directory, "test_data", "exp424", "Datafiles")
    score_ornl(filepath)
