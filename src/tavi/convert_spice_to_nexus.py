import pathlib
from typing import Optional

from tavi.data.nexus_entry import NexusEntry


def convert_spice_to_nexus(
    path_to_spice_folder: str,
    path_to_nexus_folder: Optional[str] = None,
    path_to_instrument_json: Optional[str] = None,
    path_to_sample_json: Optional[str] = None,
):
    """convert datafiles in spice folder into nexus entries

    Note:
        if path_to_nexus_folder if not provided, the nexus data entries
        will be creatred in the same parent folder as the SPICE data,
        using the defalut folder name IPTSxxxxx_INSTRU_exp0000"""

    scan_dict = NexusEntry.from_spice(
        path_to_spice_folder,
        path_to_instrument_json=path_to_instrument_json,
        path_to_sample_json=path_to_sample_json,
    )
    scan0001 = next(iter(scan_dict.values()))
    nexus_folder_name = scan0001["attrs"]["dataset_name"]

    if path_to_nexus_folder is None:
        *parent_path, _ = path_to_spice_folder.split("/")
        path_to_nexus_folder = "/".join(parent_path + [nexus_folder_name + "/"])

    # make directory
    pathlib.Path(path_to_nexus_folder).mkdir(parents=True, exist_ok=True)

    for scan_name, scan in scan_dict.items():
        scan.to_nexus(path_to_nexus_folder + scan_name + ".h5", name=scan_name)


if __name__ == "__main__":
    path_to_instrument_json = "./src/tavi/instrument/instrument_params/cg4c.json"
    path_to_sample_json = "./test_data/test_samples/nitio3.json"
    path_to_spice_folder_exp424 = "./test_data/exp424"
    # path_to_spice_folder = "./test_data/exp815"  # empty runs in exp815
    # path_to_spice_folder = "./test_data/exp813"
    convert_spice_to_nexus(
        path_to_spice_folder_exp424,
        path_to_instrument_json=path_to_instrument_json,
        path_to_sample_json=path_to_sample_json,
    )
