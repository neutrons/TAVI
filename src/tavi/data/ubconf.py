from typing import Optional

import numpy as np

from tavi.data.nxentry import NexusEntry
from tavi.utilities import Peak


class UBConf(NexusEntry):
    "Logs for UB matrix determination"

    ub_peaks: Optional[tuple[Peak, ...]] = None
    u_mat: Optional[np.ndarray] = None
    b_mat: Optional[np.ndarray] = None
    ub_mat: Optional[np.ndarray] = None
    plane_normal: Optional[np.ndarray] = None
    in_plane_ref: Optional[np.ndarray] = None


def convert_spice_ub_to_nexus(
    path_to_spice_folder: str,
    path_to_hdf5_folder: str,
    verbose: bool = False,
) -> None:
    """Convert all UBConf files into one NeXus entry named UBConf.h5"""

    if verbose:
        disp_str = (
            f"Converting SPICE UBConf files at {path_to_spice_folder} to an NeXus entry at {path_to_hdf5_folder}"
        )
        print(disp_str)

    p = Path(path_to_spice_folder)
    ub_files = sorted((p / "UBConf").glob("*.ini"))
    tmp_ub_files = sorted((p / "UBConf/tmp").glob("*.ini"))
    ub_files_all = ub_files + tmp_ub_files
    ub_conf_dicts = {ub_file.parts[-1].split(".")[0]: _read_spice_ub(ub_file) for ub_file in ub_files_all}

    with h5py.File(path_to_hdf5_folder + "UBConf.h5", "w") as root:

        for ub_name, ub_data in ub_conf_dicts.items():
            ub_entry = root.create_group(ub_name)
            ub_entry.attrs["NX_class"] = "NXcollection"
            ub_entry.attrs["X_required"] = "false"
            for k, v in ub_data.items():
                ub_entry.create_dataset(name=k, data=v, maxshape=None)
