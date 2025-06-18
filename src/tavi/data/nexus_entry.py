import os
from datetime import datetime
from typing import Optional

import h5py
import numpy as np

from tavi.data.nexus_builder import spice_data_to_nxdict


def _find_val_path(val, grp, prefix=""):
    """Find value in hdf5 groups"""
    for obj_name, obj in grp.items():
        if obj_name in ("SPICElogs", "data"):
            continue
        else:
            path = f"{prefix}/{obj_name}"
            if val == obj_name:
                if val == "detector":
                    return path + "/data"
                elif val == "monitor":
                    return path + "/monitor"
                else:
                    return path
            # test for group (go down)
            elif isinstance(obj, h5py.Group):
                gpath = _find_val_path(val, obj, path)
                if gpath:
                    return gpath


class NexusEntry(dict):
    """Read and write NeXus data files.

    Methods:
        from_nexus
        to_nexus
        get
    """

    @staticmethod
    def _getitem_recursively(obj: dict, key: str, ATTRS: bool):
        "find key in obj recursively, return None if nonexsiting"
        value = None
        if key.split("/")[0] in obj:
            for key_item in key.split("/"):
                # you could be looking in the wrong folder
                # e.g. looking for detector/data in data/detector
                # return None when this happens
                try:
                    obj = obj[key_item]
                except KeyError:
                    return None
            if ATTRS:  # get attributes
                try:
                    value = obj["attrs"]
                except KeyError:
                    print(f"Attribute of {key} does not exist.")
            else:  # get dataset
                try:
                    value = obj["dataset"]
                except KeyError:
                    print(f"Dataset of {key} does not exist.")

        for k, v in obj.items():
            if value is not None:  # value found
                break
            if k == "attrs" or k == "dataset":  # gone too far
                continue

            if isinstance(v, dict):
                value = NexusEntry._getitem_recursively(v, key, ATTRS)

        return value

    @staticmethod
    def _write_recursively(items, nexus_entry):
        """write items to nexus entry recursively
        Note:
        string encoded with utf-8
        """

        def format_dataset(value):
            """format if type is given in attributes"""
            dv = value["dataset"]

            if not (attr := value.get("attrs")):
                return dv

            attr_type = attr.get("type")
            match attr_type:
                case "NX_CHAR":
                    if isinstance(dv, list | np.ndarray):
                        dv = np.array([v.encode("utf-8") for v in dv])
                    elif isinstance(dv, str):
                        dv = dv.encode("utf-8")
                    else:
                        raise ValueError(f"Unrecogonized type for dataset={dv}")
                case "NX_FLOAT":
                    dv = np.array(dv).astype("float")
                case "NX_INT":
                    dv = np.array(dv).astype("int")
                case "NX_DATE_TIME":
                    pass
                case _:
                    if isinstance(dv, str):
                        dv = dv.encode("utf-8")
            # print(dv)
            return dv

        for key, value in items.items():
            if key == "attrs":
                for attr_key, attr_value in value.items():
                    if isinstance(attr_value, str):
                        attr_value = attr_value.encode("utf-8")
                    nexus_entry.attrs[attr_key] = attr_value
            elif isinstance(value, dict):
                if "dataset" in value.keys():
                    dv = format_dataset(value)
                    if key in nexus_entry:  # dataset exists
                        del nexus_entry[key]
                        # ds = nexus_entry[key]
                        # ds[...] = dv
                    # else:  # create dataset
                    ds = nexus_entry.create_dataset(name=key, data=dv, maxshape=None)
                    NexusEntry._write_recursively(value, ds)
                else:
                    grp = nexus_entry.require_group(name=key + "/")
                    NexusEntry._write_recursively(value, grp)

    @staticmethod
    def _read_recursively(nexus_entry, items=None):
        """read item from nexus_entry recursively"""
        if items is None:
            items = {}
        for key, value in nexus_entry.items():
            if isinstance(value, h5py.Group):
                attr_dict = {}
                for k, v in value.attrs.items():
                    attr_dict.update({k: v})
                items.update({key: {"attrs": attr_dict}})
                NexusEntry._read_recursively(value, items[key])
            elif isinstance(value, h5py.Dataset):
                attr_dict = {}
                for k, v in value.attrs.items():
                    attr_dict.update({k: v})
                try:  # unpacking as string
                    value = str(value.asstr()[...])
                except TypeError:  # arrays instead
                    value = value[...]
                    if not np.shape(value):  # single element
                        value = value.tolist()

                if not attr_dict:  # empty attributes
                    items.update({key: {"dataset": value}})
                else:
                    items.update({key: {"attrs": attr_dict, "dataset": value}})

        return items

    @staticmethod
    def _dict_to_nexus_entry(nexus_dict):
        """convert nested dict to instances of the NexusEntry class"""
        nexus_entries = {}
        for scan_num, scan_content in nexus_dict.items():
            content_list = []
            for key, val in scan_content.items():
                content_list.append((key, val))
            nexus_entries.update({scan_num: NexusEntry(content_list)})
        return nexus_entries

    # TODO read in instrument and sample configuration json files
    @classmethod
    def from_spice(
        cls,
        path_to_spice_folder: str,
        scan_num: Optional[int] = None,
        path_to_instrument_json: Optional[str] = None,
        path_to_sample_json: Optional[str] = None,
    ):
        """return a NexusEntry instance from loading a SPICE file

        Args:
            path_to_spice_folder (str): path to a SPICE folder
            scan_num (int): read all scans in folder if not None
            path_to_instrument_json: Optional[str] = None,
            path_to_sample_json: Optional[str] = None,
        """

        path_to_spice_folder = os.path.abspath(path_to_spice_folder)

        # Check that it exists
        if not os.path.exists(path_to_spice_folder):
            raise FileNotFoundError(f"Path does not exist: {path_to_spice_folder}")

        # Check that it's a directory
        if not os.path.isdir(path_to_spice_folder):
            raise NotADirectoryError(f"Path is not a directory: {path_to_spice_folder}")

        # Optional: check that "Datafiles" subdirectory exists
        datafiles_dir = os.path.join(path_to_spice_folder, "Datafiles")
        if not os.path.isdir(datafiles_dir):
            raise FileNotFoundError(f"Missing expected subdirectory: {datafiles_dir}")

        nexus_dict = spice_data_to_nxdict(
            path_to_spice_folder,
            scan_num,
            path_to_instrument_json,
            path_to_sample_json,
        )

        return cls._dict_to_nexus_entry(nexus_dict)

    @classmethod
    def from_nexus(
        cls,
        path_to_nexus: str,
        scan_num: Optional[int] = None,
        prefix: Optional[str] = None,
    ):
        """return a NexusEntry instance from loading a NeXus file

        Args:
            path_to_nexus (str): path to a NeXus file with the extension .h5
            scan_num (int): read all scans in file if not None
            prefix (str)
        """
        with h5py.File(path_to_nexus, "r") as nexus_file:
            if prefix is not None:
                nexus_file = nexus_file[prefix]
            if scan_num is None:  # read all scans in file
                nexus_dict = cls._read_recursively(nexus_file)
            else:  # read one scan only
                key = f"scan{scan_num:04}"
                nexus_dict = {key: cls._read_recursively(nexus_file[key])}

        return cls._dict_to_nexus_entry(nexus_dict)

    def to_nexus(self, path_to_nexus: str, name="scan") -> None:
        """write a NexueEntry instance to a NeXus file

        Args:
            path_to_nexus (str): path to a NeXus file with the extension .h5
        """
        with h5py.File(path_to_nexus, "a") as nexus_file:
            h5py.get_config().track_order = True
            scan_grp = nexus_file.require_group(name + "/")
            NexusEntry._write_recursively(self, scan_grp)
            if "data" in scan_grp:  # create soft link for data
                scan_grp_data = scan_grp["data"]
                def_y = scan_grp_data.attrs["signal"]
                def_x = scan_grp_data.attrs["axes"]
                path_y = _find_val_path(def_y, scan_grp)
                path_x = _find_val_path(def_x, scan_grp)

                if path_y is not None:
                    if scan_grp.get(path_y[1:]) is not None:  # remove the first "/"
                        if isinstance(scan_grp_data.get(def_y), h5py.Dataset):
                            del scan_grp_data[def_y]
                        path_y = "/" + name + path_y
                        scan_grp_data[def_y] = h5py.SoftLink(path_y)
                        scan_grp_data[def_y].attrs["target"] = path_y
                if path_x is not None:
                    if scan_grp.get(path_x[1:]) is not None:  # remove the first "/"
                        if isinstance(scan_grp_data.get(def_x), h5py.Dataset):
                            del scan_grp_data[def_x]
                        path_x = "/" + name + path_x
                        scan_grp_data[def_x] = h5py.SoftLink(path_x)
                        scan_grp_data[def_x].attrs["target"] = path_x

            # Create the ATTRIBUTES
            nexus_file.attrs["file_name"] = os.path.abspath(path_to_nexus)
            nexus_file.attrs["file_time"] = datetime.now().astimezone().isoformat()
            nexus_file.attrs["h5py_version"] = h5py.version.version
            nexus_file.attrs["HDF5_Version"] = h5py.version.hdf5_version

    def get(self, key, ATTRS=False, default=None):
        """
        Return dataset specified by key regardless of the hierarchy.
        Return attributes instead if ATTRS is True.
        Look only in the NeXus contents, ignore the DAS logs.

        Args:
            key (str): keyword or path. e.g. "s1" or "detector/data"
            ATTRS (bool): return attributes if true, dataset if false

        Note:
            Unique keys like 's1' or 'm2' can be found straight forwardly.
            To find monitor or detecor data use monitor/data or detector/data
        """
        # remove daslogs
        for log in ("SPICElogs", "DASlogs"):
            if log in self.keys():
                self.pop(log)

        # special treatment for annoying cases!!!
        if (key == "detector") and (not ATTRS):
            key = "detector/data"
        if (key == "monitor") and (not ATTRS):
            key = "monitor/monitor"
        if (key == "sample") and (not ATTRS):
            key = "sample/sample"

        value = NexusEntry._getitem_recursively(self, key, ATTRS)

        return value if value is not None else default

    def get_data_from_daslogs(self, key, default=None):
        """Look up data specified by key from das logs"""
        for log in ("SPICElogs", "DASlogs"):
            if log in self.keys():
                daslogs = self[log]
                continue
        value = daslogs.get(key)
        return value["dataset"] if value is not None else default

    def get_metadata_from_daslogs(self, key, default=None):
        """Look up metadata specified by key from das logs"""
        for log in ("SPICElogs", "DASlogs"):
            if log in self.keys():
                daslogs = self[log]
                continue
        value = daslogs["attrs"].get(key, None)
        return value if value is not None else default

    def get_dataset_names(self) -> list:
        """return the list of dataset names"""

        def _get_name(val_dict):
            for key, val in val_dict.items():
                if "dataset" not in val:
                    continue
                if type(val["dataset"]) is not np.ndarray:
                    continue
                if np.size(val["dataset"]) == num_pts:
                    dataset_names.append(key)

        dataset_names = []
        #  num_pts = len(pts) if (pts := self.get("Pt.")) is not None else 0
        num_pts = len(self.get("Pt."))

        instru = self["instrument"]
        for component, val_dict in instru.items():
            if component == "attrs":
                continue
            if component == "detector":
                if "data" in val_dict:
                    dataset_names.append("detector")
                    continue
                _get_name(val_dict)
            _get_name(val_dict)

        monitor = self["monitor"]
        if self.get("monitor/data") is not None:
            monitor.pop("data")
        _get_name(monitor)

        sample = self["sample"]
        _get_name(sample)

        return dataset_names
