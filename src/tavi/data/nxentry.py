from typing import Optional

import h5py
import numpy as np

from tavi.data.spice_reader import spice_data_to_nxdict


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
                obj = obj[key_item]
            if ATTRS:
                try:
                    value = obj["attrs"]
                except KeyError:
                    print(f"Attribute of {key} does not exist.")
            else:
                try:
                    value = obj["dataset"]
                except KeyError:
                    print(f"Dataset of {key} does not exist.")

        for k, v in obj.items():
            if value is not None:
                break
            if k == "attrs" or k == "dataset":
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
                    dv = dv.encode("utf-8")
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
                        ds = nexus_entry[key]
                        ds[...] = dv
                    else:  # create dataset
                        ds = nexus_entry.create_dataset(name=key, data=dv, maxshape=None)
                    NexusEntry._write_recursively(value, ds)
                else:
                    grp = nexus_entry.require_group(key + "/")
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

    # TODO read in instrument and sample configuratio json files
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
        # TODO validate path

        nexus_dict = spice_data_to_nxdict(
            path_to_spice_folder,
            scan_num,
            path_to_instrument_json,
            path_to_sample_json,
        )

        return NexusEntry._dict_to_nexus_entry(nexus_dict)

    @classmethod
    def from_nexus(cls, path_to_nexus: str, scan_num: Optional[int] = None):
        """return a NexusEntry instance from loading a NeXus file

        Args:
            path_to_nexus (str): path to a NeXus file with the extension .h5
            scan_num (int): read all scans in file if not None
        """
        with h5py.File(path_to_nexus, "r") as nexus_file:
            if scan_num is None:  # read all scans in file
                nexus_dict = NexusEntry._read_recursively(nexus_file)
            else:  # read one scan only
                key = f"scan{scan_num:04}"
                nexus_dict = {key: NexusEntry._read_recursively(nexus_file[key])}

        return NexusEntry._dict_to_nexus_entry(nexus_dict)

    def to_nexus(self, path_to_nexus: str, name="scan") -> None:
        """write a NexueEntry instance to a NeXus file

        Args:
            path_to_nexus (str): path to a NeXus file with the extention .h5
        """
        with h5py.File(path_to_nexus, "a") as nexus_file:
            scan_grp = nexus_file.require_group(name + "/")
            NexusEntry._write_recursively(self, scan_grp)

    def get(self, key, ATTRS=False, default=None):
        """
        Return dataset spicified by key regardless of the hierarchy.
        Return attributes instead if ATTRS is True.
        Look only in the NeXus contents, ignore the DAS logs.

        Args:
            key (str): keyword or path. e.g. "s1" or "detector/data"
            ATTRS (bool): return attributes if ture, dataset if false

        Note:
            Unique keys like 's1' or 'm2' can be found straight forwardly.
            To find monitor or detecor data use monitor/data or detector/data
        """
        # remove daslogs
        for log in ("SPICElogs", "DASlogs"):
            if log in self.keys():
                self.pop(log)
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
