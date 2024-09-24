import h5py


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

        for key, value in items.items():
            if key == "attrs":
                for attr_key, attr_value in value.items():
                    if isinstance(attr_value, str):
                        attr_value = attr_value.encode("utf-8")
                    nexus_entry.attrs[attr_key] = attr_value
            else:
                if isinstance(value, dict):
                    if "dataset" in value.keys():
                        dv = value["dataset"]
                        if isinstance(dv, str):
                            dv = dv.encode("utf-8")
                        ds = nexus_entry.create_dataset(
                            name=key,
                            data=dv,
                            maxshape=None,
                        )
                        NexusEntry._write_recursively(value, ds)
                    else:
                        grp = nexus_entry.create_group(key + "/")
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

    @classmethod
    def from_spice(cls, path_to_spice_folder: str) -> dict:
        """return a NexusEntry instance from loading a SPICE file

        Args:
            path_to_spice_folder (str): path to a SPICE folder
        """

        nexus_dict = {}
        return cls([(key, val) for key, val in nexus_dict.items()])

    @classmethod
    def from_nexus(cls, path_to_nexus: str) -> dict:
        """return a NexusEntry instance from loading a NeXus file

        Args:
            path_to_nexus (str): path to a NeXus file with the extension .h5
        """
        with h5py.File(path_to_nexus, "r") as nexus_file:
            nexus_dict = NexusEntry._read_recursively(nexus_file)
        return cls([(key, val) for key, val in nexus_dict.items()])

    def to_nexus(self, path_to_nexus: str) -> None:
        """write a NexueEntry instance to a NeXus file

        Args:
            path_to_nexus (str): path to a NeXus file with the extention .h5
        """
        with h5py.File(path_to_nexus, "w") as nexus_file:
            NexusEntry._write_recursively(self, nexus_file)

    # TODO need to catch errors when multiple IPTS are loaded in a tavi file
    def get(self, key, ATTRS=False, default=None):
        """
        Return dataset spicified by key regardless of the hierarchy.
        Return attributes instead if ATTRS is True.

        Args:
            key (str): keyword or path. e.g. "s1" or "detector/data"
            ATTRS (bool): return attributes if ture, dataset if false

        Note:
            Unique keys like 's1' or 'm2' can be found straight forwardly.
            To find monitor or detecor data use monitor/data or detector/data
        """
        value = NexusEntry._getitem_recursively(self, key, ATTRS)

        return value if value is not None else default
