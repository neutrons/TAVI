import h5py


class NexusEntry(dict):

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

    def get(self, key, ATTRS=False, default=None):
        """
        Return dataset spicified by key regardless of the hierarchy.
        Return attributes instead if ATTRS is True.

        Note:
            Unique keys like 's1' or 'm2' can be found straight forwardly.
            To find monitor or detecor data use monitor/data or detector/data
        """
        value = NexusEntry._getitem_recursively(self, key, ATTRS)
        if value is not None:
            return value
        else:
            return default

    @classmethod
    def from_nexus(cls, path_to_nexus: str) -> dict:
        with h5py.File(path_to_nexus, "r") as nexus_file:
            nexus_dict = NexusEntry._read_recursively(nexus_file)
        return cls([(key, val) for key, val in nexus_dict.items()])

    def to_nexus(self, path_to_nexus: str) -> None:
        with h5py.File(path_to_nexus, "w") as nexus_file:
            NexusEntry._write_recursively(self, nexus_file)
