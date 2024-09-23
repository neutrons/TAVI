import h5py


class NexusEntry(dict):

    @staticmethod
    def _getitem_recursively(obj, key, ATTRS):
        "find key in obj recursively, return None if nonexsiting"
        value = None
        if key in obj:
            if ATTRS:
                try:
                    value = obj[key]["attrs"]
                except KeyError:
                    print("Attribute does not exist.")
            else:
                try:
                    value = obj[key]["dataset"]
                except KeyError:
                    print("Dataset does not exist.")

        for k, v in obj.items():
            if value is not None:
                break
            if k == "attrs" or k == "dataset":
                continue

            if isinstance(v, dict):
                value = NexusEntry._getitem_recursively(v, key, ATTRS)

        return value

    def get(self, key, ATTRS=False, default=None):
        """
        Return dataset spicified by key regardless of the hierarchy.
        Return attributes instead if ATTRS is True.

        Note:
            Only works if the key is unique
        """
        value = NexusEntry._getitem_recursively(self, key, ATTRS)
        if value is not None:
            return value
        else:
            return default

    @staticmethod
    def _write_recursively(items, nexus_entry):
        """write items to nexus entry recursively"""
        for key, value in items.items():
            if key == "attrs":
                for attr_key, attr_value in value.items():
                    nexus_entry.attrs[attr_key] = attr_value
            else:
                if isinstance(value, dict):
                    if "dataset" in value.keys():
                        ds = nexus_entry.create_dataset(name=key, data=value["dataset"], maxshape=None)
                        NexusEntry._write_recursively(value, ds)
                    else:
                        grp = nexus_entry.create_group(key + "/")
                        NexusEntry._write_recursively(value, grp)

    def to_nexus(self, path_to_nexus):
        with h5py.File(path_to_nexus, "w") as nexus_file:
            NexusEntry._write_recursively(self, nexus_file)
