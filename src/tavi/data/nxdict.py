from typing import Optional

import numpy as np


class NexusDict(object):
    """Read in dictionaries from DAS logs, instrument configuration json or sample json,
    format into NeXus style: nxentry_dict = {"attrs":{"attr1":attr1, ...}, "dataset":dataset}

    Attributes:
        daslogs_dict (dict)
        instrument_dict (dict)
        sample_dict (dict)
        nxentry_dict (dict)


    Methods:
        set_attrs()
        get_attrs()
        get_dataset()
        set_dataset()

    """

    def __init__(
        self,
        daslogs_dict: Optional[dict] = None,
        instrument_dict: Optional[dict] = None,
        sample_dict: Optional[dict] = None,
    ) -> None:
        self.daslogs_dict = daslogs_dict
        self.instrument_dict = instrument_dict
        self.sample_dict = sample_dict
        self.nxentry_dict: dict = {"attrs": {}}

    def set_attrs(self, **kwargs):
        for k, v in kwargs.items():
            self.nxentry_dict["attrs"].update({k: v})

    def get_attrs(self, key: str, default=None):
        val = self.nxentry_dict["attrs"].get(key)
        return val if val is not None else default

    def get_dataset(self, key: str, default=None):
        val = self.nxentry_dict.get(key)
        return val["dataset"] if val is not None else default

    def set_dataset(self, key: str, dataset: str | dict | np.ndarray, **kwargs):
        """Set dataset with proper format if dataset is a string or an array. Take directly if dictionary."""

        attr_dict = {}
        for k, v in kwargs.items():
            attr_dict.update({k: v})
        if "type" in kwargs.keys():
            match kwargs["type"]:
                case "NX_CHAR":
                    dataset = str(dataset)
                case "NX_INT":
                    dataset = dataset
                case "NX_FLOAT":
                    dataset = dataset

        self.nxentry_dict.update({key: {"dataset": dataset}})
        self.nxentry_dict[key].update({"attrs": attr_dict})

    def set_dataset_from(
        self,
        key: str,
        source: Literal["DAS_DATA", "DAS_METADATA", "INSTRU", "SAMPLE"] = "DAS_DATA",
        daslogs_key: Optional[str] = None,
        **kwargs,
    ):
        match source:
            case "DAS_DATA":
                if self.daslogs_dict is None:
                    raise ValueError("Cannot find DAS logs.")
            case "DAS_METADATA":
                if self.daslogs_dict is None:
                    raise ValueError("Cannot find DAS logs.")
            case "INSTRU":
                if self.instrument_dict is None:
                    raise ValueError("Cannot find instrument configuration dict.")
            case "SMAPLE":
                if self.sample_config_params is None:
                    raise ValueError("Cannot find sample configuration dict.")
            case _:
                raise ValueError(f"Unrecogonized source {source}.")

        match source:
            case "DAS_DATA":
                try:
                    val = self.daslogs_dict[key] if daslogs_key is None else self.daslogs_dict[daslogs_key]
                except KeyError:
                    print(f"Variable {key} cannot be found in DAS logs.")
                    return self.nxentry_dict

                self.nxentry_dict.update({key: val})
                if not kwargs:
                    return self.nxentry_dict

                attr_dict = {}
                for k, v in kwargs.items():
                    attr_dict.update({k: v})
                self.nxentry_dict[key].update({"attrs": attr_dict})
                return self.nxentry_dict
            case "DAS_METADATA":
                pass
            case "INSTRU":
                pass
            case "SMAPLE":
                pass


def spicelogs_to_nested_dict(
    spicelogs: dict,
    instrument_dict: Optional[dict],
    sample_dict: Optional[dict],
) -> dict:

    nxsource = NexusDict()
    nxsource.set_attrs(NX_class="NXsource", EX_required="true")
    nxsource.set_dataset(key="name", dataset="HFIR", type="NX_CHAR", EX_required="true")
    nxsource.set_dataset(key="probe", dataset="neutron", type="NX_CHAR", EX_required="true")

    # Effective distance from sample Distance as seen by radiation from sample.
    # This number should be negative to signify that it is upstream of the sample.
    # nxsource.attrs["distance"] = -0.0

    nxmono = NexusDict(daslogs_dict=spicelogs)
    nxmono.set_attrs(NX_class="NXcrystal", EX_required="true")
    # nxmono.set_dataset_from(source="DAS_DATA", key="ei", type="NX_FLOAT", EX_required="true", unit="meV")

    return {}
