from typing import Optional

# def nxentry_dict(self):
#     metadata = self.daslogs_dict["attrs"]
#


def nxcoll(spicelogs, instrument_config_params):
    pass


# def nxana(spicelogs, instrument_config_params):
#     ana = {"attrs": {"NX_class": "NXcrystal", "EX_required": "true"}}
#     ana = _add_dataset_entry(ana, key="ef", daslogs=spicelogs, type="NX_FLOAT", EX_required="true", unit="meV")
#     ana.update({"type": {"dataset": spicelogs["attrs"]["analyzer"], "attrs": {"type": "NX_CHAR"}}})
#     ana.update({"sense": {"dataset": spicelogs["attrs"]["sense"][2], "attrs": {"type": "NX_CHAR"}}})
#     ana = _add_dataset_entry(ana, key="a1", daslogs=spicelogs, type="NX_FLOAT", unit="degrees")
#     ana = _add_dataset_entry(ana, key="a2", daslogs=spicelogs, type="NX_FLOAT", unit="degrees")
#     ana = _add_dataset_entry(ana, key="afocus", daslogs=spicelogs, type="NX_FLOAT")
#     for i in range(8):
#         ana = _add_dataset_entry(ana, key=f"qm{i+1}", daslogs=spicelogs, type="NX_FLOAT")
#         ana = _add_dataset_entry(ana, key=f"xm{i+1}", daslogs=spicelogs, type="NX_FLOAT")
#     return ana


# def nxdet(spicelogs, instrument_config_params):
#     det = {"attrs": {"NX_class": "NXdetector", "EX_required": "true"}}
#     det = _add_dataset_entry(
#         det, key="data", daslogs=spicelogs, daslogs_key="detector", type="NX_INT", EX_required="true", unit="counts"
#     )
#     # polar_angle
#     return det


# def nxmonitor(spicelogs, instrument_config_params):
#     monitor = {"attrs": {"NX_class": "NXmonitor", "EX_required": "true"}}
#     preset_type = spicelogs["attrs"]["preset_type"]
#     match preset_type:
#         case "countfile":  # polarization data
#             print("Countfile preset type is not supported.")

#         case "normal":
#             preset_channel = spicelogs["attrs"]["preset_channel"]
#             monitor.update({"mode": {"dataset": preset_channel, "attrs": {"type": "NX_CHAR", "EX_required": "true"}}})
#             monitor.update(
#                 {
#                     "preset": {
#                         "dataset": spicelogs["attrs"]["preset_value"],
#                         "attrs": {"type": "NX_FLOAT", "EX_required": "true"},
#                     }
#                 }
#             )
#             monitor = _add_dataset_entry(monitor, key="time", daslogs=spicelogs, type="NX_FLOAT", units="seconds")
#             monitor = _add_dataset_entry(monitor, key="monitor", daslogs=spicelogs, type="NX_INT", units="counts")
#             monitor = _add_dataset_entry(monitor, key="mcu", daslogs=spicelogs, type="NX_FLOAT")
#             monitor_data_dict = monitor[preset_channel]
#             monitor.update({"data": monitor_data_dict})

#         case _:
#             print(f"Unrecogonized preset type {preset_type}.")

#     return monitor


def nxsample(spicelogs, sample_config_params):
    return {}


def nxslit(spicelogs, instrument_config_params):
    return {}


def nxflipper(spicelogs, instrument_config_params):
    return {}


def _spicelogs_to_nexus(spicelogs, instrument_config_params, sample_config_params):

    # def nxmono(spicelogs, instrument_config_params):

    # mono = _add_dataset_entry(mono, key="ei", daslogs=spicelogs, )
    # mono.update({"type": {"dataset": spicelogs["attrs"]["monochromator"], "attrs": {"type": "NX_CHAR"}}})
    # mono.update({"sense": {"dataset": spicelogs["attrs"]["sense"][0], "attrs": {"type": "NX_CHAR"}}})
    # mono = _add_dataset_entry(mono, key="m1", daslogs=spicelogs, type="NX_FLOAT", unit="degrees")
    # mono = _add_dataset_entry(mono, key="m2", daslogs=spicelogs, type="NX_FLOAT", unit="degrees")
    # mono = _add_dataset_entry(mono, key="mfocus", daslogs=spicelogs, type="NX_FLOAT")
    # mono = _add_dataset_entry(mono, key="marc", daslogs=spicelogs, type="NX_FLOAT")
    # mono = _add_dataset_entry(mono, key="mtrans", daslogs=spicelogs, type="NX_FLOAT")
    # mono = _add_dataset_entry(mono, key="focal_length", daslogs=spicelogs, type="NX_FLOAT")

    # return mono

    scan_dict = {
        "attrs": {
            "NX_class": "NXentry",
            "EX_required": "true",
        },
        "SPICElogs": spicelogs,
        "definition": {
            "attrs": {"EX_required": "true", "type": "NX_CHAR"},
            "dataset": "NXtas",
        },
        "title": {
            "attrs": {"EX_required": "true", "type": "NX_CHAR"},
            "dataset": metadata["scan_title"],
        },
        "start_time": {
            "attrs": {"EX_required": "true", "type": "NX_DATE_TIME"},
            "dataset": start_time,
        },
        "end_time": {
            "attrs": {"type": "NX_DATE_TIME"},
            "dataset": end_time,
        },
        "instrument": {
            "attrs": {"EX_required": "true", "NX_class": "NXinstrument"},
            "name": {"attrs": {"type": "NX_CHAR"}, "dataset": metadata["instrument"]},
            "source": nxsource.entry_dict,
            # "collimator": nxcoll(spicelogs, instrument_config_params),
            # "monochromator": nxmono(spicelogs, instrument_config_params),
            # "analyser": nxana(spicelogs, instrument_config_params),
            # "detector": nxdet(spicelogs, instrument_config_params),
            # "slit": nxslit(spicelogs, instrument_config_params),
            # "flipper": nxflipper(spicelogs, instrument_config_params),
        },
        # "monitor": nxmonitor(spicelogs, instrument_config_params),
        # "sample": nxsample(spicelogs, sample_config_params),
        "data": {},
    }
    return scan_dict


def daslogs_to_nexus_dict(
    daslogs: dict,
    instrument_config_params: Optional[str],
    sample_config_params: Optional[str],
) -> dict:
    """Format DASlogs dict into NeXus dict"""

    match (das_key := tuple(daslogs.keys())[0]):
        case "SPICElogs":
            scan_dict = _spicelogs_to_nexus(
                daslogs["SPICElogs"],
                instrument_config_params,
                sample_config_params,
            )
        case _:
            raise KeyError(f"Unrecogonized DASlogs key {das_key}.")

    return scan_dict
