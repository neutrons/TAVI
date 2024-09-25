from datetime import datetime
from typing import Optional


def nxsource(spicelogs, instrument_config_params):
    source = {
        "attrs": {"NX_class": "NXsource", "EX_required": "true"},
        "name": {
            "attrs": {"type": "NX_CHAR", "EX_required": "true"},
            "dataset": "HFIR",
        },
        "probe": {
            "attrs": {"type": "NX_CHAR", "EX_required": "true"},
            "dataset": "neutron",
        },
    }
    # Effective distance from sample Distance as seen by radiation from sample.
    # This number should be negative to signify that it is upstream of the sample.
    # nxsource.attrs["distance"] = -0.0
    return source


def nxmono(spicelogs, instrument_config_params):
    metadata = spicelogs["attrs"]

    try:
        mono = {
            "attrs": {"NX_class": "NXcrystal", "EX_required": "true"},
            "ei": {
                "dataset": spicelogs["ei"],
                "attrs": {"type": "NX_FLOAT", "EX_required": "true", "units": "meV"},
            },
            "type": {"dataset": metadata["monochromator"], "attrs": {"type": "NX_CHAR"}},
            "sense": {"dataset": metadata["sense"][0], "attrs": {"type": "NX_CHAR"}},
            "m1": {
                "dataset": spicelogs["m1"],
                "attrs": {"type": "NX_FLOAT", "units": "degrees"},
            },
            "m2": {
                "dataset": spicelogs["m2"],
                "attrs": {"type": "NX_FLOAT", "units": "degrees"},
            },
            "mfocus": {"dataset": spicelogs["mfocus"], "attrs": {"type": "NX_FLOAT"}},
            "marc": {"dataset": spicelogs["marc"], "attrs": {"type": "NX_FLOAT"}},
            "mtrans": {"dataset": spicelogs["mtrans"], "attrs": {"type": "NX_FLOAT"}},
            "focal_length": {"dataset": spicelogs["focal_length"], "attrs": {"type": "NX_FLOAT"}},
        }
    except KeyError:
        pass
    return mono


def nxcoll(spicelogs, instrument_config_params):
    return {}


def nxana(spicelogs, instrument_config_params):
    return {}


def nxdet(spicelogs, instrument_config_params):
    return {}


def nxmonitor(spicelogs, instrument_config_params):
    return {}


def nxsample(spicelogs, sample_config_params):
    return {}


def nxslit(spicelogs, instrument_config_params):
    return {}


def nxflipper(spicelogs, instrument_config_params):
    return {}


def _spicelogs_to_nexus(spicelogs, instrument_config_params, sample_config_params):
    metadata = spicelogs["attrs"]
    # TODO timezone
    start_date_time = "{} {}".format(metadata["date"], metadata["time"])
    start_time = datetime.strptime(start_date_time, "%m/%d/%Y %I:%M:%S %p").isoformat()
    # if "end_time" in das_logs.attrs:  # last scan never finished
    end_date_time = metadata["end_time"]
    end_time = datetime.strptime(end_date_time, "%m/%d/%Y %I:%M:%S %p").isoformat()

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
            "source": nxsource(spicelogs, instrument_config_params),
            "collimator": nxcoll(spicelogs, instrument_config_params),
            "monochromator": nxmono(spicelogs, instrument_config_params),
            "analyser": nxana(spicelogs, instrument_config_params),
            "detector": nxdet(spicelogs, instrument_config_params),
            "slit": nxslit(spicelogs, instrument_config_params),
            "flipper": nxflipper(spicelogs, instrument_config_params),
        },
        "monitor": nxmonitor(spicelogs, instrument_config_params),
        "sample": nxsample(spicelogs, sample_config_params),
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
