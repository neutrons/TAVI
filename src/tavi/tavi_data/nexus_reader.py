import numpy as np


def nexus_to_dict(nexus_entry):
    """Reads a nexus entry, convert to dictionaries of data and meta_data

    Args:
        nexus entry

    Returns:
        meta_data (dict)
        data (dict)
    """

    def dataset_to_string(ds):
        return str(ds.asstr()[...])

    scan_info = {
        "scan": int(nexus_entry.name[5:]),
        "time": dataset_to_string(nexus_entry["start_time"]),
        "scan_title": dataset_to_string(nexus_entry["title"]),
        # "preset_type": "normal",
        "preset_channel": dataset_to_string(nexus_entry["monitor/mode"]),
        "preset_value": float(nexus_entry["monitor/preset"][...]),
        "def_y": nexus_entry["data"].attrs["signal"],
        "def_x": nexus_entry["data"].attrs["axes"],
    }

    sample_ub_info = {
        "samplename": dataset_to_string(nexus_entry["sample/name"]),
        "lattice_constants": nexus_entry["sample/unit_cell"][...],
        "ub_matrix": np.reshape(nexus_entry["sample/orientation_matrix"][...], (3, 3)),
        # "mode": # UB mode
        "palne_normal": nexus_entry["sample/plane_normal"][...],
        "ubconf": dataset_to_string(nexus_entry["sample/ub_conf"]),
    }

    instrument = dataset_to_string(nexus_entry["instrument/name"])

    instrument_info = {
        "instrument": instrument,
        # "monochromator":
        # "analyzer":
        "sense": dataset_to_string(nexus_entry["instrument/monochromator/sense"])
        + dataset_to_string(nexus_entry["sample/sense"])
        + dataset_to_string(nexus_entry["instrument/analyser/sense"]),
        "collimation": nexus_entry["instrument/collimator/divergence_x"][...],
    }

    num = np.size(nexus_entry["data/" + scan_info["def_x"]])

    data = {
        "Pt.": np.arange(start=1, stop=num + 1, step=1),
        # analyser
        "a1": nexus_entry["instrument/analyser/a1"][...],
        "a2": nexus_entry["instrument/analyser/a2"][...],
        "afocus": nexus_entry["instrument/analyser/afocus"][...],
        "ef": nexus_entry["instrument/analyser/ef"][...],
        # monochromator
        "m1": nexus_entry["instrument/monochromator/m1"][...],
        "m2": nexus_entry["instrument/monochromator/m2"][...],
        "ei": nexus_entry["instrument/monochromator/ei"][...],
        "focal_length": nexus_entry["instrument/monochromator/focal_length"][...],
        "mfocus": nexus_entry["instrument/monochromator/mfocus"][...],
        "marc": nexus_entry["instrument/monochromator/marc"][...],
        "mtrans": nexus_entry["instrument/monochromator/mtrans"][...],
        # detector
        "detector": nexus_entry["instrument/detector/data"][...],
        # monitor
        "time": nexus_entry["monitor/time"][...],
        "monitor": nexus_entry["monitor/monitor"][...],
        "mcu": nexus_entry["monitor/mcu"][...],
        # sample
        "s1": nexus_entry["sample/s1"][...],
        "s2": nexus_entry["sample/s2"][...],
        "sgl": nexus_entry["sample/sgl"][...],
        "sgu": nexus_entry["sample/sgu"][...],
        "stl": nexus_entry["sample/stl"][...],
        "stu": nexus_entry["sample/stu"][...],
        "q": nexus_entry["sample/q"][...],
        "qh": nexus_entry["sample/qh"][...],
        "qk": nexus_entry["sample/qk"][...],
        "ql": nexus_entry["sample/ql"][...],
        "en": nexus_entry["sample/en"][...],
    }

    # instrument specific motors

    match instrument:
        case "CG4C":
            # analyzer
            for i in range(8):
                data.update({f"qm{i+1}": nexus_entry[f"instrument/analyser/qm{i+1}"][...]})
                data.update({f"xm{i+1}": nexus_entry[f"instrument/analyser/xm{i+1}"][...]})
            # slits
            data.update(
                {
                    "bat": nexus_entry["instrument/slit/bat"][...],
                    "bab": nexus_entry["instrument/slit/bab"][...],
                    "bal": nexus_entry["instrument/slit/bal"][...],
                    "bar": nexus_entry["instrument/slit/bar"][...],
                    "bbt": nexus_entry["instrument/slit/bbt"][...],
                    "bbb": nexus_entry["instrument/slit/bbb"][...],
                    "bbl": nexus_entry["instrument/slit/bbl"][...],
                    "bbr": nexus_entry["instrument/slit/bbr"][...],
                }
            )
        case "HB1":
            pass
        case _:
            pass

    # temprature
    temperatue_str = (
        "coldtip",
        "tsample",
        "temp_a",
        "temp_2",
        "vti",
        "sample",
        "temp",
        "dr_tsample",
        "dr_temp",
    )
    for t in nexus_entry["sample"].keys():
        if t in temperatue_str:
            data.update({t: nexus_entry["sample/" + t][...]})

    return scan_info, sample_ub_info, instrument_info, data


def nexus_to_SPICE(nexus_entry):
    """Reads a nexus entry, convert to a SPICE scan file

    Args:
        nexus entry

    """
    pass
