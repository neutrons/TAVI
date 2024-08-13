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
        "scan": int(nexus_entry.name[-4:]),  # last 4 digits are scan number
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
        # "Pt.": nexus_entry["sample/Pt."][...],
        "Pt.": np.arange(start=1, stop=num + 1, step=1),
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

    # analyser
    ana_str = (
        ("a1", "a2", "afocus", "ef") + tuple([f"qm{i+1}" for i in range(8)]) + tuple([f"xm{i+1}" for i in range(8)])
    )
    nexus_ana_str = nexus_entry["instrument/analyser"].keys()
    for key in ana_str:
        if key in nexus_ana_str:
            data.update({key: nexus_entry["instrument/analyser/" + key][...]})

    # monochromator
    mono_str = ("m1", "m2", "ei", "focal_length", "mfocus", "marc", "mtrans")
    nexus_mono_str = nexus_entry["instrument/monochromator"].keys()
    for key in mono_str:
        if key in nexus_mono_str:
            data.update({key: nexus_entry["instrument/monochromator/" + key][...]})

    # slits
    slits_str1 = ("bat", "bab", "bal", "bar", "bbt", "bbb", "bbl", "bbr")
    slits_str2 = ("slita_lf", "slita_rt", "slita_tp", "slitb_bt", "slitb_lf", "slitb_rt", "slitb_tp")
    slits_str3 = ("slit_pre_bt", "slit_pre_lf", "slit_pre_rt", "slit_pre_tp")
    slit_str = slits_str1 + slits_str2 + slits_str3

    nexus_slit_str = nexus_entry["instrument/slit"].keys()
    for key in slit_str:
        if key in nexus_slit_str:
            data.update({key: nexus_entry["instrument/slit/" + key][...]})

    # temprature
    temperatue_str = (
        (
            "temp",
            "temp_a",
            "temp_2",
            "coldtip",
            "tsample",
            "sample",
        )
        + (
            "vti",
            "dr_tsample",
            "dr_temp",
        )
        + ("lt", "ht", "sorb_temp", "sorb", "sample_ht")
    )
    for t in nexus_entry["sample"].keys():
        if t in temperatue_str:
            data.update({t: nexus_entry["sample/" + t][...]})

    # field
    field_str = ("persistent_field",)
    for f in nexus_entry["sample"].keys():
        if f in field_str:
            data.update({f: nexus_entry["sample/" + f][...]})

    return scan_info, sample_ub_info, instrument_info, data


def nexus_to_SPICE(nexus_entry):
    """Reads a nexus entry, convert to a SPICE scan file

    Args:
        nexus entry

    """
    pass
