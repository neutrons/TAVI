import numpy as np
import h5py
import os
from pathlib import Path
from datetime import datetime

# from tavi.instrument.instrument_params.takin_test import instrument_params


def read_spice(file_name):
    """Reads an ascii generated by spice, and returns a header structure and a data table

    Args:
        file_name (str): a string containing the filename

    Returns:
        spice_data (numpy.array): an array containing all columns/rows
        headers (dict): a dictionary containing information from the commented lines.
        col_headers (list):
        unused (dict): not yet used in hdf5
    """
    with open(file_name, encoding="utf-8") as f:
        all_content = f.readlines()
        header_list = [line.strip() for line in all_content if "#" in line]
        col_name_index = header_list.index("# col_headers =") + 1
        col_names = header_list[col_name_index].strip("#").split()
        header_list.pop(col_name_index)
        header_list.pop(col_name_index - 1)

    spice_data = np.genfromtxt(file_name, comments="#")

    col_headers = col_names
    headers = {}
    unused = []

    for line in header_list:
        line = line.strip("# ")
        if "=" in line:  # empty field
            if line[-1] == "=":
                unused.append(line[:-2])  # remove  " ="
            else:
                parts = line.split("=")
                key = parts[0].strip()
                val = "=".join(parts[1:])[1:]  # remove the fisrt space charactor
                headers[key] = val
        elif "completed" in line or "stopped" in line:  # last line
            parts = line.split(" ")
            headers["end_time"] = parts[3] + " " + parts[0] + " " + parts[1]
        else:  # empty field
            unused.append(line)

    return spice_data, col_headers, headers, unused


def read_spice_ub(ub_file):
    """Reads ub info from UBConf

    Args:
        ub_file (str): a string containing the filename

    Returns:
        ubconf (dict)
    """
    ubconf = {}
    with open(ub_file, encoding="utf-8") as f:
        all_content = f.readlines()

        for idx, line in enumerate(all_content):
            if line.strip() == "":
                continue  # skip if empty
            elif line.strip()[0] == "[":
                continue  # skiplines like "[xx]"

            ub_dict = {}
            key, val = line.strip().split("=")
            if key == "Mode":
                if all_content[idx - 1].strip() == "[UBMode]":
                    ub_dict["UBMode"] = int(val)
                elif all_content[idx - 1].strip() == "[AngleMode]":
                    ub_dict["AngleMode"] = int(val)
            else:
                if "," in line:  # vector
                    ub_dict[key] = np.array([float(v) for v in val.strip('"').split(",")])
                else:  # float number
                    ub_dict[key] = float(val)

            ubconf.update(ub_dict)

    return ubconf


def spicelogs_to_nexus(nxentry):
    """Format info from SPICElogs into Nexus format"""
    spice_logs = nxentry["SPICElogs"]

    # Create the GROUPS
    nxentry.attrs["NX_class"] = "NXentry"
    nxentry.attrs["EX_required"] = "true"

    nxentry.create_group("instrument")
    nxentry["instrument"].attrs["NX_class"] = "NXinstrument"
    nxentry["instrument"].attrs["EX_required"] = "true"

    nxentry["instrument/"].create_group("source")
    nxentry["instrument/source"].attrs["NX_class"] = "NXsource"
    nxentry["instrument/source"].attrs["EX_required"] = "true"

    nxmono = nxentry["instrument/"].create_group("monochromator")
    nxmono.attrs["NX_class"] = "NXcrystal"
    nxmono.attrs["EX_required"] = "true"

    nxana = nxentry["instrument/"].create_group("analyser")
    nxana.attrs["NX_class"] = "NXcrystal"
    nxana.attrs["EX_required"] = "true"

    nxentry["instrument/"].create_group("detector")
    nxentry["instrument/detector"].attrs["NX_class"] = "NXdetector"
    nxentry["instrument/detector"].attrs["EX_required"] = "true"

    nxentry.create_group("sample")
    nxentry["sample"].attrs["NX_class"] = "NXsample"
    nxentry["sample"].attrs["EX_required"] = "true"

    nxentry.create_group("monitor")
    nxentry["monitor"].attrs["NX_class"] = "NXmonitor"
    nxentry["monitor"].attrs["EX_required"] = "true"

    nxentry.create_group("data")
    nxentry["data"].attrs["NX_class"] = "NXdata"
    nxentry["data"].attrs["EX_required"] = "true"

    nxentry["instrument/"].create_group("collimator")
    nxentry["instrument/collimator"].attrs["NX_class"] = "NXcollimator"

    nxslit = nxentry["instrument/"].create_group("slit")
    nxslit.attrs["NX_class"] = "NXslit"

    nxflipper = nxentry["instrument/"].create_group("flipper")
    nxflipper.attrs["NX_class"] = "NXflipper"

    # nxentry["instrument"].create_group("filter")
    # nxentry["instrument"].attrs["NX_class"] = "NXfilter"

    # Create the FIELDS

    nxentry.create_dataset(name="title", data=spice_logs.attrs["scan_title"], maxshape=None)
    nxentry["title"].attrs["type"] = "NX_CHAR"
    nxentry["title"].attrs["EX_required"] = "true"

    # TODO timezone
    start_date_time = "{} {}".format(spice_logs.attrs["date"], spice_logs.attrs["time"])
    start_time = datetime.strptime(start_date_time, "%m/%d/%Y %I:%M:%S %p").isoformat()

    nxentry.create_dataset(name="start_time", data=start_time, maxshape=None)
    nxentry["start_time"].attrs["type"] = "NX_DATE_TIME"
    nxentry["start_time"].attrs["EX_required"] = "true"

    if "end_time" in spice_logs.attrs:  # last scan never finished
        end_date_time = spice_logs.attrs["end_time"]
        end_time = datetime.strptime(end_date_time, "%m/%d/%Y %I:%M:%S %p").isoformat()
        nxentry.create_dataset(name="end_time", data=end_time, maxshape=None)
        nxentry["end_time"].attrs["type"] = "NX_DATE_TIME"

    # Valid enumeration values for root['/entry']['definition'] are:
    # NXtas
    nxentry.create_dataset(name="definition", data="NXtas", maxshape=None)
    nxentry["definition"].attrs["type"] = "NX_CHAR"
    nxentry["definition"].attrs["EX_required"] = "true"

    #  --------------------------- instrument ---------------------------
    nxentry["instrument"].create_dataset(name="name", data=spice_logs.attrs["instrument"], maxshape=None)
    nxentry["instrument/name"].attrs["type"] = "NX_CHAR"

    #  --------------------------- source ---------------------------
    nxentry["instrument/source"].create_dataset(name="name", data="HFIR", maxshape=None)
    nxentry["instrument/source/name"].attrs["type"] = "NX_CHAR"
    nxentry["instrument/source/name"].attrs["EX_required"] = "true"

    # Valid enumeration values for root['/entry/instrument/source']['probe'] are:
    # 	 neutron
    # 	 x-ray

    nxentry["instrument/source"].create_dataset(name="probe", data="neutron", maxshape=None)
    nxentry["instrument/source/probe"].attrs["type"] = "NX_CHAR"
    nxentry["instrument/source/probe"].attrs["EX_required"] = "true"

    # Effective distance from sample Distance as seen by radiation from sample.
    # This number should be negative to signify that it is upstream of the sample.
    # nxentry["instrument/source"].attrs["distance"] = -0.0

    #  --------------------------- collimators ---------------------------
    nxentry["instrument/collimator"].create_dataset(name="type", data="Soller", maxshape=None)
    nxentry["instrument/collimator/type"].attrs["type"] = "NX_CHAR"

    div_x = [float(v) for v in list(spice_logs.attrs["collimation"].split("-"))]
    nxentry["instrument/collimator"].create_dataset(name="divergence_x", data=div_x, maxshape=None)
    nxentry["instrument/collimator/divergence_x"].attrs["type"] = "NX_ANGLE"
    nxentry["instrument/collimator/divergence_x"].attrs["units"] = "min"

    #  --------------------------- monochromator ---------------------------
    nxmono.create_dataset(name="ei", data=spice_logs["ei"], maxshape=None)
    nxmono["ei"].attrs["type"] = "NX_FLOAT"
    nxmono["ei"].attrs["EX_required"] = "true"
    # nxmono["ei"].attrs["axis"] = "1"
    nxmono["ei"].attrs["units"] = "meV"

    nxmono.create_dataset(name="type", data=spice_logs.attrs["monochromator"], maxshape=None)
    nxmono.attrs["type"] = "NX_CHAR"

    nxmono.create_dataset(name="m1", data=spice_logs["m1"], maxshape=None)
    nxmono["m1"].attrs["type"] = "NX_FLOAT"
    nxmono["m1"].attrs["units"] = "degrees"

    nxmono.create_dataset(name="m2", data=spice_logs["m2"], maxshape=None)
    nxmono["m2"].attrs["type"] = "NX_FLOAT"
    nxmono["m2"].attrs["units"] = "degrees"

    if "mfocus" in spice_logs.keys():
        nxmono.create_dataset(name="mfocus", data=spice_logs["mfocus"], maxshape=None)
        nxmono["mfocus"].attrs["type"] = "NX_FLOAT"

    if "marc" in spice_logs.keys():
        nxmono.create_dataset(name="marc", data=spice_logs["marc"], maxshape=None)
        nxmono["marc"].attrs["type"] = "NX_FLOAT"

    if "mtrans" in spice_logs.keys():
        nxmono.create_dataset(name="mtrans", data=spice_logs["mtrans"], maxshape=None)
        nxmono["mtrans"].attrs["type"] = "NX_FLOAT"

    if "focal_length" in spice_logs.keys():
        nxmono.create_dataset(name="focal_length", data=spice_logs["focal_length"], maxshape=None)
        nxmono["focal_length"].attrs["type"] = "NX_FLOAT"

    nxmono.create_dataset(name="sense", data=spice_logs.attrs["sense"][0], maxshape=None)
    nxmono.attrs["type"] = "NX_CHAR"

    # nxmono.create_dataset(name="rotation_angle", data=1.0, maxshape=None)
    # nxmono["rotation_angle"].attrs["type"] = "NX_FLOAT"
    # nxmono["rotation_angle"].attrs["EX_required"] = "true"
    # nxmono["rotation_angle"].attrs["units"] = "NX_ANGLE"

    # --------------------------- analyzer ---------------------------

    nxana.create_dataset(name="ef", data=spice_logs["ef"], maxshape=None)
    nxana["ef"].attrs["type"] = "NX_FLOAT"
    nxana["ef"].attrs["EX_required"] = "true"
    # nxana["ef"].attrs["axis"] = "1"
    nxana["ef"].attrs["units"] = "meV"

    nxana.create_dataset(name="type", data=spice_logs.attrs["analyzer"], maxshape=None)
    nxana.attrs["type"] = "NX_CHAR"

    nxana.create_dataset(name="a1", data=spice_logs["a1"], maxshape=None)
    nxana["a1"].attrs["type"] = "NX_FLOAT"
    nxana["a1"].attrs["units"] = "degrees"

    nxana.create_dataset(name="a2", data=spice_logs["a2"], maxshape=None)
    nxana["a2"].attrs["type"] = "NX_FLOAT"
    nxana["a2"].attrs["units"] = "degrees"

    if "afocus" in spice_logs.keys():
        nxana.create_dataset(name="afocus", data=spice_logs["afocus"], maxshape=None)
        nxana["afocus"].attrs["type"] = "NX_FLOAT"

    if spice_logs.attrs["instrument"] == "CG4C":
        for i in range(8):  # qm1--qm8, xm1 -- xm8
            nxana.create_dataset(name=f"qm{i+1}", data=spice_logs[f"qm{i+1}"], maxshape=None)
            nxana.create_dataset(name=f"xm{i+1}", data=spice_logs[f"xm{i+1}"], maxshape=None)

    nxana.create_dataset(name="sense", data=spice_logs.attrs["sense"][2], maxshape=None)
    nxana.attrs["type"] = "NX_CHAR"

    # nxana.create_dataset(name="rotation_angle", data=1.0, maxshape=None)
    # nxana["rotation_angle"].attrs["type"] = "NX_FLOAT"
    # nxana["rotation_angle"].attrs["EX_required"] = "true"
    # nxana["rotation_angle"].attrs["units"] = "NX_ANGLE"

    # nxana.create_dataset(name="polar_angle", data=1.0, maxshape=None)
    # nxana["polar_angle"].attrs["type"] = "NX_FLOAT"
    # nxana["polar_angle"].attrs["EX_required"] = "true"
    # nxana["polar_angle"].attrs["units"] = "NX_ANGLE"

    # --------------------------- detector ---------------------------

    nxentry["instrument/detector"].create_dataset(name="data", data=spice_logs["detector"], maxshape=None, dtype="int")
    nxentry["instrument/detector/data"].attrs["type"] = "NX_INT"
    nxentry["instrument/detector/data"].attrs["EX_required"] = "true"
    nxentry["instrument/detector/data"].attrs["units"] = "counts"

    # TODO HB1 polarized experiment

    # nxentry["instrument/detector"].create_dataset(name="polar_angle", data=1.0, maxshape=None)
    # nxentry["instrument/detector/polar_angle"].attrs["type"] = "NX_FLOAT"
    # nxentry["instrument/detector/polar_angle"].attrs["EX_required"] = "true"
    # nxentry["instrument/detector/polar_angle"].attrs["units"] = "NX_ANGLE"

    # ---------------------------- flipper ---------------------------------
    if "fguide" in spice_logs.keys():
        nxflipper.create_dataset(name="fguide", data=spice_logs["fguide"], maxshape=None)
        nxflipper["fguide"].attrs["type"] = "NX_FLOAT"
    if "hguide" in spice_logs.keys():
        nxflipper.create_dataset(name="hguide", data=spice_logs["hguide"], maxshape=None)
        nxflipper["hguide"].attrs["type"] = "NX_FLOAT"
    if "vguide" in spice_logs.keys():
        nxflipper.create_dataset(name="vguide", data=spice_logs["vguide"], maxshape=None)
        nxflipper["vguide"].attrs["type"] = "NX_FLOAT"

    # TODO Helmohtz coils guide fields: tbguide, aguide, bguide
    #
    # --------------------------- slits ---------------------------

    slits_str1 = ("bat", "bab", "bal", "bar", "bbt", "bbb", "bbl", "bbr")
    slits_str2 = ("slita_lf", "slita_rt", "slita_tp", "slitb_bt", "slitb_lf", "slitb_rt", "slitb_tp")
    slits_str3 = ("slit_pre_bt", "slit_pre_lf", "slit_pre_rt", "slit_pre_tp")

    slits_str = (slits_str1, slits_str2, slits_str3)

    for slit_str in slits_str:
        if slit_str[0] in spice_logs.keys():
            for st in slit_str:
                nxslit.create_dataset(name=st, data=spice_logs[st])
                nxslit[st].attrs["type"] = "NX_FLOAT"
                nxslit[st].attrs["units"] = "cm"

    # --------------------------- sample ---------------------------

    nxentry["sample"].create_dataset(name="name", data=spice_logs.attrs["samplename"], maxshape=None)
    nxentry["sample/name"].attrs["type"] = "NX_CHAR"
    nxentry["sample/name"].attrs["EX_required"] = "true"

    nxentry["sample"].create_dataset(name="qh", data=spice_logs["h"], maxshape=None)
    nxentry["sample/qh"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/qh"].attrs["EX_required"] = "true"
    # nxentry["sample/qh"].attrs["axis"] = "1"
    # nxentry["sample/qh"].attrs["units"] = "NX_DIMENSIONLESS"

    nxentry["sample"].create_dataset(name="qk", data=spice_logs["k"], maxshape=None)
    nxentry["sample/qk"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/qk"].attrs["EX_required"] = "true"
    # nxentry["sample/qk"].attrs["axis"] = "1"
    # nxentry["sample/qk"].attrs["units"] = "NX_DIMENSIONLESS"

    nxentry["sample"].create_dataset(name="ql", data=spice_logs["l"], maxshape=None)
    nxentry["sample/ql"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/ql"].attrs["EX_required"] = "true"
    # nxentry["sample/ql"].attrs["axis"] = "1"
    # nxentry["sample/ql"].attrs["units"] = "NX_DIMENSIONLESS"

    nxentry["sample"].create_dataset(name="en", data=spice_logs["e"], maxshape=None)
    nxentry["sample/en"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/en"].attrs["EX_required"] = "true"
    # nxentry["sample/en"].attrs["axis"] = "1"
    nxentry["sample/en"].attrs["units"] = "meV"

    nxentry["sample"].create_dataset(name="sgu", data=spice_logs["sgu"], maxshape=None)
    nxentry["sample/sgu"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/sgu"].attrs["EX_required"] = "true"
    nxentry["sample/sgu"].attrs["units"] = "degrees"

    nxentry["sample"].create_dataset(name="sgl", data=spice_logs["sgl"], maxshape=None)
    nxentry["sample/sgl"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/sgl"].attrs["EX_required"] = "true"
    nxentry["sample/sgl"].attrs["units"] = "degrees"

    nxentry["sample"].create_dataset(name="unit_cell", data=spice_logs.attrs["latticeconstants"], maxshape=None)
    nxentry["sample/unit_cell"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/unit_cell"].attrs["EX_required"] = "true"
    # nxentry["sample/unit_cell"].attrs["units"] = "NX_LENGTH"

    nxentry["sample"].create_dataset(name="orientation_matrix", data=spice_logs.attrs["ubmatrix"], maxshape=None)
    nxentry["sample/orientation_matrix"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/orientation_matrix"].attrs["EX_required"] = "true"
    nxentry["sample/orientation_matrix"].attrs["units"] = "NX_DIMENSIONLESS"

    nxentry["sample"].create_dataset(name="ub_conf", data=spice_logs.attrs["ubconf"].split(".")[0], maxshape=None)
    nxentry["sample/ub_conf"].attrs["type"] = "NX_CHAR"

    nxentry["sample"].create_dataset(name="plane_normal", data=spice_logs.attrs["plane_normal"], maxshape=None)
    nxentry["sample/plane_normal"].attrs["type"] = "NX_FLOAT"

    nxentry["sample"].create_dataset(name="q", data=spice_logs["q"], maxshape=None)
    nxentry["sample/q"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/q"].attrs["units"] = "Angstrom^-1"

    nxentry["sample"].create_dataset(name="stu", data=spice_logs["stu"], maxshape=None)
    nxentry["sample/stu"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/stu"].attrs["units"] = "degrees"

    nxentry["sample"].create_dataset(name="stl", data=spice_logs["stl"], maxshape=None)
    nxentry["sample/stl"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/stl"].attrs["units"] = "degrees"

    nxentry["sample"].create_dataset(name="s1", data=spice_logs["s1"], maxshape=None)
    nxentry["sample/s1"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/s1"].attrs["units"] = "degrees"

    nxentry["sample"].create_dataset(name="s2", data=spice_logs["s2"], maxshape=None)
    nxentry["sample/s2"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/s2"].attrs["units"] = "degrees"

    nxentry["sample"].create_dataset(name="type", data=spice_logs.attrs["sampletype"], maxshape=None)
    nxentry["sample/type"].attrs["type"] = "NX_CHAR"

    nxentry["sample"].create_dataset(name="sense", data=spice_logs.attrs["sense"][1], maxshape=None)
    nxentry["sample"].attrs["type"] = "NX_CHAR"

    # nxentry["sample"].create_dataset(name="rotation_angle", data=1.0, maxshape=None)
    # nxentry["sample/rotation_angle"].attrs["type"] = "NX_FLOAT"
    # nxentry["sample/rotation_angle"].attrs["EX_required"] = "true"
    # nxentry["sample/rotation_angle"].attrs["units"] = "NX_ANGLE"

    # nxentry["sample"].create_dataset(name="polar_angle", data=1.0, maxshape=None)
    # nxentry["sample/polar_angle"].attrs["type"] = "NX_FLOAT"
    # nxentry["sample/polar_angle"].attrs["EX_required"] = "true"
    # nxentry["sample/polar_angle"].attrs["units"] = "NX_ANGLE"

    # --------------------------- monitor ---------------------------
    # Valid enumeration values for root['/entry/monitor']['mode'] are:
    # 	 monitor
    # 	 time
    #    mcu

    if spice_logs.attrs["preset_type"] == "normal":

        preset_channel = spice_logs.attrs["preset_channel"]

        nxentry["monitor"].create_dataset(name="mode", data=preset_channel, maxshape=None)
        nxentry["monitor/mode"].attrs["type"] = "NX_CHAR"
        nxentry["monitor/mode"].attrs["EX_required"] = "true"

        nxentry["monitor"].create_dataset(name="preset", data=spice_logs.attrs["preset_value"], maxshape=None)
        nxentry["monitor/preset"].attrs["type"] = "NX_FLOAT"
        nxentry["monitor/preset"].attrs["EX_required"] = "true"

        nxentry["monitor"].create_dataset(name="time", data=spice_logs["time"], maxshape=None)
        nxentry["monitor/time"].attrs["type"] = "NX_FLOAT"
        nxentry["monitor/time"].attrs["units"] = "seconds"

        nxentry["monitor"].create_dataset(name="monitor", data=spice_logs["monitor"], maxshape=None)
        nxentry["monitor/monitor"].attrs["type"] = "NX_INT"
        nxentry["monitor/monitor"].attrs["units"] = "counts"

        nxentry["monitor"].create_dataset(name="mcu", data=spice_logs["mcu"], maxshape=None)
        nxentry["monitor/mcu"].attrs["type"] = "NX_FLOAT"

        nxentry["monitor"].create_dataset(name="data", data=spice_logs[preset_channel], maxshape=None)
        nxentry["monitor/data"].attrs["type"] = "NX_FLOAT"
        nxentry["monitor/data"].attrs["EX_required"] = "true"
        # nxentry["monitor/data"].attrs["units"] = "counts"

    # TODO polarized exp at HB1
    elif spice_logs.attrs["preset_type"] == "countfile":
        print("Polarization data, not yet supported.")

    else:
        print("Unrecogonized preset type. ")

    # --------------------------- data links ---------------------------

    def find_val(val, grp, prefix=""):

        for obj_name, obj in grp.items():
            if obj_name in ("SPICElogs", "data"):
                continue
            else:
                path = f"{prefix}/{obj_name}"
                if val == obj_name:
                    return path
                elif isinstance(obj, h5py.Group):  # test for group (go down)
                    gpath = find_val(val, obj, path)
                    if gpath:
                        return gpath

    nexus_dict = {"h": "qh", "k": "qk", "l": "ql", "e": "en"}
    def_x = spice_logs.attrs["def_x"]
    def_y = spice_logs.attrs["def_y"]
    if def_x in nexus_dict:
        def_x = nexus_dict[def_x]

    path_x = find_val(def_x, nxentry)
    path_y = find_val(def_y, nxentry)
    if def_y == "detector" or def_y == "monitor":
        path_y += "/data"

    # Create the LINKS
    nxentry["data/" + def_y] = h5py.SoftLink(nxentry.name + path_y)
    nxentry["data/" + def_y + "/"].attrs["target"] = nxentry.name + path_y
    # if def_y == "detector" or def_y == "monitor":
    #     nxentry["data"].attrs["signal"] = "data"
    # else:
    nxentry["data"].attrs["signal"] = def_y

    # Create the LINKS
    nxentry["data/" + def_x] = h5py.SoftLink(nxentry.name + path_x)
    nxentry["data/" + def_x + "/"].attrs["target"] = nxentry.name + path_x

    if def_x in nexus_dict:
        nxentry["data"].attrs["axes"] = nexus_dict[def_x]

    else:
        nxentry["data"].attrs["axes"] = def_x

    # # Create the DOC strings
    # nxentry["definition"].attrs["EX_doc"] = "Official NeXus NXDL schema to which this file conforms "
    # nxentry["sample/name"].attrs["EX_doc"] = "Descriptive name of sample "
    # nxentry["monitor/mode"].attrs[
    #     "EX_doc"
    # ] = "Count to a preset value based on either clock time (timer) or received monitor counts (monitor). "
    # nxentry["monitor/preset"].attrs["EX_doc"] = "preset value for time or monitor "
    # nxentry["monitor/data"].attrs["EX_doc"] = "Total integral monitor counts "
    # nxentry["data"].attrs[
    #     "EX_doc"
    # ] = "One of the ei,ef,qh,qk,ql,en should get a primary=1 attribute to denote the main scan axis "

    # Create the ATTRIBUTES
    # root["/"].attrs["default"] = "entry"
    # root["/entry"].attrs["default"] = "data"
    # nxentry["data"].attrs["signal"] = "data"
    # nxentry["data/data"].attrs["signal"] = "1"

    # --------------------------------------------------------------------------------------

    # # /entry/sample
    # nxsample = nxentry.create_group("sample")
    # nxsample.attrs["name"] = headers["samplename"]
    # nxsample.attrs["type"] = headers["sampletype"]
    # nxsample.attrs["mosiac"] = headers["samplemosaic"]

    # TODO sample environment -- temperture, magnetic field, pressure
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
    for t in temperatue_str:
        if t in spice_logs.keys():
            nxentry["sample"].create_dataset(name=t, data=spice_logs[t], maxshape=None)
            nxentry["sample/" + t].attrs["type"] = "NX_FLOAT"
            nxentry["sample/" + t].attrs["units"] = "K"

    # TODO field
    # TODO pressure


# TODO
def instrument_info_to_nexus(nxentry, instrument_params):
    """Extra info missing in SPICE, for resolution calculation

    Args:
        nxentry:
        instrument_params (dict)

    Note:
        This function does NOT overwrite exsiting parameters in SPICE,
        it only adds extra info to the nexus file.
    """
    source = instrument_params["source"]
    mono = instrument_params["monochromator"]
    monitor = instrument_params["monitor"]
    ana = instrument_params["analyser"]
    detector = instrument_params["detector"]
    collimators = instrument_params["collimators"]
    #  --------------------------- source ---------------------------
    # rectangular or circular
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    nxentry["instrument/source"].create_dataset(name="shape", data=source["shape"], maxshape=None)
    nxentry["instrument/source/shape"].attrs["type"] = "NX_CHAR"

    nxentry["instrument/source"].create_dataset(name="width", data=source["width"], maxshape=None)
    nxentry["instrument/source/width"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/source/width"].attrs["units"] = "cm"

    nxentry["instrument/source"].create_dataset(name="height", data=source["height"], maxshape=None)
    nxentry["instrument/source/height"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/source/height"].attrs["units"] = "cm"

    #  --------------------------- monochromator ---------------------------
    nxentry["instrument/monochromator"].create_dataset(name="d_spacing", data=mono["d_spacing"], maxshape=None)
    nxentry["instrument/monochromator/d_spacing"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/monochromator/d_spacing"].attrs["units"] = "Angstrom"

    nxentry["instrument/monochromator"].create_dataset(name="mosaic", data=mono["mosaic"], maxshape=None)
    nxentry["instrument/monochromator/mosaic"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/monochromator/mosaic"].attrs["units"] = "min"

    nxentry["instrument/monochromator"].create_dataset(name="mosaic_v", data=mono["mosaic_v"], maxshape=None)
    nxentry["instrument/monochromator/mosaic_v"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/monochromator/mosaic_v"].attrs["units"] = "min"

    # divide by np.sqrt(12) if rectangular
    nxentry["instrument/monochromator"].create_dataset(name="width", data=mono["width"], maxshape=None)
    nxentry["instrument/monochromator/width"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/monochromator/width"].attrs["units"] = "cm"

    # divide by np.sqrt(12) if rectangular
    nxentry["instrument/monochromator"].create_dataset(name="height", data=mono["height"], maxshape=None)
    nxentry["instrument/monochromator/height"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/monochromator/height"].attrs["units"] = "cm"

    # divide by np.sqrt(12) if rectangular
    nxentry["instrument/monochromator"].create_dataset(name="depth", data=mono["depth"], maxshape=None)
    nxentry["instrument/monochromator/depth"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/monochromator/depth"].attrs["units"] = "cm"

    # horizontal focusing
    nxentry["instrument/monochromator"].create_dataset(name="curvh", data=mono["curvh"], maxshape=None)
    nxentry["instrument/monochromator/curvh"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/monochromator/curvh"].attrs["units"] = "degrees"

    # vertical focusing
    nxentry["instrument/monochromator"].create_dataset(name="curvv", data=mono["curvv"], maxshape=None)
    nxentry["instrument/monochromator/curvv"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/monochromator/curvv"].attrs["units"] = "degrees"

    #  --------------------------- monitor ---------------------------
    nxentry["monitor"].create_dataset(name="width", data=monitor["width"], maxshape=None)
    nxentry["monitor/width"].attrs["type"] = "NX_FLOAT"
    nxentry["monitor/width"].attrs["units"] = "cm"

    # divide by np.sqrt(12) if rectangular
    nxentry["monitor"].create_dataset(name="height", data=monitor["height"], maxshape=None)
    nxentry["monitor/height"].attrs["type"] = "NX_FLOAT"
    nxentry["monitor/height"].attrs["units"] = "cm"

    #  --------------------------- analyzer ---------------------------
    nxentry["instrument/analyser"].create_dataset(name="d_spacing", data=ana["d_spacing"], maxshape=None)
    nxentry["instrument/analyser/d_spacing"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/analyser/d_spacing"].attrs["units"] = "Angstrom"

    nxentry["instrument/analyser"].create_dataset(name="mosaic", data=ana["mosaic"], maxshape=None)
    nxentry["instrument/analyser/mosaic"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/analyser/mosaic"].attrs["units"] = "min"

    nxentry["instrument/analyser"].create_dataset(name="mosaic_v", data=ana["mosaic_v"], maxshape=None)
    nxentry["instrument/analyser/mosaic_v"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/analyser/mosaic_v"].attrs["units"] = "min"

    # divide by np.sqrt(12) if rectangular
    nxentry["instrument/analyser"].create_dataset(name="width", data=ana["width"], maxshape=None)
    nxentry["instrument/analyser/width"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/analyser/width"].attrs["units"] = "cm"

    # divide by np.sqrt(12) if rectangular
    nxentry["instrument/analyser"].create_dataset(name="height", data=ana["height"], maxshape=None)
    nxentry["instrument/analyser/height"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/analyser/height"].attrs["units"] = "cm"

    # divide by np.sqrt(12) if rectangular
    nxentry["instrument/analyser"].create_dataset(name="depth", data=ana["depth"], maxshape=None)
    nxentry["instrument/analyser/depth"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/analyser/depth"].attrs["units"] = "cm"

    # horizontal focusing
    nxentry["instrument/analyser"].create_dataset(name="curvh", data=ana["curvh"], maxshape=None)
    nxentry["instrument/analyser/curvh"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/analyser/curvh"].attrs["units"] = "degrees"

    # vertical focusing
    nxentry["instrument/analyser"].create_dataset(name="curvv", data=ana["curvv"], maxshape=None)
    nxentry["instrument/analyser/curvv"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/analyser/curvv"].attrs["units"] = "degrees"

    #  --------------------------- detector ---------------------------
    # rectangular or circular
    # divide by np.sqrt(12) if rectangular
    # Diameter D/4 if spherical
    nxentry["instrument/detector"].create_dataset(name="shape", data=detector["shape"], maxshape=None)
    nxentry["instrument/detector/shape"].attrs["type"] = "NX_CHAR"

    # divide by np.sqrt(12) if rectangular
    nxentry["instrument/detector"].create_dataset(name="width", data=detector["width"], maxshape=None)
    nxentry["instrument/detector/width"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/detector/width"].attrs["units"] = "cm"

    # divide by np.sqrt(12) if rectangular
    nxentry["instrument/detector"].create_dataset(name="height", data=detector["height"], maxshape=None)
    nxentry["instrument/detector/height"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/detector/height"].attrs["units"] = "cm"

    #  --------------------------- collimators ---------------------------
    nxentry["instrument/collimator"].create_dataset(name="divergence_y", data=[150, 270, 300, 600], maxshape=None)
    nxentry["instrument/collimator/divergence_y"].attrs["type"] = "NX_ANGLE"
    nxentry["instrument/collimator/divergence_y"].attrs["units"] = "min"

    #  --------------------------- distances ---------------------------
    nxentry["instrument"].create_dataset(name="dist_src_mono", data=650.0, maxshape=None)
    nxentry["instrument/dist_src_mono"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/dist_src_mono"].attrs["units"] = "cm"

    nxentry["instrument"].create_dataset(name="dist_mono_sample", data=190.0, maxshape=None)
    nxentry["instrument/dist_mono_sample"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/dist_mono_sample"].attrs["units"] = "cm"

    nxentry["instrument"].create_dataset(name="dist_sample_ana", data=160.0, maxshape=None)
    nxentry["instrument/dist_sample_ana"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/dist_sample_ana"].attrs["units"] = "cm"

    nxentry["instrument"].create_dataset(name="dist_ana_det", data=60.0, maxshape=None)
    nxentry["instrument/dist_ana_det"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/dist_ana_det"].attrs["units"] = "cm"

    nxentry["instrument"].create_dataset(name="dist_mono_monitor", data=86.0, maxshape=None)
    nxentry["instrument/dist_mono_monitor"].attrs["type"] = "NX_FLOAT"
    nxentry["instrument/dist_mono_monitor"].attrs["units"] = "cm"

    # -----------------------------------sample-------------------------------------
    nxentry["sample"].create_dataset(name="shape", data="cylindrical", maxshape=None)
    nxentry["sample/shape"].attrs["type"] = "NX_CHAR"

    nxentry["sample"].create_dataset(name="width", data=1.0, maxshape=None)
    nxentry["sample/width"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/width"].attrs["units"] = "cm"

    nxentry["sample"].create_dataset(name="height", data=1.0, maxshape=None)
    nxentry["sample/height"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/height"].attrs["units"] = "cm"

    nxentry["sample"].create_dataset(name="depth", data=1.0, maxshape=None)
    nxentry["sample/depth"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/depth"].attrs["units"] = "cm"

    nxentry["sample"].create_dataset(name="mosaic", data=1.0, maxshape=None)
    nxentry["sample/mosaic"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/mosaic"].attrs["units"] = "min"

    nxentry["sample"].create_dataset(name="mosaic_v", data=1.0, maxshape=None)
    nxentry["sample/mosaic_v"].attrs["type"] = "NX_FLOAT"
    nxentry["sample/mosaic_v"].attrs["units"] = "min"


def convert_spice_to_nexus(path_to_spice_folder, path_to_hdf5):
    """Load data from spice folder. Convert to a nexus file

    Args:
        path_to_spice_folder (str): spice folder, ends with '/'
        path_to_nexus (str): path to hdf5 data file, ends with '.h5'
        instrument_config: python file contains instrument configuration parameters
    """

    print(f"Converting {path_to_spice_folder} to {path_to_hdf5}")

    p = Path(path_to_spice_folder)

    with h5py.File(path_to_hdf5, "w") as root:

        # ----------------------------- ub info ------------------------------------
        ub_files = sorted((p / "UBConf").glob("*.ini"))
        tmp_ub_files = sorted((p / "UBConf/tmp").glob("*.ini"))
        ub_entries = root.create_group("UBConf")
        ub_entries.attrs["NX_class"] = "NXcollection"
        ub_entries.attrs["X_required"] = "false"

        for ub_file in ub_files + tmp_ub_files:
            ub_entry_name = ub_file.parts[-1].split(".")[0]
            ub_entry = ub_entries.create_group(ub_entry_name)
            ub_conf = read_spice_ub(ub_file)
            for ub_item in ub_conf.items():
                k, v = ub_item
                ub_entry.create_dataset(name=k, data=v, maxshape=None)

        # --------------------------- scan info ------------------------------------
        scans = sorted((p / "Datafiles").glob("*.dat"))
        instrument_str = scans[0].parts[-1].split("_")[0]

        # read in exp_info from the first scan and save as attibutes of the file
        _, _, headers, _ = read_spice(scans[0])
        ipts = headers["proposal"]
        exp_num = headers["experiment_number"]

        # root.attrs["name"] = headers["experiment"]
        # root.attrs["users"] = headers["users"]
        # root.attrs["local_contact"] = headers["local_contact"]

        # convert SPICE scans into nexus entries
        for scan in scans:  # ignoring unused keys
            spice_data, col_headers, headers, unused = read_spice(scan)

            # pre-processing
            scan_num = ((scan.parts[-1].split("_"))[-1]).split(".")[0]  # e.g. "scan0001"
            if "scan_title" in unused:
                headers["scan_title"] = ""

            # /entry/SPICElogs0
            nxentry = root.create_group(scan_num)
            spice_logs = nxentry.create_group("SPICElogs")
            spice_logs.attrs["NX_class"] = "NXcollection"
            spice_logs.attrs["X_required"] = "false"
            spice_logs.attrs["instrument"] = instrument_str

            # metadata to attibutes
            exp_str = ["scan_title", "users", "local_contact", "experiment"]
            for k, v in headers.items():

                if "," in v and k not in exp_str:  # vectors
                    spice_logs.attrs[k] = np.array([float(v0) for v0 in v.split(",")])
                elif v.replace(".", "").isnumeric():  # numebrs only
                    if v.isdigit():  # int
                        spice_logs.attrs[k] = int(v)
                    else:  # float
                        spice_logs.attrs[k] = float(v)
                # separate COM/FWHM and its errorbar
                elif k == "Center of Mass":
                    com, e_com = v.split("+/-")
                    spice_logs.attrs["COM"] = float(com)
                    spice_logs.attrs["COM_err"] = float(e_com)
                elif k == "Full Width Half-Maximum":
                    fwhm, e_fwhm = v.split("+/-")
                    spice_logs.attrs["FWHM"] = float(fwhm)
                    spice_logs.attrs["FWHM_err"] = float(e_fwhm)
                else:  # other crap, keep as is
                    spice_logs.attrs[k] = v

            # motor position table to datasets
            if spice_data.ndim == 1:  # empty data or 1 point only
                if len(spice_data):  # 1 point only
                    for idx, col_header in enumerate(col_headers):
                        spice_logs.create_dataset(col_header, data=spice_data[idx])
                else:  # empty
                    pass
            else:  # nomarl data
                for idx, col_header in enumerate(col_headers):
                    spice_logs.create_dataset(col_header, data=spice_data[:, idx])

            spicelogs_to_nexus(nxentry)
            # instrument_info_to_nexus(nxentry, instrument_config)

        # Create the ATTRIBUTES
        root.attrs["file_name"] = os.path.abspath(f"IPTS{ipts}_{instrument_str}_exp{exp_num}")
        root.attrs["file_time"] = datetime.now().isoformat()
        root.attrs["h5py_version"] = h5py.version.version
        root.attrs["HDF5_Version"] = h5py.version.hdf5_version


if __name__ == "__main__":
    spice_folder = "./tests/test_data_folder/exp424/"  # CG4C
    # spice_folder = "./tests/test_data_folder/exp758/" # NB3
    # spice_folder = "./tests/test_data_folder/exp932/"  # HB1 polarized
    # h5_file_name = "./tests/test_data_folder/tavi_exp758.h5"
    nexus_file_name = "./tests/test_data_folder/nexus_exp424.h5"

    convert_spice_to_nexus(spice_folder, nexus_file_name)
    # TODO update instrument parameters and sample parameters of a single scan
    #  instrument_params)