# -*- coding: utf-8 -*-
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np

from tavi.data.spice_reader import _create_spicelogs, read_spice_datafile


def _find_val(val, grp, prefix=""):
    """Find value in hdf5 groups"""
    for obj_name, obj in grp.items():
        if obj_name in ("SPICElogs", "data"):
            continue
        else:
            path = f"{prefix}/{obj_name}"
            if val == obj_name:
                return path
            # test for group (go down)
            elif isinstance(obj, h5py.Group):
                gpath = _find_val(val, obj, path)
                if gpath:
                    return gpath


def _recast_type(ds, dtype):
    if type(ds) is str:
        if "," in ds:  # vector
            if dtype == "NX_FLOAT":
                dataset = np.array([float(d) for d in ds.split(",")])
            else:
                dataset = np.array([int(d) for d in ds.split(",")])
        elif ds.replace(".", "").isnumeric():  # numebrs only
            if dtype == "NX_FLOAT":
                dataset = float(ds)
            else:
                dataset = int(ds)
    else:  # expect np.ndarray
        if dtype == "NX_FLOAT":
            dataset = np.array([float(d) for d in ds])
        else:
            dataset = np.array([int(d) for d in ds])
    return dataset


class NXdataset(dict):
    """Dataset in a format consistent with NeXus, containg attrs and dataset"""

    def __init__(self, ds, **kwargs):
        if ds is None:
            return None

        attr_dict = {}
        for k, v in kwargs.items():
            attr_dict.update({k: v})

        match kwargs:
            case {"type": "NX_CHAR"}:
                dataset = str(ds)
            case {"type": "NX_INT"} | {"type": "NX_FLOAT"}:
                dataset = _recast_type(ds, kwargs["type"])
            case {"type": "NX_DATE_TIME"}:
                dataset = datetime.strptime(ds, "%m/%d/%Y %I:%M:%S %p").isoformat()
            case _:
                dataset = ds

        return super().__init__([("dataset", dataset), ("attrs", attr_dict)])

    def get_attr(self, key: str, default=None):
        val = self["attrs"].get(key)
        return val if val is not None else default

    def get_dataset(self, default=None):
        val = self["dataset"]
        return val if val is not None else default


class NXentry(dict):
    """Entry in a format consistent with NeXus"""

    def __init__(self, **kwargs):
        attr_dict = {}
        dataset_list = []
        for k, v in kwargs.items():
            if isinstance(v, NXdataset) or isinstance(v, NXentry):
                if not v:  # ignore if empty
                    pass
                else:
                    dataset_list.append((k, v))
            else:  # must be an attribute
                attr_dict.update({k: v})
        return super().__init__([("attrs", attr_dict)] + dataset_list)

    def add_dataset(self, key: str, ds: NXdataset):
        if not ds:  # ignore if empty
            pass
        else:
            self.update({key: ds})

    def add_attribute(self, key: str, attr):
        if not attr:  # ignore if empty
            pass
        else:
            self["attrs"].update({key: attr})


def _formatted_spicelogs(spicelogs: dict) -> NXentry:
    """Format SPICE logs into NeXus dict"""
    formatted_spicelogs = NXentry(NX_class="NXcollection", EX_required="false")
    metadata = spicelogs.pop("metadata")
    for attr_key, attr_entry in metadata.items():
        formatted_spicelogs.add_attribute(attr_key, attr_entry)

    for entry_key, entry_data in spicelogs.items():
        formatted_spicelogs.add_dataset(key=entry_key, ds=NXdataset(ds=entry_data))
    return formatted_spicelogs


# TODO json support
def spice_scan_to_nxdict(
    path_to_scan_file: str,
    path_to_instrument_json: Optional[str] = None,
    path_to_sample_json: Optional[str] = None,
) -> NXentry:
    """Format SPICE data in a nested dictionary format

    Note:
        json files can overwrite the parameters in SPICE
    """

    # parse instruemnt and sample json files
    instrument_config_params = None
    if path_to_instrument_json is not None:
        instrument_config = Path(path_to_instrument_json)
        if not instrument_config.is_file():
            raise ValueError(f"Invalid instrument json file path: {path_to_instrument_json}")
        with open(instrument_config, "r", encoding="utf-8") as file:
            instrument_config_params = json.load(file)

    sample_config_params = None
    if path_to_sample_json is not None:
        sample_config = Path(path_to_sample_json)
        if not sample_config.is_file():
            raise ValueError(f"Invalid sample json file path: {path_to_instrument_json}")
        with open(sample_config, "r", encoding="utf-8") as file:
            sample_config_params = json.load(file)

    spicelogs = _create_spicelogs(path_to_scan_file)
    metadata = spicelogs["metadata"]

    # ---------------------------------------- source ----------------------------------------

    nxsource = NXentry(
        name=NXdataset(ds="HFIR", type="NX_CHAR", EX_required="true"),
        probe=NXdataset(ds="neutron", type="NX_CHAR", EX_required="true"),
        NX_class="NXsource",
        EX_required="true",
    )
    # ------------------------------------- monochromator ----------------------------------------

    nxmono = NXentry(
        ei=NXdataset(ds=spicelogs.get("ei"), type="NX_FLOAT", EX_required="true", units="meV"),
        type=NXdataset(ds=metadata.get("monochromator"), type="NX_CHAR"),
        sense=NXdataset(ds=metadata["sense"][0], type="NX_CHAR"),
        m1=NXdataset(ds=spicelogs.get("m1"), type="NX_FLOAT", units="degrees"),
        m2=NXdataset(ds=spicelogs.get("m2"), type="NX_FLOAT", units="degrees"),
        NX_class="NXcrystal",
        EX_required="true",
    )
    nxmono.add_dataset("mfocus", NXdataset(ds=spicelogs.get("mfocus"), type="NX_FLOAT"))
    nxmono.add_dataset("marc", NXdataset(ds=spicelogs.get("marc"), type="NX_FLOAT"))
    nxmono.add_dataset("mtrans", NXdataset(ds=spicelogs.get("mtrans"), type="NX_FLOAT"))
    nxmono.add_dataset("focal_length", NXdataset(ds=spicelogs.get("focal_length"), type="NX_FLOAT"))

    # ------------------------------------- analyzer ----------------------------------------

    nxana = NXentry(
        ef=NXdataset(ds=spicelogs.get("ef"), type="NX_FLOAT", EX_required="true", units="meV"),
        type=NXdataset(ds=metadata.get("analyzer"), type="NX_CHAR"),
        sense=NXdataset(ds=metadata["sense"][2], type="NX_CHAR"),
        a1=NXdataset(ds=spicelogs.get("a1"), type="NX_FLOAT", units="degrees"),
        a2=NXdataset(ds=spicelogs.get("a2"), type="NX_FLOAT", units="degrees"),
        afocus=NXdataset(ds=spicelogs.get("afocus"), type="NX_FLOAT"),
        NX_class="NXcrystal",
        EX_required="true",
    )
    for i in range(8):  # CG4C horizontal focusing
        nxana.add_dataset(key=f"qm{i+1}", ds=NXdataset(ds=spicelogs.get(f"qm{i+1}"), type="NX_FLOAT"))
        nxana.add_dataset(key=f"xm{i+1}", ds=NXdataset(ds=spicelogs.get(f"xm{i+1}"), type="NX_FLOAT"))

    # ------------------------------------- detector ----------------------------------------

    nxdet = NXentry(
        data=NXdataset(ds=spicelogs.get("detector"), type="NX_INT", EX_required="true", units="counts"),
        NX_class="NXdetector",
        EX_required="true",
    )

    # ------------------------------------- collimators ----------------------------------------

    div_x = [float(v) for v in list(metadata["collimation"].split("-"))]
    nxcoll = NXentry(
        type=NXdataset(ds="Soller", type="NX_CHAR"),
        NX_class="NXcollimator",
    )
    nxcoll.add_dataset(key="divergence_x", ds=NXdataset(ds=div_x, type="NX_ANGLE", units="minutes of arc"))

    # ------------------------------------- slits ----------------------------------------

    nxslits = NXentry(NX_class="NXslit")
    slits_str1 = tuple([f"b{idx}{loc}" for idx in ("a", "b") for loc in ("t", "b", "l", "r")])
    slits_str2 = tuple([f"slit{idx}_{loc}" for idx in ("a", "b") for loc in ("lf", "rt", "tp", "bt")])
    slits_str3 = tuple([f"slit_{idx}_{loc}" for idx in ("pre",) for loc in ("lf", "rt", "tp", "bt")])
    slits_str = (slits_str1, slits_str2, slits_str3)
    for slit_str in slits_str:
        for st in slit_str:
            nxslits.add_dataset(key=st, ds=NXdataset(ds=spicelogs.get(st), type="NX_FLOAT", units="cm"))
    # ------------------------------------- flipper ----------------------------------------

    nxflipper = NXentry(NX_class="NXflipper")
    nxflipper.add_dataset(key="fguide", ds=NXdataset(ds=spicelogs.get("fguide"), type="NX_FLOAT"))
    nxflipper.add_dataset(key="hguide", ds=NXdataset(ds=spicelogs.get("hguide"), type="NX_FLOAT"))
    nxflipper.add_dataset(key="vguide", ds=NXdataset(ds=spicelogs.get("vguide"), type="NX_FLOAT"))

    # ---------------------------------------- instrument ---------------------------------------------

    nxinstrument = NXentry(
        source=nxsource,
        monochromator=nxmono,
        collimator=nxcoll,
        analyzer=nxana,
        detector=nxdet,
        slits=nxslits,
        flipper=nxflipper,
        name=NXdataset(ds=metadata.get("instrument"), type="NX_CHAR"),
        NX_class="NXinstrument",
        EX_required="true",
    )
    # ---------------------------------------- monitor ---------------------------------------------

    preset_type = metadata.get("preset_type")
    if preset_type == "normal":
        preset_channel = metadata.get("preset_channel")
        nxmonitor = NXentry(
            mode=NXdataset(ds=preset_channel, type="NX_CHAR", EX_required="true"),
            preset=NXdataset(ds=metadata.get("preset_value"), type="NX_FLOAT", EX_required="true"),
            time=NXdataset(ds=spicelogs.get("time"), type="NX_FLOAT", units="seconds"),
            monitor=NXdataset(ds=spicelogs.get("monitor"), type="NX_INT", units="counts"),
            mcu=NXdataset(ds=spicelogs.get("mcu"), type="NX_FLOAT", units="mcu"),
            # TODO link to the right channel
            # data=NXdataset(),
            NX_class="NXmonitor",
            EX_required="true",
        )

    # TODO polarized exp at HB1
    elif preset_type == "countfile":
        print("Polarization data, not yet supported.")
        nxmonitor = NXentry(NX_class="NXmonitor", EX_required="true")
    else:
        print(f"Unrecogonized preset type {preset_type}.")
        nxmonitor = NXentry(NX_class="NXmonitor", EX_required="true")

    # ---------------------------------------- sample ---------------------------------------------

    nxsample = NXentry(
        name=NXdataset(ds=metadata.get("samplename"), type="NX_CHAR", EX_required="true"),
        type=NXdataset(ds=metadata.get("sampletype"), type="NX_CHAR", EX_required="true"),
        unit_cell=NXdataset(ds=metadata.get("latticeconstants"), type="NX_FLOAT", EX_required="true"),
        qh=NXdataset(ds=spicelogs.get("h"), type="NX_FLOAT", EX_required="true"),
        qk=NXdataset(ds=spicelogs.get("k"), type="NX_FLOAT", EX_required="true"),
        ql=NXdataset(ds=spicelogs.get("l"), type="NX_FLOAT", EX_required="true"),
        en=NXdataset(ds=spicelogs.get("e"), type="NX_FLOAT", EX_required="true", units="meV"),
        q=NXdataset(ds=spicelogs.get("q"), type="NX_FLOAT"),
        sense=NXdataset(ds=metadata["sense"][1], type="NX_CHAR"),
        NX_class="NXsample",
        EX_required="true",
    )
    nxsample.add_dataset(key="Pt.", ds=NXdataset(ds=spicelogs.get("Pt."), type="NX_INT"))

    # Motor angles
    nxsample.add_dataset(key="s1", ds=NXdataset(ds=spicelogs.get("s1"), type="NX_FLOAT", units="degrees"))
    nxsample.add_dataset(key="s2", ds=NXdataset(ds=spicelogs.get("s2"), type="NX_FLOAT", units="degrees"))
    nxsample.add_dataset(key="sgu", ds=NXdataset(ds=spicelogs.get("sgu"), type="NX_FLOAT", units="degrees"))
    nxsample.add_dataset(key="sgl", ds=NXdataset(ds=spicelogs.get("sgl"), type="NX_FLOAT", units="degrees"))
    nxsample.add_dataset(key="stu", ds=NXdataset(ds=spicelogs.get("stu"), type="NX_FLOAT", units="degrees"))
    nxsample.add_dataset(key="stl", ds=NXdataset(ds=spicelogs.get("stl"), type="NX_FLOAT", units="degrees"))

    # UB info
    nxsample.add_dataset(
        key="orientation_matrix",
        ds=NXdataset(ds=metadata.get("ubmatrix"), type="NX_FLOAT", EX_required="true", units="NX_DIMENSIONLESS"),
    )
    nxsample.add_dataset(key="ub_conf", ds=NXdataset(ds=metadata["ubconf"].split(".")[0], type="NX_CHAR"))
    nxsample.add_dataset(key="plane_normal", ds=NXdataset(ds=metadata.get("plane_normal"), type="NX_FLOAT"))

    # ---------------------------------------- sample environment ---------------------------------------------

    # TODO all sample environment variable names needed!
    temperatue_str = (
        ("temp", "temp_a", "temp_2", "coldtip", "tsample", "sample")
        + ("vti", "dr_tsample", "dr_temp")
        + ("lt", "ht", "sorb_temp", "sorb", "sample_ht")
    )
    for te in temperatue_str:
        nxsample.add_dataset(key=te, ds=NXdataset(ds=spicelogs.get(te), type="NX_FLOAT", units="Kelvin"))

    field_str = ("persistent_field", "mag_i")
    for fi in field_str:
        nxsample.add_dataset(key=fi, ds=NXdataset(ds=spicelogs.get(fi), type="NX_FLOAT", units="Tesla"))
    # ---------------------------------------- data ---------------------------------------------
    nexus_keywork_conversion_dict = {"h": "qh", "k": "qk", "l": "ql", "e": "en"}
    def_x = metadata.get("def_x")
    def_y = metadata.get("def_y")
    if def_x in nexus_keywork_conversion_dict:
        def_x = nexus_keywork_conversion_dict[def_x]

    nxdata = NXentry(NX_class="NXdata", EX_required="true", signal=def_y, axes=def_x)
    # ---------------------------------------- scan ---------------------------------------------

    # TODO timezone
    start_date_time = "{} {}".format(metadata.get("date"), metadata.get("time"))
    # TODO what is last scan never finished?
    # if "end_time" in das_logs.attrs:
    end_date_time = metadata.get("end_time")

    nxscan = NXentry(
        SPICElogs=_formatted_spicelogs(spicelogs),
        data=nxdata,
        definition=NXdataset(ds="NXtas", type="NX_CHAR", EX_required="true"),
        title=NXdataset(ds=metadata.get("scan_title"), type="NX_CHAR", EX_required="true"),
        start_time=NXdataset(ds=start_date_time, type="NX_DATE_TIME", EX_required="true"),
        end_time=NXdataset(ds=end_date_time, type="NX_DATE_TIME"),
        instrument=nxinstrument,
        monitor=nxmonitor,
        sample=nxsample,
        NX_class="NXentry",
        EX_required="true",
    )
    return nxscan


def spice_data_to_nxdict(
    path_to_spice_folder: str,
    scan_num: Optional[int] = None,
    path_to_instrument_json: Optional[str] = None,
    path_to_sample_json: Optional[str] = None,
) -> dict:
    """Read SPICE data files into nested dictionary, ready to be converted to NeXus format

    Args:
           path_to_spice_folder (str): path to a SPICE folder
           scan_num (int): read all scans in folder if not None, otherwise read one scan only
           path_to_instrument_json: Optional[str] = None,
           path_to_sample_json: Optional[str] = None,
    """
    if path_to_spice_folder[-1] != "/":
        path_to_spice_folder += "/"

    scan_list = os.listdir(path_to_spice_folder + "Datafiles/")
    if scan_num is None:  # read all scans in folder
        filter_keyword = ".dat"
    else:  # read one scan only
        filter_keyword = f"scan{scan_num:04}.dat"

    scan_list = [path_to_spice_folder + "Datafiles/" + scan for scan in scan_list if scan.endswith(filter_keyword)]
    scan_list.sort()

    # get IPTS number and instrument string
    first_scan = scan_list[0]
    (_, _, headers, _, _) = read_spice_datafile(first_scan)
    ipts = headers["proposal"]
    spice_file_name = first_scan.split("/")[-1]  # e.g. CG4C_exp0424_scan0001.dat
    instrument_str, exp_num, _ = spice_file_name.split("_")
    dataset_name = f"IPTS{ipts}_{instrument_str}_{exp_num}"

    nexus_dict = {}
    for path_to_scan_file in scan_list:
        pass

        scan_name, scan_dict = spice_scan_to_nxdict(
            path_to_scan_file,
            path_to_instrument_json,
            path_to_sample_json,
        )

        nexus_dict.update({scan_name: scan_dict})
        nexus_dict[scan_name]["attrs"].update({"dataset_name": dataset_name})

    return nexus_dict