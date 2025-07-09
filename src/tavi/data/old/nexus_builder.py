# -*- coding: utf-8 -*-
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Optional
from zoneinfo import ZoneInfo

import h5py
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
    if isinstance(ds, str):
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
    elif isinstance(ds, float | int):
        if dtype == "NX_FLOAT":
            dataset = float(ds)
        else:
            dataset = int(ds)
    elif isinstance(ds, np.ndarray) or isinstance(ds, list):
        ds = np.array(ds)
        sz = np.shape(ds)
        if dtype == "NX_FLOAT":
            dataset = np.array([float(d) for d in ds.flat])
        else:
            dataset = np.array([int(d) for d in ds.flat])
        dataset = dataset.reshape(sz)
    else:
        raise TypeError(f"dataset={ds} has unrecogonized type {type(ds)}")
    return dataset


class NXdataset(dict):
    """Dataset in a format consistent with NeXus, containg attrs and dataset"""

    def __init__(self, ds, **kwargs):
        if ds is None:
            return None
        elif isinstance(ds, str):
            if not ds:
                return None

        attr_dict = {}
        for k, v in kwargs.items():
            attr_dict.update({k: v})

        match kwargs:
            case {"type": "NX_CHAR"}:
                dataset = np.array(ds) if isinstance(ds, list) else str(ds)

            case {"type": "NX_INT"} | {"type": "NX_FLOAT"}:
                dataset = _recast_type(ds, kwargs["type"])
            case {"type": "NX_DATE_TIME"}:
                dt = datetime.strptime(ds, "%m/%d/%Y %I:%M:%S %p")
                date = datetime(
                    dt.year,
                    dt.month,
                    dt.day,
                    dt.hour,
                    dt.minute,
                    dt.second,
                    tzinfo=ZoneInfo("America/New_York"),
                )
                dataset = date.isoformat()
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
    if spicelogs.get("ub_conf") is not None:
        ub_conf = spicelogs.pop("ub_conf")
        ub_file_path = ub_conf.pop("file_path")
        formatted_ub_conf = NXentry(file_path=ub_file_path, NX_class="NXcollection", EX_required="false")
        for entry_key, entry_data in ub_conf.items():
            formatted_ub_conf.add_dataset(key=entry_key, ds=NXdataset(ds=entry_data))

    metadata = spicelogs.pop("metadata")
    for attr_key, attr_entry in metadata.items():
        formatted_spicelogs.add_attribute(attr_key, attr_entry)

    for entry_key, entry_data in spicelogs.items():
        formatted_spicelogs.add_dataset(key=entry_key, ds=NXdataset(ds=entry_data))

    return formatted_spicelogs


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
    # print(f"scan#{metadata['scan']}")
    # ---------------------------------------- user ----------------------------------------

    nxuser = NXentry(NX_class="NXuser")
    nxuser.add_dataset("name", NXdataset(ds=metadata["users"].split(","), type="NX_CHAR"))

    # ---------------------------------------- source ----------------------------------------

    nxsource = NXentry(
        name=NXdataset(ds="HFIR", type="NX_CHAR", EX_required="true"),
        probe=NXdataset(ds="neutron", type="NX_CHAR", EX_required="true"),
        NX_class="NXsource",
        EX_required="true",
    )

    if instrument_config_params is not None:
        if (source_params := instrument_config_params.get("source")) is not None:
            for param, value in source_params.items():
                match param:
                    case "shape":
                        nxsource.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_CHAR"))
                    case "width" | "height":
                        nxsource.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_FLOAT", units="cm"))
                    case _:
                        nxsource.add_dataset(key=param, ds=NXdataset(ds=value))

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

    if instrument_config_params is not None:
        if (mono_params := instrument_config_params.get("monochromator")) is not None:
            for param, value in mono_params.items():
                match param:
                    case "shape":
                        nxmono.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_CHAR"))
                    case "width" | "height" | "depth":
                        nxmono.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_FLOAT", units="cm"))
                    case "mosaic_h" | "mosaic_v":
                        nxmono.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_FLOAT", units="mins of arc"))
                    case "curved_h" | "optimally_curved_h" | "curved_v" | "optimally_curved_v":
                        nxmono.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_CHAR"))
                    case "curvh" | "curvv":
                        nxmono.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_FLOAT", units="cm^-1"))
                    case _:
                        nxmono.add_dataset(key=param, ds=NXdataset(ds=value))

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
        nxana.add_dataset(key=f"qm{i + 1}", ds=NXdataset(ds=spicelogs.get(f"qm{i + 1}"), type="NX_FLOAT"))
        nxana.add_dataset(key=f"xm{i + 1}", ds=NXdataset(ds=spicelogs.get(f"xm{i + 1}"), type="NX_FLOAT"))

    if instrument_config_params is not None:
        if (ana_params := instrument_config_params.get("analyzer")) is not None:
            for param, value in ana_params.items():
                match param:
                    case "shape":
                        nxana.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_CHAR"))
                    case "width" | "height" | "depth":
                        nxana.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_FLOAT", units="cm"))
                    case "mosaic_h" | "mosaic_v":
                        nxana.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_ANGLE", units="minutes of arc"))
                    case "curved_h" | "optimally_curved_h" | "curved_v" | "optimally_curved_v":
                        nxana.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_CHAR"))
                    case "curvh" | "curvv":
                        nxana.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_FLOAT", units="cm^-1"))
                    case _:
                        nxana.add_dataset(key=param, ds=NXdataset(ds=value))

    # ------------------------------------- detector ----------------------------------------

    nxdet = NXentry(NX_class="NXdetector", EX_required="true")

    if (ds := spicelogs.get("detector")) is not None:  # single detector
        nxdet.add_dataset(key="data", ds=NXdataset(ds=ds, type="NX_INT", EX_required="true", units="counts"))

    if instrument_config_params is not None:
        if (det_params := instrument_config_params.get("detector")) is not None:
            for param, value in det_params.items():
                match param:
                    case "shape":
                        nxdet.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_CHAR"))
                    case "width" | "height":
                        nxdet.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_FLOAT", units="cm"))
                    case _:
                        nxdet.add_dataset(key=param, ds=NXdataset(ds=value))

    # for HB1 in polarized mode
    idx = 1
    while f"detector_{idx}" in spicelogs.keys():
        nxdet.add_dataset(
            key=f"detector_{idx}", ds=NXdataset(ds=spicelogs.get(f"detector_{idx}"), type="NX_INT", units="counts")
        )
        idx += 1

    # ------------------------------------- collimators ----------------------------------------

    div_x = [float(v) for v in list(metadata["collimation"].split("-"))]
    nxcoll = NXentry(
        type=NXdataset(ds="Soller", type="NX_CHAR"),
        NX_class="NXcollimator",
    )
    nxcoll.add_dataset(key="divergence_x", ds=NXdataset(ds=div_x, type="NX_ANGLE", units="minutes of arc"))

    if instrument_config_params is not None:
        if (coll_params := instrument_config_params.get("collimators")) is not None:
            try:
                div_x = [
                    coll_params["h_pre_mono"],
                    coll_params["h_pre_sample"],
                    coll_params["h_post_sample"],
                    coll_params["h_post_ana"],
                ]
                nxcoll.add_dataset(key="divergence_x", ds=NXdataset(ds=div_x, type="NX_ANGLE", units="minutes of arc"))
            except KeyError:
                pass

            try:
                div_y = [
                    coll_params["v_pre_mono"],
                    coll_params["v_pre_sample"],
                    coll_params["v_post_sample"],
                    coll_params["v_post_ana"],
                ]
                nxcoll.add_dataset(key="divergence_y", ds=NXdataset(ds=div_y, type="NX_ANGLE", units="minutes of arc"))
            except KeyError:
                pass

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

    # ------------------------------------- goniometer ----------------------------------------

    nxgoni = NXentry(
        sense=NXdataset(ds=metadata["sense"][1], type="NX_CHAR"),
        NX_class="NXtransformations",
    )
    # Motor angles
    nxgoni.add_dataset(key="s1", ds=NXdataset(ds=spicelogs.get("s1"), type="NX_FLOAT", units="degrees"))
    nxgoni.add_dataset(key="s2", ds=NXdataset(ds=spicelogs.get("s2"), type="NX_FLOAT", units="degrees"))
    nxgoni.add_dataset(key="sgu", ds=NXdataset(ds=spicelogs.get("sgu"), type="NX_FLOAT", units="degrees"))
    nxgoni.add_dataset(key="sgl", ds=NXdataset(ds=spicelogs.get("sgl"), type="NX_FLOAT", units="degrees"))
    nxgoni.add_dataset(key="stu", ds=NXdataset(ds=spicelogs.get("stu"), type="NX_FLOAT", units="degrees"))
    nxgoni.add_dataset(key="stl", ds=NXdataset(ds=spicelogs.get("stl"), type="NX_FLOAT", units="degrees"))
    nxgoni.add_dataset(key="omega", ds=NXdataset(ds=spicelogs.get("omega"), type="NX_FLOAT", units="degrees"))
    nxgoni.add_dataset(key="two_theta", ds=NXdataset(ds=spicelogs.get("two_theta"), type="NX_FLOAT", units="degrees"))
    nxgoni.add_dataset(key="2theta", ds=NXdataset(ds=spicelogs.get("2theta"), type="NX_FLOAT", units="degrees"))
    nxgoni.add_dataset(key="chi", ds=NXdataset(ds=spicelogs.get("chi"), type="NX_FLOAT", units="degrees"))
    nxgoni.add_dataset(key="phi", ds=NXdataset(ds=spicelogs.get("phi"), type="NX_FLOAT", units="degrees"))

    if instrument_config_params is not None:
        if (goni_params := instrument_config_params.get("goniometer")) is not None:
            for param, value in goni_params.items():
                match param:
                    case "type":
                        nxgoni.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_CHAR"))

    # ---------------------------------------- monitor ---------------------------------------------

    preset_type = metadata.get("preset_type")
    if preset_type == "normal":
        preset_channel = metadata.get("preset_channel")
        nxmonitor = NXentry(
            preset_type=NXdataset(ds="normal", type="NX_CHAR"),
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

    # polarized exp at HB1
    elif preset_type == "countfile":
        countfile = metadata.get("countfile")
        # lines = countfile.split(",")
        # preset_list = []
        # for line in lines:
        #     if line.startswith(" count preset"):
        #         *_, mode, preset = line.split(" ")
        #         preset_list.append((mode, preset))

        nxmonitor = NXentry(
            preset_type=NXdataset(ds="countfile", type="NX_CHAR"),
            countfile=NXdataset(ds=countfile.split(","), type="NX_CHAR"),
            NX_class="NXmonitor",
            EX_required="true",
        )
        idx = 1
        while f"detector_{idx}" in spicelogs.keys():
            mode = metadata.get("norm_channels")[idx - 1]
            preset = metadata.get("norm_vals")[idx - 1]
            nxmonitor.add_dataset(key=f"mode_{idx}", ds=NXdataset(ds=mode, type="NX_CHAR"))
            nxmonitor.add_dataset(key=f"preset_{idx}", ds=NXdataset(ds=preset, type="NX_FLOAT"))
            time = spicelogs.get(f"time_{idx}")
            nxmonitor.add_dataset(key=f"time_{idx}", ds=NXdataset(ds=time, type="NX_FLOAT", units="seconds"))
            monitor = spicelogs.get(f"monitor_{idx}")
            nxmonitor.add_dataset(key=f"monitor_{idx}", ds=NXdataset(ds=monitor, type="NX_INT", units="counts"))
            mcu = spicelogs.get(f"mcu_{idx}")
            nxmonitor.add_dataset(key=f"mcu_{idx}", ds=NXdataset(ds=mcu, type="NX_FLOAT", units="mcu"))
            label = metadata.get("labels")[idx - 1]
            nxmonitor.add_dataset(key=f"label_{idx}", ds=NXdataset(ds=label, type="NX_CHAR"))
            idx += 1
    else:
        print(f"Unrecogonized preset type {preset_type} for scan #{metadata.get('scan')}.")
        nxmonitor = NXentry(
            time=NXdataset(ds=spicelogs.get("time"), type="NX_FLOAT", units="seconds"),
            monitor=NXdataset(ds=spicelogs.get("monitor"), type="NX_INT", units="counts"),
            mcu=NXdataset(ds=spicelogs.get("mcu"), type="NX_FLOAT", units="mcu"),
            NX_class="NXmonitor",
            EX_required="true",
        )

    if instrument_config_params is not None:
        if (monitor_params := instrument_config_params.get("monitor")) is not None:
            for param, value in monitor_params.items():
                match param:
                    case "shape":
                        nxmonitor.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_CHAR"))

    # ---------------------------------------- instrument ---------------------------------------------
    nxinstrument = NXentry(
        source=nxsource,
        monochromator=nxmono,
        collimator=nxcoll,
        analyser=nxana,
        detector=nxdet,
        slits=nxslits,
        flipper=nxflipper,
        goniometer=nxgoni,
        name=NXdataset(ds=metadata.get("instrument"), type="NX_CHAR"),
        NX_class="NXinstrument",
        EX_required="true",
    )

    # ---------------------------------------- sample ---------------------------------------------

    nxsample = NXentry(
        name=NXdataset(ds=metadata.get("samplename"), type="NX_CHAR", EX_required="true"),
        type=NXdataset(ds=metadata.get("sampletype"), type="NX_CHAR", EX_required="true"),
        unit_cell=NXdataset(ds=metadata.get("latticeconstants"), type="NX_FLOAT", EX_required="true"),
        qh=NXdataset(ds=spicelogs.get("h"), type="NX_FLOAT", EX_required="true", units="r.l.u."),
        qk=NXdataset(ds=spicelogs.get("k"), type="NX_FLOAT", EX_required="true", units="r.l.u."),
        ql=NXdataset(ds=spicelogs.get("l"), type="NX_FLOAT", EX_required="true", units="r.l.u."),
        en=NXdataset(ds=spicelogs.get("e"), type="NX_FLOAT", EX_required="true", units="meV"),
        q=NXdataset(ds=spicelogs.get("q"), type="NX_FLOAT", units="Angstrom^-1"),
        NX_class="NXsample",
        EX_required="true",
    )
    nxsample.add_dataset(key="Pt.", ds=NXdataset(ds=spicelogs.get("Pt."), type="NX_INT"))
    #  UB info

    nxsample.add_dataset(
        key="orientation_matrix",
        ds=NXdataset(ds=metadata.get("ubmatrix"), type="NX_FLOAT", EX_required="true", units="NX_DIMENSIONLESS"),
    )
    # nxsample.add_dataset(key="ub_conf", ds=NXdataset(ds=metadata["ubconf"].split(".")[0], type="NX_CHAR"))
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

    if sample_config_params is not None:
        for param, value in sample_config_params.items():
            match param:
                case "a" | "b" | "c" | "alpha" | "beta" | "gamma":
                    pass
                case "shape":
                    nxsample.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_CHAR"))
                case "width" | "height" | "depth" | "plane_normal" | "in_plane_ref":
                    nxsample.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_FLOAT", units="cm"))
                case "mosaic_h" | "mosaic_v":
                    nxsample.add_dataset(key=param, ds=NXdataset(ds=value, type="NX_ANGLE", units="minutes of arc"))
                case "ub_matrix":
                    nxsample.add_dataset(
                        key="orientation_matrix", ds=NXdataset(ds=np.reshape(np.array(value), (3, 3)), type="NX_FLOAT")
                    )
                case _:
                    nxsample.add_dataset(key=param, ds=NXdataset(ds=value))

    # ---------------------------------------- data ---------------------------------------------
    nexus_keywork_conversion_dict = {"h": "qh", "k": "qk", "l": "ql", "e": "en"}
    def_x = metadata.get("def_x")
    def_y = metadata.get("def_y")
    if def_x in nexus_keywork_conversion_dict:
        def_x = nexus_keywork_conversion_dict[def_x]

    nxdata = NXentry(NX_class="NXdata", EX_required="true", signal=def_y, axes=def_x)
    # ---------------------------------------- scan ---------------------------------------------

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
        user=nxuser,
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

    # Ensure trailing separator is not required
    datafiles_dir = os.path.join(path_to_spice_folder, "Datafiles")

    # Get all files in the Datafiles directory
    scan_list = os.listdir(datafiles_dir)

    # Set filter based on whether a scan number is given
    if scan_num is None:
        filter_keyword = ".dat"
    else:
        filter_keyword = f"scan{scan_num:04}.dat"

    # Filter and join the full paths
    scan_list = [os.path.join(datafiles_dir, scan) for scan in scan_list if scan.endswith(filter_keyword)]
    scan_list.sort()

    # get IPTS number and instrument string
    first_scan = scan_list[0]
    _, _, headers, _, _ = read_spice_datafile(first_scan)
    ipts = headers.get("proposal", "UNKNOWN")

    spice_file_name = os.path.basename(first_scan)  # safe way to get filename

    try:
        instrument_str, exp_num, _ = spice_file_name.split("_")
    except ValueError:
        raise ValueError(f"Unexpected filename format: {spice_file_name}")

    dataset_name = f"IPTS{ipts}_{instrument_str}_{exp_num}"

    nexus_dict = {}
    for path_to_scan_file in scan_list:
        file_name = Path(path_to_scan_file).name  # 'INSTR_expXXXX_scanXXXX.dat'

        try:
            # Extract 'scanXXXX' part (3rd item)
            scan_part = file_name.split("_")[2]  # 'scanXXXX.dat'
            scan_name = scan_part.split(".")[0]  # 'scanXXXX'
        except IndexError:
            raise ValueError(f"Unexpected filename format: {file_name}")

        nxentry = spice_scan_to_nxdict(
            path_to_scan_file,
            path_to_instrument_json,
            path_to_sample_json,
        )
        nxentry["attrs"].update({"dataset_name": dataset_name})
        nexus_dict.update({scan_name: nxentry})

    return nexus_dict
