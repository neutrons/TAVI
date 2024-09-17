# -*- coding: utf-8 -*-
from dataclasses import dataclass
from typing import Optional

import h5py
import numpy as np


@dataclass
class ScanInfo:
    """Metadata containing scan information"""

    scan_num: Optional[int] = None
    start_time: Optional[str] = None
    end_time: Optional[str] = None
    scan_title: str = ""
    # preset_type: str = "normal"
    preset_channel: str = "time"
    preset_value: float = 1.0
    def_y: str = "detector"
    def_x: str = "s1"


@dataclass
class SampleUBInfo:
    """Metadata about sample and UB matrix"""

    sample_name: Optional[str] = None
    lattice_constants: tuple = (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)
    ub_matrix: Optional[np.ndarray] = None
    # ub_mode: int = 0  # mode for UB determination in SPICE
    # angle_mode: int = 0  # mode for goni angle calculation in SPICE
    plane_normal: Optional[np.ndarray] = None
    # in_plane_ref: Optional[np.ndarray] = None
    ubconf: Optional[str] = None  # path to UB configration file


# TODO
@dataclass
class InstrumentInfo:
    """Metadata about instrument configuration"""

    instrument_name: str = ""
    # monochromator:
    # analyzer:

    sense: int = 1
    collimation: tuple[float, float, float, float] = (60, 60, 60, 60)
    # TODO vertical collimation??


@dataclass
class ScanData:
    """Data points in a measured scan"""

    pt: Optional[tuple[int, ...]] = None
    detector: Optional[tuple[int, ...]] = None
    # monitor
    time: Optional[tuple[float, ...]] = None
    monitor: Optional[tuple[int, ...]] = None
    mcu: Optional[tuple[float, ...]] = None
    # monochromator
    m1: Optional[tuple[float, ...]] = None
    m2: Optional[tuple[float, ...]] = None
    ei: Optional[tuple[float, ...]] = None
    focal_length: Optional[tuple[float, ...]] = None
    mfocus: Optional[tuple[float, ...]] = None
    marc: Optional[tuple[float, ...]] = None
    mtrans: Optional[tuple[float, ...]] = None
    # analyzer
    ef: Optional[tuple[float, ...]] = None
    a1: Optional[tuple[float, ...]] = None
    a2: Optional[tuple[float, ...]] = None
    afocus: Optional[tuple[float, ...]] = None
    # ctax double-focused analyzer
    # qm: Optional[tuple] = None  # 8 blades
    # xm: Optional[tuple] = None  # 8 blades
    # goiometer motor angles
    s1: Optional[tuple[float, ...]] = None
    s2: Optional[tuple[float, ...]] = None
    sgl: Optional[tuple[float, ...]] = None
    sgu: Optional[tuple[float, ...]] = None
    stl: Optional[tuple[float, ...]] = None
    stu: Optional[tuple[float, ...]] = None
    chi: Optional[tuple[float, ...]] = None
    phi: Optional[tuple[float, ...]] = None
    # slits
    bat: Optional[tuple[float, ...]] = None
    bab: Optional[tuple[float, ...]] = None
    bal: Optional[tuple[float, ...]] = None
    bar: Optional[tuple[float, ...]] = None
    bbt: Optional[tuple[float, ...]] = None
    bbb: Optional[tuple[float, ...]] = None
    bbl: Optional[tuple[float, ...]] = None
    bbr: Optional[tuple[float, ...]] = None
    slita_bt: Optional[tuple[float, ...]] = None
    slita_tp: Optional[tuple[float, ...]] = None
    slita_lf: Optional[tuple[float, ...]] = None
    slita_rt: Optional[tuple[float, ...]] = None
    slitb_bt: Optional[tuple[float, ...]] = None
    slitb_tp: Optional[tuple[float, ...]] = None
    slitb_lf: Optional[tuple[float, ...]] = None
    slitb_rt: Optional[tuple[float, ...]] = None
    slit_pre_bt: Optional[tuple[float, ...]] = None
    slit_pre_tp: Optional[tuple[float, ...]] = None
    slit_pre_lf: Optional[tuple[float, ...]] = None
    slit_pre_rt: Optional[tuple[float, ...]] = None
    # Q-E space
    q: Optional[tuple[float, ...]] = None
    qh: Optional[tuple[float, ...]] = None
    qk: Optional[tuple[float, ...]] = None
    ql: Optional[tuple[float, ...]] = None
    en: Optional[tuple[float, ...]] = None
    # temperature
    temp: Optional[tuple[float, ...]] = None
    temp_a: Optional[tuple[float, ...]] = None
    temp_2: Optional[tuple[float, ...]] = None
    coldtip: Optional[tuple[float, ...]] = None
    tsample: Optional[tuple[float, ...]] = None
    sample: Optional[tuple[float, ...]] = None
    vti: Optional[tuple[float, ...]] = None
    dr_tsample: Optional[tuple[float, ...]] = None
    dr_temp: Optional[tuple[float, ...]] = None
    lt: Optional[tuple[float, ...]] = None
    ht: Optional[tuple[float, ...]] = None
    sorb_temp: Optional[tuple[float, ...]] = None
    sorb: Optional[tuple[float, ...]] = None
    sample_ht: Optional[tuple[float, ...]] = None
    # field
    persistent_field: Optional[tuple[float, ...]] = None

    def _visit_func(self, name, node):
        if isinstance(node, h5py.Dataset):
            attr_name = name.split("/")[-1]
            if attr_name not in (
                "sense",
                "type",
                "Pt.",
                "name",
                "orientation_matrix",
                "plane_normal",
                "ub_conf",
                "unit_cell",
            ):
                setattr(self, attr_name, node[...])

    def get_data_from(self, nexus_entry):
        nexus_entry["instrument/monochromator"].visititems(self._visit_func)
        nexus_entry["instrument/analyser"].visititems(self._visit_func)
        nexus_entry["instrument/slit"].visititems(self._visit_func)
        nexus_entry["sample"].visititems(self._visit_func)
        # Pt.
        self.pt = nexus_entry.get("sample/Pt.")[...]
        # detector
        self.detector = nexus_entry.get("instrument/detector/data")[...]
        # monitor
        self.time = nexus_entry.get("monitor/time")[...]
        self.monitor = nexus_entry.get("monitor/monitor")[...]
        self.mcu = nexus_entry.get("monitor/mcu")[...]
        return self


def nexus_entry_to_scan(nexus_entry):

    def get_string(entry_string):
        metadata = nexus_entry.get(entry_string)
        if metadata is not None:
            return str(metadata.asstr()[...])

    def get_data(entry_string):
        data = nexus_entry.get(entry_string)
        if data is not None:
            return data[...]

    scan_index = nexus_entry.name.split("/")[-1]  # e.g. "scan0001"

    dataset_name = nexus_entry.attrs["dataset_name"]  # e.g. "IPTS32124_CG4C_exp0424"
    scan_info = ScanInfo(
        scan_num=int(scan_index[4:]),
        start_time=get_string("start_time"),
        end_time=get_string("end_time"),
        scan_title=get_string("title"),
        # preset_type="normal",
        preset_channel=get_string("monitor/mode"),
        preset_value=float(nexus_entry["monitor/preset"][...]),
        def_y=nexus_entry["data"].attrs["signal"],
        def_x=nexus_entry["data"].attrs["axes"],
    )
    sample_ub_info = SampleUBInfo(
        sample_name=get_string("sample/name"),
        lattice_constants=np.array(get_data("sample/unit_cell")),
        ub_matrix=np.array(get_data("sample/orientation_matrix")).reshape(3, 3),
        # ub_mode=0,
        # angle_mode=0,
        plane_normal=get_data("sample/plane_normal"),
        # in_plane_ref=None,
        ubconf=get_string("sample/ub_conf"),
    )
    instrument_info = InstrumentInfo()

    scan_data = ScanData().get_data_from(nexus_entry)

    return (dataset_name, scan_info, sample_ub_info, instrument_info, scan_data)


# TODO
def spice_entry_to_scan(spice_entry):
    scan_info = ScanInfo()
    sample_ub_info = SampleUBInfo()
    instrument_info = InstrumentInfo()
    scan_data = ScanData()
    return (scan_info, sample_ub_info, instrument_info, scan_data)
