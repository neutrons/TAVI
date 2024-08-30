# -*- coding: utf-8 -*-
from typing import NamedTuple, Optional

import numpy as np


class ScanInfo(NamedTuple):
    """Metadata containing scan information"""

    scan_num: Optional[int] = None
    time: Optional[str] = None
    scan_title: str = ""
    preset_type: str = "normal"
    preset_channel: str = "time"
    preset_value: float = 1.0
    def_y: str = "detector"
    def_x: str = "s1"


class SampleUBInfo(NamedTuple):
    """Metadata about sample and UB matrix"""

    sample_name: Optional[str] = None
    lattice_constants: tuple[
        float,
        float,
        float,
        float,
        float,
        float,
    ] = (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)
    ub_matrix: Optional[np.ndarray] = None
    # mode: int = 0 # mode for UB determination in SPICE
    plane_normal: Optional[np.ndarray] = None
    in_plane_ref: Optional[np.ndarray] = None
    ubconf: Optional[str] = None  # path to UB configration file


class InstrumentInfo(NamedTuple):
    """Metadata about instrument configuration"""

    instrument_name: str = ""
    # monochromator:
    # analyzer:
    # TODO
    sense: int = 1
    collimation: tuple[float, float, float, float] = (60, 60, 60, 60)
    # TODO vertical collimation??


class ScanData(NamedTuple):
    """Data points in a measured scan"""

    # Pt.
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
    qm: Optional[
        tuple[
            tuple[float, ...],
            tuple[float, ...],
            tuple[float, ...],
            tuple[float, ...],
            tuple[float, ...],
            tuple[float, ...],
            tuple[float, ...],
            tuple[float, ...],
        ]
    ] = None
    xm: Optional[
        tuple[
            tuple[float, ...],
            tuple[float, ...],
            tuple[float, ...],
            tuple[float, ...],
            tuple[float, ...],
            tuple[float, ...],
            tuple[float, ...],
            tuple[float, ...],
        ]
    ] = None
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
    slit_pre_bt: Optional[tuple[float, ...]] = None
    slit_pre_tp: Optional[tuple[float, ...]] = None
    slit_pre_lf: Optional[tuple[float, ...]] = None
    slit_pre_rt: Optional[tuple[float, ...]] = None
    slit_aft_bt: Optional[tuple[float, ...]] = None
    slit_aft_tp: Optional[tuple[float, ...]] = None
    slit_aft_lf: Optional[tuple[float, ...]] = None
    slit_aft_rt: Optional[tuple[float, ...]] = None
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


def dataset_to_string(ds):
    return str(ds.asstr()[...])


def unpack_nexus(nexus_entry):
    scan_name = nexus_entry.filename.split("/")[-1]  # e.g. "scan0001"
    scan_info = ScanInfo(
        scan_num=int(scan_name[4:-3]),
        time=dataset_to_string(nexus_entry["start_time"]),
        scan_title=dataset_to_string(nexus_entry["title"]),
        preset_type="normal",
        preset_channel=dataset_to_string(nexus_entry["monitor/mode"]),
        preset_value=float(nexus_entry["monitor/preset"][...]),
        def_y=nexus_entry["data"].attrs["signal"],
        def_x=nexus_entry["data"].attrs["axes"],
    )
    sample_ub_info = SampleUBInfo()
    instrument_info = InstrumentInfo()
    scan_data = ScanData()

    return (scan_info, sample_ub_info, instrument_info, scan_data)
