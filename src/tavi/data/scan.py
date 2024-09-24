# -*- coding: utf-8 -*-
from dataclasses import dataclass
from typing import Optional

import numpy as np

from tavi.data.nxentry import NexusEntry


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

    sample_name: str = ""
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

    sense: str = "+-+"
    collimation: tuple[float, float, float, float] = (60.0, 60.0, 60.0, 60.0)


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
    ei: Optional[tuple[float, ...]] = None
    # analyzer
    ef: Optional[tuple[float, ...]] = None
    # goiometer motor angles
    s1: Optional[tuple[float, ...]] = None
    s2: Optional[tuple[float, ...]] = None
    sgl: Optional[tuple[float, ...]] = None
    sgu: Optional[tuple[float, ...]] = None
    stl: Optional[tuple[float, ...]] = None
    stu: Optional[tuple[float, ...]] = None
    chi: Optional[tuple[float, ...]] = None
    phi: Optional[tuple[float, ...]] = None
    # Q-E space
    q: Optional[tuple[float, ...]] = None
    qh: Optional[tuple[float, ...]] = None
    qk: Optional[tuple[float, ...]] = None
    ql: Optional[tuple[float, ...]] = None
    en: Optional[tuple[float, ...]] = None


class Scan(NexusEntry):
    """
    Manage a single measued scan

    Attributes:
        scan_info (ScanInfo):
        sample_ub_info (SampleUBInfo):
        instrument_info (InstrumentInfo):
        data (ScanData): dictionary contains lists of scan data

    Methods:
        load_scan
        generate_curve
        plot_curve

    """

    def __init__(self, scan_info, sample_ub_info, instrument_info, data) -> None:
        self.scan_info: ScanInfo = scan_info
        self.sample_ub_info: SampleUBInfo = sample_ub_info
        self.instrument_info: InstrumentInfo = instrument_info
        self.data: ScanData = data
