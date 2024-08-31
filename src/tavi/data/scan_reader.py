# -*- coding: utf-8 -*-
from typing import NamedTuple, Optional

import numpy as np


class ScanInfo(NamedTuple):
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


class SampleUBInfo(NamedTuple):
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
class InstrumentInfo(NamedTuple):
    """Metadata about instrument configuration"""

    instrument_name: str = ""
    # monochromator:
    # analyzer:

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
    qm: Optional[tuple] = None  # 8 blades
    xm: Optional[tuple] = None  # 8 blades
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


def get_string(nexus_entry, entry_string):
    try:
        return str(nexus_entry[entry_string].asstr()[...])
    except KeyError:
        return None


def get_data(nexus_entry, entry_string):
    try:
        return nexus_entry[entry_string][...]
    except KeyError:
        return None


def nexus_entry_to_scan(nexus_entry):
    scan_name = nexus_entry.filename.split("/")[-1]  # e.g. "scan0001"
    scan_info = ScanInfo(
        scan_num=int(scan_name[4:-3]),
        start_time=get_string(nexus_entry, "start_time"),
        end_time=get_string(nexus_entry, "end_time"),
        scan_title=get_string(nexus_entry, "title"),
        # preset_type="normal",
        preset_channel=get_string(nexus_entry, "monitor/mode"),
        preset_value=float(nexus_entry["monitor/preset"][...]),
        def_y=nexus_entry["data"].attrs["signal"],
        def_x=nexus_entry["data"].attrs["axes"],
    )
    sample_ub_info = SampleUBInfo(
        sample_name=get_string(nexus_entry, "sample/name"),
        lattice_constants=np.array(get_data(nexus_entry, "sample/unit_cell")),
        ub_matrix=np.array(get_data(nexus_entry, "sample/orientation_matrix")).reshape(3, 3),
        # ub_mode=0,
        # angle_mode=0,
        plane_normal=get_data(nexus_entry, "sample/plane_normal"),
        # in_plane_ref=None,
        ubconf=get_string(nexus_entry, "sample/ub_conf"),
    )
    instrument_info = InstrumentInfo()

    scan_data = ScanData(
        detector=get_data(nexus_entry, "instrument/detector/data"),
        # monitor
        time=get_data(nexus_entry, "monitor/time"),
        monitor=get_data(nexus_entry, "monitor/monitor"),
        mcu=get_data(nexus_entry, "monitor/mcu"),
        # monochromator
        m1=get_data(nexus_entry, "instrument/monochromator/m1"),
        m2=get_data(nexus_entry, "instrument/monochromator/m2"),
        ei=get_data(nexus_entry, "instrument/monochromator/ei"),
        focal_length=get_data(nexus_entry, "instrument/monochromator/focal_length"),
        mfocus=get_data(nexus_entry, "instrument/monochromator/mfocus"),
        marc=get_data(nexus_entry, "instrument/monochromator/marc"),
        mtrans=get_data(nexus_entry, "instrument/monochromator/mtrans"),
        # analyzer
        ef=get_data(nexus_entry, "instrument/analyser/ef"),
        a1=get_data(nexus_entry, "instrument/analyser/a1"),
        a2=get_data(nexus_entry, "instrument/analyser/a2"),
        afocus=get_data(nexus_entry, "instrument/analyser/afocus"),
        # ctax double-focused analyzer
        qm=np.array([get_data(nexus_entry, f"instrument/analyser/qm{i}") for i in range(1, 9, 1)]),
        xm=np.array([get_data(nexus_entry, f"instrument/analyser/xm{i}") for i in range(1, 9, 1)]),
        # goiometer motor angles
        s1=get_data(nexus_entry, "sample/s1"),
        s2=get_data(nexus_entry, "sample/s2"),
        sgl=get_data(nexus_entry, "sample/sgl"),
        sgu=get_data(nexus_entry, "sample/sgu"),
        stl=get_data(nexus_entry, "sample/stl"),
        stu=get_data(nexus_entry, "sample/stu"),
        chi=get_data(nexus_entry, "sample/chi"),
        phi=get_data(nexus_entry, "sample/phi"),
        # slits
        bat=get_data(nexus_entry, "instrument/slit/bat"),
        bab=get_data(nexus_entry, "instrument/slit/bab"),
        bal=get_data(nexus_entry, "instrument/slit/bal"),
        bar=get_data(nexus_entry, "instrument/slit/bar"),
        bbt=get_data(nexus_entry, "instrument/slit/bbt"),
        bbb=get_data(nexus_entry, "instrument/slit/bbb"),
        bbl=get_data(nexus_entry, "instrument/slit/bbl"),
        bbr=get_data(nexus_entry, "instrument/slit/bbr"),
        slita_bt=get_data(nexus_entry, "instrument/slit/slita_bt"),
        slita_tp=get_data(nexus_entry, "instrument/slit/slita_tp"),
        slita_lf=get_data(nexus_entry, "instrument/slit/slita_lf"),
        slita_rt=get_data(nexus_entry, "instrument/slit/slita_rt"),
        slitb_bt=get_data(nexus_entry, "instrument/slit/slitb_bt"),
        slitb_tp=get_data(nexus_entry, "instrument/slit/slitb_tp"),
        slitb_lf=get_data(nexus_entry, "instrument/slit/slitb_lf"),
        slitb_rt=get_data(nexus_entry, "instrument/slit/slitb_rt"),
        slit_pre_bt=get_data(nexus_entry, "instrument/slit/slit_pre_bt"),
        slit_pre_tp=get_data(nexus_entry, "instrument/slit/slit_pre_tp"),
        slit_pre_lf=get_data(nexus_entry, "instrument/slit/slit_pre_lf"),
        slit_pre_rt=get_data(nexus_entry, "instrument/slit/slit_pre_rt"),
        # Q-E space
        q=get_data(nexus_entry, "sample/q"),
        qh=get_data(nexus_entry, "sample/qh"),
        qk=get_data(nexus_entry, "sample/qk"),
        ql=get_data(nexus_entry, "sample/ql"),
        en=get_data(nexus_entry, "sample/en"),
        # temperature
        temp=get_data(nexus_entry, "sample/temp"),
        temp_a=get_data(nexus_entry, "sample/temp_a"),
        temp_2=get_data(nexus_entry, "sample/temp_2"),
        coldtip=get_data(nexus_entry, "sample/coldtip"),
        tsample=get_data(nexus_entry, "sample/tsample"),
        sample=get_data(nexus_entry, "sample/sample"),
        vti=get_data(nexus_entry, "sample/vti"),
        dr_tsample=get_data(nexus_entry, "sample/dr_tsample"),
        dr_temp=get_data(nexus_entry, "sample/dr_temp"),
        lt=get_data(nexus_entry, "sample/lt"),
        ht=get_data(nexus_entry, "sample/ht"),
        sorb_temp=get_data(nexus_entry, "sample/sorb_temp"),
        sorb=get_data(nexus_entry, "sample/sorb"),
        sample_ht=get_data(nexus_entry, "sample/sample_ht"),
        # field
        persistent_field=get_data(nexus_entry, "sample/persistent_field"),
    )

    return (scan_info, sample_ub_info, instrument_info, scan_data)


# TODO
def spice_entry_to_scan(spice_entry):
    scan_info = ScanInfo()
    sample_ub_info = SampleUBInfo()
    instrument_info = InstrumentInfo()
    scan_data = ScanData()
    return (scan_info, sample_ub_info, instrument_info, scan_data)
