# -*- coding: utf-8 -*-
from dataclasses import dataclass
from typing import Optional

import numpy as np

from tavi.data.nxentry import NexusEntry
from tavi.sample.xtal import Xtal


@dataclass
class ScanInfo:
    """Metadata containing scan information"""

    scan_num: int
    scan_title: str = ""
    def_y: str = "detector"
    def_x: str = "s1"
    start_time: Optional[str] = None
    end_time: Optional[str] = None
    preset_type: str = "normal"
    preset_channel: str = "time"
    preset_value: float = 1.0


@dataclass
class SampleUBInfo:
    """Metadata about sample and UB matrix"""

    sample_name: str = ""
    lattice_constants: tuple = (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)
    ub_matrix: Optional[np.ndarray] = None
    type: str = "crystal"
    u: Optional[np.ndarray] = None
    v: Optional[np.ndarray] = None
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
    monochromator: str = "PG002"
    analyzer: str = "Pg002"
    sense: str = "+-+"
    collimation: tuple[float, float, float, float] = (60.0, 60.0, 60.0, 60.0)


class Scan(object):
    """
    Manage a single measued scan

    Attributes:
        data: dictionary contains lists of scan data

    Methods:
        load_scan
        generate_curve
        plot_curve

    """

    def __init__(self, name: str, nexus_dict: NexusEntry) -> None:
        self.name: str = name
        self.nexus_dict: NexusEntry = nexus_dict
        # always using spice convention
        self.SPICE_CONVENTION = True

    @classmethod
    def from_spice(
        cls,
        path_to_spice_folder: str,
        scan_num: int,
        path_to_instrument_json: Optional[str] = None,
        path_to_sample_json: Optional[str] = None,
    ):
        ((name, nexus_dict),) = NexusEntry.from_spice(
            path_to_spice_folder,
            scan_num,
            path_to_instrument_json,
            path_to_sample_json,
        ).items()
        return cls(name, nexus_dict)

    @classmethod
    def from_nexus(cls, path_to_nexus: str, scan_num: Optional[int] = None):
        ((name, nexus_dict),) = NexusEntry.from_nexus(
            path_to_nexus,
            scan_num,
        ).items()
        return cls(name, nexus_dict)

    @property
    def scan_info(self):
        scan_info = ScanInfo(
            scan_num=int(self.name[-4:]),
            scan_title=self.nexus_dict.get("title"),
            def_y=self.nexus_dict.get("data", ATTRS=True)["signal"],
            def_x=self.nexus_dict.get("data", ATTRS=True)["axes"],
            start_time=self.nexus_dict.get("start_time"),
            end_time=self.nexus_dict.get("end_time"),
            preset_type="normal",
            preset_channel=self.nexus_dict.get("monitor/mode"),
            preset_value=self.nexus_dict.get("monitor/preset"),
        )
        return scan_info

    @property
    def sample_ub_info(self):
        sample_type = self.nexus_dict.get("sample/type")
        ub_matrix = self.nexus_dict.get("sample/orientation_matrix")
        lattice_constants = self.nexus_dict.get("sample/unit_cell")
        if sample_type == "crystal" and (ub_matrix is not None):
            xtal = Xtal(lattice_constants)
            (u, v) = xtal.spice_ub_matrix_to_uv(ub_matrix.reshape(3, 3))
        else:
            u = None
            v = None

        sample_ub_info = SampleUBInfo(
            sample_name=self.nexus_dict.get("sample/name"),
            lattice_constants=lattice_constants,
            ub_matrix=ub_matrix,
            type=sample_type,
            u=u,
            v=v,
            # ub_mode: int = 0  # mode for UB determination in SPICE
            # angle_mode: int = 0  # mode for goni angle calculation in SPICE
            plane_normal=self.nexus_dict.get("sample/plane_normal"),
            # in_plane_ref: Optional[np.ndarray] = None
            ubconf=None,  # path to UB configration file)
        )

        return sample_ub_info

    @property
    def instrument_info(self):
        instru_info = InstrumentInfo(
            monochromator=self.nexus_dict.get("monochromator/type"),
            analyzer=self.nexus_dict.get("analyser/type"),
            sense=(
                self.nexus_dict.get("monochromator/sense")
                + self.nexus_dict.get("sample/sense")
                + self.nexus_dict.get("analyser/sense")
            ),
            collimation=self.nexus_dict.get("divergence_x"),
        )
        return instru_info

    @property
    def data(self):
        data_dict = {}
        names = self.nexus_dict.get_dataset_names()
        for name in names:
            data_dict.update({name: self.nexus_dict.get(name)})
        return data_dict
