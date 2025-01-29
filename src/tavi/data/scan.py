# -*- coding: utf-8 -*-
from dataclasses import dataclass
from typing import Literal, Optional, Union

import matplotlib.pyplot as plt
import numpy as np

from tavi.data.nexus_entry import NexusEntry
from tavi.data.scan_data import ScanData1D
from tavi.plotter import Plot1D
from tavi.sample.xtal import Xtal
from tavi.utilities import spice_to_mantid


@dataclass
class ScanInfo:
    """Metadata containing scan information"""

    exp_id: str
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
        name (str): scan name
        _nexus_dict (NexusEntry)
        data (dict): dictionary contains lists of scan data

    Methods:
        load_scan
        generate_curve
        plot_curve

    """

    def __init__(self, name: str, nexus_dict: NexusEntry) -> None:

        self.name: str = name
        self._nexus_dict: NexusEntry = nexus_dict
        self.data: dict = self.get_data_columns()

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
        ((name, nexus_dict),) = NexusEntry.from_nexus(path_to_nexus, scan_num).items()
        return cls(name, nexus_dict)

    @property
    def scan_info(self):
        preset_type = self._nexus_dict.get("preset_type")

        scan_info = ScanInfo(
            exp_id=self._nexus_dict["attrs"].get("dataset_name"),
            scan_num=int(self.name[-4:]),
            scan_title=self._nexus_dict.get("title"),
            def_y=self._nexus_dict.get("data", ATTRS=True)["signal"],
            def_x=self._nexus_dict.get("data", ATTRS=True)["axes"],
            start_time=self._nexus_dict.get("start_time"),
            end_time=self._nexus_dict.get("end_time"),
            preset_type=preset_type,
            preset_channel=(
                self._nexus_dict.get("monitor/mode")
                if preset_type == "normal"
                else self._nexus_dict.get("monitor/mode_1")
            ),
            preset_value=(
                self._nexus_dict.get("monitor/preset")
                if preset_type == "normal"
                else self._nexus_dict.get("monitor/preset_1")
            ),
        )
        return scan_info

    @property
    def sample_ub_info(self):
        sample_type = self._nexus_dict.get("sample/type")
        ub_matrix = self._nexus_dict.get("sample/orientation_matrix").reshape(3, 3)
        lattice_constants = self._nexus_dict.get("sample/unit_cell")
        if sample_type == "crystal" and (ub_matrix is not None):
            xtal = Xtal(lattice_constants)
            (u, v) = xtal.ub_matrix_to_uv(spice_to_mantid(ub_matrix))
        else:
            (u, v) = (None, None)

        sample_ub_info = SampleUBInfo(
            sample_name=self._nexus_dict.get("sample/name"),
            lattice_constants=lattice_constants,
            ub_matrix=ub_matrix,
            type=sample_type,
            u=u,
            v=v,
            # ub_mode: int = 0  # mode for UB determination in SPICE
            # angle_mode: int = 0  # mode for goni angle calculation in SPICE
            plane_normal=self._nexus_dict.get("sample/plane_normal"),
            # in_plane_ref: Optional[np.ndarray] = None
            ubconf=None,  # path to UB configration file
        )

        return sample_ub_info

    @property
    def instrument_info(self):
        instru_info = InstrumentInfo(
            monochromator=self._nexus_dict.get("monochromator/type"),
            analyzer=self._nexus_dict.get("analyser/type"),
            sense=(
                self._nexus_dict.get("monochromator/sense")
                + self._nexus_dict.get("goniometer/sense")
                + self._nexus_dict.get("analyser/sense")
            ),
            collimation=self._nexus_dict.get("divergence_x"),
        )
        return instru_info

    def get_data_columns(self) -> dict[str, np.ndarray]:
        """Get scan data as a dictionary"""
        data_dict = {}
        names = self._nexus_dict.get_dataset_names()
        for name in names:
            data_dict.update({name: self._nexus_dict.get(name)})
        return data_dict

    @staticmethod
    def validate_rebin_params(rebin_params: float | int | tuple) -> tuple:
        if isinstance(rebin_params, tuple):
            if len(rebin_params) != 3:
                raise ValueError("Rebin parameters should have the form (min, max, step)")
            rebin_min, rebin_max, rebin_step = rebin_params
            if (rebin_min >= rebin_max) or (rebin_step < 0):
                raise ValueError(f"Nonsensical rebin parameters {rebin_params}")

        elif isinstance(rebin_params, float | int):
            if rebin_params < 0:
                raise ValueError("Rebin step needs to be greater than zero.")
            rebin_params = (None, None, float(rebin_params))

        else:
            raise ValueError(f"Unrecogonized rebin parameters {rebin_params}")
        return rebin_params

    def _get_del_q(self):
        """Calculate del_q for aeither a s1 scan or a th2th scan"""
        qs = self.data["q"]
        q_diff = np.max(qs) - np.min(qs)
        mid_idx = int((len(qs) - 1) / 2)
        if q_diff > 1.1e-4:  # q changing, must be a th2th scan
            return qs - qs[mid_idx]
        else:  # q not changing, must be a s1 scan
            q_abs = np.mean(qs)
            s1s = self.data["s1"]
            return np.deg2rad(s1s - s1s[mid_idx]) * q_abs

    def get_data(
        self,
        axes: tuple[Optional[str], Optional[str]] = (None, None),
        norm_to: Optional[tuple[float, Literal["time", "monitor", "mcu"]]] = None,
        **rebin_params_dict: Union[float, tuple],
    ) -> ScanData1D:
        """Generate a curve from a single scan to plot, with the options
        to normalize the y-axis and rebin x-axis.

        Args:
            axes (x_str, y_str): x-axis and y-axis variables
            norm_to (norm_val (float), norm_channel(str)): value and channel for normalization
                norm_channel should be "time", "monitor" or"mcu".
            rebin_type (str | None): "tol" or "grid"
            rebin_params (float | tuple(flot, float, float)): take as step size if a numer is given,
                take as (min, max, step) if a tuple of size 3 is given
        """

        x_str, y_str = axes
        x_str = self.scan_info.def_x if x_str is None else x_str
        y_str = self.scan_info.def_y if y_str is None else y_str

        if x_str == "del_q":
            scan_data_1d = ScanData1D(x=self._get_del_q(), y=self.data[y_str])
        else:
            scan_data_1d = ScanData1D(x=self.data[x_str], y=self.data[y_str])

        label = "scan " + str(self.scan_info.scan_num)
        title = f"{label}: {self.scan_info.scan_title}"

        for rebin_type in ["grid", "tol"]:
            rebin_params = rebin_params_dict.get(rebin_type)
            if rebin_params is not None:
                break

        if not rebin_params:  # no rebin
            if norm_to is not None:  # normalize y-axis without rebining along x-axis
                norm_val, norm_channel = norm_to
                scan_data_1d.renorm(norm_col=self.data[norm_channel] / norm_val)
            else:  # equivalent to normalizing to preset
                norm_to = (self.scan_info.preset_value, self.scan_info.preset_channel)

            scan_data_1d.make_labels((x_str, y_str), norm_to, label, title)
            return scan_data_1d

        # Rebin, first validate rebin params
        rebin_params_tuple = Scan.validate_rebin_params(rebin_params)

        match rebin_type:
            case "tol":
                if norm_to is None:  # x weighted by preset channel
                    weight_channel = self.scan_info.preset_channel
                    scan_data_1d.rebin_tol(rebin_params_tuple, weight_col=self.data[weight_channel])
                    norm_to = (self.scan_info.preset_value, self.scan_info.preset_channel)
                else:  # x weighted by normalization channel
                    norm_val, norm_channel = norm_to
                    scan_data_1d.rebin_tol_renorm(
                        rebin_params_tuple,
                        norm_col=self.data[norm_channel],
                        norm_val=norm_val,
                    )
            case "grid":
                if norm_to is None:
                    scan_data_1d.rebin_grid(rebin_params_tuple)
                    norm_to = (self.scan_info.preset_value, self.scan_info.preset_channel)
                else:
                    norm_val, norm_channel = norm_to
                    scan_data_1d.rebin_grid_renorm(
                        rebin_params_tuple,
                        norm_col=self.data[norm_channel],
                        norm_val=norm_val,
                    )
            case _:
                raise ValueError('Unrecogonized rebin_params_dict keyword. Needs to be "tol" or "grid".')

        scan_data_1d.make_labels((x_str, y_str), norm_to, label, title)
        return scan_data_1d

    def get_coherent_fwhm(self, tas, projection: tuple = ((1, 0, 0), (0, 1, 0), (0, 0, 1)), axis: int = 0):

        qh = self.data.get("qh")
        qk = self.data.get("qk")
        ql = self.data.get("ql")
        ei = self.data.get("ei")
        ef = self.data.get("ef")
        rez_params_list = []

        for i in range(len(qh)):
            rez = tas.rez(hkl_list=(qh[i], qk[i], ql[i]), ei=ei[i], ef=ef[i], R0=True, projection=projection)
            coh_fwhm = rez.coh_fwhms(axis=axis)
            prefactor = rez.r0
            rez_params_list.append((prefactor, coh_fwhm))
        return tuple(rez_params_list)

    def plot(
        self,
        axes: tuple[Optional[str], Optional[str]] = (None, None),
        norm_to: Optional[tuple[float, Literal["time", "monitor", "mcu"]]] = None,
        **rebin_params_dict: Union[float, tuple],
    ):
        """Plot a 1D curve gnerated from a singal scan in a new window"""

        scan_data_1d = self.get_data(axes, norm_to, **rebin_params_dict)

        fig, ax = plt.subplots()
        plot1d = Plot1D()
        plot1d.add_scan(scan_data_1d, c="C0", fmt="o")
        plot1d.plot(ax)
        fig.show()
