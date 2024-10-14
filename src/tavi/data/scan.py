# -*- coding: utf-8 -*-
from dataclasses import dataclass
from typing import Literal, Optional, Union

import matplotlib.pyplot as plt
import numpy as np

from tavi.data.nxentry import NexusEntry
from tavi.data.plotter import Plot1D
from tavi.sample.xtal import Xtal
from tavi.utilities import spice_to_mantid


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
        self.data: dict = self.get_data()

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
            scan_title=self._nexus_dict.get("title"),
            def_y=self._nexus_dict.get("data", ATTRS=True)["signal"],
            def_x=self._nexus_dict.get("data", ATTRS=True)["axes"],
            start_time=self._nexus_dict.get("start_time"),
            end_time=self._nexus_dict.get("end_time"),
            preset_type="normal",
            preset_channel=self._nexus_dict.get("monitor/mode"),
            preset_value=self._nexus_dict.get("monitor/preset"),
        )
        return scan_info

    @property
    def sample_ub_info(self):
        sample_type = self._nexus_dict.get("sample/type")
        ub_matrix = self._nexus_dict.get("sample/orientation_matrix")
        lattice_constants = self._nexus_dict.get("sample/unit_cell")
        if sample_type == "crystal" and (ub_matrix is not None):
            xtal = Xtal(lattice_constants)
            (u, v) = xtal.ub_matrix_to_uv(spice_to_mantid(ub_matrix.reshape(3, 3)))
        else:
            u = None
            v = None

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
                + self._nexus_dict.get("sample/sense")
                + self._nexus_dict.get("analyser/sense")
            ),
            collimation=self._nexus_dict.get("divergence_x"),
        )
        return instru_info

    def get_data(self) -> dict:
        """Get scan data as a dictionary"""
        data_dict = {}
        names = self._nexus_dict.get_dataset_names()
        for name in names:
            data_dict.update({name: self._nexus_dict.get(name)})
        return data_dict

    def _rebin_tol(
        self,
        x_raw: np.ndarray,
        y_raw: np.ndarray,
        y_str: str,
        rebin_params: tuple,
        norm_channel: Literal["time", "monitor", "mcu", None],
        norm_val: float,
    ):
        """Rebin with tolerance"""

        rebin_min, rebin_max, rebin_step = rebin_params
        if rebin_min is None:
            rebin_min = np.min(x_raw)
        if rebin_max is None:
            rebin_max = np.max(x_raw)

        x_grid = np.arange(rebin_min + rebin_step / 2, rebin_max + rebin_step / 2, rebin_step)
        x = np.zeros_like(x_grid)
        y = np.zeros_like(x_grid)
        counts = np.zeros_like(x_grid)
        weights = np.zeros_like(x_grid)
        yerr = None

        if norm_channel is None:  # rebin, no renorm
            weight_channel = self.scan_info.preset_channel
            weight = self.data[weight_channel]
            for i, x0 in enumerate(x_raw):
                idx = np.nanargmax(x_grid + rebin_step / 2 >= x0)
                y[idx] += y_raw[i]
                x[idx] += x_raw[i] * weight[i]
                weights[idx] += weight[i]
                counts[idx] += 1

            # errror bars for detector only
            if "detector" in y_str:
                yerr = np.sqrt(y) / counts
            y = y / counts
            x = x / weights
            return (x, y, yerr)

        # rebin and renorm
        norm = self.data[norm_channel]
        for i, x0 in enumerate(x_raw):
            idx = np.nanargmax(x_grid + rebin_step / 2 >= x0)
            y[idx] += y_raw[i]
            x[idx] += x_raw[i] * norm[i]
            counts[idx] += norm[i]

        # errror bars for detector only
        if "detector" in y_str:
            yerr = np.sqrt(y) / counts * norm_val
        y = y / counts * norm_val
        x = x / counts
        return (x, y, yerr)

    def _rebin_grid(
        self,
        x_raw: np.ndarray,
        y_raw: np.ndarray,
        y_str: str,
        rebin_params: tuple,
        norm_channel: Literal["time", "monitor", "mcu", None],
        norm_val: float,
    ):
        """Rebin with a regular grid"""

        rebin_min, rebin_max, rebin_step = rebin_params
        if rebin_min is None:
            rebin_min = np.min(x_raw)
        if rebin_max is None:
            rebin_max = np.max(x_raw)

        x = np.arange(rebin_min + rebin_step / 2, rebin_max + rebin_step / 2, rebin_step)
        y = np.zeros_like(x)
        cts = np.zeros_like(x)
        yerr = None
        # rebin, no renorm
        if norm_channel is None:
            for i, x0 in enumerate(x_raw):
                idx = np.nanargmax(x + rebin_step / 2 >= x0)
                y[idx] += y_raw[i]
                cts[idx] += 1

            # errror bars for detector only
            if "detector" in y_str:
                yerr = np.sqrt(y) / cts
                y = y / cts
            return (x, y, yerr)

        # rebin and renorm
        norm = self.data[norm_channel]
        for i, x0 in enumerate(x_raw):
            idx = np.nanargmax(x + rebin_step / 2 >= x0)
            y[idx] += y_raw[i]
            cts[idx] += norm[i]

        # errror bars for detector only
        if "detector" in y_str:
            yerr = np.sqrt(y) / cts * norm_val
        y = y / cts * norm_val
        return (x, y, yerr)

    def generate_curve(
        self,
        x_str: Optional[str] = None,
        y_str: Optional[str] = None,
        norm_channel: Literal["time", "monitor", "mcu", None] = None,
        norm_val: float = 1.0,
        rebin_type: Literal["tol", "grid", None] = None,
        rebin_params: Union[float, tuple] = 0.0,
    ) -> Plot1D:
        """Generate a curve from a single scan to plot, with the options
        to normalize the y-axis and rebin x-axis.

        Args:
            x_str (str): x-axis variable
            y_str (str): y-axis variable
            norm_channel (str | None): choose from "time", "monitor", "mcu"
            norm_val (float): value to normalized to
            rebin_type (str | None): "tol" or "grid"
            rebin_params (float | tuple(flot, float, float)): take as step size if a numer is given,
                take as (min, max, step) if a tuple of size 3 is given
        """

        if x_str is None:
            x_str = self.scan_info.def_x

        if y_str is None:
            y_str = self.scan_info.def_y

        x_raw = self.data[x_str]
        y_raw = self.data[y_str]

        yerr = None

        if rebin_type is None:  # no rebin
            x = x_raw
            y = y_raw
            # errror bars for detector only
            if "detector" in y_str:
                yerr = np.sqrt(y)
            # normalize y-axis without rebining along x-axis
            if norm_channel is not None:
                norm = self.data[norm_channel] / norm_val
                y = y / norm
                if yerr is not None:
                    yerr = yerr / norm

            plot1d = Plot1D(x=x, y=y, yerr=yerr)
            plot1d.make_labels(x_str, y_str, norm_channel, norm_val, self.scan_info)

            return plot1d

        # validate rebin params
        if isinstance(rebin_params, tuple):
            if len(rebin_params) != 3:
                raise ValueError("Rebin parameters should have the form (min, max, step)")
            rebin_min, rebin_max, rebin_step = rebin_params
            if (rebin_min >= rebin_max) or (rebin_step < 0):
                raise ValueError(f"Nonsensical rebin parameters {rebin_params}")

        elif isinstance(rebin_params, float | int):
            if rebin_params < 0:
                raise ValueError("Rebin step needs to be greater than zero.")
            rebin_params = (None, float(rebin_params), None)
        else:
            raise ValueError(f"Unrecogonized rebin parameters {rebin_params}")

        match rebin_type:
            case "tol":  # x weighted by normalization channel
                x, y, yerr = self._rebin_tol(x_raw, y_raw, y_str, rebin_params, norm_channel, norm_val)
            case "grid":
                x, y, yerr = self._rebin_grid(x_raw, y_raw, y_str, rebin_params, norm_channel, norm_val)
            case _:
                raise ValueError('Unrecogonized rebin type. Needs to be "tol" or "grid".')

        plot1d = Plot1D(x=x, y=y, yerr=yerr)
        plot1d.make_labels(x_str, y_str, norm_channel, norm_val, self.scan_info)
        return plot1d

    def plot_curve(
        self,
        x_str: Optional[str] = None,
        y_str: Optional[str] = None,
        norm_channel: Literal["time", "monitor", "mcu", None] = None,
        norm_val: float = 1.0,
        rebin_type: Literal["tol", "grid", None] = None,
        rebin_step: float = 0.0,
    ):
        """Plot a 1D curve gnerated from a singal scan in a new window"""

        plot1d = self.generate_curve(x_str, y_str, norm_channel, norm_val, rebin_type, rebin_step)

        fig, ax = plt.subplots()
        plot1d.plot_curve(ax)
        fig.show()
