# -*- coding: utf-8 -*-
from typing import Literal, Optional

import h5py
import matplotlib.pyplot as plt
import numpy as np

from tavi.data.plotter import Plot1D
from tavi.data.scan_reader import (
    InstrumentInfo,
    SampleUBInfo,
    ScanData,
    ScanInfo,
    nexus_entry_to_scan,
    spice_entry_to_scan,
)


class Scan(object):
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

    @classmethod
    def from_nexus_entry(cls, nexus_entry):
        """Load a single scan from a NeXus entry"""
        (dataset_name, *scan_data) = nexus_entry_to_scan(nexus_entry)
        (scan_info, sample_ub_info, instrument_info, data) = scan_data
        return dataset_name, cls(scan_info, sample_ub_info, instrument_info, data)

    @classmethod
    def from_nexus_file(cls, nexus_file_name: str, scan_num: Optional[int] = None):
        """Load a signle scan from NeXus file.

        Note:
            If scan_name is None, nexus_file_name must ends with
            scan number e.g. scan0001.h5
        """

        if scan_num is None:
            scan_name = (nexus_file_name.split("/")[-1]).split(".")[0]
            if scan_name[0:4] != "scan":
                raise ValueError("Unrecogonized NeXus file name. Needs to end with scan number such as scan0001.h5")
        else:
            scan_name = f"scan{scan_num:04}"

        with h5py.File(nexus_file_name, "r") as nexus_file:
            (dataset_name, *scan_data) = nexus_entry_to_scan(nexus_file[scan_name])

        (scan_info, sample_ub_info, instrument_info, data) = scan_data
        return dataset_name, cls(scan_info, sample_ub_info, instrument_info, data)

    # TODO not yet implemented
    @classmethod
    def from_spice(cls, spice_path):
        spice_entry = spice_path
        (dataset_name, *scan_data) = spice_entry_to_scan(spice_entry)
        scan_info, sample_ub_info, instrument_info, data = scan_data
        return dataset_name, cls(scan_info, sample_ub_info, instrument_info, data)

    def get_scan_info(self):
        """Return scan_info in metadata."""
        return self.scan_info

    def get_sample_ub_info(self):
        """Return sample_UB_info in metadata."""
        return self.sample_ub_info

    def get_instrument_info(self):
        """Return instrument_info in metadata."""
        return self.instrument_info

    # def save_metadata(self, metadata_entry):
    #     """Save metadata_entry into file
    #     pass

    def get_data_entry(self, entry_name):
        """Return data entry based on entry_name"""
        return self.data[entry_name]

    def _make_labels(
        self,
        plot1d: Plot1D,
        x_str: str,
        y_str: str,
        norm_channel: Literal["time", "monitor", "mcu", None],
        norm_val: float,
    ) -> Plot1D:
        """generate labels and title"""
        if norm_channel is not None:
            if norm_channel == "time":
                norm_channel_str = "seconds"
            else:
                norm_channel_str = norm_channel
            if norm_val == 1:
                plot1d.ylabel = y_str + "/ " + norm_channel_str
            else:
                plot1d.ylabel = y_str + f" / {norm_val} " + norm_channel_str
        else:
            preset_val = self.scan_info.preset_value
            plot1d.ylabel = y_str + f" / {preset_val} " + self.scan_info.preset_channel

        plot1d.xlabel = x_str
        plot1d.label = "scan " + str(self.scan_info.scan_num)
        plot1d.title = plot1d.label + ": " + self.scan_info.scan_title

        return plot1d

    def _rebin_tol(
        self,
        x_raw: np.ndarray,
        y_raw: np.ndarray,
        y_str: str,
        rebin_step: float,
        norm_channel: Literal["time", "monitor", "mcu", None],
        norm_val: float,
    ):
        x_grid = np.arange(np.min(x_raw) + rebin_step / 2, np.max(x_raw) + rebin_step / 2, rebin_step)
        x = np.zeros_like(x_grid)
        y = np.zeros_like(x_grid)
        counts = np.zeros_like(x_grid)
        weights = np.zeros_like(x_grid)
        yerr = None

        if norm_channel is None:  # rebin, no renorm
            weight_channel = self.scan_info.preset_channel
            weight = getattr(self.data, weight_channel)
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
        norm = getattr(self.data, norm_channel)
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
        rebin_step: float,
        norm_channel: Literal["time", "monitor", "mcu", None],
        norm_val: float,
    ):
        x = np.arange(
            np.min(x_raw) + rebin_step / 2,
            np.max(x_raw) + rebin_step / 2,
            rebin_step,
        )
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
        norm = getattr(self.data, norm_channel)
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
        rebin_step: float = 0.0,
    ) -> Plot1D:
        """Generate a curve from a single scan to plot,
        with the options to normalize the y-axis and rebin x-axis.
        """

        if x_str is None:
            x_str = self.scan_info.def_x

        if y_str is None:
            y_str = self.scan_info.def_y

        x_raw = getattr(self.data, x_str)
        y_raw = getattr(self.data, y_str)

        yerr = None

        if rebin_type is None:  # no rebin
            x = x_raw
            y = y_raw
            # errror bars for detector only
            if "detector" in y_str:
                yerr = np.sqrt(y)
            # normalize y-axis without rebining along x-axis
            if norm_channel is not None:
                norm = getattr(self.data, norm_channel) / norm_val
                y = y / norm
                if yerr is not None:
                    yerr = yerr / norm

            plot1d = Plot1D(x=x, y=y, yerr=yerr)
            plot1d = self._make_labels(plot1d, x_str, y_str, norm_channel, norm_val)

            return plot1d

        if not rebin_step > 0:
            raise ValueError("Rebin step needs to be greater than zero.")

        match rebin_type:
            case "tol":  # x weighted by normalization channel
                x, y, yerr = self._rebin_tol(x_raw, y_raw, y_str, rebin_step, norm_channel, norm_val)
                plot1d = Plot1D(x=x, y=y, yerr=yerr)
                plot1d = self._make_labels(plot1d, x_str, y_str, norm_channel, norm_val)
                return plot1d
            case "grid":
                x, y, yerr = self._rebin_grid(x_raw, y_raw, y_str, rebin_step, norm_channel, norm_val)
                plot1d = Plot1D(x=x, y=y, yerr=yerr)
                plot1d = self._make_labels(plot1d, x_str, y_str, norm_channel, norm_val)
                return plot1d
            case _:
                raise ValueError('Unrecogonized rebin type. Needs to be "tol" or "grid".')

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
