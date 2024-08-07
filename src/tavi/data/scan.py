import matplotlib.pyplot as plt
import numpy as np
from tavi.data.nexus_reader import nexus_to_dict


class Scan(object):
    """
    Manage a single measued scan

    Attributes:
        scan_info (dict):
        sample_ub_info (dict):
        instrument_info (dict):
        data (dict): dictionary contains lists of scan data

    Methods:
        load_scan
        generate_curve
        plot_curve

    """

    def __init__(self, nexus_entry=None):
        """Initialze an empty scan if nexus entry not provided"""

        self.scan_info = None
        self.sample_ub_info = None
        self.instrument_info = None
        self.data = None

        if nexus_entry is not None:
            self.load_scan(nexus_entry)

    def load_scan(self, nexus_entry):
        """Unpack metadata and data from nexus_entry

        Args:
            nexus_entry:
        """

        scan_info, sample_ub_info, instrument_info, data = nexus_to_dict(nexus_entry)
        self.scan_info = scan_info
        self.sample_ub_info = sample_ub_info
        self.instrument_info = instrument_info
        self.data = data

    def get_scan_info(self):
        """Return scan_info in metadata.

        Returns:
            dict: dictionay of scan_info metadata.
        """
        return self.scan_info

    def get_sample_ub_info(self):
        """Return sample_UB_info in metadata.

        Returns:
            dict: dictionay of sample_UB_info metadata.
        """
        return self.sample_ub_info

    def get_instrument_info(self):
        """Return instrument_info in metadata.

        Returns:
            dict: dictionay of instrument_info metadata.
        """
        return self.instrument_info

    # def save_metadata(self, metadata_entry):
    #     """Save metadata_entry into file

    #     Args:
    #        metadata_entry (dict): {key: value}
    #     """
    #     pass

    def get_data_entry(self, entry_name):
        """Return data entry based on entry_name

        Args:
            entry_name (str): a key of the dictionay in data

        Returns:
            tuple: data entry
        """
        return self.data[entry_name]

    def generate_curve(
        self,
        x_str=None,
        y_str=None,
        norm_channel=None,
        norm_val=1,
        rebin_type=None,
        rebin_step=0,
    ):
        """Generate a curve from a single scan to plot,
        with the options to
            normalize the y-axis and rebin x-axis.

        Args:
            x_str (str): string of x axis
            y_str (str): string of x axis
            norm_channel (str):  None, "time", "monitor" or "mcu"
            norm_val (float):
            rebin_type (str): None, "tol", or "grid"
            rebin_size (float):

        Returns:
            x:
            y:
            xerr: if rebin
            yerr: if "detector"
            xlabel:
            ylabel:
            title: "scan number: scan title"
            label: run number as legend if overplotting multiple curves
        """

        if x_str is None:
            x_str = self.scan_info["def_x"]

        if y_str is None:
            y_str = self.scan_info["def_y"]

        x_raw = self.data[x_str]
        y_raw = self.data[y_str]

        # xerr NOT used
        xerr = None
        yerr = None

        if rebin_type is None:  # no rebin
            x = x_raw
            y = y_raw
            # normalize y-axis without rebining along x-axis
            if norm_channel is not None:
                norm = self.data[norm_channel] / norm_val
                y = y / norm
                if yerr is not None:
                    yerr = yerr / norm
            # errror bars for detector only
            if "det" in y_str:
                yerr = np.sqrt(y)
        else:
            if rebin_step > 0:
                match rebin_type:
                    case "tol":  # x weighted by normalization channel
                        x_grid = np.arange(
                            np.min(x_raw) + rebin_step / 2,
                            np.max(x_raw) + rebin_step / 2,
                            rebin_step,
                        )
                        x = np.zeros_like(x_grid)
                        y = np.zeros_like(x_grid)
                        cts = np.zeros_like(x_grid)
                        wts = np.zeros_like(x_grid)

                        if norm_channel is not None:  # rebin and renorm
                            norm = self.data[norm_channel]
                            for i, x0 in enumerate(x_raw):
                                idx = np.nanargmax(x_grid + rebin_step / 2 >= x0)
                                y[idx] += y_raw[i]
                                x[idx] += x_raw[i] * norm[i]
                                cts[idx] += norm[i]

                            # errror bars for detector only
                            if "det" in y_str:
                                yerr = np.sqrt(y) / cts * norm_val
                            y = y / cts * norm_val
                            x = x / cts

                        else:  # rebin, no renorm
                            weight_channel = self.scan_info["preset_channel"]
                            weight = self.data[weight_channel]
                            for i, x0 in enumerate(x_raw):
                                idx = np.nanargmax(x_grid + rebin_step / 2 >= x0)
                                y[idx] += y_raw[i]
                                x[idx] += x_raw[i] * weight[i]
                                wts[idx] += weight[i]
                                cts[idx] += 1

                            # errror bars for detector only
                            if "det" in y_str:
                                yerr = np.sqrt(y) / cts
                            y = y / cts
                            x = x / wts
                    case "grid":
                        x = np.arange(
                            np.min(x_raw) + rebin_step / 2,
                            np.max(x_raw) + rebin_step / 2,
                            rebin_step,
                        )
                        y = np.zeros_like(x)
                        cts = np.zeros_like(x)

                        if norm_channel is not None:  # rebin and renorm
                            norm = self.data[norm_channel]
                            for i, x0 in enumerate(x_raw):
                                idx = np.nanargmax(x + rebin_step / 2 >= x0)
                                y[idx] += y_raw[i]
                                cts[idx] += norm[i]

                            # errror bars for detector only
                            if "det" in y_str:
                                yerr = np.sqrt(y) / cts * norm_val
                            y = y / cts * norm_val

                        else:  # rebin, no renorm
                            for i, x0 in enumerate(x_raw):
                                idx = np.nanargmax(x + rebin_step / 2 >= x0)
                                y[idx] += y_raw[i]
                                cts[idx] += 1

                            # errror bars for detector only
                            if "det" in y_str:
                                yerr = np.sqrt(y) / cts
                            y = y / cts

                    case _:
                        print('Unrecogonized rebin type. Needs to be "tol" or "grid".')
            else:
                print("Rebin step needs to be greater than zero.")

        # generate labels and title
        if norm_channel is not None:
            if norm_channel == "time":
                norm_channel = "seconds"
            if norm_val == 1:
                ylabel = y_str + "/ " + norm_channel
            else:
                ylabel = y_str + f" / {norm_val} " + norm_channel
        else:
            preset_val = self.scan_info["preset_value"]
            ylabel = y_str + f" / {preset_val} " + self.scan_info["preset_channel"]

        xlabel = x_str
        label = "scan " + str(self.scan_info["scan"])
        title = label + ": " + self.scan_info["scan_title"]

        return (x, y, xerr, yerr, xlabel, ylabel, title, label)

    def plot_curve(
        self,
        x_str=None,
        y_str=None,
        norm_channel=None,
        norm_val=1,
        rebin_type=None,
        rebin_step=0,
    ):
        """Plot a 1D curve gnerated from a singal scan in a new window"""

        x, y, xerr, yerr, xlabel, ylabel, title, _ = self.generate_curve(
            x_str,
            y_str,
            norm_channel,
            norm_val,
            rebin_type,
            rebin_step,
        )

        fig, ax = plt.subplots()
        ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt="o")
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.6)

        fig.show()
