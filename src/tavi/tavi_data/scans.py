import numpy as np
from tavi.tavi_data.nexus_reader import nexus_to_dict


class Scan(object):
    """
    Manage a single measued scan

    Attributes:
        scan_info (dict):
        sample_ub_info (dict):
        instrument_info (dict):
        data (dict): dictionary contains lists of scan data

    """

    def __init__(self, nexus_entry=None):
        """Initialze an empty scan"""
        self.scan_info = None
        self.sample_ub_info = None
        self.instrument_info = None
        self.data = None

        if nexus_entry is not None:
            self.load_scan(nexus_entry)

    def load_scan(self, nexus_entry):
        """Unpack metadata and data from scan_data

        Args:
            nexus_entry:
        """

        scan_info, sample_ub_info, instrument_info, data = nexus_to_dict(nexus_entry)
        self.scan_info = scan_info
        self.sample_ub_info = sample_ub_info
        self.instrument_info = instrument_info
        self.data = data

    # def set_metadata(self, meta_data):
    #     """Set metadata"""
    #     self.meta_data = meta_data

    # def set_data(self, data):
    #     """Set metadata"""
    #     self.data = data

    def get_scan_info(self):
        """Return scan_info in metadata.

        Returns:
            dict: dictionay of scan_info metadata.
        """
        return None

    def get_sample_UB_info(self):
        """Return sample_UB_info in metadata.

        Returns:
            dict: dictionay of sample_UB_info metadata.
        """
        return None

    def get_instrument_info(self):
        """Return instrument_info in metadata.

        Returns:
            dict: dictionay of instrument_info metadata.
        """
        return None

    def save_metadata(self, metadata_entry):
        """Save metadata_entry into file

        Args:
           metadata_entry (dict): {key: value}
        """
        pass

    def get_data_entry(self, entry_name):
        """Return data entry based on entry_name

        Args:
            entry_name (str): a key of the dictionay in data

        Returns:
            tuple: data entry
        """

    def curve_gen(
        self,
        x_str=None,
        y_str=None,
        norm=None,
    ):
        """Generate a curve to plot

        Args:
            x_str (str): string of x axis
            y_str (str): string of x axis
            norm (dict): None, "time", "monitor" or "mcu"

        Returns:
            tuple: data entry"""

        if x_str == None:
            x_str = self.scan_info["def_x"]

        if y_str == None:
            y_str = self.scan_info["def_y"]

        x = self.data[x_str]
        y = self.data[y_str]
        xerr = None
        yerr = np.sqrt(y)
        xlabel = x_str
        ylabel = y_str
        label = "scan " + str(self.scan_info["scan"])
        title = label + ": " + self.scan_info["scan_title"]

        return (x, y, xerr, yerr, xlabel, ylabel, title, label)
