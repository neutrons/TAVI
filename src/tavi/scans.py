import numpy as np


class Scan(object):
    """
    Manage a single measued scan

    Attributes:
        metadata (dict): scan properties, includes scan_info,
            sample_UB_info, instrument_info
        data (dict): dictionary contains lists of scan data

    """

    def __init__(self, scan_id=None):
        """Initialze an empty scan if scan_id is None, otherwise load_scan"""
        self.metadata = {}
        self.data = {}
        if scan_id is not None:
            self.load_scan(scan_id)

    def load_scan(self, scan_id):
        """Load metadata and data in scan

        Args:
            scan_id (str):
        """
        pass

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
            norm (dict): key="time" or "monitor"

        Returns:
            tuple: data entry"""

        if x_str == None:
            x_str = self.metadata["def_x"]

        if y_str == None:
            y_str = self.metadata["def_y"]

        x = self.data[x_str]
        y = self.data[y_str]
        xerr = None
        yerr = None
        xlabel = None
        ylabel = None
        title = None

        return (x, y, xerr, yerr, xlabel, ylabel, title)
