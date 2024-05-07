import numpy as np


class Scan(object):
    """
    Scan class
    """

    def __init__(self):
        self.metadata = {}
        self.data = {}

    def get_metadata(self):
        return self.metadata

    def get_data(self):
        return self.data

    def plot_gen(
        self,
        x_str=None,
        y_str=None,
        xerr_str=None,
        yerr_str=None,
        norm_str=None,
        norm_val=None,
    ):
        if x_str == None:
            x_str = self.get_metadata()["def_x"]

        if y_str == None:
            y_str = self.get_metadata()["def_y"]

        x = self.data[x_str]
        y = self.data[y_str]
        xerr = None
        yerr = None

        return x, y, xerr, yerr
