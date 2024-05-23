import numpy as np

# from tavi.resolution.hb3 import config_params
from tavi.resolution.takin_test import instrument_params, sample_params


class TAS(object):
    """Triple-axis instrument class

    Attibutes:


    Methods:

    """

    def __init__(self):
        self.source = {}
        self.collimators = {}
        self.guide = {}
        self.monochromator = {}
        self.monitor = {}
        self.sample = {}
        self.analyzer = {}
        self.detector = {}
        self.distances = {}
        self.sample = {}

    def load_config(self, config_params, sample_params=None):
        """Load a dictornary of instrument configuration"""
        self.source = config_params["source"]
        self.collimators = config_params["collimators"]
        self.guide = config_params["guide"]
        self.monochromator = config_params["monochromator"]
        self.monitor = config_params["monitor"]
        self.analyzer = config_params["analyzer"]
        self.detector = config_params["detector"]
        self.distances = config_params["distances"]

        if sample_params is not None:
            self.sample = sample_params

    def save_config(self):
        """Save configuration into a dictionary"""
        pass

    # def load_sample(self, sample):
    #     self.sample = sample

    def cooper_nathans(self, ei, ef, q, R0=False):
        """
        Calculate resolution using Cooper-Nathans method

        Args:
            ei (float): incident energy, in units of meV
            ef (float): final energy, in units of meV
            q (float): momentum transfer, in units of inverse angstrom
            R0 (bool): calculate normalization factor if True
        """

        res = 0
        return res


if __name__ == "__main__":

    takin = TAS()
    takin.load_config(instrument_params, sample_params)

    print(takin.analyzer["type"])
    print(takin.sample["sense"])
    print(takin.sample["xtal"].a)

    # ki = kf = 1.4
    takin.cooper_nathans(ei=1.4**2 * 2.072124855, ef=1.4**2 * 2.072124855, q=1.777, R0=False)

    # describe and plot ellipses
    # ellipses = reso.calc_ellipses(res["reso"], verbose)
    # reso.plot_ellipses(ellipses, verbose)
