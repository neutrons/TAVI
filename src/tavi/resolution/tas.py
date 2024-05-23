import numpy as np

# from tavi.resolution.hb3 import config_params
from tavi.resolution.takin_test import config_params
from tavi.sample import Sample


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
        self.analyzer = {}
        self.detector = {}
        self.distances = {}

    def load_config(self, config_params):
        """Load a dictornary of instrument configuration"""
        self.source = config_params["source"]
        self.collimators = config_params["collimators"]
        self.guide = config_params["guide"]
        self.monochromator = config_params["monochromator"]
        self.monitor = config_params["monitor"]
        self.analyzer = config_params["analyzer"]
        self.detector = config_params["detector"]
        self.distances = config_params["distances"]

    def save_config(self):
        """Save configuration into a dictionary"""
        pass

    def cooper_nathans(self, ei, ef, q):
        """
        Calculate resolution using Cooper-Nathans method
        Args:
           ei
           ef
            q
        """
        res = 0
        return res


if __name__ == "__main__":

    hb3 = TAS()
    hb3.load_config(config_params)
    print(hb3.analyzer["type"])

    hb3_sample = Sample(lattice_params=(5.3995, 5.64, 11.75, 90, 90, 90))

    # hb3.cooper_nathans(ki=1.4, kf=1.4, en, q=1.777))

    # describe and plot ellipses
    # ellipses = reso.calc_ellipses(res["reso"], verbose)
    # reso.plot_ellipses(ellipses, verbose)
