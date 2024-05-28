# from tavi.resolution.hb3 import config_params
from tavi.instrument_params.takin_test import instrument_params, sample_params


class TAS(object):
    """Triple-axis instrument class. Manage instrument congigutarion parameters

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


if __name__ == "__main__":

    takin_instru = TAS()
    takin_instru.load_config(instrument_params, sample_params)

    print(takin_instru.analyzer["type"])
    print(takin_instru.sample["sense"])
    print(takin_instru.sample["xtal"].a)

    # describe and plot ellipses
    # ellipses = reso.calc_ellipses(res["reso"], verbose)
    # reso.plot_ellipses(ellipses, verbose)
