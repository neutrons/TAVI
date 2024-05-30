# from tavi.resolution.hb3 import config_params
from tavi.instrument_params.takin_test import instrument_params
from tavi.sample.sample_test import test_xtal


class TAS(object):
    """Triple-axis instrument class. Manage instrument congigutarion parameters

    Attibutes:
        source (dict):
        collimators (dict):
        guide (dict):
        monochromator (dict):
        monitor (dict):
        goniometer (dict):
        analyzer (dict):
        detector (dict):
        distances (dict):
        sample (class):


    Methods:
        load_instrument(config_params):
        load_sample(sample_params):


    """

    def __init__(self):
        self.source = {}
        self.collimators = {}
        self.guide = {}
        self.monochromator = {}
        self.monitor = {}
        self.goniometer = {}
        self.analyzer = {}
        self.detector = {}
        self.distances = {}
        self.sample = None

    def load_instrument(self, config_params):
        """Load a dictornary of instrument configuration"""
        self.source = config_params["source"]
        self.collimators = config_params["collimators"]
        self.guide = config_params["guide"]
        self.monochromator = config_params["monochromator"]
        self.monitor = config_params["monitor"]
        self.goniometer = config_params["goniometer"]
        self.analyzer = config_params["analyzer"]
        self.detector = config_params["detector"]
        self.distances = config_params["distances"]

    def save_instrument(self):
        """Save configuration into a dictionary"""
        pass

    def load_sample(self, sample):
        """Load sample info"""
        self.sample = sample


if __name__ == "__main__":

    takin_instru = TAS()
    takin_instru.load_instrument(instrument_params)
    takin_instru.load_sample(test_xtal)

    print(takin_instru.analyzer["type"])
    print(takin_instru.sample.a)

    # describe and plot ellipses
    # ellipses = reso.calc_ellipses(res["reso"], verbose)
    # reso.plot_ellipses(ellipses, verbose)
