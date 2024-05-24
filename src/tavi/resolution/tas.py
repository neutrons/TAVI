import numpy as np
from tavi.utilities import *

# from tavi.resolution.hb3 import config_params
from tavi.resolution.takin_test import instrument_params, sample_params


class TAS(object):
    """Triple-axis instrument class

    Attibutes:


    Methods:

    """

    # 4 soller slits collimators
    NUM_COLLS = 4
    IDX_COLL0_H = 0
    IDX_COLL0_V = 2
    IDX_COLL1_H = 1
    IDX_COLL1_V = 3
    IDX_COLL2_H = 4
    IDX_COLL2_V = 6
    IDX_COLL3_H = 5
    IDX_COLL3_V = 7
    # 1 monochromator and 1 analyzer
    NUM_MONOS = 1
    NUM_ANAS = 1
    IDX_MONO0_H = 0
    IDX_MONO0_V = 1
    IDX_ANA0_H = 2
    IDX_ANA0_V = 3

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

    def _mono_ana_transfer_matrix(self):
        """Transfer matrix for monochromator or anlyzer"""
        pass

    def cooper_nathans(self, ei, ef, q, R0=False):
        """
        Calculate resolution using Cooper-Nathans method

        Args:
            ei (float): incident energy, in units of meV
            ef (float): final energy, in units of meV
            q (float): momentum transfer, in units of inverse angstrom
            R0 (bool): calculate normalization factor if True
        """

        ki = np.sqrt(ei / ksq2E)
        kf = np.sqrt(ef / ksq2E)
        en = ei - ef

        two_theta = get_angle(ki, kf, q) * self.sample["sense"]
        theta_s = two_theta / 2
        theta_m = get_angle_Bragg(ki, self.monochromator["d_spacing"]) * self.monochromator["sense"]
        theta_a = get_angle_Bragg(kf, self.analyzer["d_spacing"]) * self.analyzer["sense"]

        # TODO
        # curved monochromator and analyzer

        # TODO
        # reflection efficiency

        # matrix G, divergence of collimators, [pop75] Appendix 1
        mat_g = np.zeros((TAS.NUM_COLLS * 2, TAS.NUM_COLLS * 2))
        mat_g[TAS.IDX_COLL0_H, TAS.IDX_COLL0_H] = 1.0 / self.collimators["h_pre_mono"] ** 2
        mat_g[TAS.IDX_COLL0_V, TAS.IDX_COLL0_V] = 1.0 / self.collimators["v_pre_mono"] ** 2
        mat_g[TAS.IDX_COLL1_H, TAS.IDX_COLL1_H] = 1.0 / self.collimators["h_pre_sample"] ** 2
        mat_g[TAS.IDX_COLL1_V, TAS.IDX_COLL1_V] = 1.0 / self.collimators["v_pre_sample"] ** 2
        mat_g[TAS.IDX_COLL2_H, TAS.IDX_COLL2_H] = 1.0 / self.collimators["h_post_sample"] ** 2
        mat_g[TAS.IDX_COLL2_V, TAS.IDX_COLL2_V] = 1.0 / self.collimators["v_post_sample"] ** 2
        mat_g[TAS.IDX_COLL3_H, TAS.IDX_COLL3_H] = 1.0 / self.collimators["h_post_ana"] ** 2
        mat_g[TAS.IDX_COLL3_V, TAS.IDX_COLL3_V] = 1.0 / self.collimators["v_post_ana"] ** 2

        # matrix F, divergence of monochromator and analyzer, [pop75] Appendix 1
        mat_f = np.zeros(((TAS.NUM_MONOS + TAS.NUM_ANAS) * 2, (TAS.NUM_MONOS + TAS.NUM_ANAS) * 2))
        mat_f[TAS.IDX_MONO0_H, TAS.IDX_MONO0_H] = 1.0 / self.monochromator["mosaic"] ** 2
        mat_f[TAS.IDX_MONO0_V, TAS.IDX_MONO0_V] = 1.0 / self.monochromator["mosaic_v"] ** 2
        mat_f[TAS.IDX_ANA0_H, TAS.IDX_ANA0_H] = 1.0 / self.analyzer["mosaic"] ** 2
        mat_f[TAS.IDX_ANA0_V, TAS.IDX_ANA0_V] = 1.0 / self.analyzer["mosaic_v"] ** 2

        # matrix C, constrinat between mono/ana mosaic and collimator divergence
        mat_c = np.zeros(((TAS.NUM_MONOS + TAS.NUM_ANAS) * 2, TAS.NUM_COLLS * 2))
        mat_c[TAS.IDX_MONO0_H, TAS.IDX_COLL0_H] = 0.5
        mat_c[TAS.IDX_MONO0_H, TAS.IDX_COLL1_H] = 0.5
        # mat_c[TAS.IDX_MONO0_V, TAS.IDX_COLL0_V] = 0.5
        # mat_c[TAS.IDX_MONO0_V, TAS.IDX_COLL1_V] = 0.5

        mat_c[TAS.IDX_ANA0_H, TAS.IDX_COLL2_H] = 0.5
        mat_c[TAS.IDX_ANA0_H, TAS.IDX_COLL3_H] = 0.5

        # matrix A, ki kf distribution to collimators

        # matrix B, lab frame to local q-frame
        mat_b = np.zeros((4, 6))

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
