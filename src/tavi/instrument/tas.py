# from tavi.resolution.hb3 import config_params
from tavi.instrument.instrument_params.takin_test import instrument_params
from tavi.sample.sample_test import test_xtal
from tavi.utilities import *
from tavi.instrument.tas_cmponents import *


class TAS(object):
    """Triple-axis instrument class. Manage instrument congigutarion parameters

    Attibutes:

    Methods:
        load_instrument(config_params):
        load_sample(sample_params):


    """

    def __init__(self):
        """load default instrument configuration"""
        # TODO
        self.source = None
        self.collimators = None
        self.guide = None
        self.monochromator = None
        self.monitor = None
        self.goniometer = None
        self.analyzer = None
        self.detector = None
        self.arms = None
        self.sample = None

    def load_instrument(self, config_params):
        """Load a dictornary of instrument configuration"""
        self.source = Source(config_params["source"])
        self.collimators = Collimators(config_params["collimators"])
        # self.guide = Guide(config_params["guide"])
        self.monochromator = Monochromator(config_params["monochromator"])
        self.monitor = Monitor(config_params["monitor"])
        self.goniometer = Goniometer(config_params["goniometer"])
        self.analyzer = Analyzer(config_params["analyzer"])
        self.detector = Detector(config_params["detector"])
        self.arms = Arms(config_params["distances"])

    def save_instrument(self):
        """Save configuration into a dictionary"""
        pass

    def load_sample(self, sample):
        """Load sample info"""
        self.sample = sample

    def _find_u_from_2peaks(self):
        """Calculate UB matrix from two peaks

        Args:
            peaks (list)
            angles (list)
            ei (float): incident neutron energy, in meV"""
        q_hkl1 = b_mat @ peak1
        q_hkl2 = b_mat @ peak2
        q_hkl3 = np.cross(q_hkl1, q_hkl2)
        q_hkl_2p = np.cross(q_hkl3, q_hkl1)

        q_hkl_mat = np.array(
            [
                q_hkl1 / np.linalg.norm(q_hkl1),
                q_hkl_2p / np.linalg.norm(q_hkl_2p),
                q_hkl3 / np.linalg.norm(q_hkl3),
            ]
        ).T

        # find r_inv

        two_theta1, omega1, sgl1, sgu1 = angles1
        two_theta2, omega2, sgl2, sgu2 = angles2

        q_lab1 = k * np.array(
            [
                -np.sin(two_theta1 / rad2deg),
                0,
                1 - np.cos(two_theta1 / rad2deg),
            ]
        )
        q_lab2 = k * np.array(
            [
                -np.sin(two_theta2 / rad2deg),
                0,
                1 - np.cos(two_theta2 / rad2deg),
            ]
        )

        q_sample1 = r_inv_mat1 @ q_lab1
        q_sample2 = r_inv_mat2 @ q_lab2

        q_sample3 = np.cross(q_sample1, q_sample2)
        q_sample2p = np.cross(q_sample3, q_sample1)

        q_sample_mat = np.array(
            [
                q_sample1 / np.linalg.norm(q_sample1),
                q_sample2p / np.linalg.norm(q_sample2p),
                q_sample3 / np.linalg.norm(q_sample3),
            ]
        ).T

        u_mat = q_sample_mat @ np.linalg.inv(q_hkl_mat)
        return u_mat

    # TODO check goniometer order, sign convention,
    def find_ub(self, peaks, angles, ei=13.5):
        """calculate UB matrix

        Args:
            peaks (list)
            angles (list)
            ei (float): incident neutron energy, in meV

        """

        k = np.sqrt(ei / ksq2eng)

        if not len(peaks) == len(angles):
            print("Number of peaks and angles provided do not match.")

        if len(peaks) == 2:
            b_mat = self.sample.b_mat()
            u_mat = self._find_u_from_2peaks(peaks, angles, k)

        elif len(peaks) == 3:
            pass
        elif len(peaks) > 3:
            pass
        else:
            print("I don't even know what you're doing.")

        ub_matrix = u_mat @ b_mat
        self.sample.ub_matrix = u_mat @ b_mat
        return ub_matrix


if __name__ == "__main__":

    takin_instru = TAS()
    takin_instru.load_instrument(instrument_params)
    takin_instru.load_sample(test_xtal)

    print(takin_instru.analyzer.type)
    print(takin_instru.sample.a)

    peak1 = (0, 0, 2)
    peak2 = (0, 2, 0)
    angles1 = (-51.530388, -45.220125, -0.000500, -2.501000)
    angles2 = (-105.358735, 17.790125, -0.000500, -2.501000)

    peak_list = [peak1, peak2]
    angles_list = [angles1, angles2]

    takin_instru.find_ub(peaks=peak_list, angles=angles_list, ei=13.5)

    # describe and plot ellipses
    # ellipses = reso.calc_ellipses(res["reso"], verbose)
    # reso.plot_ellipses(ellipses, verbose)
