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

    @staticmethod
    def q_lab(two_theta, k):
        """Momentum transfer q in lab frame.

        Note:
            Only for a single detector in the scattering plane.
            Depends on wavevector k and two theta, phi is zero.
        """
        return k * np.array(
            [
                -np.sin(two_theta / rad2deg),
                0,
                1 - np.cos(two_theta / rad2deg),
            ]
        )

    def _find_u_from_2peaks(self, peaks, angles, k):
        """Calculate UB matrix from two peaks

        Args:
            peaks (list)
            angles (list)
            ei (float): incident neutron energy, in meV
        """

        b_mat = self.sample.b_mat()
        peak1, peak2 = peaks
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
        angles1, angles2 = angles
        two_theta1, _, _, _ = angles1
        two_theta2, _, _, _ = angles2

        q_lab1 = TAS.q_lab(two_theta1, k)
        q_lab2 = TAS.q_lab(two_theta2, k)

        q_sample1 = self.goniometer.r_mat_inv(angles1) @ q_lab1
        q_sample2 = self.goniometer.r_mat_inv(angles2) @ q_lab2
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
        """calculate UB matrix from peaks and motor positions

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
            ub_matrix = u_mat @ b_mat

        elif len(peaks) == 3:
            pass
        elif len(peaks) > 3:
            pass
        else:
            print("I don't even know what you're doing.")

        self.sample.ub_matrix = ub_matrix
        print(ub_matrix)
        return ub_matrix


if __name__ == "__main__":

    takin_instru = TAS()
    takin_instru.load_instrument(instrument_params)
    takin_instru.load_sample(test_xtal)

    print(takin_instru.analyzer.type)
    print(takin_instru.sample.a)

    peak_list = [
        (0, 0, 2),
        (0, 2, 0),
    ]
    angles_list = [
        (-51.530388, -45.220125, -0.000500, -2.501000),
        (-105.358735, 17.790125, -0.000500, -2.501000),
    ]

    takin_instru.find_ub(peaks=peak_list, angles=angles_list, ei=13.5)
    pass

    # describe and plot ellipses
    # ellipses = reso.calc_ellipses(res["reso"], verbose)
    # reso.plot_ellipses(ellipses, verbose)
