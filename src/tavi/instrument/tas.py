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
    def q_lab(two_theta, ki, kf):
        """Momentum transfer q in lab frame.

        Note:
            Only for a single detector in the scattering plane.
            Depends on inciendnt and final wavevector ki, kf
            and two theta, phi is zero.
        """
        return np.array(
            [
                -kf * np.sin(two_theta / rad2deg),
                0,
                ki - kf * np.cos(two_theta / rad2deg),
            ]
        )

    def _find_u_from_2peaks(self, peaks, angles, ki, kf):
        """Calculate UB matrix from two peaks

        Args:
            peaks (list)
            angles (list)
            ki (float): incident neutron wavelength, in inv Ang
            kf (float): final neutron wavelength, in inv Ang
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

        q_lab1 = TAS.q_lab(two_theta1, ki=ki, kf=kf)
        q_lab2 = TAS.q_lab(two_theta2, ki=ki, kf=kf)

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
    def find_ub(self, peaks, angles, ei=13.5, ef=None):
        """calculate UB matrix from peaks and motor positions

        Args:
            peaks (list)
            angles (list)
            eng (float): incident neutron energy, in meV

        """

        ki = np.sqrt(ei / ksq2eng)
        if ef is None:
            kf = ki
        else:
            kf = np.sqrt(ef / ksq2eng)

        if not len(peaks) == len(angles):
            print("Number of peaks and angles provided do not match.")

        if len(peaks) == 2:
            b_mat = self.sample.b_mat()
            u_mat = self._find_u_from_2peaks(peaks, angles, ki, kf)
            ub_matrix = u_mat @ b_mat

        elif len(peaks) == 3:  # find_ub_from_3peaks
            pass
        elif len(peaks) > 3:  # find_ub_from_mulitple_peaks
            pass
        else:
            print("I don't even know what you're doing.")

        self.sample.ub_peaks = peaks
        self.sample.ub_angles = angles
        self.sample.ub_matrix = ub_matrix
        self.sample.inv_ub_matrix = np.linalg.inv(ub_matrix)
        # print(np.round(ub_matrix, 6))
        return ub_matrix

    def find_angles(self, peak, ei=13.5, ef=None):
        """calculate motor positions for a given peak if UB matrix has been determined

        Args:
            peak (tuple): (h, k, l) of a peak
            ei (float): incident neutron energy, in meV
            ef (float): final neutron energy, in meV

        Note:
            two_theta = 2 * acrsin(|q|/(2k))

        """
        hkl = np.array(peak)
        ki = np.sqrt(ei / ksq2eng)
        if ef is None:
            kf = ki
        else:
            kf = np.sqrt(ef / ksq2eng)

        b_mat = self.sample.b_mat()

        q_sq = 4 * np.pi**2 * hkl.T @ b_mat.T @ b_mat @ hkl
        q_norm = np.sqrt(q_sq)
        two_theta = np.arccos((ki**2 + kf**2 - q_sq) / (2 * ki * kf)) * self.goniometer.sense

        q = self.sample.ub_matrix @ hkl
        t1 = q / np.linalg.norm(q)
        t2p = np.array([t1[2], 0, -t1[0]])
        t3 = np.cross(t1, t2p)
        t2 = np.cross(t3, t1)

        # t2 = [-0.332928, 0.014534, -0.942840]
        # t3 = [-0.043637, 0.999047, 0.000009]

        t_mat = np.array(
            [t1, t2 / np.linalg.norm(t2), t3 / np.linalg.norm(t3)],
        ).T

        t_mat_inv = np.linalg.inv(t_mat)

        q_lab1 = TAS.q_lab(two_theta * rad2deg, ki, kf) / q_norm
        q_lab2 = [q_lab1[2], 0, -q_lab1[0]]
        q_lab3 = [0, 1, 0]

        q_lab_mat = np.array([q_lab1, q_lab2, q_lab3]).T
        r_mat = q_lab_mat @ t_mat_inv

        angles = self.goniometer.angles_from_r_mat(r_mat)
        angles = (two_theta * rad2deg,) + angles

        print(angles)
        return angles


if __name__ == "__main__":

    takin_instru = TAS()
    takin_instru.load_instrument(instrument_params)
    takin_instru.load_sample(test_xtal)

    # print(takin_instru.analyzer.type)
    # print(takin_instru.sample.a)

    peak_list = [
        (0, 0, 2),
        (0, 2, 0),
    ]
    angles_list = [
        (-51.530388, -45.220125, -0.000500, -2.501000),
        (-105.358735, 17.790125, -0.000500, -2.501000),
    ]

    takin_instru.find_ub(peaks=peak_list, angles=angles_list, ei=13.500172, ef=13.505137)
    takin_instru.find_angles(peak=(0, 0, 2), ei=13.50172, ef=13.505137)

    # describe and plot ellipses
    # ellipses = reso.calc_ellipses(res["reso"], verbose)
    # reso.plot_ellipses(ellipses, verbose)
