import json

import numpy as np
from tavi.instrument.tas_cmponents import *
from tavi.sample.powder import Powder
from tavi.sample.sample import Sample
from tavi.sample.xtal import Xtal
from tavi.utilities import *


class TAS(object):
    """Triple-axis instrument class. Manage instrument congigutarion parameters

    Attibutes:

    Methods:
        load_instrument(config_params):
        load_sample(sample_params):


    """

    def __init__(self, path_to_json=None):
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

        if path_to_json is not None:
            self.load_instrument_from_json(path_to_json)

    def load_instrument_from_dicts(self, config_params):
        """Load instrument configuration from a json file"""

        self.source = Source(config_params["source"])
        self.collimators = Collimators(config_params["collimators"])
        self.guide = Guide(config_params["guide"])
        self.monochromator = Monochromator(config_params["monochromator"])
        self.monitor = Monitor(config_params["monitor"])
        self.goniometer = Goniometer(config_params["goniometer"])
        self.analyzer = Analyzer(config_params["analyzer"])
        self.detector = Detector(config_params["detector"])
        self.arms = Arms(config_params["distances"])

    def load_instrument_from_json(self, path_to_json):
        """Load instrument configuration from a json file"""

        with open(path_to_json, "r", encoding="utf-8") as file:
            config_params = json.load(file)

        self.load_instrument_from_dicts(config_params)

    # def save_instrument(self):
    #     """Save configuration into a dictionary"""
    #     # convert python dictionary to json file
    #   with open("./src/tavi/instrument/instrument_params/takin_test.json", "w") as file:
    #     json.dump(instrument_config, file)

    def load_sample(self, sample):
        """Load sample info"""
        self.sample = sample

    def load_sample_from_json(self, path_to_json):
        """Load sample info"""

        with open(path_to_json, "r", encoding="utf-8") as file:
            sample_params = json.load(file)

        if sample_params["type"] == "xtal":
            sample = Xtal.from_json(sample_params)
        elif sample_params["type"] == "powder":
            sample = Powder.from_json(sample_params)
        else:
            sample = Sample.from_json(sample_params)

        self.load_sample(sample)

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
                -kf * np.sin(np.deg2rad(two_theta)),
                0,
                ki - kf * np.cos(np.deg2rad(two_theta)),
            ]
        )

    def _find_u_from_2peaks(self, peaks, angles, ki, kf):
        """Calculate UB matrix from two peaks

        Args:
            peaks (list): lists of two tuples, [(h1,k1,l1), (h2,k2,l2)]
            angles (list): lists of goniometer angles
            ki (float): incident neutron wavelength, in inv Ang
            kf (float): final neutron wavelength, in inv Ang

        Returns:
            u_mat: 3 by 3 unitary matrix
            in_plane_ref: vector of peak1 in q_sample frame (zero goniometer angles)
            plane_normal: normal vector of the scattering plane in q_sample frame,
                            always pointing up (positive y-axis)
        """

        b_mat = self.sample.b_mat()
        peak1, peak2 = peaks
        q_hkl1 = b_mat @ peak1
        q_hkl2p = b_mat @ peak2
        q_hkl3 = np.cross(q_hkl1, q_hkl2p)
        q_hkl_2 = np.cross(q_hkl3, q_hkl1)

        q_hkl_mat = np.array(
            [
                q_hkl1 / np.linalg.norm(q_hkl1),
                q_hkl_2 / np.linalg.norm(q_hkl_2),
                q_hkl3 / np.linalg.norm(q_hkl3),
            ]
        ).T

        # find r_inv
        angles1, angles2 = angles
        two_theta1, _, _, _ = angles1
        two_theta2, _, _, _ = angles2

        q_lab1 = TAS.q_lab(two_theta1, ki=ki, kf=kf)
        q_lab2 = TAS.q_lab(two_theta2, ki=ki, kf=kf)

        # Goniometer angles all zeros in q_sample frame
        q_sample1 = self.goniometer.r_mat_inv(angles1[1:]) @ q_lab1
        q_sample2p = self.goniometer.r_mat_inv(angles2[1:]) @ q_lab2
        q_sample3 = np.cross(q_sample1, q_sample2p)
        q_sample2 = np.cross(q_sample3, q_sample1)

        q_sample1 = q_sample1 / np.linalg.norm(q_sample1)
        q_sample2 = q_sample2 / np.linalg.norm(q_sample2)
        q_sample3 = q_sample3 / np.linalg.norm(q_sample3)

        q_sample_mat = np.array([q_sample1, q_sample2, q_sample3]).T

        u_mat = q_sample_mat @ np.linalg.inv(q_hkl_mat)

        plane_normal = q_sample3
        if plane_normal[1] < 0:  # plane normal always up along +Y
            plane_normal = -plane_normal

        in_plane_ref = q_sample1

        return u_mat, in_plane_ref, plane_normal

    def find_ub(self, peaks, angles, ei=13.5, ef=None):
        """calculate UB matrix from peaks and motor positions

        Args:
            peaks (list)
            angles (list)
            ei (float): incident neutron energy, in meV
            ef (float): final neutron energy, in meV

        """

        self.sample.ub_peaks = peaks
        self.sample.ub_angles = angles

        ki = np.sqrt(ei / ksq2eng)
        if ef is None:
            kf = ki
        else:
            kf = np.sqrt(ef / ksq2eng)

        if not len(peaks) == len(angles):
            print("Number of peaks and angles provided do not match.")

        if len(peaks) == 2:
            b_mat = self.sample.b_mat()
            (
                u_mat,
                in_plane_ref,
                plane_normal,
            ) = self._find_u_from_2peaks(peaks, angles, ki, kf)
            ub_matrix = u_mat @ b_mat

            self.sample.in_plane_ref = in_plane_ref
            self.sample.plane_normal = plane_normal

        # TODO
        elif len(peaks) == 3:  # find_ub_from_3peaks
            pass
        # TODO
        elif len(peaks) > 3:  # find_ub_from_mulitple_peaks
            pass
        else:
            print("I don't even know what you're doing.")

        self.sample.ub_matrix = ub_matrix
        inv_ub_matrix = np.linalg.inv(ub_matrix)
        self.sample.inv_ub_matrix = inv_ub_matrix

        # print(np.round(ub_matrix, 6))

    def find_two_theta(self, peak, ei=13.5, ef=None):
        """calculate two theta angle of a given peak

        Args:
            peak (tuple): (h, k, l) of a peak
            ei (float): incident neutron energy, in meV
            ef (float): final neutron energy, in meV
        Return
            two_theta: in degrees
        """
        S2_MIN_DEG = 1

        hkl = np.array(peak)
        ki = np.sqrt(ei / ksq2eng)
        if ef is None:
            kf = ki
        else:
            kf = np.sqrt(ef / ksq2eng)

        b_mat = self.sample.b_mat()
        q_sq = 4 * np.pi**2 * hkl.T @ b_mat.T @ b_mat @ hkl
        q_norm = np.sqrt(q_sq)

        # two_theta = np.arccos((ki**2 + kf**2 - q_sq) / (2 * ki * kf))
        two_theta = get_angle(ki, kf, q_norm)
        if two_theta is None:
            print(f"Triangle cannot be closed at q={hkl}, en={ei-ef} meV.")
            return None
        elif np.rad2deg(two_theta) < S2_MIN_DEG:
            print(f"s2 is smaller than {S2_MIN_DEG} deg at q={hkl}.")
            return None
        else:
            return np.rad2deg(two_theta) * self.goniometer.sense

    def find_angles(self, peak, ei=13.5, ef=None):
        """calculate motor positions for a given peak if UB matrix has been determined

        Args:
            peak (tuple): (h, k, l) of a peak
            ei (float): incident neutron energy, in meV
            ef (float): final neutron energy, in meV

        """
        S2_MIN_DEG = 1

        hkl = np.array(peak)
        ki = np.sqrt(ei / ksq2eng)
        if ef is None:
            kf = ki
        else:
            kf = np.sqrt(ef / ksq2eng)

        b_mat = self.sample.b_mat()

        q_sq = 4 * np.pi**2 * hkl.T @ b_mat.T @ b_mat @ hkl
        q_norm = np.sqrt(q_sq)

        # two_theta = np.arccos((ki**2 + kf**2 - q_sq) / (2 * ki * kf)) * self.goniometer.sense
        two_theta = get_angle(ki, kf, q_norm)

        if two_theta is None:
            print(f"Triangle cannot be closed at q={hkl}, en={ei-ef} meV.")
            angles = None
        elif np.rad2deg(two_theta) < S2_MIN_DEG:
            print(f"s2 is smaller than {S2_MIN_DEG} deg at q={hkl}.")
            angles = None
        else:
            two_theta = np.rad2deg(two_theta) * self.goniometer.sense
            q = self.sample.ub_matrix @ hkl
            t1 = q / np.linalg.norm(q)

            eps = 1e-8  # zero
            plane_normal = self.sample.plane_normal
            in_plane_ref = self.sample.in_plane_ref

            # Minimal tilts
            if np.dot(t1, plane_normal) < eps:  # t1 in plane
                t3 = plane_normal
                t2 = np.cross(t3, t1)
            else:  # t1 not in plane, need to change tilts
                if np.linalg.norm(np.cross(plane_normal, t1)) < eps:
                    # oops, t1 along plane_normal
                    t2 = in_plane_ref
                    t3 = np.cross(t1, t2)
                else:
                    t2p = np.cross(plane_normal, t1)
                    t3 = np.cross(t1, t2p)
                    t2 = np.cross(t3, t1)

            t_mat = np.array(
                [t1, t2 / np.linalg.norm(t2), t3 / np.linalg.norm(t3)],
            ).T

            t_mat_inv = np.linalg.inv(t_mat)

            q_lab1 = TAS.q_lab(two_theta, ki, kf) / q_norm
            q_lab2 = np.array([q_lab1[2], 0, -q_lab1[0]])
            q_lab3 = np.array([0, 1, 0])

            q_lab_mat = np.array([q_lab1, q_lab2, q_lab3]).T
            r_mat = q_lab_mat @ t_mat_inv

            angles = np.round((two_theta,) + self.goniometer.angles_from_r_mat(r_mat), 6)

        return angles
