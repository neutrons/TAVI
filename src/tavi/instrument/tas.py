# -*- coding: utf-8 -*-
from typing import Optional

import numpy as np

from tavi.instrument.tas_base import TASBase
from tavi.sample.xtal import Xtal
from tavi.utilities import MotorAngles, Peak, UBConf, en2q, get_angle_from_triangle, mantid_to_spice, spice_to_mantid


class TAS(TASBase):
    """
    Triple-axis instrument class. Handles angle and UB calculations

    Methods:
        calculate_two_theta
        calculate_ub_matrix

    """

    def __init__(self, SPICE_CONVENTION: bool = True):
        super().__init__()
        self.SPICE_CONVENTION = SPICE_CONVENTION  # use coordination system defined in SPICE

    @staticmethod
    def q_lab(
        two_theta_deg: float,
        ei: float,
        ef: float,
    ):
        """
        Reutrn momentum transfer q in lab frame, using Mantid convention

        Note:
            Only for a single detector in the scattering plane.
        """

        ki = en2q(ei)
        kf = en2q(ef)
        two_theta = np.deg2rad(two_theta_deg)
        q = np.array([-kf * np.sin(two_theta), 0, ki - kf * np.cos(two_theta)])
        return q

    def calculate_two_theta(
        self,
        hkl: tuple[float, float, float],
        ei: float,
        ef: Optional[float] = None,
    ) -> Optional[float]:
        """find two theta angle for a given peak

        Args:
            hkl (tuple): miller indice of a peak
            ei (float): incident neutron energy, in emV
            ef (float): final neutron energy, in meV
        Returns:
            two_theta (float): two theta angle, in degrees

        Note:
            if ef is None, assume ef equals ei
        """

        if ei is None:
            raise ValueError("Cannot calculate two thetha without incident energy ei.")

        ki = en2q(ei)
        if ef is None:
            ef = ei
            kf = ki
        else:
            kf = en2q(ef)

        b_mat = self.sample.b_mat_from_lattice()
        qh, qk, ql = hkl
        hkl_array = np.array([qh, qk, ql])
        q_sq = 4 * np.pi**2 * hkl_array.T @ b_mat.T @ b_mat @ hkl_array
        q_norm = np.sqrt(q_sq)
        two_theta_radian = get_angle_from_triangle(ki, kf, q_norm)
        if two_theta_radian is None:
            print(f"Triangle cannot be closed at q=({qh}, {qk}, {ql}), en={ei-ef} meV.")
            return None
        # elif np.rad2deg(two_theta_radian) < S2_MIN_DEG:
        # pass

        return np.rad2deg(two_theta_radian) * self.goniometer.sense

    # TODO complete finding UB from 3 or more peaks
    def calculate_ub_matrix(self, peaks: tuple[Peak, ...]):
        """Find UB matrix from a list of observed peaks"""
        if not isinstance(self.sample, Xtal):
            raise ValueError("sample needs to be Xtal class for UB calculation.")

        match (num_of_peaks := len(peaks)):
            case 2:
                ubconf = self._find_u_from_two_peaks(peaks)
                self.sample.set_orientation(ubconf)
            case 3:
                ubconf = self._find_ub_from_three_peaks(peaks)
                self.sample.set_orientation(ubconf)
            case _ if num_of_peaks > 3:
                ubconf = self._find_ub_from_multiple_peaks(peaks)
                self.sample.set_orientation(ubconf)
            case _:
                raise ValueError("Not enough peaks for UB matrix determination.")

    def _find_u_from_two_peaks(
        self,
        peaks: tuple[Peak, Peak],
    ) -> UBConf:
        """Calculate UB matrix from two peaks for a given goniometer"""

        peak1, peak2 = peaks
        b_mat = self.sample.b_mat_from_lattice()
        q_hkl1 = np.matmul(b_mat, np.array(peak1.hkl))
        q_hkl2p = np.matmul(b_mat, np.array(peak2.hkl))
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
        r_mat_inv = self.goniometer.r_mat_inv
        q_lab1 = TAS.q_lab(peak1.angles.two_theta, ei=peak1.ei, ef=peak1.ef)
        q_lab2 = TAS.q_lab(peak2.angles.two_theta, ei=peak2.ei, ef=peak2.ef)

        # Goniometer angles all zeros in q_sample frame
        q_sample1 = np.matmul(r_mat_inv(peak1.angles), q_lab1)
        q_sample2p = np.matmul(r_mat_inv(peak2.angles), q_lab2)
        q_sample3 = np.cross(q_sample1, q_sample2p)
        q_sample2 = np.cross(q_sample3, q_sample1)

        q_sample1 = q_sample1 / np.linalg.norm(q_sample1)
        q_sample2 = q_sample2 / np.linalg.norm(q_sample2)
        q_sample3 = q_sample3 / np.linalg.norm(q_sample3)

        q_sample_mat = np.array([q_sample1, q_sample2, q_sample3]).T

        u_mat = np.matmul(q_sample_mat, np.linalg.inv(q_hkl_mat))

        plane_normal = q_sample3
        if plane_normal[1] < 0:  # plane normal always up along +Y
            plane_normal = -plane_normal
        # suffle the order following SPICE convention
        in_plane_ref = q_sample1

        ub_mat = np.matmul(u_mat, b_mat)

        if self.SPICE_CONVENTION:
            plane_normal = mantid_to_spice(plane_normal)
            in_plane_ref = mantid_to_spice(in_plane_ref)
            ub_mat = mantid_to_spice(ub_mat)

        return UBConf(
            peaks,
            u_mat,
            None,
            ub_mat,
            plane_normal,
            in_plane_ref,
        )

    # TODO
    def _find_ub_from_three_peaks(
        self,
        peaks: tuple[Peak, Peak, Peak],
    ) -> UBConf:
        """Find UB matrix from three observed peaks for a given goniomete"""
        ubconf = UBConf()
        return ubconf

    # TODO
    def _find_ub_from_multiple_peaks(
        self,
        peaks: tuple[Peak, ...],
    ) -> UBConf:
        """Find UB matrix from more than three observed peaks for a given goniomete"""
        ubconf = UBConf()
        return ubconf

    @staticmethod
    def norm_mat(t1, t2, t3):
        mat = np.array([t1 / np.linalg.norm(t1), t2 / np.linalg.norm(t2), t3 / np.linalg.norm(t3)]).T
        return mat

    def _t_mat_minimal_tilt(self, hkl: np.ndarray):
        """Build matrix T assuming minimal goniometer tilt angles"""

        if not isinstance(self.sample, Xtal):
            raise ValueError("Sample needs to be Xtal class for UB calculation.")
        if self.sample.ub_mat is None:
            raise ValueError("UB matrix is unknown.")
        if self.sample.plane_normal is None:
            raise ValueError("Plane normal vector is not known.")
        if self.sample.in_plane_ref is None:
            raise ValueError("In-plnae reference vector is not known.")

        EPS = 1e-8  # zero

        plane_normal = np.array(self.sample.plane_normal)
        in_plane_ref = np.array(self.sample.in_plane_ref)
        ub_mat = self.sample.ub_mat

        if self.SPICE_CONVENTION:  # suffle the order following SPICE convention
            plane_normal = spice_to_mantid(plane_normal)
            in_plane_ref = spice_to_mantid(in_plane_ref)
            ub_mat = spice_to_mantid(ub_mat)

        q = ub_mat @ hkl
        t1 = q / np.linalg.norm(q)

        if np.dot(t1, plane_normal) < EPS:  # t1 in plane
            t3 = plane_normal
            t2 = np.cross(t3, t1)
            return TAS.norm_mat(t1, t2, t3)

        # t1 not in plane, need to change tilts
        if np.linalg.norm(np.cross(plane_normal, t1)) < EPS:
            # oops, t1 along plane_normal
            t2 = in_plane_ref
            t3 = np.cross(t1, t2)
            return TAS.norm_mat(t1, t2, t3)

        else:
            t2p = np.cross(plane_normal, t1)
            t3 = np.cross(t1, t2p)
            t2 = np.cross(t3, t1)
            return TAS.norm_mat(t1, t2, t3)

    def calculate_motor_angles(self, peak: tuple, ei: float, ef: Optional[float] = None):
        """calculate motor positions for a given peak if UB matrix has been determined

        Args:
            peak (tuple): (h, k, l) of a peak
            ei (float): incident neutron energy, in meV
            ef (float): final neutron energy, in meV

        """
        S2_MIN_DEG = 1

        try:
            h, k, l = peak
        except ValueError:
            print("hkl should have the format (h, k, l).")
        hkl = np.array((h, k, l))

        ki = en2q(ei)
        if ef is None:
            ef = ei
            kf = ki
        else:
            kf = en2q(ef)

        if self.sample.b_mat is None:
            b_mat = self.sample.b_mat_from_lattice()
        else:
            b_mat = self.sample.b_mat
        q_norm = 2 * np.pi * np.sqrt(hkl.T @ b_mat.T @ b_mat @ hkl)

        two_theta = get_angle_from_triangle(ki, kf, q_norm)

        if two_theta is None:
            print(f"Triangle cannot be closed at q={hkl}, en={ei-ef} meV.")
            return None
        if np.rad2deg(two_theta) < S2_MIN_DEG:
            print(f"s2 is smaller than {S2_MIN_DEG} deg at q={hkl}.")
            return None

        two_theta_deg = np.rad2deg(two_theta) * self.goniometer.sense
        t_mat = self._t_mat_minimal_tilt(hkl)
        t_mat_inv = np.linalg.inv(t_mat)

        q_lab1 = TAS.q_lab(two_theta_deg, ei, ef) / q_norm
        q_lab2 = np.array([q_lab1[2], 0, -q_lab1[0]])
        q_lab3 = np.array([0, 1, 0])

        q_lab_mat = TAS.norm_mat(q_lab1, q_lab2, q_lab3)
        r_mat = q_lab_mat @ t_mat_inv

        _, omega, sgl, sgu, chi, phi = self.goniometer.angles_from_r_mat(r_mat)

        return MotorAngles(two_theta_deg, omega, sgl, sgu, chi, phi)
