# -*- coding: utf-8 -*-
from typing import Optional

import numpy as np

from tavi.instrument.tas_base import TASBase
from tavi.utilities import Peak, UBConf, en2q, get_angle_from_triangle


class TAS(TASBase):
    """
    Triple-axis instrument class. Handles angle and UB calculations

    Methods:
        calculate_two_theta
        calculate_ub_matrix

    """

    @staticmethod
    def q_lab(
        two_theta: float,
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

        q = np.array(
            [
                -kf * np.sin(np.deg2rad(two_theta)),
                0,
                ki - kf * np.cos(np.deg2rad(two_theta)),
            ]
        )
        return q

    def calculate_two_theta(
        self,
        hkl: tuple[float],
        ei: Optional[float] = None,
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

    # TODO
    def calculate_ub_matrix(
        self,
        peaks: tuple[Peak],
    ) -> np.ndarray:
        """Find UB matrix from a list of observed peaks"""

        num_of_peaks = len(peaks)
        if num_of_peaks == 2:
            ubconf = self._find_u_from_two_peaks(peaks)
            self.sample.set_orientation(ubconf)

        elif num_of_peaks == 3:
            ubconf = self._find_ub_from_three_peaks(peaks)
        elif num_of_peaks > 3:
            ubconf = self._find_ub_from_multiple_peaks(peaks)
        else:
            raise ValueError("Not enough peaks for UB matrix determination.")

    def _find_u_from_two_peaks(
        self,
        peaks: tuple[Peak],
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
        q_lab1 = TAS.q_lab(
            peak1.angles.two_theta,
            ei=peak1.ei,
            ef=peak1.ef,
        )
        q_lab2 = TAS.q_lab(
            peak2.angles.two_theta,
            ei=peak2.ei,
            ef=peak2.ef,
        )

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

        in_plane_ref = q_sample1

        return UBConf(
            peaks,
            u_mat,
            None,
            np.matmul(u_mat, b_mat),
            plane_normal,
            in_plane_ref,
        )

    # TODO
    def _find_ub_from_three_peaks(
        self,
        peaks: tuple[Peak],
    ) -> UBConf:
        """Find UB matrix from three observed peaks for a given goniomete"""
        ubconf = None
        return ubconf

    # TODO
    def _find_ub_from_multiple_peaks(
        self,
        peaks: tuple[Peak],
    ) -> UBConf:
        """Find UB matrix from more than three observed peaks for a given goniomete"""
        ubconf = None
        return ubconf
