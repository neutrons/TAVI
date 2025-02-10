# -*- coding: utf-8 -*-
from typing import Optional

import numpy as np

from tavi.instrument.tas_base import TASBase
from tavi.instrument.ub_algorithm import find_u_from_two_peaks
from tavi.sample.xtal import Xtal
from tavi.utilities import MotorAngles, en2q, get_angle_from_triangle


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
    def q_lab(two_theta_deg: float, ei: float, ef: float):
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

    def get_two_theta(
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

        return np.rad2deg(two_theta_radian) * self.goniometer._sense

    def calculate_ub_matrix(self, peaks):
        """Find UB matrix from a list of observed peaks"""
        if not isinstance(self.sample, Xtal):
            raise ValueError("sample needs to be Xtal class for UB calculation.")

        match (num_of_peaks := len(peaks)):
            case 2:
                b_mat = self.sample.b_mat_from_lattice()
                ubconf = find_u_from_two_peaks(peaks, b_mat)
                self.sample.set_orientation(ubconf)

            case 3:
                pass
                # ubconf = self._find_ub_from_three_peaks(peaks)
                # self.sample.set_orientation(ubconf)
            case _ if num_of_peaks > 3:
                pass
                # ubconf = self._find_ub_from_multiple_peaks(peaks)
                # self.sample.set_orientation(ubconf)
            case _:
                raise ValueError("Not enough peaks for UB matrix determination.")

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

        two_theta_deg = np.rad2deg(two_theta) * self.goniometer._sense
        t_mat = self._t_mat_minimal_tilt(hkl)
        t_mat_inv = np.linalg.inv(t_mat)

        q_lab1 = TAS.q_lab(two_theta_deg, ei, ef) / q_norm
        q_lab2 = np.array([q_lab1[2], 0, -q_lab1[0]])
        q_lab3 = np.array([0, 1, 0])

        q_lab_mat = TAS.norm_mat(q_lab1, q_lab2, q_lab3)
        r_mat = q_lab_mat @ t_mat_inv

        _, omega, sgl, sgu, chi, phi = self.goniometer.angles_from_r_mat(r_mat)

        return MotorAngles(two_theta_deg, omega, sgl, sgu, chi, phi)
