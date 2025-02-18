# -*- coding: utf-8 -*-
from typing import Optional

import numpy as np

from tavi.instrument.tas_base import TASBase

# from tavi.sample.xtal import Xtal
from tavi.ub_algorithm import find_u_from_two_peaks, two_theta_from_hkl
from tavi.utilities import MotorAngles, Peak, en2q, get_angle_from_triangle, mantid_to_spice


class TAS(TASBase):
    """
    Triple-axis instrument class. Handles angle and UB calculations

    Attributes:
        s2_limit (tuple | None): (min, max) of s2 angle
        sense ("+" | "-"): "+" if s2 is right-hand
        fixed_ei (float | None): set if ei is fixed
        fixed_ef (float | None): set if ef is fixed

    Methods:
        calculate_two_theta
        calculate_ub_matrix

    """

    def __init__(
        self,
        spice_convention: bool = True,
        fixed_ei: Optional[float] = None,
        fixed_ef: Optional[float] = None,
    ):
        super().__init__()
        self.spice_convention = spice_convention  # use coordination system defined in SPICE
        self.fixed_ei = fixed_ei
        self.fixed_ef = fixed_ef

    def __repr__(self):
        cls = self.__class__.__name__
        cls_str = (
            f"{cls}(fixed_ei={self.fixed_ei!r}, fixed_ef={self.fixed_ef!r}, spice_convention={self.spice_convention})"
        )
        return cls_str

    def get_two_theta(
        self,
        hkl: tuple[float, float, float],
        ei: Optional[float] = None,
        ef: Optional[float] = None,
    ) -> Optional[float]:
        """find two theta angle for a given peak

        Args:
            hkl (tuple): miller indice of a peak
            ei (float): incident neutron energy, in emV
            ef (float | None): final neutron energy, in meV. Use ei if not given
        Returns:
            two_theta (float | None): two theta angle, in degrees. Reutrn None if can't reach.
        """

        if ef is None:
            if self.fixed_ef is None:
                raise ValueError("Cannot calculate two thetha without knowing final energy ef.")
            else:
                ef = self.fixed_ef
        if ei is None:
            if self.fixed_ei is None:
                ei = ef  # assuming ei and ef are the same
            else:
                ei = self.fixed_ei

        two_theta_radian = two_theta_from_hkl(hkl, ei, ef, self.sample.b_mat)

        if two_theta_radian is None:
            print(f"Triangle cannot be closed at q={hkl}, en={ei-ef} meV.")
            return None
        # elif np.rad2deg(two_theta_radian) < S2_MIN_DEG:
        # pass
        else:
            two_theta = np.rad2deg(two_theta_radian) * self.goniometer._sense
            return two_theta

    def calculate_ub_matrix(self, peaks: tuple[Peak, ...]):
        """Find UB matrix from a list of observed peaks"""

        if (ef := self.fixed_ef) is None:
            if (ei := self.fixed_ei) is None:
                raise ValueError("Ei and Ef needs to be provided to calculate q_lab.")
            else:  # fixed Ei
                ef = ei
        else:  # fixed Ef
            ei = ef

        match (num_of_peaks := len(peaks)):
            case 2:
                b_mat = self.sample.b_mat
                u_mat = find_u_from_two_peaks(peaks, b_mat, self.goniometer.r_mat_inv, ei, ef)
                ub_mat = np.matmul(u_mat, b_mat)

            case 3:
                pass

            case _ if num_of_peaks > 3:
                pass

            case _:
                raise ValueError("Not enough peaks for UB matrix determination.")

        # suffle the order following SPICE convention
        if self.spice_convention:
            # plane_normal = mantid_to_spice(plane_normal)
            # in_plane_ref = mantid_to_spice(in_plane_ref)
            ub_mat = mantid_to_spice(ub_mat)

        # ubconf = UBConf(peaks, u_mat, None, ub_mat, plane_normal, in_plane_ref)
        # self.sample.set_orientation(ubconf)

    # TODO
    def calculate_motor_angles(self, peak: tuple[float, float, float], en: float):
        """calculate motor positions for a given peak if UB matrix has been determined

        Args:
            peak (tuple): Miller indice (h, k, l) of a peak
            en (float): energy transfer, in meV. en = ei - ef
        """

        try:
            h, k, l = peak
        except ValueError:
            print(f"hkl ={peak} should have the format (h,k,l).")
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
