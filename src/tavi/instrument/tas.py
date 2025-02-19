# -*- coding: utf-8 -*-
from typing import Optional

import numpy as np

from tavi.instrument.tas_base import TASBase
from tavi.ub_algorithm import (
    find_u_from_two_peaks,
    find_ub_from_multiple_peaks,
    find_ub_from_three_peaks,
    r_matrix_with_minimal_tilt,
    two_theta_from_hkle,
)
from tavi.utilities import MotorAngles, Peak, UBConf, mantid_to_spice


class TAS(TASBase):
    """
    Triple-axis instrument class. Handles angle and UB calculations

    Attributes:
        spice_convention (bool): False if using Mandid/International Crystallography Table convention
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

    def _get_ei_ef(
        self, ei: Optional[float] = None, ef: Optional[float] = None, en: float = 0.0
    ) -> tuple[float, float]:
        """Determine Ei and Ef based on if Ei or Ef is fixed"""
        if self.fixed_ef is not None:
            ef = self.fixed_ef
            ei = en + ef
        elif self.fixed_ei is not None:
            ei = self.fixed_ei
            ef = ei - en
        else:
            raise ValueError(f"{self} should has either Ei or Ef fixed.")
        return (ei, ef)

    def get_two_theta(self, hkl: tuple[float, float, float], en: float = 0.0) -> Optional[float]:
        """find two theta angle for a given peak

        Args:
            hkl (tuple): miller indice of a peak
            ei (float): incident neutron energy, in emV
            ef (float | None): final neutron energy, in meV. Use ei if not given
        Returns:
            two_theta (float | None): two theta angle, in degrees. Reutrn None if can't reach.
        """

        ei, ef = self._get_ei_ef(en=en)
        two_theta_radian = two_theta_from_hkle(hkl, ei, ef, self.sample.b_mat)

        if two_theta_radian is None:
            print(f"Triangle cannot be closed at q={hkl}, en={en} meV.")
            return None
        # elif np.rad2deg(two_theta_radian) < S2_MIN_DEG:
        #   pass
        else:
            two_theta = np.rad2deg(two_theta_radian) * self.goniometer._sense
            return two_theta

    def calculate_ub_matrix(self, peaks: tuple[Peak]) -> UBConf:
        """Find UB matrix from a list of observed peaks"""

        ei, ef = self._get_ei_ef()

        match (num_of_peaks := len(peaks)):
            case 2:
                b_mat = self.sample.b_mat
                u_mat, plane_normal, in_plane_ref = find_u_from_two_peaks(
                    peaks, b_mat, self.goniometer.r_mat_inv, ei, ef
                )
                ub_mat = np.matmul(u_mat, b_mat)

            case 3:
                (u_mat, b_mat, ub_mat, plane_normal, in_plane_ref) = find_ub_from_three_peaks(
                    peaks, self.goniometer.r_mat_inv, ei, ef
                )
                self.sample.update_lattice_parametres_from_b_mat(b_mat)

            case _ if num_of_peaks > 3:
                (u_mat, b_mat, ub_mat, plane_normal, in_plane_ref) = find_ub_from_multiple_peaks(
                    peaks, self.goniometer.r_mat_inv, ei, ef
                )
                self.sample.update_lattice_parametres_from_b_mat(b_mat)

            case _:
                raise ValueError("Not enough peaks for UB matrix determination.")

        # suffle the order following SPICE convention
        if self.spice_convention:
            plane_normal = mantid_to_spice(plane_normal)
            in_plane_ref = mantid_to_spice(in_plane_ref)
            ub_mat = mantid_to_spice(ub_mat)

        ub_conf = UBConf(ub_mat, plane_normal, in_plane_ref, u_mat, b_mat, peaks)
        self.sample.ub_conf = ub_conf
        return ub_conf

    def calculate_motor_angles(
        self,
        hkl: tuple[float, float, float],
        en: float = 0.0,
    ) -> Optional[MotorAngles]:
        """calculate motor positions for a given peak if UB matrix has been determined

        Args:
            peak (tuple): Miller indice (h, k, l) of a peak
            en (float): energy transfer, in meV. en = ei - ef
        """

        if len(hkl) != 3:
            raise ValueError(f"hkl ={hkl} should have the format (h,k,l).")

        ei, ef = self._get_ei_ef(en=en)
        b_mat = self.sample.b_mat

        two_theta_radian = two_theta_from_hkle(hkl, ei, ef, b_mat)

        if two_theta_radian is None:
            print(f"Position hkl={hkl}, en={en} can't be reached.")
            return None

        two_theta = np.rad2deg(two_theta_radian * self.goniometer._sense)
        ub_conf = self.sample.ub_conf
        if ub_conf is None:
            raise ValueError(f"UB info not found. ub_conf={ub_conf}.")
        if (
            ((mat := ub_conf.ub_mat) is None)
            or ((n := ub_conf.plane_normal) is None)
            or ((i := ub_conf.in_plane_ref) is None)
        ):
            raise ValueError(f"Missing UB info. ub_mat={mat}, plane_normal={n}, in_plne_ref={i}.")

        r_mat = r_matrix_with_minimal_tilt(hkl, ei, ef, two_theta, ub_conf)
        angles = self.goniometer.angles_from_r_mat(r_mat, two_theta)

        return angles
