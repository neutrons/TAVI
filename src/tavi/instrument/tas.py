# -*- coding: utf-8 -*-
from typing import Optional

import numpy as np

from tavi.instrument.tas_base import TASBase
from tavi.ub_algorithm import (
    find_u_from_one_peak_and_scattering_plane,
    find_u_from_two_peaks,
    find_ub_from_multiple_peaks,
    find_ub_from_three_peaks,
    plane_normal_from_one_peak,
    plane_normal_from_two_peaks,
    psi_from_hkle,
    q_lab,
    r_matrix_with_minimal_tilt,
    two_theta_from_hkle,
)
from tavi.utilities import MotorAngles, Peak, UBConf, spice_to_mantid


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
            if self.fixed_ei is not None:
                ei = self.fixed_ei  # fixed both Ef and Ei
            else:  # fixed Ef only
                ei = en + ef
        elif self.fixed_ei is not None:  # fixed Ei only
            ei = self.fixed_ei
            ef = ei - en
        else:
            raise ValueError(f"{self} should has either Ei or Ef fixed.")
        return (ei, ef)

    def get_two_theta(self, hkl: tuple[float, float, float], en: float = 0.0) -> Optional[float]:
        """find two theta angle for a given peak.

        Note:
            two theta is the angle between ki and kf
            the sign is determined by the sense of gomiometer

        Args:
            hkl (tuple): miller indice of a peak
            en (float): energy trnasfer en = ei - ef, in emV

        Returns:
            two_theta (float | None): two theta angle, in degrees. Reutrn None if can't reach.
        """

        ei, ef = self._get_ei_ef(en=en)
        two_theta_radian = two_theta_from_hkle(hkl, ei, ef, self.sample.b_mat)

        if two_theta_radian is None:
            print(f"Triangle cannot be closed at q={hkl}, en={en} meV.")
            return None
        else:
            two_theta = np.rad2deg(two_theta_radian) * self.goniometer._sense
            return two_theta

    def get_psi(self, hkl: tuple[float, float, float], en: float = 0.0) -> Optional[float]:
        """find psi angle for a given peak

        Note:
            psi is the angle between ki and Q = ki - kf
            psi will always has the opposite sign of two_theta

        Args:
            hkl (tuple): miller indice of a peak
            en (float): energy transfer en = ei -ef, in emV
        Returns:
            psi (float | None): psi angle, in degrees. Reutrn None if can't reach.
        """

        ei, ef = self._get_ei_ef(en=en)
        psi_radian = psi_from_hkle(hkl, ei, ef, self.sample.b_mat)

        if psi_radian is None:
            print(f"Triangle cannot be closed at q={hkl}, en={en} meV.")
            return None
        else:
            psi = np.rad2deg(psi_radian) * self.goniometer._sense * (-1)
            return psi

    def calculate_ub_matrix(self, peaks: tuple[Peak, ...], scattering_plane=None) -> Optional[UBConf]:
        """Find UB matrix from a list of observed peaks"""

        ei, ef = self._get_ei_ef()

        if isinstance(peaks, Peak):
            peaks = (peaks,)
            if scattering_plane is None:
                raise ValueError("Scattering plane cannot be None.")

        match (num_of_peaks := len(peaks)):
            case 1:
                b_mat = self.sample.b_mat
                u_mat = find_u_from_one_peak_and_scattering_plane(
                    peaks[0], scattering_plane, b_mat, self.goniometer.r_mat_inv, ei, ef
                )
                ub_mat = np.matmul(u_mat, b_mat)
                plane_normal, in_plane_ref = plane_normal_from_one_peak(
                    peaks[0].hkl,
                    peaks[0].angles,
                    self.goniometer.r_mat_inv,
                    ub_mat,
                )

                ub_conf = UBConf(
                    spice_convention=self.spice_convention,
                    ub_mat=ub_mat,
                    plane_normal=plane_normal,
                    in_plane_ref=in_plane_ref,
                    ub_peaks=peaks,
                )

            case 2:
                peak1, peak2 = peaks
                b_mat = self.sample.b_mat
                u_mat = find_u_from_two_peaks((peak1, peak2), b_mat, self.goniometer.r_mat_inv, ei, ef)
                plane_normal, in_plane_ref = plane_normal_from_two_peaks(u_mat, b_mat, peak1.hkl, peak2.hkl)
                ub_mat = np.matmul(u_mat, b_mat)

                ub_conf = UBConf(
                    spice_convention=self.spice_convention,
                    ub_mat=ub_mat,
                    plane_normal=plane_normal,
                    in_plane_ref=in_plane_ref,
                    ub_peaks=peaks,
                )

            case 3:
                peak1, peak2, peak3 = peaks
                ZERO = 1e-6
                if np.dot(peak1.hkl, np.cross(peak2.hkl, peak3.hkl)) < ZERO:
                    print("Cannot use three coplanar peaks to determine UB matrix.")
                    return None
                ub_mat = find_ub_from_three_peaks((peak1, peak2, peak3), self.goniometer.r_mat_inv, ei, ef)
                g_star_mat = np.matmul(ub_mat.T, ub_mat)
                self.sample.update_lattice_parametres_from_g_star_mat(g_star_mat)

                ub_conf = UBConf(spice_convention=self.spice_convention, ub_mat=ub_mat, ub_peaks=peaks)

            case _ if num_of_peaks > 3:
                ub_mat = find_ub_from_multiple_peaks(peaks, self.goniometer.r_mat_inv, ei, ef)
                g_star_mat = np.matmul(ub_mat.T, ub_mat)
                self.sample.update_lattice_parametres_from_g_star_mat(g_star_mat)

                ub_conf = UBConf(spice_convention=self.spice_convention, ub_mat=ub_mat, ub_peaks=peaks)

            case _:
                raise ValueError("Not enough peaks for UB matrix determination.")

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

        Return:
            Return MotorAngles
            Return None if the intented position is out of reach

        Note"
            Convert UB matrix, plane_normal, in_plane_ref to Mantind/International Table
            converntion before performing the calculation
        """

        if len(hkl) != 3:
            raise ValueError(f"hkl ={hkl} should have the format (h,k,l).")

        two_theta = self.get_two_theta(hkl=hkl, en=en)
        psi = self.get_psi(hkl=hkl, en=en)
        if (two_theta is None) or (psi is None):
            return None

        ei, ef = self._get_ei_ef(en=en)
        # checking if UB configuration info exists
        ub_conf = self.sample.ub_conf
        if ub_conf is None:
            raise ValueError(f"UB info not found. ub_conf={ub_conf}.")
        if (mat := ub_conf.ub_mat) is None:
            raise ValueError("UB matrix is not unknown.")

        goni_mode = self.goniometer.mode
        if goni_mode == "bisect":

            # convert to Mantid convention if needed
            if self.spice_convention:
                ub_conf_mantid = UBConf(ub_mat=spice_to_mantid(mat))

            angles = self.goniometer.angles_in_bisect_mode(hkl, two_theta, psi, ub_conf_mantid)

        elif goni_mode is None:  # default is minumal tilt
            if ((n := ub_conf.plane_normal) is None) or ((i := ub_conf.in_plane_ref) is None):
                raise ValueError(f"Missing UB info. plane_normal={n}, in_plne_ref={i}.")

            # convert to Mantid convention if needed
            if self.spice_convention:
                ub_conf_mantid = UBConf(
                    ub_mat=spice_to_mantid(mat),
                    plane_normal=spice_to_mantid(n),
                    in_plane_ref=spice_to_mantid(i),
                    ub_peaks=ub_conf.ub_peaks,
                )
            else:
                ub_conf_mantid = ub_conf

            r_mat = r_matrix_with_minimal_tilt(hkl, ei, ef, two_theta, ub_conf_mantid)
            angles = self.goniometer.angles_from_r_mat(r_mat, two_theta)

        return angles

    def calcvulate_hkl_from_angles(self, angles: MotorAngles) -> Optional[np.ndarray]:
        ei, ef = self._get_ei_ef()
        qlab = q_lab(ei=ei, ef=ef, theta=angles.two_theta)
        r = self.goniometer.r_mat(angles=angles)
        if (ub_conf := self.sample.ub_conf) is None:
            print("Cannot get hkl from motor angles without knowing the UB matrix.")
            return None

        ub = spice_to_mantid(ub_conf.ub_mat) if self.spice_convention else ub_conf.ub_mat
        rub = np.matmul(r, ub)
        rub_inv = np.linalg.inv(rub)
        hkl = rub_inv.dot(qlab) / (2 * np.pi)
        return hkl
