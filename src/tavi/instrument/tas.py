# -*- coding: utf-8 -*-
from typing import Literal, Optional, Union

import numpy as np

from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.instrument.resolution.ellipsoid import ResoEllipsoid
from tavi.instrument.tas_base import TASBase
from tavi.ub_algorithm import (
    UBConf,
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
from tavi.utilities import MotorAngles, Peak


class TAS(TASBase):
    """
    Triple-axis instrument class. Handles angle and UB calculations

    Attributes:
        convention (str): "Mantid" or "Spice".
        fixed_ei (float | None): set if ei is fixed
        fixed_ef (float | None): set if ef is fixed

    Methods:
        calculate_two_theta
        calculate_ub_matrix

    Note:
        Definition of different conventions can be found in the UBConf class

    """

    def __init__(
        self,
        convention: Literal["Mantid", "Spice"] = "Spice",
        fixed_ei: Optional[float] = None,
        fixed_ef: Optional[float] = None,
    ):
        super().__init__()
        self.convention = convention  # use coordination system defined in SPICE
        self.fixed_ei = fixed_ei
        self.fixed_ef = fixed_ef

    def __repr__(self):
        cls = self.__class__.__name__
        cls_str = f"{cls}(fixed_ei={self.fixed_ei!r}, fixed_ef={self.fixed_ef!r}, convention={self.convention})"
        return cls_str

    def _get_ei_ef(
        self, ei: Optional[float] = None, ef: Optional[float] = None, en: float = 0.0
    ) -> tuple[float, float]:
        """Determine Ei and Ef based on if Ei or Ef is fixed"""
        if self.fixed_ef is not None:
            ef = self.fixed_ef
            if self.fixed_ei is not None:
                ei = self.fixed_ei  # fixed both Ef and Ei
                if not np.allclose(en, 0.0, atol=1e-6):
                    raise ValueError(f"{self} has both Ei and Ef fixed, No energy trnsfer allowed.")
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
                    convention=self.convention,
                    b_mat=b_mat,
                    _ub_mat=ub_mat,
                    _plane_normal=plane_normal,
                    _in_plane_ref=in_plane_ref,
                    ub_peaks=peaks,
                )

            case 2:
                peak1, peak2 = peaks
                b_mat = self.sample.b_mat
                u_mat = find_u_from_two_peaks((peak1, peak2), b_mat, self.goniometer.r_mat_inv, ei, ef)
                plane_normal, in_plane_ref = plane_normal_from_two_peaks(u_mat, b_mat, peak1.hkl, peak2.hkl)
                ub_mat = np.matmul(u_mat, b_mat)

                ub_conf = UBConf(
                    convention=self.convention,
                    b_mat=b_mat,
                    _ub_mat=ub_mat,
                    _plane_normal=plane_normal,
                    _in_plane_ref=in_plane_ref,
                    ub_peaks=peaks,
                )

            case 3:  # Three non-colinear peaks
                peak1, peak2, peak3 = peaks
                ZERO = 1e-6
                if np.dot(peak1.hkl, np.cross(peak2.hkl, peak3.hkl)) < ZERO:
                    print("Cannot use three coplanar peaks to determine UB matrix.")
                    return None
                ub_mat = find_ub_from_three_peaks((peak1, peak2, peak3), self.goniometer.r_mat_inv, ei, ef)
                g_star_mat = np.matmul(ub_mat.T, ub_mat)
                self.sample.update_lattice_parametres_from_g_star_mat(g_star_mat)

                ub_conf = UBConf(
                    convention=self.convention,
                    b_mat=self.sample.b_mat,
                    _ub_mat=ub_mat,
                    ub_peaks=peaks,
                )

            case _ if num_of_peaks > 3:
                ub_mat = find_ub_from_multiple_peaks(peaks, self.goniometer.r_mat_inv, ei, ef)
                g_star_mat = np.matmul(ub_mat.T, ub_mat)
                self.sample.update_lattice_parametres_from_g_star_mat(g_star_mat)

                ub_conf = UBConf(
                    convention=self.convention,
                    b_mat=self.sample.b_mat,
                    _ub_mat=ub_mat,
                    ub_peaks=peaks,
                )

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
            angles = self.goniometer.angles_in_bisect_mode(hkl, two_theta, psi, ub_conf)

        elif goni_mode is None:  # default is minumal tilt
            if ((n := ub_conf.plane_normal) is None) or ((i := ub_conf.in_plane_ref) is None):
                raise ValueError(f"Missing UB info. plane_normal={n}, in_plne_ref={i}.")

            r_mat = r_matrix_with_minimal_tilt(hkl, ei, ef, two_theta, ub_conf)
            angles = self.goniometer.angles_from_r_mat(r_mat, two_theta)

        return angles

    def calculate_hkl_from_angles(self, angles: MotorAngles) -> Optional[np.ndarray]:
        ei, ef = self._get_ei_ef()
        qlab = q_lab(ei=ei, ef=ef, theta=angles.two_theta)
        r = self.goniometer.r_mat(angles=angles)
        if (ub_conf := self.sample.ub_conf) is None:
            print("Cannot get hkl from motor angles without knowing the UB matrix.")
            return None

        rub = np.matmul(r, ub_conf._ub_mat)
        rub_inv = np.linalg.inv(rub)
        hkl = rub_inv.dot(qlab) / (2 * np.pi)
        return hkl

    def cooper_nathans(
        self,
        hkl: Union[tuple[float, float, float], list[tuple[float, float, float]]],
        en: Union[float, list[float]],
        projection: tuple = ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        R0: bool = False,
    ) -> Union[ResoEllipsoid, list[ResoEllipsoid]]:
        """Calculated resolution ellipsoid at given (h,k,l,e) position for given projection

        Args:
            hkl: momentum transfer, miller indices in reciprocal lattice
            ei: incident energy, in units of meV
            ef: final energy, in units of meV
            R0: calculate normalization factor if True
            projection (tuple): three non-coplaner vectors. If projection is None, the calculation is done in local Q frame
        Return
            A ResoEllipdois instance or a list of ResoEllipdois instances
        """
        cn = CooperNathans(instrument=self)
        return cn.calculate(hkl=hkl, en=en, projection=projection, R0=R0)
