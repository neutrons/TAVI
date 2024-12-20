from typing import Optional, Union

import numpy as np

from tavi.instrument.components.collimators import Collimators
from tavi.instrument.components.mono_ana import MonoAna
from tavi.instrument.resolution.ellipsoid import ResoEllipsoid
from tavi.instrument.tas import TAS
from tavi.utilities import (
    get_angle_bragg,
    get_angle_from_triangle,
    ksq2eng,
    rotation_matrix_2d,
    sig2fwhm,
    spice_to_mantid,
)


class CN(TAS):
    """Copper-Nathans method

    Methods:
        validate_instrument_parameters
        cooper_nathans

    """

    # 4 soller slits collimators
    NUM_COLLS = 4
    IDX_COLL0_H, IDX_COLL0_V = 0, 2
    IDX_COLL1_H, IDX_COLL1_V = 1, 3
    IDX_COLL2_H, IDX_COLL2_V = 4, 6
    IDX_COLL3_H, IDX_COLL3_V = 5, 7
    # 1 monochromator and 1 analyzer
    NUM_MONOS, NUM_ANAS = 1, 1
    IDX_MONO0_H, IDX_MONO0_V = 0, 1
    IDX_ANA0_H, IDX_ANA0_V = 2, 3

    def __init__(self, SPICE_CONVENTION: bool = True) -> None:
        """Load instrument configuration from json if provided"""
        super().__init__(SPICE_CONVENTION=SPICE_CONVENTION)

        # constants independent of q and eng
        self._mat_f: Optional[np.ndarray] = None
        self._mat_g: Optional[np.ndarray] = None

    @staticmethod
    def calc_mat_f(mono: MonoAna, ana: MonoAna) -> np.ndarray:
        # matrix F, divergence of monochromator and analyzer, [pop75] Appendix 1
        mat_f = np.zeros(((CN.NUM_MONOS + CN.NUM_ANAS) * 2, (CN.NUM_MONOS + CN.NUM_ANAS) * 2))
        mat_f[CN.IDX_MONO0_H, CN.IDX_MONO0_H] = 1.0 / mono._mosaic_h**2
        mat_f[CN.IDX_MONO0_V, CN.IDX_MONO0_V] = 1.0 / mono._mosaic_v**2
        mat_f[CN.IDX_ANA0_H, CN.IDX_ANA0_H] = 1.0 / ana._mosaic_h**2
        mat_f[CN.IDX_ANA0_V, CN.IDX_ANA0_V] = 1.0 / ana._mosaic_v**2
        return mat_f

    @staticmethod
    def calc_mat_g(coll: Collimators):
        (
            coll_h_pre_mono,
            coll_h_pre_sample,
            coll_h_post_sample,
            coll_h_post_ana,
        ) = coll._horizontal_divergence

        (
            coll_v_pre_mono,
            coll_v_pre_sample,
            coll_v_post_sample,
            coll_v_post_ana,
        ) = coll._vertical_divergence

        # matrix G, divergence of collimators, [pop75] Appendix 1
        mat_g = np.zeros((CN.NUM_COLLS * 2, CN.NUM_COLLS * 2))
        mat_g[CN.IDX_COLL0_H, CN.IDX_COLL0_H] = 1.0 / coll_h_pre_mono**2
        mat_g[CN.IDX_COLL0_V, CN.IDX_COLL0_V] = 1.0 / coll_v_pre_mono**2
        mat_g[CN.IDX_COLL1_H, CN.IDX_COLL1_H] = 1.0 / coll_h_pre_sample**2
        mat_g[CN.IDX_COLL1_V, CN.IDX_COLL1_V] = 1.0 / coll_v_pre_sample**2
        mat_g[CN.IDX_COLL2_H, CN.IDX_COLL2_H] = 1.0 / coll_h_post_sample**2
        mat_g[CN.IDX_COLL2_V, CN.IDX_COLL2_V] = 1.0 / coll_v_post_sample**2
        mat_g[CN.IDX_COLL3_H, CN.IDX_COLL3_H] = 1.0 / coll_h_post_ana**2
        mat_g[CN.IDX_COLL3_V, CN.IDX_COLL3_V] = 1.0 / coll_v_post_ana**2

        return mat_g

    @staticmethod
    def calc_mat_a(ki, kf, theta_m, theta_a):
        """matrix A,Y=AU, tranform from collimators angular divergence to  ki-kf frame"""
        mat_a = np.zeros((6, 2 * CN.NUM_COLLS))
        mat_a[0, CN.IDX_COLL0_H] = 0.5 * ki / np.tan(theta_m)
        mat_a[0, CN.IDX_COLL1_H] = -0.5 * ki / np.tan(theta_m)
        mat_a[1, CN.IDX_COLL1_H] = ki
        mat_a[2, CN.IDX_COLL1_V] = -ki

        mat_a[3, CN.IDX_COLL2_H] = 0.5 * kf / np.tan(theta_a)
        mat_a[3, CN.IDX_COLL3_H] = -0.5 * kf / np.tan(theta_a)
        mat_a[4, CN.IDX_COLL2_H] = kf
        mat_a[5, CN.IDX_COLL2_V] = kf
        return mat_a

    @staticmethod
    def calc_mat_b(ki, kf, phi, two_theta):
        """matrix B, X=BY, transform from ki-kf frame to momentum transfer q-frame"""
        mat_b = np.zeros((4, 6))
        mat_b[0:3, 0:3] = rotation_matrix_2d(phi)
        mat_b[0:3, 3:6] = rotation_matrix_2d(phi - two_theta) * (-1)
        mat_b[3, 0] = 2 * ksq2eng * ki
        mat_b[3, 3] = -2 * ksq2eng * kf
        return mat_b

    @staticmethod
    def calc_mat_c(theta_m, theta_a):
        """matrix C, constrinat between mono/ana mosaic and collimator divergence"""
        mat_c = np.zeros(((CN.NUM_MONOS + CN.NUM_ANAS) * 2, CN.NUM_COLLS * 2))
        mat_c[CN.IDX_MONO0_H, CN.IDX_COLL0_H] = 0.5
        mat_c[CN.IDX_MONO0_H, CN.IDX_COLL1_H] = 0.5
        mat_c[CN.IDX_MONO0_V, CN.IDX_COLL0_V] = 0.5 / np.sin(theta_m)
        mat_c[CN.IDX_MONO0_V, CN.IDX_COLL1_V] = -0.5 / np.sin(theta_m)

        mat_c[CN.IDX_ANA0_H, CN.IDX_COLL2_H] = 0.5
        mat_c[CN.IDX_ANA0_H, CN.IDX_COLL3_H] = 0.5
        mat_c[CN.IDX_ANA0_V, CN.IDX_COLL2_V] = 0.5 / np.sin(theta_a)
        mat_c[CN.IDX_ANA0_V, CN.IDX_COLL3_V] = -0.5 / np.sin(theta_a)
        return mat_c

    @staticmethod
    def _generate_hkle_list(
        hkl_list: Union[tuple, list[tuple]],
        ei: Union[float, list[float]],
        ef: Union[float, list[float]],
    ) -> list[tuple[tuple[float, float, float], float, float]]:
        """Generate a list containing tuple ((h, k, l), ei, ef)"""
        hkle_list = []
        if not isinstance(ei, list):
            ei = [ei]
        if not isinstance(ef, list):
            ef = [ef]
        if not isinstance(hkl_list, list):
            hkl_list = [hkl_list]

        if isinstance(hkl_list, list):
            for one_ei in ei:
                for one_ef in ef:
                    for hkl in hkl_list:
                        if isinstance(hkl, tuple) and len(hkl) == 3:
                            hkle_list.append((hkl, one_ei, one_ef))
                        else:
                            raise ValueError(f"hkl={hkl} is not a tuple of length 3.")
        return hkle_list

    def validate_instrument_parameters(self):
        """Check if enough instrument parameters are provided for Cooper-Nathans mehtod"""

        try:  # monochromator
            mono = self.monochromator
        except AttributeError:
            print("Monochromator info are missing.")

        if None in (mono_mosaic := (mono.mosaic_h, mono.mosaic_v)):
            raise ValueError("Mosaic of monochromator is missing.")
        elif not all(val > 0 for val in mono_mosaic):
            raise ValueError("Mosaic of monochromator cannot be negative.")

        try:  # analyzer
            ana = self.analyzer
        except AttributeError:
            print("Analyzer info are missing.")

        if None in (ana_mosaic := (ana.mosaic_h, ana.mosaic_v)):
            raise ValueError("Mosaic of analyzer is missing.")
        elif not all(val > 0 for val in ana_mosaic):
            raise ValueError("Mosaic of analyzer cannot be negative.")

        # collimators
        if (coll := self.collimators) is None:
            raise ValueError("Collimators info are missing.")
        elif not all(val > 0 for val in coll.horizontal_divergence):
            raise ValueError("Horizontal divergence of collimators cannot be negative.")
        elif not all(val > 0 for val in coll.vertical_divergence):
            raise ValueError("Vertical divergence of collimators cannot be negative.")

        # sample
        # if self.sample is None:
        #     raise ValueError("Sample info are missing.")

    def cooper_nathans(
        self,
        hkl_list: Union[tuple[float, float, float], list[tuple[float, float, float]]],
        ei: Union[float, list[float]],
        ef: Union[float, list[float]],
        projection: tuple = ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        R0: bool = False,
    ) -> Union[ResoEllipsoid, list[ResoEllipsoid]]:
        """Calculate resolution using Cooper-Nathans method

        Args:
            hkl: momentum transfer, miller indices in reciprocal lattice
            ei: incident energy, in units of meV
            ef: final energy, in units of meV
            R0: calculate normalization factor if True
            projection (tuple): three non-coplaner vectors. If projection is None, the calculation is done in local Q frame

        """
        self.validate_instrument_parameters()

        self._mat_f = CN.calc_mat_f(self.monochromator, self.analyzer)
        self._mat_g = CN.calc_mat_g(self.collimators)

        hkle_list = self._generate_hkle_list(hkl_list, ei, ef)
        rez_list = []
        for hkl, ei, ef in hkle_list:

            # q_lab = conv_mat @ hkl
            # q_mod = np.linalg.norm(q_lab)
            q_mod = np.linalg.norm(self.sample.b_mat @ hkl) * 2 * np.pi

            ki = np.sqrt(ei / ksq2eng)
            kf = np.sqrt(ef / ksq2eng)

            rez = ResoEllipsoid(hkle=hkl + (ei - ef,), projection=projection, sample=self.sample)

            try:
                two_theta = get_angle_from_triangle(ki, kf, q_mod) * self.goniometer._sense
            except TypeError:
                rez.STATUS = False
                print(f"Cannot close triangle for ei={ei}, ef={ef}, hkl={np.round(hkl,3)}.")

                rez._set_labels()
                rez_list.append(rez)
                continue

            # phi = <ki to q>, always has the oppositie sign of s2
            phi = get_angle_from_triangle(ki, q_mod, kf) * self.goniometer._sense * (-1)

            theta_m = get_angle_bragg(ki, self.monochromator.d_spacing) * self.monochromator._sense
            theta_a = get_angle_bragg(kf, self.analyzer.d_spacing) * self.analyzer._sense

            # TODO
            # curved monochromator and analyzer

            # TODO
            # reflection efficiency

            mat_a = CN.calc_mat_a(ki, kf, theta_m, theta_a)
            mat_b = CN.calc_mat_b(ki, kf, phi, two_theta)
            mat_c = CN.calc_mat_c(theta_m, theta_a)

            mat_h = mat_c.T @ self._mat_f @ mat_c + self._mat_g
            mat_h_inv = np.linalg.inv(mat_h)
            mat_ba = mat_b @ mat_a
            mat_cov = mat_ba @ mat_h_inv @ mat_ba.T

            # TODO how to add smaple mosaic in cooper-nathans?
            mat_cov[1, 1] += q_mod**2 * self.sample._mosaic_h**2
            mat_cov[2, 2] += q_mod**2 * self.sample._mosaic_v**2

            mat_reso = np.linalg.inv(mat_cov) * sig2fwhm**2

            motor_angles = self.calculate_motor_angles(peak=hkl, ei=ei, ef=ef)
            if motor_angles is None:
                rez.STATUS = False
                rez._set_labels()
                rez_list.append(rez)
                continue

            r_mat = self.goniometer.r_mat(motor_angles)
            ub_mat = self.sample.ub_mat

            if self.SPICE_CONVENTION:
                ub_mat = spice_to_mantid(ub_mat)

            conv_mat = 2 * np.pi * np.matmul(r_mat, ub_mat)
            rez._project_to_frame(mat_reso, phi, conv_mat)

            # TODO check normalization factor
            # -------------------------------------------------------------------------
            # - if the instruments works in kf=const mode and the scans are counted for
            #   or normalised to monitor counts no ki^3 or kf^3 factor is needed.
            # - if the instrument works in ki=const mode the kf^3 factor is needed.

            if R0:  # calculate
                r0 = np.pi**2 / 4 / np.sin(theta_m) / np.sin(theta_a)
                r0 *= np.sqrt(np.linalg.det(self._mat_f) / np.linalg.det(mat_h))
            else:
                r0 = 0

            rez.r0 = r0

            if np.isnan(rez.r0) or np.isinf(rez.r0) or np.isnan(rez.mat.any()) or np.isinf(rez.mat.any()):
                rez.STATUS = False
            else:
                rez.STATUS = True

            rez._set_labels()
            rez_list.append(rez)

        if len(rez_list) == 1:
            return rez_list[0]
        return rez_list
