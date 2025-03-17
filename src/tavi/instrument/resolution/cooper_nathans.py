from typing import Union

import numpy as np

from tavi.instrument.resolution.ellipsoid import ResoEllipsoid
from tavi.instrument.resolution.resolution_calculator import ResolutionCalculator
from tavi.utilities import get_angle_bragg, get_angle_from_triangle, ksq2eng, rotation_matrix_2d, sig2fwhm


class CooperNathans(ResolutionCalculator):
    """Cooper-Nathans method

    Note:
        [pop75] Popovici, Acta Cryst. (1975). A31, 507
    """

    # 1 monochromator and 1 analyzer
    NUM_MONOS = 1
    IDX_MONO0_H, IDX_MONO0_V = 0, 1
    NUM_ANAS = 1
    IDX_ANA0_H, IDX_ANA0_V = 2, 3
    # 4 soller slits collimators
    NUM_COLLS = 4
    IDX_COLL0_H, IDX_COLL0_V = 0, 2
    IDX_COLL1_H, IDX_COLL1_V = 1, 3
    IDX_COLL2_H, IDX_COLL2_V = 4, 6
    IDX_COLL3_H, IDX_COLL3_V = 5, 7

    def mat_f(self) -> np.ndarray:
        """matrix F, divergence of monochromator and analyzer, [pop75] Appendix 1
        Note: No conversion between sigma and FWHM
        """
        cls = type(self)
        monochromator = self.instrument.monochromator
        analyzer = self.instrument.analyzer
        i, j = cls.NUM_MONOS, cls.NUM_ANAS
        mat_f = np.zeros(((i + j) * 2, (i + j) * 2))
        mat_f[cls.IDX_MONO0_H, cls.IDX_MONO0_H] = 1.0 / monochromator._mosaic_h**2
        mat_f[cls.IDX_MONO0_V, cls.IDX_MONO0_V] = 1.0 / monochromator._mosaic_v**2
        mat_f[cls.IDX_ANA0_H, cls.IDX_ANA0_H] = 1.0 / analyzer._mosaic_h**2
        mat_f[cls.IDX_ANA0_V, cls.IDX_ANA0_V] = 1.0 / analyzer._mosaic_v**2
        return mat_f

    def mat_g(self) -> np.ndarray:
        """matrix G, divergence of monochromator and analyzer, [pop75] Appendix 1
        Note: No conversion between sigma and FWHM
        """
        cls = type(self)
        collimators = self.instrument.collimators
        (h_pre_mono, h_pre_sample, h_post_sample, h_post_ana) = collimators._horizontal_divergence
        (v_pre_mono, v_pre_sample, v_post_sample, v_post_ana) = collimators._vertical_divergence

        # matrix G, divergence of collimators, [pop75] Appendix 1
        mat_g = np.zeros((cls.NUM_COLLS * 2, cls.NUM_COLLS * 2))
        mat_g[cls.IDX_COLL0_H, cls.IDX_COLL0_H] = 1.0 / h_pre_mono**2
        mat_g[cls.IDX_COLL0_V, cls.IDX_COLL0_V] = 1.0 / v_pre_mono**2
        mat_g[cls.IDX_COLL1_H, cls.IDX_COLL1_H] = 1.0 / h_pre_sample**2
        mat_g[cls.IDX_COLL1_V, cls.IDX_COLL1_V] = 1.0 / v_pre_sample**2
        mat_g[cls.IDX_COLL2_H, cls.IDX_COLL2_H] = 1.0 / h_post_sample**2
        mat_g[cls.IDX_COLL2_V, cls.IDX_COLL2_V] = 1.0 / v_post_sample**2
        mat_g[cls.IDX_COLL3_H, cls.IDX_COLL3_H] = 1.0 / h_post_ana**2
        mat_g[cls.IDX_COLL3_V, cls.IDX_COLL3_V] = 1.0 / v_post_ana**2

        return mat_g

    @classmethod
    def calc_mat_a(cls, ki, kf, theta_m, theta_a):
        """matrix A,Y=AU, tranform from collimators angular divergence to  ki-kf frame"""
        mat_a = np.zeros((6, 2 * cls.NUM_COLLS))
        mat_a[0, cls.IDX_COLL0_H] = 0.5 * ki / np.tan(theta_m)
        mat_a[0, cls.IDX_COLL1_H] = -0.5 * ki / np.tan(theta_m)
        mat_a[1, cls.IDX_COLL1_H] = ki
        mat_a[2, cls.IDX_COLL1_V] = ki  # negative in Takin

        mat_a[3, cls.IDX_COLL2_H] = 0.5 * kf / np.tan(theta_a)
        mat_a[3, cls.IDX_COLL3_H] = -0.5 * kf / np.tan(theta_a)
        mat_a[4, cls.IDX_COLL2_H] = kf
        mat_a[5, cls.IDX_COLL2_V] = kf
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

    @classmethod
    def calc_mat_c(cls, theta_m, theta_a):
        """matrix C, constrinat between mono/ana mosaic and collimator divergence"""
        mat_c = np.zeros(((cls.NUM_MONOS + cls.NUM_ANAS) * 2, cls.NUM_COLLS * 2))
        mat_c[cls.IDX_MONO0_H, cls.IDX_COLL0_H] = 0.5
        mat_c[cls.IDX_MONO0_H, cls.IDX_COLL1_H] = 0.5
        mat_c[cls.IDX_MONO0_V, cls.IDX_COLL0_V] = 0.5 / np.sin(theta_m)
        mat_c[cls.IDX_MONO0_V, cls.IDX_COLL1_V] = -0.5 / np.sin(theta_m)

        mat_c[cls.IDX_ANA0_H, cls.IDX_COLL2_H] = 0.5
        mat_c[cls.IDX_ANA0_H, cls.IDX_COLL3_H] = 0.5
        mat_c[cls.IDX_ANA0_V, cls.IDX_COLL2_V] = 0.5 / np.sin(theta_a)
        mat_c[cls.IDX_ANA0_V, cls.IDX_COLL3_V] = -0.5 / np.sin(theta_a)
        return mat_c

    def validate_instrument_parameters(self):
        base = self.__class__.__base__
        base.validate_instrument_parameters(self)

    def calculate(
        self,
        hkl: Union[tuple[float, float, float], list[tuple[float, float, float]]],
        en: Union[float, list[float]],
        projection: tuple = ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        R0: bool = False,
    ):
        """Calculate resolution using Cooper-Nathans method

        Args:
            hkl (tuple | list(tuple)): momentum transfer, miller indices in reciprocal lattice
            en (float | list(float)): en = ei - ef is energy trnsfer, in units of meV
            projection (tuple): three non-coplaner vectors. If projection is None,
                                the calculation is done in local Q frame (Q_para, Q_perp, Q_up)
            R0: calculate normalization factor if True
        """

        cls = type(self)

        self.validate_instrument_parameters()
        hkle_list = self.generate_hkle_list(hkl, en)

        instru = self.instrument

        rez_list = []
        for hkl, ei, ef in hkle_list:
            # q_lab = conv_mat @ hkl
            # q_mod = np.linalg.norm(q_lab)
            q_mod = np.linalg.norm(instru.sample.b_mat @ hkl) * 2 * np.pi

            ki = np.sqrt(ei / ksq2eng)
            kf = np.sqrt(ef / ksq2eng)

            try:
                two_theta = get_angle_from_triangle(ki, kf, q_mod) * instru.goniometer._sense
            except TypeError:
                print(f"Cannot close triangle for ei={ei}, ef={ef}, hkl={np.round(hkl,3)}.")
                rez = ResoEllipsoid(hkle=hkl + (ei - ef,), projection=projection, sample=instru.sample)
                rez.STATUS = False
                rez._set_labels()
                rez_list.append(rez)
                continue

            # phi = <ki to q>, always has the oppositie sign of s2
            phi = get_angle_from_triangle(ki, q_mod, kf) * instru.goniometer._sense * (-1)

            theta_m = get_angle_bragg(ki, instru.monochromator.d_spacing) * instru.monochromator._sense
            theta_a = get_angle_bragg(kf, instru.analyzer.d_spacing) * instru.analyzer._sense

            # TODO
            # curved monochromator and analyzer

            # TODO
            # reflection efficiency

            mat_a = cls.calc_mat_a(ki, kf, theta_m, theta_a)
            mat_b = cls.calc_mat_b(ki, kf, phi, two_theta)
            mat_c = cls.calc_mat_c(theta_m, theta_a)

            mat_h = mat_c.T @ self.mat_f() @ mat_c + self.mat_g()
            mat_h_inv = np.linalg.inv(mat_h)
            mat_ba = mat_b @ mat_a
            mat_cov = mat_ba @ mat_h_inv @ mat_ba.T

            # TODO how to add smaple mosaic in cooper-nathans?
            mat_cov[1, 1] += q_mod**2 * instru.sample._mosaic_h**2
            mat_cov[2, 2] += q_mod**2 * instru.sample._mosaic_v**2

            mat_reso = np.linalg.inv(mat_cov) * sig2fwhm**2

            motor_angles = instru.calculate_motor_angles(hkl=hkl, en=ei - ef)
            if motor_angles is None:
                print(f"Cannot reach hkl={np.round(hkl,3)}, ei={ei}, ef={ef}.")
                rez = ResoEllipsoid(hkle=hkl + (ei - ef,), projection=projection, sample=instru.sample)
                rez.STATUS = False
                rez._set_labels()
                rez_list.append(rez)
                continue

            r_mat = instru.goniometer.r_mat(motor_angles)
            ub_mat = instru.sample.ub_conf._ub_mat

            conv_mat = 2 * np.pi * np.matmul(r_mat, ub_mat)
            rez = ResoEllipsoid(hkle=hkl + (ei - ef,), projection=projection, sample=instru.sample)
            rez._project_to_frame(mat_reso, phi, conv_mat)

            # TODO check normalization factor
            # -------------------------------------------------------------------------
            # - if the instruments works in kf=const mode and the scans are counted for
            #   or normalised to monitor counts no ki^3 or kf^3 factor is needed.
            # - if the instrument works in ki=const mode the kf^3 factor is needed.

            if R0:  # calculate
                # r0 = np.pi**2 / 4 / np.sin(theta_m) / np.sin(theta_a)
                # r0 *= np.sqrt(np.linalg.det(mat_f) / np.linalg.det(mat_h))
                r0 = 1
            else:
                r0 = 0

            rez.r0 = r0

            if np.isnan(rez.r0) or np.isinf(rez.r0) or np.isnan(rez.mat.any()) or np.isinf(rez.mat.any()):
                rez.STATUS = False
            else:
                rez.STATUS = True

            rez._set_labels()
            rez_list.append(rez)

        return rez_list[0] if len(rez_list) == 1 else rez_list
