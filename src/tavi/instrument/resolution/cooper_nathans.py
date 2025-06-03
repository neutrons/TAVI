from functools import partial
from typing import Union

import numpy as np

from tavi.instrument.resolution.ellipsoid import ResoEllipsoid
from tavi.instrument.resolution.resolution_calculator import ResolutionCalculator
from tavi.utilities import en2q, ksq2eng, rotation_matrix_2d, sig2fwhm


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

    @classmethod
    def mat_f(cls, monochromator, analyzer) -> np.ndarray:
        """matrix F, divergence of monochromator and analyzer, [pop75] Appendix 1
        Note: No conversion between sigma and FWHM
        """
        i, j = cls.NUM_MONOS, cls.NUM_ANAS
        mat_f = np.zeros(((i + j) * 2, (i + j) * 2))
        mat_f[cls.IDX_MONO0_H, cls.IDX_MONO0_H] = 1.0 / monochromator._mosaic_h**2
        mat_f[cls.IDX_MONO0_V, cls.IDX_MONO0_V] = 1.0 / monochromator._mosaic_v**2
        mat_f[cls.IDX_ANA0_H, cls.IDX_ANA0_H] = 1.0 / analyzer._mosaic_h**2
        mat_f[cls.IDX_ANA0_V, cls.IDX_ANA0_V] = 1.0 / analyzer._mosaic_v**2
        return mat_f

    @classmethod
    def mat_g(cls, collimators) -> np.ndarray:
        """matrix G, divergence of monochromator and analyzer, [pop75] Appendix 1
        Note: No conversion between sigma and FWHM
        """

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
        """matrix A,Y=AU, transform from collimators angular divergence to  ki-kf frame"""
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
        self.__class__.__base__.validate_instrument_parameters(self)
        # TODO
        # add more method-specif validation here

    @staticmethod
    def _calculate_at_hkle(hkl, ei, ef, instrument):
        """Calculate the resolution matrix and R0 factor in the local Q frame"""

        # q_mod = np.linalg.norm(tas.sample.b_mat @ hkl) * 2 * np.pi
        q_norm = instrument.sample.get_q_norm(hkl)
        ki, kf = en2q(ei), en2q(ef)
        motor_angles = instrument.calculate_motor_angles(hkl=hkl, en=ei - ef)

        # phi = <ki to q>, always has the oppositie sign of s2
        psi = np.radians(instrument.get_psi(hkl, en=ei - ef))
        two_theta = np.radians(motor_angles.two_theta)
        theta_m = np.radians(instrument.get_theta_m(ei))
        theta_a = np.radians(instrument.get_theta_a(ef))

        # TODO
        # curved monochromator and analyzer

        # TODO
        # reflection efficiency

        mat_a = CooperNathans.calc_mat_a(ki, kf, theta_m, theta_a)
        mat_b = CooperNathans.calc_mat_b(ki, kf, psi, two_theta)
        mat_c = CooperNathans.calc_mat_c(theta_m, theta_a)
        mat_f = CooperNathans.mat_f(instrument.monochromator, instrument.analyzer)
        mat_g = CooperNathans.mat_g(instrument.collimators)

        mat_h = mat_c.T @ mat_f @ mat_c + mat_g
        mat_h_inv = np.linalg.inv(mat_h)
        mat_ba = mat_b @ mat_a
        mat_cov = mat_ba @ mat_h_inv @ mat_ba.T

        # TODO how to add sample mosaic in cooper-nathans?
        mat_cov[1, 1] += q_norm**2 * instrument.sample._mosaic_h**2
        mat_cov[2, 2] += q_norm**2 * instrument.sample._mosaic_v**2

        mat_reso = np.linalg.inv(mat_cov) * sig2fwhm**2

        # TODO check normalization factor
        # -------------------------------------------------------------------------
        # - if the instruments works in kf=const mode and the scans are counted for
        #   or normalised to monitor counts no ki^3 or kf^3 factor is needed.
        # - if the instrument works in ki=const mode the kf^3 factor is needed.

        # monochromator and analyzer reflectivity taken to be 1
        # Note that the ki**3 or kf**3 factors are not considered
        rm_theta = 1
        ra_theta = 1
        rm_k = rm_theta * ki**3 / np.tan(theta_m)
        ra_k = ra_theta * kf**3 / np.tan(theta_a)

        r0 = rm_k * ra_k
        r0 *= np.pi**2 / 4 / np.sin(theta_m) / np.sin(theta_a)
        r0 *= np.sqrt(np.linalg.det(mat_f) / np.linalg.det(mat_h))

        return (mat_reso, r0)

    def calculate(
        self,
        hkl: Union[tuple[float, float, float], list[tuple[float, float, float]]],
        en: Union[float, list[float]],
        projection: tuple = ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    ) -> Union[None, ResoEllipsoid, tuple[ResoEllipsoid, ...]]:
        """Calculate resolution using Cooper-Nathans method

        Args:
            hkl (tuple | list(tuple)): momentum transfer, miller indices in reciprocal lattice
            en (float | list(float)): en = ei - ef is energy transfer, in units of meV
            projection (tuple): three non-coplaner vectors. If projection is None,
                                the calculation is done in local Q frame (Q_para, Q_perp, Q_up)

        Return:
            tuple (ResoEllipsoid)
        """

        self.validate_instrument_parameters()
        calculate_at_hkle = partial(CooperNathans._calculate_at_hkle, instrument=self.instrument)
        hkle_list = self.generate_hkle_list(hkl, en)

        rez_list = []
        for q_hkl, ei, ef in hkle_list:
            # check if the (Q,E) position can be reached
            motor_angles = self.instrument.calculate_motor_angles(hkl=tuple(q_hkl), en=ei - ef)
            if motor_angles is None:
                rez_list.append(None)
                continue

            reso_mat, r0 = calculate_at_hkle(q_hkl, ei, ef)
            rez = ResoEllipsoid(
                instrument=self.instrument,
                hkle=tuple(q_hkl) + (ei - ef,),
                projection=projection,
                reso_mat=reso_mat,
                r0=r0,
            )

            # if np.isnan(rez.r0) or np.isinf(rez.r0) or np.isnan(rez.mat.any()) or np.isinf(rez.mat.any()):
            #     rez.STATUS = False
            # else:
            #     rez.STATUS = True

            rez_list.append(rez)

        match len(rez_list):
            case 0:
                return None
            case 1:
                return rez_list[0]
            case _:
                return tuple(rez_list)
