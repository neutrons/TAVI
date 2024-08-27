import numpy as np
import numpy.linalg as la

from tavi.instrument.resolution.reso_ellipses import ResoEllipsoid
from tavi.instrument.tas import TAS
from tavi.utilities import *


class CN(TAS):
    """Copper-Nathans method

    Attibutes:
        _mat_f
        _mat_g

    Methods:
        cooper_nathans

    """

    g_esp = 1e-8

    # 4 soller slits collimators
    NUM_COLLS = 4
    IDX_COLL0_H = 0
    IDX_COLL0_V = 2
    IDX_COLL1_H = 1
    IDX_COLL1_V = 3
    IDX_COLL2_H = 4
    IDX_COLL2_V = 6
    IDX_COLL3_H = 5
    IDX_COLL3_V = 7
    # 1 monochromator and 1 analyzer
    NUM_MONOS = 1
    NUM_ANAS = 1
    IDX_MONO0_H = 0
    IDX_MONO0_V = 1
    IDX_ANA0_H = 2
    IDX_ANA0_V = 3

    def __init__(self) -> None:
        """Load instrument configuration from json if provided"""
        super().__init__()

        # constants independent of q and eng
        self._mat_f = None
        self._mat_g = None

    def cooper_nathans(
        self,
        ei: float,
        ef: float,
        hkl: tuple[float],
        projection: tuple[tuple] = ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        r0: bool = False,
    ):
        """Calculate resolution using Cooper-Nathans method

        Args:
            ei (float): incident energy, in units of meV
            ef (float): final energy, in units of meV
            hkl (tuple of floats): momentum transfer, miller indices in reciprocal lattice
            R0 (bool): calculate normalization factor if True
            projection (tuple): three non-coplaner vectors. If projection is None, the calculation is done in local Q frame

        """

        rez = ResoEllipsoid()
        rez.projection = projection

        if self._mat_f is None:
            # matrix F, divergence of monochromator and analyzer, [pop75] Appendix 1
            mat_f = np.zeros(((CN.NUM_MONOS + CN.NUM_ANAS) * 2, (CN.NUM_MONOS + CN.NUM_ANAS) * 2))
            mat_f[CN.IDX_MONO0_H, CN.IDX_MONO0_H] = 1.0 / self.monochromator._mosaic_h**2
            mat_f[CN.IDX_MONO0_V, CN.IDX_MONO0_V] = 1.0 / self.monochromator._mosaic_v**2
            mat_f[CN.IDX_ANA0_H, CN.IDX_ANA0_H] = 1.0 / self.analyzer._mosaic_h**2
            mat_f[CN.IDX_ANA0_V, CN.IDX_ANA0_V] = 1.0 / self.analyzer._mosaic_v**2
            self._mat_f = mat_f

        if self._mat_g is None:
            (
                coll_h_pre_mono,
                coll_h_pre_sample,
                coll_h_post_sample,
                coll_h_post_ana,
            ) = self.collimators._horizontal_divergence

            (
                coll_v_pre_mono,
                coll_v_pre_sample,
                coll_v_post_sample,
                coll_v_post_ana,
            ) = self.collimators._vertical_divergence

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

            self._mat_g = mat_g

        # determine frame
        if isinstance(hkl, tuple | list) and len(hkl) == 3:
            if projection is None:  # Local Q frame
                rez.frame = "q"
                rez.angles = (90, 90, 90)
                rez.STATUS = True

            elif projection == ((1, 0, 0), (0, 1, 0), (0, 0, 1)):  # HKL
                rez.frame = "hkl"
                rez.q = hkl
                rez.hkl = hkl
                rez.angles = (
                    self.sample.gamma_star,
                    self.sample.alpha_star,
                    self.sample.beta_star,
                )
                rez.STATUS = True

            else:  # customized projection
                p1, p2, p3 = projection
                reciprocal_vecs = [
                    self.sample.a_star_vec,
                    self.sample.b_star_vec,
                    self.sample.c_star_vec,
                ]
                v1 = np.sum([p1[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)
                v2 = np.sum([p2[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)
                v3 = np.sum([p3[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)

                if np.dot(v1, np.cross(v2, v3)) < CN.g_esp:
                    # TODO
                    print("Left handed!")
                if np.abs(np.dot(v1, np.cross(v2, v3))) < CN.g_esp:
                    print("Projection vectors need to be non-coplanar.")
                else:
                    mat_w = np.array([p1, p2, p3]).T
                    mat_w_inv = np.array(
                        [
                            np.cross(p2, p3),
                            np.cross(p3, p1),
                            np.cross(p1, p2),
                        ]
                    ) / np.dot(p1, np.cross(p2, p3))

                    hkl_prime = mat_w_inv @ hkl
                    rez.frame = "proj"
                    rez.q = hkl_prime
                    rez.hkl = hkl

                    rez.angles = (
                        get_angle_vec(v1, v2),
                        get_angle_vec(v2, v3),
                        get_angle_vec(v3, v1),
                    )

            angles = self.find_angles(hkl, ei, ef)  # s2, s1, sgl, sgu
            if angles is not None:
                r_mat = self.goniometer.r_mat(angles[1:])  # s1, sgl, sgu
                ub_mat = self.sample.ub_matrix
                conv_mat = 2 * np.pi * r_mat @ ub_mat
                q_lab = conv_mat @ hkl
                q_mod = np.linalg.norm(q_lab)
                rez.STATUS = True
            else:
                rez.STATUS = False

        else:
            print("q needs to be a tupe of size 3.")
            rez.STATUS = False

        if rez.STATUS:
            ki = np.sqrt(ei / ksq2eng)
            kf = np.sqrt(ef / ksq2eng)
            en = ei - ef
            rez.en = en

            two_theta = get_angle(ki, kf, q_mod) * self.goniometer.sense

            # phi = <ki, q>, always has the oppositie sign of s2
            phi = get_angle(ki, q_mod, kf) * self.goniometer.sense * (-1)

            theta_m = get_angle_bragg(ki, self.monochromator.d_spacing) * self.monochromator.sense
            theta_a = get_angle_bragg(kf, self.analyzer.d_spacing) * self.analyzer.sense

            # TODO
            # curved monochromator and analyzer

            # TODO
            # reflection efficiency

            # TODO
            # matrix A,Y=AU, tranform from collimators angular divergence to  ki-kf frame
            mat_a = np.zeros((6, 2 * CN.NUM_COLLS))
            mat_a[0, CN.IDX_COLL0_H] = 0.5 * ki / np.tan(theta_m)
            mat_a[0, CN.IDX_COLL1_H] = -0.5 * ki / np.tan(theta_m)
            mat_a[1, CN.IDX_COLL1_H] = ki
            mat_a[2, CN.IDX_COLL1_V] = -ki

            mat_a[3, CN.IDX_COLL2_H] = 0.5 * kf / np.tan(theta_a)
            mat_a[3, CN.IDX_COLL3_H] = -0.5 * kf / np.tan(theta_a)
            mat_a[4, CN.IDX_COLL2_H] = kf
            mat_a[5, CN.IDX_COLL2_V] = kf

            # matrix B, X=BY, transform from ki-kf frame to momentum transfer q-frame
            mat_b = np.zeros((4, 6))
            mat_b[0:3, 0:3] = rotation_matrix_2d(phi)
            mat_b[0:3, 3:6] = rotation_matrix_2d(phi - two_theta) * (-1)
            mat_b[3, 0] = 2 * ksq2eng * ki
            mat_b[3, 3] = -2 * ksq2eng * kf

            # matrix C, constrinat between mono/ana mosaic and collimator divergence
            mat_c = np.zeros(((CN.NUM_MONOS + CN.NUM_ANAS) * 2, CN.NUM_COLLS * 2))
            mat_c[CN.IDX_MONO0_H, CN.IDX_COLL0_H] = 0.5
            mat_c[CN.IDX_MONO0_H, CN.IDX_COLL1_H] = 0.5
            mat_c[CN.IDX_MONO0_V, CN.IDX_COLL0_V] = 0.5 / np.sin(theta_m)
            mat_c[CN.IDX_MONO0_V, CN.IDX_COLL1_V] = -0.5 / np.sin(theta_m)

            mat_c[CN.IDX_ANA0_H, CN.IDX_COLL2_H] = 0.5
            mat_c[CN.IDX_ANA0_H, CN.IDX_COLL3_H] = 0.5
            mat_c[CN.IDX_ANA0_V, CN.IDX_COLL2_V] = 0.5 / np.sin(theta_a)
            mat_c[CN.IDX_ANA0_V, CN.IDX_COLL3_V] = -0.5 / np.sin(theta_a)

            mat_h = mat_c.T @ self._mat_f @ mat_c + self._mat_g
            mat_h_inv = la.inv(mat_h)
            mat_ba = mat_b @ mat_a
            mat_cov = mat_ba @ mat_h_inv @ mat_ba.T

            # TODO hwo to add smaple mosaic i cooper-nathans?
            mat_cov[1, 1] += q_mod**2 * self.sample._mosaic_h**2
            mat_cov[2, 2] += q_mod**2 * self.sample._mosaic_v**2

            mat_reso = la.inv(mat_cov) * sig2fwhm**2

            if rez.frame == "q":
                rez.mat = mat_reso
                rez.q = (q_mod, 0, 0)

            elif rez.frame == "hkl":
                conv_mat_4d = np.eye(4)
                conv_mat_4d[0:3, 0:3] = (
                    np.array(
                        [
                            [np.sin(phi), 0, np.cos(phi)],
                            [np.cos(phi), 0, -np.sin(phi)],
                            [0, 1, 0],
                        ]
                    )
                    @ conv_mat
                )
                rez.mat = conv_mat_4d.T @ mat_reso @ conv_mat_4d

            elif rez.frame == "proj":
                conv_mat_4d = np.eye(4)
                conv_mat_4d[0:3, 0:3] = (
                    np.array(
                        [
                            [np.sin(phi), 0, np.cos(phi)],
                            [np.cos(phi), 0, -np.sin(phi)],
                            [0, 1, 0],
                        ]
                    )
                    @ conv_mat
                    @ mat_w
                )
                rez.mat = conv_mat_4d.T @ mat_reso @ conv_mat_4d

            # TODO check normalization factor
            # -------------------------------------------------------------------------
            # - if the instruments works in kf=const mode and the scans are counted for
            #   or normalised to monitor counts no ki^3 or kf^3 factor is needed.
            # - if the instrument works in ki=const mode the kf^3 factor is needed.

            if r0:  # calculate
                r0 = np.pi**2 / 4 / np.sin(theta_m) / np.sin(theta_a)
                r0 *= np.sqrt(np.linalg.det(self._mat_f) / np.linalg.det(mat_h))
            else:
                r0 = 0

            rez.r0 = r0

            if np.isnan(rez.r0) or np.isinf(rez.r0) or np.isnan(rez.mat.any()) or np.isinf(rez.mat.any()):
                rez.STATUS = False
            else:
                rez.STATUS = True

        rez.set_labels()
        return rez
