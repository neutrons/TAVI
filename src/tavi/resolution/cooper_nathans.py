import numpy as np
import numpy.linalg as la
from tavi.utilities import *
from tavi.tas import TAS
from tavi.resolution.reso_ellipses import Reso_Ellipsoid
from tavi.instrument_params.takin_test import instrument_params
from tavi.sample.sample_test import test_xtal


class CN(TAS):
    """Copper-Nathans method

    Attibutes:
        _mat_f
        _mat_g


    Methods:
        cooper_nathans

    """

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

    def __init__(self):
        super().__init__()

        # constants independent of q and eng
        self._mat_f = None
        self._mat_g = None

    def cooper_nathans(self, ei, ef, q, R0=False, frame="local"):
        """
        Calculate resolution using Cooper-Nathans method

        Args:
            ei (float): incident energy, in units of meV
            ef (float): final energy, in units of meV
            q (float | tuple of floats): momentum transfer, in units of inverse angstrom
            R0 (bool): calculate normalization factor if True
            frame (str): "local", "mod_q" for powder or "hkle" for single crystal
        """
        if self._mat_f is None:
            # matrix F, divergence of monochromator and analyzer, [pop75] Appendix 1
            mat_f = np.zeros(((CN.NUM_MONOS + CN.NUM_ANAS) * 2, (CN.NUM_MONOS + CN.NUM_ANAS) * 2))
            mat_f[CN.IDX_MONO0_H, CN.IDX_MONO0_H] = 1.0 / self.monochromator["mosaic"] ** 2
            mat_f[CN.IDX_MONO0_V, CN.IDX_MONO0_V] = 1.0 / self.monochromator["mosaic_v"] ** 2
            mat_f[CN.IDX_ANA0_H, CN.IDX_ANA0_H] = 1.0 / self.analyzer["mosaic"] ** 2
            mat_f[CN.IDX_ANA0_V, CN.IDX_ANA0_V] = 1.0 / self.analyzer["mosaic_v"] ** 2
            self._mat_f = mat_f

        if self._mat_g is None:
            # matrix G, divergence of collimators, [pop75] Appendix 1
            mat_g = np.zeros((CN.NUM_COLLS * 2, CN.NUM_COLLS * 2))
            mat_g[CN.IDX_COLL0_H, CN.IDX_COLL0_H] = 1.0 / self.collimators["h_pre_mono"] ** 2
            mat_g[CN.IDX_COLL0_V, CN.IDX_COLL0_V] = 1.0 / self.collimators["v_pre_mono"] ** 2
            mat_g[CN.IDX_COLL1_H, CN.IDX_COLL1_H] = 1.0 / self.collimators["h_pre_sample"] ** 2
            mat_g[CN.IDX_COLL1_V, CN.IDX_COLL1_V] = 1.0 / self.collimators["v_pre_sample"] ** 2
            mat_g[CN.IDX_COLL2_H, CN.IDX_COLL2_H] = 1.0 / self.collimators["h_post_sample"] ** 2
            mat_g[CN.IDX_COLL2_V, CN.IDX_COLL2_V] = 1.0 / self.collimators["v_post_sample"] ** 2
            mat_g[CN.IDX_COLL3_H, CN.IDX_COLL3_H] = 1.0 / self.collimators["h_post_ana"] ** 2
            mat_g[CN.IDX_COLL3_V, CN.IDX_COLL3_V] = 1.0 / self.collimators["v_post_ana"] ** 2

            self._mat_g = mat_g

        if isinstance(q, tuple | list):
            if len(q) == 3:
                h, k, l = q
                q = self.sample.hkl2q((h, k, l))
            else:
                print("Length of q is not 3.")
        elif isinstance(q, float):
            pass
        else:
            print("q is not float or tupe or list.")

        ki = np.sqrt(ei / ksq2E)
        kf = np.sqrt(ef / ksq2E)
        en = ei - ef

        two_theta = get_angle(ki, kf, q) * self.goniometer["sense"]
        # theta_s = two_theta / 2

        # phi = <ki, q>, always has the oppositie sign of theta_s
        phi = get_angle(ki, q, kf) * self.goniometer["sense"] * (-1)

        theta_m = get_angle_bragg(ki, self.monochromator["d_spacing"]) * self.monochromator["sense"]
        theta_a = get_angle_bragg(kf, self.analyzer["d_spacing"]) * self.analyzer["sense"]

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
        mat_b[3, 0] = 2 * ksq2E * ki
        mat_b[3, 3] = -2 * ksq2E * kf

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

        mat_h = mat_c.T @ self._mat_f @ mat_c
        mat_hg_inv = la.inv(mat_h + self._mat_g)
        mat_ba = mat_b @ mat_a
        mat_cov = mat_ba @ mat_hg_inv @ mat_ba.T

        # TODO smaple mosaic????
        # mat_cov[1, 1] += q**2 * self.sample["mosaic"] ** 2
        # mat_cov[2, 2] += q**2 * self.sample["mosaic_v"] ** 2

        mat_reso = la.inv(mat_cov) * sig2fwhm**2

        # TODO normalization factor
        if R0:  # calculate
            r0 = 1
        else:
            r0 = 0

        rez = Reso_Ellipsoid(mat_reso, r0)

        return rez


if __name__ == "__main__":

    takin_instru = CN()
    takin_instru.load_instrument(instrument_params)
    takin_instru.load_sample(test_xtal)

    print(takin_instru.analyzer["type"])
    print(takin_instru.sample.a)

    # ki = kf = 1.4
    # rez = takin_instru.cooper_nathans(
    #     ei=1.4**2 * 2.072124855,
    #     ef=1.4**2 * 2.072124855,
    #     q=1.777,
    #     R0=False,
    #     frame="local",
    # )

    rez = takin_instru.cooper_nathans(
        ei=1.4**2 * 2.072124855,
        ef=1.4**2 * 2.072124855,
        q=(1, 1, 0),
        R0=False,
        frame="hkle",
    )

    print(rez.STATUS)
    print(rez.mat)
    # describe and plot ellipses
    # ellipses = reso.calc_ellipses(res["reso"], verbose)
    # reso.plot_ellipses(ellipses, verbose)
