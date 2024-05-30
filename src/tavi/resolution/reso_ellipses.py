import numpy as np
import numpy.linalg as la
import mpl_toolkits.mplot3d as mplot3d
import matplotlib
import matplotlib.pyplot as plot
from tavi.utilities import *


np.set_printoptions(floatmode="fixed", precision=4)


class Reso_Ellipsoid(object):
    """Manage the resolution ellipoid

    Attributs:
        STATUS (bool): True if resolution calculation is successful
        mat (float): 4 by 4 resolution matrix
        r0 (float | None): Normalization factor

    Methods:

    """

    g_eps = 1e-8

    def __init__(self, mat, r0):

        self.mat = mat
        self.r0 = r0

        if np.isnan(r0) or np.isinf(r0) or np.isnan(mat.any()) or np.isinf(mat.any()):
            self.STATUS = False
        else:
            self.STATUS = True

    def ellipsoid_volume(self):
        """volume of the ellipsoid"""
        pass

    @staticmethod
    def quadric_proj(quadric, idx):
        """projects along one axis of the quadric"""

        if np.abs(quadric[idx, idx]) < Reso_Ellipsoid.g_eps:
            return np.delete(np.delete(quadric, idx, axis=0), idx, axis=1)

        # row/column along which to perform the orthogonal projection
        vec = 0.5 * (quadric[idx, :] + quadric[:, idx])  # symmetrise if not symmetric
        vec /= np.sqrt(quadric[idx, idx])  # normalise to indexed component
        proj_op = np.outer(vec, vec)  # projection operator
        ortho_proj = quadric - proj_op  # projected quadric

        # comparing with simple projection
        # rank = len(quadric)
        # vec /= np.sqrt(np.dot(vec, vec))
        # proj_op = np.outer(vec, vec)
        # ortho_proj = np.dot((np.identity(rank) - proj_op), quadric)

        # remove row/column that was projected out
        # print("\nProjected row/column %d:\n%s\n->\n%s.\n" % (idx, str(quadric), str(ortho_proj)))
        return np.delete(np.delete(ortho_proj, idx, axis=0), idx, axis=1)

    @property
    def coh_fwhms(self):
        """Coherent FWHMs"""
        reso = self.mat
        fwhms = []

        for i in range(len(reso)):
            fwhms.append(sig2fwhm / np.sqrt(reso[i, i]))
        fwhms = np.array(fwhms)

        print(f"Coherent-elastic fwhms:{fwhms}")

    @property
    def incoh_fwhms(self):
        """Incoherent FWHMs"""

        reso = self.mat
        proj_q_para = Reso_Ellipsoid.quadric_proj(reso, 3)
        proj_q_para = Reso_Ellipsoid.quadric_proj(proj_q_para, 2)
        proj_q_para = Reso_Ellipsoid.quadric_proj(proj_q_para, 1)

        proj_q_perp = Reso_Ellipsoid.quadric_proj(reso, 3)
        proj_q_perp = Reso_Ellipsoid.quadric_proj(proj_q_perp, 2)
        proj_q_perp = Reso_Ellipsoid.quadric_proj(proj_q_perp, 0)

        proj_q_up = Reso_Ellipsoid.quadric_proj(reso, 3)
        proj_q_up = Reso_Ellipsoid.quadric_proj(proj_q_up, 1)
        proj_q_up = Reso_Ellipsoid.quadric_proj(proj_q_up, 0)

        proj_en = Reso_Ellipsoid.quadric_proj(reso, 2)
        proj_en = Reso_Ellipsoid.quadric_proj(proj_en, 1)
        proj_en = Reso_Ellipsoid.quadric_proj(proj_en, 0)

        fwhms = np.array(
            [
                1.0 / np.sqrt(np.abs(proj_q_para[0, 0])) * sig2fwhm,
                1.0 / np.sqrt(np.abs(proj_q_perp[0, 0])) * sig2fwhm,
                1.0 / np.sqrt(np.abs(proj_q_up[0, 0])) * sig2fwhm,
                1.0 / np.sqrt(np.abs(proj_en[0, 0])) * sig2fwhm,
            ]
        )

        print(f"Incoherent-elastic fwhms: {fwhms}")

    @property
    def principal_fwhms(self):
        """FWHMs are principal axes"""

        evals, _ = la.eig(self.mat)
        fwhms = np.array(1.0 / np.sqrt(np.abs(evals)) * sig2fwhm)
        print(f"Principal axes fwhms: {fwhms}")

    # def descr_ellipse(self):
    #     """calculates the characteristics of a given ellipse by principal axis trafo"""
    #     [evals, evecs] = la.eig(self.mat)

    #     fwhms = 1.0 / np.sqrt(np.abs(evals)) * sig2fwhm

    #     angles = np.array([])
    #     if len(quadric) == 2:
    #         angles = np.array([np.arctan2(evecs[1][0], evecs[0][0])])

    #     return [fwhms, angles / np.pi * 180.0, evecs]

    def plot_ellipses(
        self,
        ellis,
        verbose=True,
        plot_results=True,
        file="",
        dpi=600,
        ellipse_points=128,
        use_tex=False,
    ):
        """Plot 2D ellipses"""
        matplotlib.rc("text", usetex=use_tex)
