import numpy as np
import numpy.linalg as la
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

    def coh_fwhms(self):
        """Coherent FWHM"""
        reso = self.mat
        vecFwhms = []

        for i in range(len(reso)):
            vecFwhms.append(sig2fwhm / np.sqrt(reso[i, i]))

        return np.array(vecFwhms)

    def incoh_fwhms(self):
        """Incoherent FWHM"""

        reso = self.mat
        Qres_proj_Qpara = Reso_Ellipsoid.quadric_proj(reso, 3)
        Qres_proj_Qpara = Reso_Ellipsoid.quadric_proj(Qres_proj_Qpara, 2)
        Qres_proj_Qpara = Reso_Ellipsoid.quadric_proj(Qres_proj_Qpara, 1)

        Qres_proj_Qperp = Reso_Ellipsoid.quadric_proj(reso, 3)
        Qres_proj_Qperp = Reso_Ellipsoid.quadric_proj(Qres_proj_Qperp, 2)
        Qres_proj_Qperp = Reso_Ellipsoid.quadric_proj(Qres_proj_Qperp, 0)

        Qres_proj_Qup = Reso_Ellipsoid.quadric_proj(reso, 3)
        Qres_proj_Qup = Reso_Ellipsoid.quadric_proj(Qres_proj_Qup, 1)
        Qres_proj_Qup = Reso_Ellipsoid.quadric_proj(Qres_proj_Qup, 0)

        Qres_proj_E = Reso_Ellipsoid.quadric_proj(reso, 2)
        Qres_proj_E = Reso_Ellipsoid.quadric_proj(Qres_proj_E, 1)
        Qres_proj_E = Reso_Ellipsoid.quadric_proj(Qres_proj_E, 0)

        return np.array(
            [
                1.0 / np.sqrt(np.abs(Qres_proj_Qpara[0, 0])) * sig2fwhm,
                1.0 / np.sqrt(np.abs(Qres_proj_Qperp[0, 0])) * sig2fwhm,
                1.0 / np.sqrt(np.abs(Qres_proj_Qup[0, 0])) * sig2fwhm,
                1.0 / np.sqrt(np.abs(Qres_proj_E[0, 0])) * sig2fwhm,
            ]
        )

    def descr_ellipse(self):
        """calculates the characteristics of a given ellipse by principal axis trafo"""
        [evals, evecs] = la.eig(self.mat)

        fwhms = 1.0 / np.sqrt(np.abs(evals)) * sig2fwhm

        angles = np.array([])
        if len(quadric) == 2:
            angles = np.array([np.arctan2(evecs[1][0], evecs[0][0])])

        return [fwhms, angles / np.pi * 180.0, evecs]
