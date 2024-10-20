from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
from mpl_toolkits.axisartist import Axes

from tavi.instrument.resolution.ellipse_curve import ResoCurve, ResoEllipse
from tavi.sample.xtal import Xtal
from tavi.utilities import get_angle_vec, sig2fwhm

np.set_printoptions(floatmode="fixed", precision=4)


class ResoEllipsoid(object):
    """Manage the 4D resolution ellipoid

    Attributs:
        STATUS (None | bool): True if resolution calculation is successful
        frame (str): "q", "hkl", or "proj"
        projection (tuple): three non-coplanar vectors
        q (tuple): momentum transfer (h', k', l') in the coordinate system specified by projection
        hkl (tuple): momentum transfer (h, k, l)
        en (float): energy transfer
        agnles (tuple): angles between plotting axes

        mat (np.ndarray): 4 by 4 resolution matrix
        r0 (float | None): Normalization factor

    Methods:
        generate_ellipse(axes=(0, 1), PROJECTION=False, ORIGIN=True)
        set_labels()
        plot(): plot all six combination of axes

    """

    EPS = 1e-8

    def __init__(
        self,
        hkl: tuple[float, float, float],
        en: float,
        sample: Xtal,
        projection: Optional[tuple] = ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    ) -> None:

        self.STATUS: Optional[bool] = None
        self.q: Optional[tuple[float, float, float]] = None
        self.hkl: tuple[float, float, float] = hkl
        self.en: float = en

        self.projection: Optional[tuple] = projection
        self.angles: tuple[float, float, float] = (90.0, 90.0, 90.0)
        self.axes_labels = None

        self.mat = None
        self.r0 = None

        match self.projection:
            case None:  # Local Q frame
                self.frame = "q"
                self.angles = (90.0, 90.0, 90.0)
                self.q = (np.linalg.norm(sample.b_mat @ hkl) * 2 * np.pi, 0.0, 0.0)

            case ((1, 0, 0), (0, 1, 0), (0, 0, 1)):  # HKL
                self.frame = "hkl"
                self.q = hkl
                self.angles = (sample.gamma_star, sample.alpha_star, sample.beta_star)
            case _:  # customized projection
                p1, p2, p3 = self.projection
                reciprocal_vecs = [sample.a_star_vec, sample.b_star_vec, sample.c_star_vec]
                v1 = np.sum([p1[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)
                v2 = np.sum([p2[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)
                v3 = np.sum([p3[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)

                if np.dot(v1, np.cross(v2, v3)) < ResoEllipsoid.EPS:
                    raise ValueError("Projection is left handed! Please use right-handed projection")
                if np.abs(np.dot(v1, np.cross(v2, v3))) < ResoEllipsoid.EPS:
                    raise ValueError("Projection vectors need to be non-coplanar.")

                mat_w_inv = np.array([np.cross(p2, p3), np.cross(p3, p1), np.cross(p1, p2)]) / np.dot(
                    p1, np.cross(p2, p3)
                )
                hkl_prime = mat_w_inv @ self.hkl
                self.frame = "proj"
                self.q = hkl_prime
                self.angles = (get_angle_vec(v1, v2), get_angle_vec(v2, v3), get_angle_vec(v3, v1))

        self.set_labels()

    def project_to_frame(self, mat_reso, phi, conv_mat):
        """determinate the frame from the projection vectors"""

        mat_lab_to_local = np.array(
            [
                [np.sin(phi), 0, np.cos(phi)],
                [np.cos(phi), 0, -np.sin(phi)],
                [0, 1, 0],
            ]
        )
        match self.frame:
            case "q":
                self.mat = mat_reso

            case "hkl":
                conv_mat_4d = np.eye(4)
                conv_mat_4d[0:3, 0:3] = mat_lab_to_local @ conv_mat
                self.mat = conv_mat_4d.T @ mat_reso @ conv_mat_4d

            case "proj":
                p1, p2, p3 = self.projection
                mat_w = np.array([p1, p2, p3]).T

                conv_mat_4d = np.eye(4)
                conv_mat_4d[0:3, 0:3] = mat_lab_to_local @ conv_mat @ mat_w
                self.mat = conv_mat_4d.T @ mat_reso @ conv_mat_4d

    def volume(self):
        """volume of the ellipsoid"""
        pass

    def coh_fwhms(self, axis=None):
        """Coherent FWHM"""

        curve = ResoCurve()
        idx = int(axis)
        curve.fwhm = np.array(sig2fwhm / np.sqrt(self.mat[idx, idx]))
        if idx == 3:
            curve.cen = self.en
        else:
            curve.cen = self.q[idx]

        curve.r0 = self.r0
        curve.xlabel = self.axes_labels[idx]
        curve.title = f"q={np.round(self.hkl,3)}, en={np.round(self.en,3)}"
        curve.legend = f"coherent FWHM={np.round(curve.fwhm, 3)}"
        return curve

    def incoh_fwhms(self, axis=None):
        """Incoherent FWHMs"""

        curve = ResoCurve()
        idx = int(axis)

        reso = self.mat
        for i in (3, 2, 1, 0):
            if not i == idx:
                reso = ResoEllipsoid.quadric_proj(reso, i)

        curve.fwhm = (1.0 / np.sqrt(np.abs(reso[0, 0])) * sig2fwhm,)
        if idx == 3:
            curve.cen = self.en
        else:
            curve.cen = self.q[idx]

        curve.r0 = self.r0
        curve.xlabel = self.axes_labels[idx]
        curve.title = f"q={np.round(self.hkl,3)}, en={np.round(self.en,3)}"
        curve.legend = f"incoherent FWHM={np.round(curve.fwhm , 3)}"
        return curve

    # def principal_fwhms_calc(self):
    #     """FWHMs of principal axes in 4D"""

    #     evals, _ = la.eig(self.mat)
    #     fwhms = np.array(1.0 / np.sqrt(np.abs(evals)) * sig2fwhm)
    #     self.principal_fwhms = fwhms

    @staticmethod
    def quadric_proj(quadric, idx):
        """projects along one axis of the quadric"""

        if np.abs(quadric[idx, idx]) < ResoEllipsoid.EPS:
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

    def generate_ellipse(self, axes=(0, 1), PROJECTION=False, ORIGIN=True):
        """Gnerate a 2D ellipse by either making a cut or projection

        Arguments:
            axes(tuple)
            PROJECTION: slice if False
            ORIGIN: shift the center if True

        """
        ellipse = ResoEllipse()
        ellipse.ORIGIN = ORIGIN
        qe_list = np.concatenate((self.q, self.en), axis=None)
        # axes = np.sort(axes)
        # match tuple(np.sort(axes)):
        match tuple(axes):
            case (0, 1) | (1, 0):
                ellipse.angle = np.round(self.angles[0], 2)
            case (1, 2) | (2, 1):
                ellipse.angle = np.round(self.angles[1], 2)
            case (0, 2) | (2, 0):
                ellipse.angle = np.round(self.angles[2], 2)
            case _:
                ellipse.angle = 90.0

        if PROJECTION:  # make a projection
            q_res = self.mat
            ellipse.centers = (qe_list[axes[0]], qe_list[axes[1]])
            ellipse.axes_labels = (self.axes_labels[axes[0]], self.axes_labels[axes[1]])
            # Qres_QxE_proj = np.delete(np.delete(self.mat, 2, axis=0), 2, axis=1)
            for i in (3, 2, 1, 0):
                if i not in axes:
                    q_res = ResoEllipsoid.quadric_proj(q_res, i)
            ellipse.mat = q_res

        else:  # make a slice
            q_res = self.mat
            ellipse.centers = (qe_list[axes[0]], qe_list[axes[1]])
            ellipse.axes_labels = (self.axes_labels[axes[0]], self.axes_labels[axes[1]])
            for i in (3, 2, 1, 0):
                if i not in axes:
                    q_res = np.delete(np.delete(q_res, i, axis=0), i, axis=1)
            ellipse.mat = q_res

        evals, evecs = la.eig(ellipse.mat)  # evecs normalized already
        fwhms = 1.0 / np.sqrt(np.abs(evals)) * sig2fwhm
        # angles = np.array([np.arctan2(evecs[1][0], evecs[0][0])])
        ellipse.fwhms = fwhms
        ellipse.vecs = evecs
        ellipse.generate_axes()

        return ellipse

    def set_labels(self):
        """Set axes labels based on the frame"""
        match self.frame:
            case "q":
                self.axes_labels = ("Q_para (1/A)", "Q_perp (1/A)", "Q_up (1/A)", "E (meV)")
            case "hkl":
                self.axes_labels = ("H (r.l.u.)", "K (r.l.u.)", "L (r.l.u.)", "E (meV)")
            case "proj":
                self.axes_labels = (
                    f"{self.projection[0]}",
                    f"{self.projection[1]}",
                    f"{self.projection[2]}",
                    "E (meV)",
                )

    # TODO
    def calc_ellipses(self, verbose=True):
        """Calculate FWHMs and rotation angles"""

        # 2d sliced ellipses

        elps_qx_en = self.generate_ellipse(axes=(0, 3), PROJECTION=False)
        # elps_qx_en.plot()

        elps_qy_en = self.generate_ellipse(axes=(1, 3), PROJECTION=False)
        elps_qz_en = self.generate_ellipse(axes=(2, 3), PROJECTION=False)

        elps_qx_qy = self.generate_ellipse(axes=(0, 1), PROJECTION=False)
        elps_qx_qy.plot()
        elps_qy_qz = self.generate_ellipse(axes=(1, 2), PROJECTION=False)
        elps_qx_qz = self.generate_ellipse(axes=(0, 2), PROJECTION=False)
        elps_qx_qz.plot()
        plt.show()

        # 2d projected ellipses
        # Qres_QxE_proj = np.delete(np.delete(self.mat, 2, axis=0), 2, axis=1)

        elps_proj_qx_en = self.generate_ellipse(axes=(0, 3), PROJECTION=True)
        elps_proj_qy_en = self.generate_ellipse(axes=(1, 3), PROJECTION=True)
        elps_proj_qz_en = self.generate_ellipse(axes=(2, 3), PROJECTION=True)

        elps_proj_qx_qy = self.generate_ellipse(axes=(0, 1), PROJECTION=True)
        elps_proj_qy_qz = self.generate_ellipse(axes=(1, 2), PROJECTION=True)
        elps_proj_qx_qz = self.generate_ellipse(axes=(0, 2), PROJECTION=True)

        return None

    def plot(self):
        """Plot all 2D ellipses"""

        # fig = plt.figure()
        fig = plt.figure(figsize=(10, 6))
        elps_qx_en = self.generate_ellipse(axes=(0, 3), PROJECTION=False)
        ax = fig.add_subplot(231, axes_class=Axes, grid_helper=elps_qx_en.grid_helper)
        elps_qx_en.generate_plot(ax, c="black", linestyle="solid")
        elps_proj_qx_en = self.generate_ellipse(axes=(0, 3), PROJECTION=True)
        elps_proj_qx_en.generate_plot(ax, c="black", linestyle="dashed")

        elps_qy_en = self.generate_ellipse(axes=(1, 3), PROJECTION=False)
        ax = fig.add_subplot(232, axes_class=Axes, grid_helper=elps_qy_en.grid_helper)
        elps_qy_en.generate_plot(ax, c="black", linestyle="solid")
        elps_proj_qy_en = self.generate_ellipse(axes=(1, 3), PROJECTION=True)
        elps_proj_qy_en.generate_plot(ax, c="black", linestyle="dashed")

        elps_qz_en = self.generate_ellipse(axes=(2, 3), PROJECTION=False)
        ax = fig.add_subplot(233, axes_class=Axes, grid_helper=elps_qz_en.grid_helper)
        elps_qz_en.generate_plot(ax, c="black", linestyle="solid")
        elps_proj_qz_en = self.generate_ellipse(axes=(2, 3), PROJECTION=True)
        elps_proj_qz_en.generate_plot(ax, c="black", linestyle="dashed")

        elps_qx_qy = self.generate_ellipse(axes=(0, 1), PROJECTION=False)
        ax = fig.add_subplot(234, axes_class=Axes, grid_helper=elps_qx_qy.grid_helper)
        elps_qx_qy.generate_plot(ax, c="black", linestyle="solid")
        elps_proj_qx_qy = self.generate_ellipse(axes=(0, 1), PROJECTION=True)
        elps_proj_qx_qy.generate_plot(ax, c="black", linestyle="dashed")

        elps_qy_qz = self.generate_ellipse(axes=(1, 2), PROJECTION=False)
        ax = fig.add_subplot(235, axes_class=Axes, grid_helper=elps_qy_qz.grid_helper)
        elps_qy_qz.generate_plot(ax, c="black", linestyle="solid")
        elps_proj_qy_qz = self.generate_ellipse(axes=(1, 2), PROJECTION=True)
        elps_proj_qy_qz.generate_plot(ax, c="black", linestyle="dashed")

        elps_qx_qz = self.generate_ellipse(axes=(0, 2), PROJECTION=False)
        ax = fig.add_subplot(236, axes_class=Axes, grid_helper=elps_qx_qz.grid_helper)
        elps_qx_qz.generate_plot(ax, c="black", linestyle="solid")
        elps_proj_qx_qz = self.generate_ellipse(axes=(0, 2), PROJECTION=True)
        elps_proj_qx_qz.generate_plot(ax, c="black", linestyle="dashed")

        fig.tight_layout(pad=2)
