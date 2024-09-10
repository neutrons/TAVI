import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
from mpl_toolkits.axisartist import Axes, Subplot
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

from tavi.utilities import get_angle_vec, sig2fwhm

np.set_printoptions(floatmode="fixed", precision=4)


class ResoCurve(object):
    """1D Gaussian curve

    Attributes:
        cen (float)
        fwhm (float)
        xlabel
        ylabel
        title
        legend

    Methods:
        generate_curve
        generate_plot
    """

    def __init__(self):
        self.center: float = None
        self.fwhm = None
        self.r0 = None
        self.x = None
        self.y = None
        self.xlabel = None
        self.ylabel = None
        self.title = None
        self.legend = None

    def generate_curve(self, num_of_sigmas=3, points=100):
        """Generate points on the Gaussian curve"""

        sigma = self.fwhm / sig2fwhm
        cen = self.cen
        self.x = np.linspace(
            -num_of_sigmas * sigma + cen,
            num_of_sigmas * sigma + cen,
            points,
        )
        amp = self.r0 / np.sqrt(2 * np.pi) / sigma
        self.y = amp * np.exp(-((self.x - cen) ** 2) / sigma**2 / 2)

    def generate_plot(self, ax, c="black", linestyle="solid"):
        """Generate the plot"""
        self.generate_curve()

        s = ax.plot(
            self.x,
            self.y,
            c=c,
            linestyle=linestyle,
            label=self.legend,
        )
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.set_title(self.title)
        ax.grid(alpha=0.6)
        ax.legend()


class ResoEllipse(object):
    """2D ellipses

    Attributes:
        mat(2 by 2 matrix)
        centers (tuple): 2 by 1
        fwhms (tuple): 2 by 1
        vecs: eigen-vectors
        angle: angle between the two plotting axes, in degrees
        axes_labels (str)
        ORIGIN (bool|None): shift the origin if True

    Methods:
        generate_ellipse(ellipse_points=128)
        generate_axes()
        generate_plot(ax, c="black", linestyle="solid")
        plot()

    """

    g_esp = 1e-8

    def __init__(self):
        self.mat = None
        self.centers = None
        self.fwhms = None
        self.vecs = None
        self.angle = None
        self.axes_labels = None

        self.ORIGIN = None
        self.grid_helper = None

    def generate_ellipse(self, ellipse_points=128):
        """Generate points on a ellipse"""

        phi = np.linspace(0, 2.0 * np.pi, ellipse_points)

        pts = np.dot(
            self.vecs,
            np.array(
                [
                    self.fwhms[0] / 2 * np.cos(phi),
                    self.fwhms[1] / 2 * np.sin(phi),
                ],
            ),
        )
        if self.ORIGIN:
            pts[0] += self.centers[0]
            pts[1] += self.centers[1]
        return pts

    def _tr(self, x, y):
        x, y = np.asarray(x), np.asarray(y)
        return x + y / np.tan(self.angle / 180 * np.pi), y

    def _inv_tr(self, x, y):
        x, y = np.asarray(x), np.asarray(y)
        return x - y / np.tan(self.angle / 180 * np.pi), y

    def generate_axes(self):
        """Generate grid helper"""

        if not np.abs(self.angle - 90) < 1e-2:  # regular axes
            self.grid_helper = GridHelperCurveLinear(
                (self._tr, self._inv_tr),
                grid_locator1=MaxNLocator(integer=True, steps=[1]),
                grid_locator2=MaxNLocator(integer=True, steps=[1]),
            )

    def generate_plot(self, ax, c="black", linestyle="solid"):
        """Gnerate the ellipse for plotting"""

        pts = self.generate_ellipse()

        if self.grid_helper is None:

            s = ax.plot(
                pts[0],
                pts[1],
                c=c,
                linestyle=linestyle,
            )
        else:  # askew axes
            s = ax.plot(
                *self._tr(pts[0], pts[1]),
                c=c,
                linestyle=linestyle,
            )

        ax.set_xlabel(self.axes_labels[0])
        ax.set_ylabel(self.axes_labels[1])
        ax.grid(alpha=0.6)

        return None

    def plot(self):
        """Plot the ellipsis."""

        fig = plt.figure()
        ax = Subplot(fig, 1, 1, 1, grid_helper=self.grid_helper)
        fig.add_subplot(ax)
        self.generate_plot(ax)
        fig.show()


class ResoEllipsoid(object):
    """Manage the 4D resolution ellipoid

    Attributs:
        STATUS (None | bool): True if resolution calculation is successful
        frame (str): "q", "hkl", or "proj"
        projection (tuple): three non-coplanar vectors
        angles (tuple): angles between the projection vectors
        q (tuple): momentum transfer (h', k', l') in the coordinate system specified by projection
        hkl (tuple): momentum transfer (h, k, l)
        en (float): energy transfer

        mat (float): 4 by 4 resolution matrix
        r0 (float | None): Normalization factor

    Methods:
        generate_ellipse(axes=(0, 1), PROJECTION=False, ORIGIN=True)
        set_labels()
        plot(): plot all six combination of axes

    """

    g_eps = 1e-8

    def __init__(self):

        self.STATUS = None
        self.q = None
        self.hkl = None
        self.en = None
        self.frame = None
        self.projection = None
        self.angles = None
        self.axes_labels = None

        self.mat = None
        self.r0 = None

    def determinate_frame(self, sample):
        """determinate the frame from the projection vectors"""
        match self.projection:
            case None:  # Local Q frame
                self.frame = "q"
                self.angles = (90, 90, 90)
                self.STATUS = True

            case ((1, 0, 0), (0, 1, 0), (0, 0, 1)):  # HKL
                self.frame = "hkl"
                self.q = self.hkl
                self.angles = (sample.gamma_star, sample.alpha_star, sample.beta_star)
                self.STATUS = True

            case _:  # customized projection
                p1, p2, p3 = self.projection
                reciprocal_vecs = [sample.a_star_vec, sample.b_star_vec, sample.c_star_vec]
                v1 = np.sum([p1[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)
                v2 = np.sum([p2[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)
                v3 = np.sum([p3[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)

                if np.dot(v1, np.cross(v2, v3)) < ResoEllipsoid.g_esp:
                    raise ValueError("Projection is left handed! Please use right-handed projection")
                if np.abs(np.dot(v1, np.cross(v2, v3))) < ResoEllipsoid.g_esp:
                    raise ValueError("Projection vectors need to be non-coplanar.")

                # mat_w = np.array([p1, p2, p3]).T
                mat_w_inv = np.array([np.cross(p2, p3), np.cross(p3, p1), np.cross(p1, p2)]) / np.dot(
                    p1, np.cross(p2, p3)
                )
                hkl_prime = mat_w_inv @ self.hkl
                self.frame = "proj"
                self.q = hkl_prime
                self.angles = (get_angle_vec(v1, v2), get_angle_vec(v2, v3), get_angle_vec(v3, v1))
                self.STATUS = True

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

        if np.abs(quadric[idx, idx]) < ResoEllipsoid.g_eps:
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
                ellipse.angle = 90

        if PROJECTION:
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
        # determine frame
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
