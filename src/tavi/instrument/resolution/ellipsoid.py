from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axisartist import Axes

from tavi.instrument.resolution.ellipse import ResoEllipse
from tavi.plotter import Plot2D

# from tavi.sample.xtal import Xtal
from tavi.utilities import get_angle_vec, sig2fwhm


class ResoEllipsoid(object):
    """Manage the 4D resolution ellipoid

    Attributs:
        frame (str): "q", "hkl", or "proj"
        projection (tuple): three non-coplanar vectors
        q (tuple): momentum transfer (h', k', l') in the coordinate system specified by projection
        hkle (tuple): momentum transfer (h, k, l) and energy transfer
        agnles (tuple): angles between plotting axes

        mat (np.ndarray): 4 by 4 resolution matrix
        r0 (float | None): Normalization factor

    Methods:
        generate_ellipse(axes=(0, 1), PROJECTION=False, ORIGIN=True)
        plot(): plot all six combination of axes

    """

    def __init__(
        self,
        instrument,
        hkle: tuple[float, float, float, float],
        projection: Optional[tuple] = ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        reso_mat: Optional[np.ndarray] = None,
        r0: Optional[float] = None,
    ) -> None:
        *hkl, en = hkle
        self.hkl: tuple = tuple(hkl)
        self.en: float = en

        self.projection: Optional[tuple] = projection
        self.mat = reso_mat
        self.r0 = r0

        if projection is None:
            self.frame = "q"
            self.angles = (90.0, 90.0, 90.0)
            self.q = (np.linalg.norm(instrument.sample.b_mat @ self.hkl) * 2 * np.pi, 0.0, 0.0)
            self.axes_labels = ResoEllipsoid.labels_from_projection(projection)
        else:
            self._project_to_frame(instrument)

    def __repr__(self):
        h, k, l = self.hkl
        return f"ResoEllipsoid at hkl=({h:.3f}, {k:.3f}, {l:.3f}), en={self.en:.3f} meV"

    @staticmethod
    def labels_from_projection(projection: Optional[tuple] = ((1, 0, 0), (0, 1, 0), (0, 0, 1))):
        if projection is None:
            return ("Q_para (1/A)", "Q_perp (1/A)", "Q_up (1/A)", "E (meV)")
        elif projection == ((1, 0, 0), (0, 1, 0), (0, 0, 1)):  # HKL
            return (
                "(H, 0, 0) (r.l.u.)",
                "(0, K, 0) (r.l.u.)",
                "(0, 0, L) (r.l.u.)",
                "E (meV)",
            )
        hkl_str = ("H", "K", "L")
        labels = []
        for i, p in enumerate(projection):
            miller_str = hkl_str[i]
            label_str = "("
            for j in range(3):
                # Check if the value is an integer or a float that is close to an integer
                val = p[j]
                if isinstance(val, int) or (isinstance(val, float) and val.is_integer()):
                    num = int(val)
                    if num == 0:
                        label_str += "0"
                    elif num == 1:
                        label_str += miller_str
                    elif num == -1:
                        label_str += "-" + miller_str
                    else:
                        label_str += f"{num}" + miller_str
                else:
                    label_str += f"{p[j]:.3f}" + miller_str
                if j == 2:
                    label_str += ") (r.l.u.)"
                else:
                    label_str += ", "
            labels.append(label_str)

        labels.append("E (meV)")
        return tuple(labels)

    def _project_to_frame(self, instrument):
        """determinate the frame from the projection vectors"""

        ZERO = 1e-8

        motor_angles = instrument.calculate_motor_angles(hkl=self.hkl, en=self.en)
        r_mat = instrument.goniometer.r_mat(motor_angles)
        ub_mat = instrument.sample.ub_conf._ub_mat
        conv_mat = 2 * np.pi * np.matmul(r_mat, ub_mat)
        psi = np.radians(instrument.get_psi(self.hkl, en=self.en))

        mat_lab_to_local = np.array(
            [
                [np.sin(psi), 0, np.cos(psi)],
                [np.cos(psi), 0, -np.sin(psi)],
                [0, 1, 0],
            ]
        )
        if self.projection == ((1, 0, 0), (0, 1, 0), (0, 0, 1)):  # HKL
            self.frame = "hkl"
            self.q = self.hkl
            *_, alpha_star, beta_star, gamma_star = instrument.sample.reciprocal_latt_params
            self.angles = (gamma_star, alpha_star, beta_star)

            conv_mat_4d = np.eye(4)
            conv_mat_4d[0:3, 0:3] = mat_lab_to_local @ conv_mat
            self.mat = conv_mat_4d.T @ self.mat @ conv_mat_4d
            self.axes_labels = ResoEllipsoid.labels_from_projection()

        else:  # customized projection
            p1, p2, p3 = self.projection
            a_star_vec, b_star_vec, c_star_vec = instrument.sample.reciprocal_space_vectors
            reciprocal_vecs = [a_star_vec, b_star_vec, c_star_vec]
            v1 = np.sum([p1[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)
            v2 = np.sum([p2[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)
            v3 = np.sum([p3[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)

            if np.dot(v1, np.cross(v2, v3)) < ZERO:
                raise ValueError("Projection is left handed! Please use right-handed projection")
            if np.abs(np.dot(v1, np.cross(v2, v3))) < ZERO:
                raise ValueError("Projection vectors need to be non-coplanar.")

            mat_w_inv = np.array([np.cross(p2, p3), np.cross(p3, p1), np.cross(p1, p2)]) / np.dot(p1, np.cross(p2, p3))
            hkl_prime = mat_w_inv @ self.hkl
            self.frame = "proj"
            self.q = hkl_prime
            self.angles = (get_angle_vec(v1, v2), get_angle_vec(v2, v3), get_angle_vec(v3, v1))

            mat_w = np.array([p1, p2, p3]).T

            conv_mat_4d = np.eye(4)
            conv_mat_4d[0:3, 0:3] = mat_lab_to_local @ conv_mat @ mat_w
            self.mat = conv_mat_4d.T @ self.mat @ conv_mat_4d

            self.axes_labels = ResoEllipsoid.labels_from_projection(self.projection)

    def volume(self):
        """volume of the ellipsoid"""
        pass

    def coh_fwhms(self, axis=None):
        """Coherent FWHM"""
        idx = int(axis)
        return sig2fwhm / np.sqrt(self.mat[idx, idx])

    def incoh_fwhms(self, axis=None):
        """Incoherent FWHMs"""
        idx = int(axis)

        reso = self.mat
        for i in (3, 2, 1, 0):
            if not i == idx:
                reso = ResoEllipsoid.quadric_proj(reso, i)

        return sig2fwhm / np.sqrt(np.abs(reso[0, 0]))

    # def coh_fwhms(self, axis=None):
    #     """Coherent FWHM"""

    #     curve = ResoCurve()
    #     idx = int(axis)
    #     curve.fwhm = np.array(sig2fwhm / np.sqrt(self.mat[idx, idx]))
    #     if idx == 3:
    #         curve.cen = self.en
    #     else:
    #         curve.cen = self.q[idx]

    #     curve.r0 = self.r0
    #     curve.xlabel = self.axes_labels[idx]
    #     curve.title = f"q={np.round(self.hkl,3)}, en={np.round(self.en,3)}"
    #     curve.legend = f"coherent FWHM={np.round(curve.fwhm, 3)}"
    #     return curve

    # def incoh_fwhms(self, axis=None):
    #     """Incoherent FWHMs"""

    #     curve = ResoCurve()
    #     idx = int(axis)

    #     reso = self.mat
    #     for i in (3, 2, 1, 0):
    #         if not i == idx:
    #             reso = ResoEllipsoid.quadric_proj(reso, i)

    #     curve.fwhm = (1.0 / np.sqrt(np.abs(reso[0, 0])) * sig2fwhm,)
    #     if idx == 3:
    #         curve.cen = self.en
    #     else:
    #         curve.cen = self.q[idx]

    #     curve.r0 = self.r0
    #     curve.xlabel = self.axes_labels[idx]
    #     curve.title = f"q={np.round(self.hkl,3)}, en={np.round(self.en,3)}"
    #     curve.legend = f"incoherent FWHM={np.round(curve.fwhm , 3)}"
    #     return curve

    # def principal_fwhms_calc(self):
    #     """FWHMs of principal axes in 4D"""

    #     evals, _ = la.eig(self.mat)
    #     fwhms = np.array(1.0 / np.sqrt(np.abs(evals)) * sig2fwhm)
    #     self.principal_fwhms = fwhms

    @staticmethod
    def quadric_proj(quadric, idx):
        """projects along one axis of the quadric"""
        ZERO = 1e-8

        if np.abs(quadric[idx, idx]) < ZERO:
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

    def get_ellipse(
        self,
        axes: tuple[int, int] = (0, 1),
        PROJECTION: bool = False,
        ORIGIN: bool = True,
    ) -> ResoEllipse:
        """Gnerate a 2D ellipse by either making a cut or projection

        Arguments:
            axes(tuple)
            PROJECTION: slice if False
            ORIGIN: shift the center if True

        """

        x_axis, y_axis = axes
        qe_list = np.concatenate((self.q, self.en), axis=None)
        # axes = np.sort(axes)
        # match tuple(np.sort(axes)):
        match tuple(axes):
            case (0, 1) | (1, 0):
                angle = np.round(self.angles[0], 2)
            case (1, 2) | (2, 1):
                angle = np.round(self.angles[1], 2)
            case (0, 2) | (2, 0):
                angle = np.round(self.angles[2], 2)
            case _:
                angle = 90.0
        axes_labels = (self.axes_labels[x_axis], self.axes_labels[y_axis])

        if ORIGIN:
            centers = (qe_list[x_axis], qe_list[y_axis])
        else:
            centers = (0.0, 0.0)

        q_res = self.mat
        if PROJECTION:  # make a projection
            # Qres_QxE_proj = np.delete(np.delete(self.mat, 2, axis=0), 2, axis=1)
            for i in (3, 2, 1, 0):
                if i not in axes:
                    q_res = ResoEllipsoid.quadric_proj(q_res, i)
            mat = q_res

        else:  # make a slice
            for i in (3, 2, 1, 0):
                if i not in axes:
                    q_res = np.delete(np.delete(q_res, i, axis=0), i, axis=1)
            mat = q_res

        return ResoEllipse(mat, centers, angle, axes_labels)

    def plot_ellipses(self):
        """Plot all 2D ellipses"""

        # fig = plt.figure()
        fig = plt.figure(figsize=(12, 8), constrained_layout=True)

        for i, indices in enumerate([(0, 3), (1, 3), (2, 3), (0, 1), (1, 2), (0, 2)]):
            ellipse_co = self.get_ellipse(axes=indices, PROJECTION=False)
            ellipse_inco = self.get_ellipse(axes=indices, PROJECTION=True)

            p = Plot2D()
            if indices == (2, 3):
                p.add_reso(ellipse_co, c="k", linestyle="solid", label="Coherent")
                p.add_reso(ellipse_inco, c="k", linestyle="dashed", label="Incoherent")

            else:
                p.add_reso(ellipse_co, c="k", linestyle="solid")
                p.add_reso(ellipse_inco, c="k", linestyle="dashed")

            ax = fig.add_subplot(
                int(f"23{i + 1}"),
                axes_class=Axes,
                grid_helper=p.grid_helper(ellipse_co.angle),
            )
            p.plot(ax)

        fig.suptitle(f"Q={self.hkl}, En={self.en} meV")
        # fig.tight_layout(pad=2)
