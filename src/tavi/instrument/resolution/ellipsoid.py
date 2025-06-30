import re
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axisartist import Axes

from tavi.instrument.resolution.ellipse import ResoEllipse
from tavi.plotter import Plot2D
from tavi.utilities import get_angle_vec, labels_from_projection, sig2fwhm


class ResoEllipsoid(object):
    """Manage the 4D resolution ellipoid

    Attributs:
        frame (str): "q", "hkl", or "proj"
        projection (tuple): three non-coplanar vectors
        q_vec (tuple): momentum transfer (h', k', l') in the coordinate system specified by projection
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
        axes: Optional[tuple] = ((1, 0, 0), (0, 1, 0), (0, 0, 1), "en"),
        reso_mat: Optional[np.ndarray] = None,
        r0: Optional[float] = None,
    ) -> None:
        *hkl, en = hkle
        self.hkl: tuple = tuple(hkl)
        self.en: float = en
        self.q = np.linalg.norm(instrument.sample.b_mat @ self.hkl) * 2 * np.pi

        self.axes: Optional[tuple] = axes
        self.mat = reso_mat
        self.r0 = r0
        self.method: Optional[str] = None
        self.instrument_params: Optional[str] = None

        if axes is None:
            self.frame = "q"
            self.angles = (90.0, 90.0, 90.0)
            self.q_vec = (self.q, 0.0, 0.0)
            self.axes_labels = labels_from_projection(axes)
        else:
            self._project_to_frame(instrument)

    def __repr__(self):
        h, k, l = self.hkl
        return f"ResoEllipsoid at hkl=({h:.4g}, {k:.4g}, {l:.4g}), en={self.en:.4g} meV"

    def __str__(self):
        h, k, l = self.hkl
        r_mat_str = pd.DataFrame(
            self.mat,
            index=[re.sub(r"\([^()]*\)$", "", label) for label in self.axes_labels],
            columns=[re.sub(r"\([^()]*\)$", "", label) for label in self.axes_labels],
        ).to_string(float_format="{:.4g}".format)
        summary_str = [
            f"Resolution ellipsoid centered at (h,k,l)=({h:.4g}, {k:.4g}, {l:.4g}) (r.l.u.), en={self.en:.4g} meV. (Q = {self.q:.4g} A^-1)",
            f"Calculated using {self.method} method.",
            "Projection axes are " + ", ".join(self.axes_labels),
            "resolution matrix R=",
            r_mat_str,
            f"prefactor R0={self.r0:.4g}",
            "Coherent FWHMs are: ",
            f"{self.coh_fwhms(0):.4g} in {self.axes_labels[0]}, ",
            f"{self.coh_fwhms(1):.4g} in {self.axes_labels[1]}, ",
            f"{self.coh_fwhms(2):.4g} in {self.axes_labels[2]}, ",
            f"{self.coh_fwhms(3):.4g} in {self.axes_labels[3]}, ",
            "Incoherent FWHMs are: ",
            f"{self.incoh_fwhms(0):.4g} in {self.axes_labels[0]}, ",
            f"{self.incoh_fwhms(1):.4g} in {self.axes_labels[1]}, ",
            f"{self.incoh_fwhms(2):.4g} in {self.axes_labels[2]}, ",
            f"{self.incoh_fwhms(3):.4g} in {self.axes_labels[3]}, ",
            "",
        ]
        return "\n".join(summary_str) + self.instrument_params

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
        if self.axes == ((1, 0, 0), (0, 1, 0), (0, 0, 1), "en"):  # HKL
            self.frame = "hkl"
            self.q_vec = self.hkl
            *_, alpha_star, beta_star, gamma_star = instrument.sample.reciprocal_latt_params
            self.angles = (gamma_star, alpha_star, beta_star)

            conv_mat_4d = np.eye(4)
            conv_mat_4d[0:3, 0:3] = mat_lab_to_local @ conv_mat
            self.mat = conv_mat_4d.T @ self.mat @ conv_mat_4d
            self.axes_labels = labels_from_projection()

        else:  # customized projection
            p1, p2, p3 = np.array([item for item in self.axes if isinstance(item, tuple) and len(item) == 3])
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
            self.q_vec = hkl_prime
            self.angles = (get_angle_vec(v1, v2), get_angle_vec(v2, v3), get_angle_vec(v3, v1))

            mat_w = np.array([p1, p2, p3]).T

            conv_mat_4d = np.eye(4)
            conv_mat_4d[0:3, 0:3] = mat_lab_to_local @ conv_mat @ mat_w
            self.mat = conv_mat_4d.T @ self.mat @ conv_mat_4d

            # roll axes
            en_axis = self.axes.index("en")
            new_idx = [0, 1, 2]
            new_idx.insert(en_axis, 3)
            self.mat = self.mat[new_idx, :][:, new_idx]

            self.axes_labels = labels_from_projection(self.axes)

    # TODO
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
        try:
            en_idx = self.axes.index("en")
        except AttributeError:
            en_idx = 3
        qe_list = np.insert(self.q_vec, en_idx, self.en)
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

    def plot_ellipses(self, fig):
        index_pairs = [(0, 3), (1, 3), (2, 3), (0, 1), (1, 2), (0, 2)]

        for i, indices in enumerate(index_pairs):
            ellipse_inco = self.get_ellipse(axes=indices, PROJECTION=True)
            ellipse_co = self.get_ellipse(axes=indices, PROJECTION=False)

            p = Plot2D()
            if indices == (2, 3):
                p.add_reso(ellipse_inco, c="k", linestyle="dashed", label="Incoherent")
                p.add_reso(ellipse_co, c="k", linestyle="solid", label="Coherent")

            else:
                p.add_reso(ellipse_inco, c="k", linestyle="dashed")
                p.add_reso(ellipse_co, c="k", linestyle="solid")

            ax = fig.add_subplot(
                2,
                3,
                i + 1,
                axes_class=Axes,
                grid_helper=p.grid_helper(ellipse_co.angle, nbins=(5, 5)),
            )
            x_fwhm, y_fwhm = [self.incoh_fwhms(idx) * 0.6 for idx in indices]
            x_cen, y_cen = ellipse_co.centers
            p.xlim = (x_cen - x_fwhm, x_cen + x_fwhm)
            p.ylim = (y_cen - y_fwhm, y_cen + y_fwhm)
            p.plot(ax)
            # ax.axis["bottom"].major_ticklabels.set_text("$\\mathdefault{2}$")  # optional rotation

        h, k, l = self.hkl
        fig.suptitle(f"Q=({h:.4g}, {k:.4g}, {l:.4g}), En={self.en} meV")
        fig.tight_layout(pad=2)

    def plot(self):
        """Plot all 2D ellipses"""

        fig = plt.figure(figsize=(12, 8), constrained_layout=True)
        self.plot_ellipses(fig)
