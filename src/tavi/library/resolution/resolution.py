"""Resolution manager."""

from logging import raiseExceptions
from typing import Literal, Optional, Tuple

import numpy as np

from tavi.library.experiment.experiment import Experiment
from tavi.library.experiment.utilities import SE2K, get_angle_vec, quadric_proj
from tavi.library.geometry.sample import Sample
from tavi.library.Instrument.instrument import Instrument

MODEL_CHOICES = Literal["Cooper-Nathans"]


class Resolution:
    """Resolution manager class."""

    def __init__(
        self,
        model: MODEL_CHOICES,
        instrument: Instrument,
        sample: Sample,
        experiment: Experiment,
        scan_idx: Optional[int] = 0,
        pt_idx: Optional[int] = 0,
        axes: Optional[Tuple] = ((1, 0, 0), (0, 1, 0), (0, 0, 1), "e"),
    ) -> None:
        """
        Init components.

        Args:
            model: Resolution model to use. Currently supports ``"Cooper-Nathans"``.
            instrument: Triple-axis instrument configuration.
            sample: Sample associated with the experiment.
            experiment: Experiment context providing scan geometry.
            scan_idx: Scan index within the experiment.
            pt_idx: Point index within the scan.
            axes: Projection axes for the resolution ellipsoid; last entry is ``"e"`` for energy.

        """
        if model == "Cooper-Nathans":
            from tavi.library.resolution.cooper_nathans import CooperNathans

            self.model = CooperNathans()
        else:
            raise ValueError(f"Unknown resolution model '{model}'. Choose from: 'CooperNathans'.")
        self.instrument = instrument
        self.sample = sample
        self.experiment = experiment
        self.scan_idx = scan_idx
        self.pt_idx = pt_idx
        self.axes = axes

    def get_resolution(self, hkl: Tuple[float, float, float], ei: float, ef: float) -> Tuple[np.ndarray, float]:
        """Get resolution matrix and r0 from a selected model at hkl."""
        q_norm = self.sample.ol.q_norm_from_hkl(hkl)
        ki, kf = SE2K(ei), SE2K(ef)
        psi = self.experiment.get_psi(q_norm, ei, ef) * (
            -self.instrument.goni.sense
        )  # negative sign ensures opposite sign of s2
        two_theta = self.experiment.get_two_theta(q_norm, ei, ef) * (self.instrument.goni.sense)
        theta_m = self.instrument.monochromater.theta(ei) * (self.instrument.monochromater.sense)
        theta_a = self.instrument.analyzer.theta(ef) * (self.instrument.analyzer.sense)  # set sense
        res_q = self.model.resolution_matrix(
            self.instrument, self.sample, q_norm, ki, kf, psi, two_theta, theta_m, theta_a
        )
        if not self.axes:
            return res_q
        else:
            res_ellipsoid = ResoEllipsoid(res_q[0], res_q[1], self.axes)
            r_mat = self.r_matrix_with_minimal_tilt(hkl, ei, ef)
            res_proj = res_ellipsoid.project_to_frame(r_mat, psi, self.sample.ol.UB)
            return res_proj

    def get_axes_angles(self) -> Tuple[float, float, float]:
        """Compute the angles between the basis of axes."""
        tol = 1e-8
        self.angles = (90, 90, 90)
        if self.axes is None:
            return self.angles
        elif self.axes == ((1, 0, 0), (0, 1, 0), (0, 0, 1), "e"):  # HKL
            alpha_star, beta_star, gamma_star = (
                self.sample.ol.alpha_star,
                self.sample.ol.beta_star,
                self.sample.ol.gamma_star,
            )
            self.angles = (np.degrees(gamma_star), np.degrees(alpha_star), np.degrees(beta_star))
        else:
            p1, p2, p3 = np.array([item for item in self.axes if isinstance(item, tuple) and len(item) == 3])
            reciprocal_vecs = self.sample.ol.reciprocal_vectors()
            # Eq. 64
            v1 = np.sum([p1[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)
            v2 = np.sum([p2[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)
            v3 = np.sum([p3[i] * vec for (i, vec) in enumerate(reciprocal_vecs)], axis=0)

            if np.dot(v1, np.cross(v2, v3)) < tol:
                raise ValueError("Projection is left handed! Please use right-handed projection")
            if np.abs(np.dot(v1, np.cross(v2, v3))) < tol:
                raise ValueError("Projection vectors need to be non-coplanar.")
            self.angles = (get_angle_vec(v1, v2), get_angle_vec(v2, v3), get_angle_vec(v3, v1))
        return self.angles

    def get_ellipse(
        self, res_mat: np.ndarray, ellipse_axes: tuple[int, int] = (0, 1), PROJECTION: bool = False, ORIGIN: bool = True
    ) -> np.ndarray:
        """Generate 2D resolution ellipse."""
        res_2d = res_mat

        for idx in (3, 2, 1, 0):
            # if idx not selected, handle it
            if idx not in ellipse_axes:
                # make an integral
                if PROJECTION:
                    res_2d = quadric_proj(res_2d, idx)
                # make a slice/cut
                else:
                    res_2d = np.delete(np.delete(res_2d, idx, axis=0), idx, axis=1)
        axes_angles = self.get_axes_angles()
        match tuple(ellipse_axes):
            case (0, 1) | (1, 0):
                axes_angle = np.round(axes_angles[0], 2)
            case (1, 2) | (2, 1):
                axes_angle = np.round(axes_angles[1], 2)
            case (0, 2) | (2, 0):
                axes_angle = np.round(axes_angles[2], 2)
            case _:
                axes_angle = 90.0
        return res_2d, axes_angle

    def r_matrix_with_minimal_tilt(
        self,
        hkl: Tuple[float, float, float],
        ei: float,
        ef: float,
    ) -> np.ndarray:
        """
        Calculate R matrix when the tilt from the scattering plane is minimal.

        Args:
            hkl: Miller indices of the reflection.
            ei: incident energy, in meV.
            ef: final energy, in meV.

        """
        ki = SE2K(ei)
        kf = SE2K(ef)
        ub_mat = self.sample.ol.UB
        q_norm = self.sample.ol.q_norm_from_hkl(hkl)
        plane_normal = self.sample.ol.plane_normal
        in_plane_ref = self.sample.ol.in_plane_ref
        if not plane_normal or not in_plane_ref:
            raise ValueError("Both plane_normal and in_plane_ref must be set.")
        if plane_normal is None or in_plane_ref is None:
            raiseExceptions("plane_normal and in_plane_ref must be set.")

        two_theta = self.experiment.get_two_theta(q_norm, ei, ef) * (self.instrument.goni.sense)
        # with minimal tilt, we are considering scenario described above Eq.114
        q_lab1 = np.array([-kf * np.sin(two_theta), 0, ki - kf * np.cos(two_theta)]) / q_norm
        q_lab2 = np.array([ki - kf * np.cos(two_theta), 0, kf * np.sin(two_theta)]) / q_norm
        q_lab3 = np.array([0, 1, 0])

        Q_lab = np.array([q_lab1, q_lab2, q_lab3]).T
        tol = 1e-5

        # Eq.106
        t1 = ub_mat @ np.array(hkl)
        if np.abs(t1 @ plane_normal) < tol:
            # t1 in plane
            t3 = plane_normal
            t2 = np.cross(t3, t1)
        elif np.linalg.norm(np.cross(plane_normal, t1)) < tol:
            # calculate R when t1 is along plane_normal. Easiest way to construct three perp
            # axes are use t2 as in-plane axis. This already considers a 90 degree rotation
            # to bring t1 in-plane.
            t2 = in_plane_ref  # set t2 in plane
            t3 = np.cross(t1, t2)  # t3 along y
        else:
            # t1 not in plane. np.cross(plane_normal, t1) already incorporate a minimal rotation along
            # t2 to bring t1 in-plane
            t2 = np.cross(plane_normal, t1)
            t3 = np.cross(t1, t2)

        T_mat = np.array(
            [
                t1 / np.linalg.norm(t1),
                t2 / np.linalg.norm(t2),
                t3 / np.linalg.norm(t3),
            ]
        ).T
        r_mat = Q_lab @ np.linalg.inv(T_mat)
        return r_mat
