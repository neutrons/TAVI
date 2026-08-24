"""
Calculate UB matrix and related quantities.

Follows Andrei Savici's UB Matrix Formalism used in mantid at
https://github.com/mantidproject/documents/blob/main/Design/UBMatriximplementationnotes.pdf version:March 06, 2011.
Equations listed in the comments refer to the document above.
"""

from logging import raiseExceptions
from typing import Tuple

import numpy as np


class OrientedLattice:
    """Oritented lattice class."""

    def __init__(
        self,
        a: float = 1,
        b: float = 1,
        c: float = 1,
        alpha: float = 90,
        beta: float = 90,
        gamma: float = 90,
        u_mat: np.ndarray = np.eye(3),
        plane_normal: np.ndarray | None = None,
        in_plane_ref: np.ndarray | None = None,
        powder: bool = False,
    ) -> None:
        """Initialize OrientedLattice class."""
        self._a = a
        self._b = b
        self._c = c
        self._alpha = alpha
        self._beta = beta
        self._gamma = gamma
        self._powder = powder
        self._u_mat = u_mat
        self._B = None
        self.calculate_B()
        self._UB = self._u_mat @ self._B
        self._plane_normal = plane_normal
        self._in_plane_ref = in_plane_ref

    @property
    def a(self) -> float:
        """Get lattice parameters."""
        return self._a

    @a.setter
    def a(self, value: float) -> None:
        """Return updated lattice parameters value."""
        self._a = value
        self.calculate_B()
        self._UB = self._u_mat @ self._B

    @property
    def b(self) -> float:
        """Get lattice parameters."""
        return self._b

    @b.setter
    def b(self, value: float) -> None:
        """Return updated lattice parameters value."""
        self._b = value
        self.calculate_B()
        self._UB = self._u_mat @ self._B

    @property
    def c(self) -> float:
        """Get lattice parameters."""
        return self._c

    @c.setter
    def c(self, value: float) -> None:
        """Return updated lattice parameters value."""
        self._c = value
        self.calculate_B()
        self._UB = self._u_mat @ self._B

    @property
    def alpha(self) -> float:
        """Get lattice parameters."""
        return self._alpha

    @alpha.setter
    def alpha(self, value: float) -> None:
        """Return updated lattice parameters value."""
        self._alpha = value
        self.calculate_B()
        self._UB = self._u_mat @ self._B

    @property
    def beta(self) -> float:
        """Get lattice parameters."""
        return self._beta

    @beta.setter
    def beta(self, value: float) -> None:
        """Return updated lattice parameters value."""
        self._beta = value
        self.calculate_B()
        self._UB = self._u_mat @ self._B

    @property
    def gamma(self) -> float:
        """Get lattice parameters."""
        return self._gamma

    @gamma.setter
    def gamma(self, value: float) -> None:
        """Return updated lattice parameters value."""
        self._gamma = value
        self.calculate_B()
        self._UB = self._u_mat @ self._B

    @property
    def B(self) -> np.ndarray:
        """Get B matrix."""
        return self._B

    @B.setter
    def B(self, mat: np.ndarray) -> None:
        """Update B matrix."""
        self._B = mat
        self._UB = self._u_mat @ self._B
        self.update_lattice_parameters_from_B()

    @property
    def U(self) -> np.ndarray:
        """Get U matrix."""
        return self._u_mat

    @U.setter
    def U(self, value: np.ndarray) -> None:
        """Set U matrix."""
        self._u_mat = value
        self._UB = self._u_mat @ self._B

    @property
    def UB(self) -> np.ndarray:
        """Get UB matrix."""
        return self._UB

    @UB.setter
    def UB(self, mat: np.ndarray) -> None:
        """Set UB, recalculate everything."""
        self._UB = mat
        self.update_B_from_UB()
        self._u_mat = self._UB @ np.linalg.inv(self._B)

    @property
    def a_star(self) -> np.ndarray:
        """Calculate a* from ub."""
        G_star = self._UB.T @ self._UB
        return np.sqrt(G_star[0, 0])

    @property
    def b_star(self) -> np.ndarray:
        """Calculate a* from ub."""
        G_star = self._UB.T @ self._UB
        return np.sqrt(G_star[1, 1])

    @property
    def c_star(self) -> np.ndarray:
        """Calculate a* from ub."""
        G_star = self._UB.T @ self._UB
        return np.sqrt(G_star[2, 2])

    @property
    def plane_normal(self) -> np.ndarray:
        """Return plane normal vector."""
        return self._plane_normal

    @plane_normal.setter
    def plane_normal(self, val: np.ndarray) -> None:
        """Set plane normal vector."""
        self._plane_normal = val

    @property
    def in_plane_ref(self) -> np.ndarray:
        """Return plane normal vector."""
        return self._in_plane_ref

    @in_plane_ref.setter
    def in_plane_ref(self, val: np.ndarray) -> None:
        """Set plane normal vector."""
        self._in_plane_ref = val

    def update_B_from_UB(self) -> None:
        """Calculate B matrix from G*. Also update lattice parameters."""
        G_star = self._UB.T @ self._UB
        a_star = np.sqrt(G_star[0, 0])
        b_star = np.sqrt(G_star[1, 1])
        c_star = np.sqrt(G_star[2, 2])
        self.alpha_star = np.arccos(G_star[1, 2] / (b_star * c_star))
        self.beta_star = np.arccos(G_star[0, 2] / (a_star * c_star))
        self.gamma_star = np.arccos(G_star[0, 1] / (a_star * b_star))

        G = np.linalg.inv(G_star)
        self._a = np.sqrt(G[0, 0])
        self._b = np.sqrt(G[1, 1])
        self._c = np.sqrt(G[2, 2])
        self._alpha = np.degrees(np.arccos(G[1, 2] / (self._b * self._c)))
        self._beta = np.degrees(np.arccos(G[0, 2] / (self._a * self._c)))
        self._gamma = np.degrees(np.arccos(G[0, 1] / (self._a * self._b)))

        self._B = np.array(
            [
                [a_star, b_star * np.cos(self.gamma_star), c_star * np.cos(self.beta_star)],
                [
                    0,
                    b_star * np.sin(self.gamma_star),
                    -c_star * np.sin(self.beta_star) * np.cos(np.radians(self._alpha)),
                ],
                [0, 0, 1.0 / self._c],
            ]
        )

    def set_u_from_uv(self, uv: tuple[np.ndarray, np.ndarray]) -> None:
        """Compute UB matrix from u,v and lattice parameters. Eq.76-81."""
        u, v = uv
        t1_c = self._B @ u
        t3_c = np.cross(self._B @ u, self._B @ v)
        t2_c = np.cross(t3_c, t1_c)
        T_c = np.array(
            [
                t1_c / np.linalg.norm(t1_c),
                t2_c / np.linalg.norm(t2_c),
                t3_c / np.linalg.norm(t3_c),
            ]
        ).T
        T_epsilon = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]).T
        T_c = np.array(
            [
                t1_c / np.linalg.norm(t1_c),
                t2_c / np.linalg.norm(t2_c),
                t3_c / np.linalg.norm(t3_c),
            ]
        ).T
        T_epsilon = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]).T
        self._u_mat = T_epsilon @ T_c.T
        self._UB = self._u_mat @ self._B

    def calculate_B(self) -> None:
        """Calculate b matrix from lattice parameters.Eq.52."""
        cos_alpha = np.cos(np.radians(self._alpha))
        sin_alpha = np.sin(np.radians(self._alpha))

        cos_beta = np.cos(np.radians(self._beta))
        sin_beta = np.sin(np.radians(self._beta))

        cos_gamma = np.cos(np.radians(self._gamma))
        sin_gamma = np.sin(np.radians(self._gamma))

        Vabg = np.sqrt(1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 2 * cos_alpha * cos_beta * cos_gamma)

        a_star = sin_alpha / (self._a * Vabg)
        b_star = sin_beta / (self._b * Vabg)
        c_star = sin_gamma / (self._c * Vabg)

        cos_alpha_star = (cos_beta * cos_gamma - cos_alpha) / (sin_beta * sin_gamma)
        cos_beta_star = (cos_gamma * cos_alpha - cos_beta) / (sin_gamma * sin_alpha)
        sin_beta_star = np.sqrt(1 - cos_beta_star**2)
        cos_gamma_star = (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta)
        sin_gamma_star = np.sqrt(1 - cos_gamma_star**2)

        self.alpha_star = np.arccos(cos_alpha_star)
        self.beta_star = np.arccos(cos_beta_star)
        self.gamma_star = np.arccos(cos_gamma_star)

        cos_gamma_star = (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta)
        sin_gamma_star = np.sqrt(1 - cos_gamma_star**2)
        self._B = np.array(
            [
                [a_star, b_star * cos_gamma_star, c_star * cos_beta_star],
                [0, b_star * sin_gamma_star, -c_star * sin_beta_star * cos_alpha],
                [0, 0, 1.0 / self._c],
            ]
        )

    def update_lattice_parameters_from_B(self) -> None:
        """Calculate lattice parameters from B matrix."""
        G_star = self._B.T @ self._B
        G = np.linalg.inv(G_star)
        self._a = np.sqrt(G[0, 0])
        self._b = np.sqrt(G[1, 1])
        self._c = np.sqrt(G[2, 2])
        self._alpha = np.degrees(np.arccos(G[1, 2] / (self._b * self._c)))
        self._beta = np.degrees(np.arccos(G[0, 2] / (self._a * self._c)))
        self._gamma = np.degrees(np.arccos(G[0, 1] / (self._a * self._b)))

    def get_uv(self) -> tuple[np.ndarray, np.ndarray]:
        """
        If ubmatrix is given, return uv vector.

        Given in TAS and mantid case, u is along [0, 0, 1], v is along [1, 0, 0].
        """
        u = np.linalg.inv(self._u_mat @ self._B) @ np.array([0, 0, 1]).T
        v = np.linalg.inv(self._u_mat @ self._B) @ np.array([1, 0, 0]).T
        return (u, v)

    def calculate_d_star(self, hkl: Tuple[float, float, float]) -> float:
        """Calculate d* from hkl."""
        return np.abs(self._B @ np.array(hkl).T)

    def calculate_d(self, hkl: Tuple[float, float, float]) -> float:
        """Calculate d from hkl."""
        return 1 / self.calculate_d_star(hkl)

    def get_angle_between_two_hkl(self, hkl1: Tuple[float, float, float], hkl2: Tuple[float, float, float]) -> float:
        """Calculate the angle between two hkl positions. Good to have for moving from one peak to another."""
        p1, p2 = np.array(hkl1), np.array(hkl2)
        cos_angle = (p1 @ p2) / (np.linalg.norm(p1) * np.linalg.norm(p2))
        return np.degrees(cos_angle)

    def q_norm_from_hkl(self, hkl: tuple[float, float, float]) -> float:
        """
        Return norm of q for given (h,k,l).

        Note:
            Either b_mat or ub_mat would work, since U^T.U=U^-1.U=1

        """
        q_norm = 2 * np.pi * np.linalg.norm(self._B @ np.array(hkl))
        return q_norm

    def real_space_vectors(self) -> tuple:
        """Compute real lattice vector based on Eq.9-21."""
        cos_alpha = np.cos(np.radians(self._alpha))
        sin_alpha = np.sin(np.radians(self._alpha))

        cos_beta = np.cos(np.radians(self._beta))
        sin_beta = np.sin(np.radians(self._beta))

        cos_gamma = np.cos(np.radians(self._gamma))
        sin_gamma = np.sin(np.radians(self._gamma))

        self.Vabg = np.sqrt(1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 2 * cos_alpha * cos_beta * cos_gamma)

        a_vec = np.array([self._a, 0, 0])
        b_vec = np.array([self._b * cos_gamma, self._b * sin_gamma, 0])

        c_vec = np.array(
            [
                self._c * cos_beta,
                self._c * (cos_alpha - cos_gamma * cos_beta) / sin_gamma,
                self._c * self.Vabg / sin_gamma,
            ]
        )
        return a_vec, b_vec, c_vec

    def reciprocal_vectors(self)->tuple:
        """Compute reciprocal lattice vectors."""
        a_vec, b_vec, c_vec = self.real_space_vectors()
        V = self.Vabg*self._a*self._b*self._c
        prefactor = 2 * np.pi / V
        a_star_vec = np.cross(b_vec, c_vec) * prefactor
        b_star_vec = np.cross(c_vec, a_vec) * prefactor
        c_star_vec = np.cross(a_vec, b_vec) * prefactor

        return (a_star_vec, b_star_vec, c_star_vec)


    def rot_matrix_with_minimal_tilt(
        self,
        hkl: Tuple[float, float, float],
        ki: float,
        kf: float,
        two_theta: float,
    ) -> np.ndarray:
        """
        Calculate R matrix when the tilt from the scattering plane is minimal.

        Args:
            hkl: Miller indices of the reflection.
            ki: incident wavevector, in inverse Angstrom.
            kf: final wavevector, in inverse Angstrom.
            two_theta: scattering angle with instrument.goni.sense considered.

        """
        ub_mat = self.UB
        q_norm = self.q_norm_from_hkl(hkl)
        plane_normal = self.plane_normal
        in_plane_ref = self.in_plane_ref
        if plane_normal is None or in_plane_ref is None:
            raiseExceptions("plane_normal and in_plane_ref must be set.")

        # two_theta = self.experiment.get_two_theta(q_norm, ei, ef) * (self.instrument.goni.sense)
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
        rot_mat = Q_lab @ np.linalg.inv(T_mat)
        return rot_mat
