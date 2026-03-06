"""
Calculate UB matrix and related quantities.

Follows Andrei Savici's UB Matrix Formlism used in mantid at
https://github.com/mantidproject/documents/blob/main/Design/UBMatriximplementationnotes.pdf version:March 06, 2011.
Equations listed in the comments refer to the document above.
"""

from typing import Tuple

import numpy as np
from pydantic import BaseModel

from tavi.library.experiment.utilities import SE2K, get_angle_from_triangle, q_norm_from_hkl


class OrientedLattice(BaseModel):
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

    @property
    def a(self) -> float:
        """Get lattice parameters."""
        return self._a

    @a.setter
    def a(self, value: float) -> None:
        """Return updated lattice parameters value."""
        self._a = value
        self.calculate_B()
        self.update_lattice_parameters_from_B()
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
        self.update_lattice_parameters_from_B()
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
        self.update_lattice_parameters_from_B()
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
        self.update_lattice_parameters_from_B()
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
        self.update_lattice_parameters_from_B()
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
        self.update_lattice_parameters_from_B()
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
    def UB(self, mat: np.ndarray) -> np.ndarray:
        """Set UB, recalculate everything."""
        self._UB = mat
        self.update_B_from_UB()
        self._u_mat = self._UB @ np.linalg.inv(self._B)

    @property
    def a_star(self) -> np.ndarray:
        """Calculate a* from ub."""
        return self._UB @ np.ndarray([1, 0, 0]).T

    @property
    def b_star(self) -> np.ndarray:
        """Calculate a* from ub."""
        return self._UB @ np.ndarray([0, 1, 0]).T

    @property
    def c_star(self) -> np.ndarray:
        """Calculate a* from ub."""
        return self._UB @ np.ndarray([0, 1, 0]).T

    def update_B_from_UB(self) -> None:
        """Calculate B matri from G*. Also update lattice parameters."""
        G_star = self._UB.T @ self._UB
        a_star = np.sqrt(G_star[0, 0])
        b_star = np.sqrt(G_star[1, 1])
        c_star = np.sqrt(G_star[2, 2])
        alpha_star = np.arccos(G_star[1, 2] / (b_star * c_star))
        beta_star = np.arccos(G_star[0, 2] / (a_star * c_star))
        gamma_star = np.arccos(G_star[0, 1] / (a_star * b_star))

        G = np.linalg.inv(G_star)
        self._a = np.sqrt(G[0, 0])
        self._b = np.sqrt(G[1, 1])
        self._c = np.sqrt(G[2, 2])
        self._alpha = np.radians(np.arccos(G[1, 2] / (self._b * self._c)))
        self._beta = np.radians(np.arccos(G[0, 2] / (self._a * self._c)))
        self._gamma = np.radians(np.arccos(G[0, 1] / (self._a * self._b)))

        self._B = np.array(
            [
                [a_star, b_star * np.cos(gamma_star), c_star * np.cos(beta_star)],
                [0, b_star * np.sin(gamma_star), -c_star * np.sin(beta_star) * np.cos(self._alpha)],
                [0, 0, 1.0 / self._c],
            ]
        )

    def get_ub_from_uv(self, uv: tuple[np.ndarray, np.ndarray]) -> np.ndarray:
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
        return self._u_mat @ self._B

    def getUB(self) -> np.ndarray:
        """Get UB matrix."""
        return self._u_mat @ self._B

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

        cos_gamma_star = (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta)
        sin_gamma_star = np.sqrt(1 - cos_gamma_star**2)
        self._B = np.array(
            [
                [a_star, b_star * cos_gamma_star, c_star * cos_beta_star],
                [0, b_star * sin_gamma_star, -c_star * sin_beta_star * cos_alpha],
                [0, 0, 1.0 / self._c],
            ]
        )
        print(self._u_mat, self._B)

    def update_lattice_parameters_from_B(self) -> None:
        """Calculate lattice parameters from B matrix."""
        G_star = self._B.T @ self._B
        G = np.linalg.inv(G_star)
        self._a = np.sqrt(G[0, 0])
        self._b = np.sqrt(G[1, 1])
        self._c = np.sqrt(G[2, 2])
        self._alpha = np.radians(np.arccos(G[1, 2] / (self._b * self._c)))
        self._beta = np.radians(np.arccos(G[0, 2] / (self._a * self._c)))
        self._gamma = np.radians(np.arccos(G[0, 1] / (self._a * self._b)))

    def get_uv(self) -> tuple[np.ndarray, np.ndarray]:
        """
        If ubmatrix is given, return uv vector.

        Given in TAS and mantid case, u is along [0, 0, 1], v is along [1, 0, 0].
        """
        u = np.linalg.inv(self._u_mat @ self._B) @ np.array([0, 0, 1]).T
        v = np.linalg.inv(self._u_mat @ self._B) @ np.array([1, 0, 0]).T
        return (u, v)

    def two_theta_from_hkl(self, hkl: Tuple[float, float, float], ei: float, ef: float) -> float:
        """Compute two_theta from hkl."""
        ki = SE2K(ei)
        kf = SE2K(ef)
        q_norm = q_norm_from_hkl(hkl, self._B)
        two_theta = get_angle_from_triangle(ki, kf, q_norm)  # in radians
        return np.degrees(two_theta)

    def calculate_d_star(self, hkl: Tuple[float, float, float]) -> float:
        """Calculate d* from hkl."""
        return np.abs(self._B @ np.array(hkl).T)

    def calculate_d(self, hkl: Tuple[float, float, float]) -> float:
        """Calculate d from hkl."""
        return 1 / self.calculate_d_star(hkl)
