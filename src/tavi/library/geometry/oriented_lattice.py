"""
Calculate UB matrix and related quantities.

Follows Andrei Savici's UB Matrix Formlism used in mantid at
https://github.com/mantidproject/documents/blob/main/Design/UBMatriximplementationnotes.pdf version:March 06, 2011.
Equations listed in the comments refer to the document above.
"""

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
        powder: bool = False,
    ) -> None:
        """Initialize OrientedLattice class."""
        self._a = a
        self._b = b
        self._c = c
        self._alpha = alpha
        self._gamma = gamma
        self._beta = beta
        self._u_mat = u_mat
        self._B = self.calculate_B(self._a, self._b, self._c, self._alpha, self._beta, self._gamma)

    @property
    def a(self) -> None:
        """Get lattice parameters."""
        return self._a

    @a.setter
    def a(self, value: float) -> None:
        """Set lattice parameters and recalculate B matrix."""
        self._B = self.calculate_B(value, self._b, self._c, self._alpha, self._beta, self._gamma)
        self._a

    @property
    def b(self) -> float:
        """Get lattice parameters."""
        return self._b

    @b.setter
    def b(self, value: float) -> None:
        """Set lattice parameters and recalculate B matrix."""
        self._B = self.calculate_B(self._a, value, self._c, self._alpha, self._beta, self._gamma)
        self._b

    @property
    def c(self) -> float:
        """Get lattice parameters."""
        return self._c

    @c.setter
    def c(self, value: float) -> None:
        """Set lattice parameters and recalculate B matrix."""
        self._B = self.calculate_B(self._a, self._b, value, self._alpha, self._beta, self._gamma)
        self._a

    @property
    def alpha(self) -> float:
        """Get lattice parameters."""
        return self._alpha

    @alpha.setter
    def alpha(self, value: float) -> None:
        """Set lattice parameters and recalculate B matrix."""
        self._B = self.calculate_B(self._a, self._b, self._c, value, self._beta, self._gamma)
        self._alpha

    @property
    def beta(self) -> float:
        """Get lattice parameters."""
        return self._beta

    @beta.setter
    def beta(self, value: float) -> None:
        """Set lattice parameters and recalculate B matrix."""
        self._B = self.calculate_B(self._a, self._b, self._c, self._alpha, value, self._gamma)
        self._beta

    @property
    def gamma(self) -> float:
        """Get lattice parameters."""
        return self._gamma

    @gamma.setter
    def gamma(self, value: float) -> None:
        """Set lattice parameters and recalculate B matrix."""
        self._B = self.calculate_B(self._a, self._b, self._c, self._alpha, self._beta, value)
        self._gamma

    @property
    def B(self) -> np.ndarray:
        """Get B matrix."""
        return self._B

    @B.setter
    def B(self, mat: np.ndarray) -> None:
        """Set B matrix. Should also update lattice parameters."""
        self._B = mat
        # TODO
        # implement refine lattice parameter algo

    @property
    def U(self) -> np.ndarray:
        """Get U matrix."""
        return self._u_mat

    @U.setter
    def U(self, uv: tuple[np.ndarray, np.ndarray]) -> np.ndarray:
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
        u_mat = T_epsilon @ T_c.T
        return u_mat @ self._B

    @property
    def getUB(self) -> np.ndarray:
        """Get UB matrix."""
        return self._u_mat @ self._B

    def overwrite_U(self, u_mat: np.ndarray) -> None:
        """Manual overwrite U."""
        self._u_mat = u_mat

    def calculate_B(self, a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
        """Calculate b matrix from lattice parameters.Eq.52."""
        cos_alpha = np.cos(np.radians(alpha))
        sin_alpha = np.sin(np.radians(alpha))

        cos_beta = np.cos(np.radians(beta))
        sin_beta = np.sin(np.radians(beta))

        cos_gamma = np.cos(np.radians(gamma))
        sin_gamma = np.sin(np.radians(gamma))

        Vabg = np.sqrt(1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 2 * cos_alpha * cos_beta * cos_gamma)

        a_star = sin_alpha / (a * Vabg)
        b_star = sin_beta / (b * Vabg)
        c_star = sin_gamma / (c * Vabg)

        cos_alpha_star = (cos_beta * cos_gamma - cos_alpha) / (sin_beta * sin_gamma)
        cos_beta_star = (cos_gamma * cos_alpha - cos_beta) / (sin_gamma * sin_alpha)
        sin_beta_star = np.sqrt(1 - cos_beta_star**2)
        cos_gamma_star = (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta)
        sin_gamma_star = np.sqrt(1 - cos_gamma_star**2)

        cos_gamma_star = (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta)
        sin_gamma_star = np.sqrt(1 - cos_gamma_star**2)

        B = np.array(
            [
                [a_star, b_star * cos_gamma_star, c_star * cos_beta_star],
                [0, b_star * sin_gamma_star, -c_star * sin_beta_star * cos_alpha],
                [0, 0, 1.0 / c],
            ]
        )
        return B

    def get_uv(self) -> tuple[np.ndarray, np.ndarray]:
        """
        If ubmatrix is given, return uv vector.

        Given in TAS and mantid case, u is along [0, 0, 1], v is along [1, 0, 0].
        """
        u = np.matmul(np.linalg.inv(self._u_mat), np.array([0, 0, 1]))
        v = np.matmul(np.linalg.inv(self._u_mat), np.array([1, 0, 0]))
        return (u, v)
