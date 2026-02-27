"""
Calculate UB matrix and related quantities. 
Follows Andrei Savici's UB Matrix Formlism used in mantid at 
https://github.com/mantidproject/documents/blob/main/Design/UBMatriximplementationnotes.pdf.
Equations listed in the comments refer to the document above.
"""

import numpy as np
from tavi.library.utilities import SE2K


def q_lab(ei: float, ef: float, theta: float, phi: float) -> np.ndarray:
    """
    Return qlab matrix. Follow's mantid convention. Eq.8.

    Args:
        ei: incident energy, in meV.
        ef: final energy, in meV.
        theta: angle between scattered beam and z axis, indegree.
        phi: the angle between scattered beam's projection onto the (x,y) plane and x axis.

    """
    ki = SE2K(ei)
    kf = SE2K(ef)
    return np.array([
        -kf*np.sin(np.radians(theta))*np.cos(np.radians(phi)),
        -kf*np.sin(np.radians(theta))*np.sin(np.radians(phi)),
        ki - kf*np.cos(np.radians(theta))
    ])

def ub_to_uv(ub_mat:np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """If ubmatrix is given, return uv vector. Given in TAS and mantid case, u is along [0, 0, 1], v is along [1, 0, 0]"""
    u = np.matmul(np.linalg.inv(ub_mat), np.array([0, 0, 1]))
    v = np.matmul(np.linalg.inv(ub_mat), np.array([1, 0, 0]))
    return (u, v)


def b_mat(lattice_params: tuple[float, float, float, float, float, float]) -> np.ndarray:
    """Calculate b matrix from lattice parameters.Eq.52."""
    a, b, c, alpha, beta, gamma = lattice_params
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

    B = np.array(
        [
            [a_star, b_star * cos_gamma_star, c_star * cos_beta_star],
            [0, b_star * sin_gamma_star, -c_star * sin_beta_star * cos_alpha],
            [0, 0, 1.0 / c],
        ]
    )
    return B


def uv_to_ub(
    u: np.ndarray, v: np.ndarray, lattice_params: tuple[float, float, float, float, float, float]
) -> np.ndarray:
    """Compute UB matrix from u,v and lattice parameters. Eq.76-81."""
    B = b_mat(lattice_params)
    t1_c = B @ u
    t3_c = np.cross(B @ u, B @ v)
    t2_c = np.cross(t3_c, t1_c)
    T_c = np.array(
        [
            t1_c / np.linalg.norm(t1_c),
            t2_c / np.linalg.norm(t2_c),
            t3_c / np.linalg.norm(t3_c),
        ]
    ).T
    T_epsilon = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]).T
    u_mat = T_epsilon @ T_c.T
    return u_mat @ B

#-----------Calculate UB Matrix----------------------------
