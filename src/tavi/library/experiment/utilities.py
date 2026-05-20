"""General utilities for tas related functions and classes."""

import numpy as np
from scipy.constants import e, hbar, m_n

# import scipy.constants as co
# ksq2E = (co.Planck / co.elementary_charge / 2.0 / np.pi) ** 2.0 * co.elementary_charge / 2.0 / co.neutron_mass * 1e23
# calculated with scipy.constants using the formula above, it converts SI unit to angstrom, meV unit
ksq2eng = 2.072124855

# convert sigma to FWHM, FWHM = 2*sqrt(2*ln(2))*sigma
sig2fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0))


def SE2K(E: float) -> float:
    """Convert energy E to momentum transfer k."""
    E2K = np.sqrt(2e-3 * e * m_n) * 1e-10 / hbar
    return np.sqrt(E) * E2K


def q_lab(ei: float, ef: float, theta: float, phi: float = 0) -> np.ndarray:
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
    return np.array(
        [
            -kf * np.sin(np.radians(theta)) * np.cos(np.radians(phi)),
            -kf * np.sin(np.radians(theta)) * np.sin(np.radians(phi)),
            ki - kf * np.cos(np.radians(theta)),
        ]
    )


def get_angle_from_triangle(a: float, b: float, c: float) -> float:
    """
    In a triangle with sides a,b and c, get angle between a and b in radian.

    Note:
        return value in [0,pi]

    """
    zero = 1e-6
    if (np.abs(a) < zero) or (np.abs(b) < zero):
        raise ValueError("Triangle cannot be closed.")
    arc_cos = (a**2 + b**2 - c**2) / (2 * a * b)
    if arc_cos > 1 or arc_cos < -1:
        raise ValueError("Triangle cannot be closed.")
    return np.arccos(arc_cos)


def rot_x(theta: float) -> np.ndarray:
    """Rotation matrix for rotation about x axis by theta degree."""
    t = np.radians(theta)
    return np.array([[1, 0, 0], [0, np.cos(t), -np.sin(t)], [0, np.sin(t), np.cos(t)]])


def rot_y(theta: float) -> np.ndarray:
    """Rotation matrix for rotation about y axis by theta degree."""
    t = np.radians(theta)
    return np.array([[np.cos(t), 0, np.sin(t)], [0, 1, 0], [-np.sin(t), 0, np.cos(t)]])


def rot_z(theta: float) -> np.ndarray:
    """Rotation matrix for rotation about z axis by theta degree."""
    t = np.radians(theta)
    return np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])


def rotation_matrix_2d(phi: float) -> np.ndarray:
    """
    Rotate the coordination system by angle of phi about z-axis instead of a vector.

    Notice the negative sign from rot_z
    """
    return np.array([[np.cos(phi), np.sin(phi), 0], [-np.sin(phi), np.cos(phi), 0], [0, 0, 1]])


def quadric_proj(quadric: np.ndarray, idx: int) -> np.ndarray:
    """
    Project out a specific dimension.

    Implementation is equivalent of R.A.Robinson et al. Analytical Techniques for Instrument
    Design - Matrix Methods Eq. 13 (demonimeter should not be squared). Or Eck[14] 57. Calculation follows Takin's ellipse.h example.

    We can verify with:
        k = 0
        test = np.zeros((4,4))
        for i in range(4):
            for j in range(4):
                if i == k or j == k:
                    continue
                test[i][j] = quadric[i][j] - quadric[i][k]*quadric[j][k]/quadric[k][k].

    """
    atol = 1e-8

    if np.abs(quadric[idx, idx]) < atol:
        return np.delete(np.delete(quadric, idx, axis=0), idx, axis=1)

    # row/column along which to perform the orthogonal projection
    vec = 0.5 * (quadric[idx, :] + quadric[:, idx])  # symmetrise if not symmetric
    vec /= np.sqrt(quadric[idx, idx])  # normalise to indexed component
    proj_op = np.outer(vec, vec)  # projection operator
    ortho_proj = quadric - proj_op  # projected quadric
    return np.delete(np.delete(ortho_proj, idx, axis=0), idx, axis=1)
