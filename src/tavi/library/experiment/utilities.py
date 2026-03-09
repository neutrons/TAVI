"""General utilities for tas related functions and classes."""

import numpy as np
from scipy.constants import e, hbar, m_n


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
