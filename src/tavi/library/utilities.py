"""General utilities for tas related functions and classes."""

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.constants import e, hbar, m_n


def SE2K(E: float) -> float:
    """Convert energy E to momentum transfer k."""
    E2K = np.sqrt(2e-3 * e * m_n) * 1e-10 / hbar
    return np.sqrt(E) * E2K


@dataclass(frozen=True)
class Peak:
    """
    Description of a peak and motor angles.

    Phyical/virtual monitor positions

    hkl: miller indice (h,k,l)
    angles:
    """

    hkl: tuple[float, float, float]
    angles: Optional[MotorAngles] = None


@dataclass(frozen=True)
class MotorAngles:
    """
    Moter anlges.

    two_theta: s2 angle, in degree
    omega: s1 angle, in degree
    sgl: sample goniometer lower, in degree
    sgu: sample goniometer upper, in degree
    chi: chi angle for a four-circle goniometer, in degree
    phi: phi angle for a four-circle goniometer, in degree

    Note:
        use angles = (two_theta, omega, sgl, sgu) for a Huber table,
        angles = (two_theta, omega, chi, phi) for a four-circle in the bisect mode

    """

    two_theta: float
    omega: Optional[float] = None
    sgl: Optional[float] = None
    sgu: Optional[float] = None
    chi: Optional[float] = None
    phi: Optional[float] = None


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
