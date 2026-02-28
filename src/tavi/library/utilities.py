"""General utilities for tas related functions and classes."""

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.constants import hbar, m_n


def SE2K(e: float) -> float:
    """Convert energy E to momentum transfer k."""
    E2K = np.sqrt(2e-3 * e * m_n) * 1e-10 / hbar
    return np.sqrt(e) * E2K


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
