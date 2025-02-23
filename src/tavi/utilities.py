# -*- coding: utf-8 -*-
from typing import NamedTuple, Optional

import numpy as np

# --------------------------------------------------------------------------
# constants
# --------------------------------------------------------------------------

# import scipy.constants as co
# ksq2E = (co.Planck / co.elementary_charge / 2.0 / np.pi) ** 2.0 * co.elementary_charge / 2.0 / co.neutron_mass * 1e23
ksq2eng = 2.072124855  # calculated with scipy.constants using the formula above

sig2fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0))
cm2angstrom = 1e8
min2rad = 1.0 / 60.0 / 180.0 * np.pi
rad2deg = 180.0 / np.pi

# --------------------------------------------------------------------------
# Named tuples
# --------------------------------------------------------------------------


class MotorAngles(NamedTuple):
    """Moter anlges

    two_theta: s2 angle, in degree
    omega: s1 angle, in degree
    sgl: sample goniometer lower, in degree
    sgu: sample goniometer upper, in degree
    chi: chi angle for a four-circle goniometer, in degree
    phi: phi angle for a four-circle goniometer, in degree

    Note:
        use angles = (two_theta, omega, sgl, sgu) for a Huber table,
        angles = (two_theta, omega, chi, phi) for a four-circle in the bisect mode"""

    two_theta: float
    omega: Optional[float] = None
    sgl: Optional[float] = None
    sgu: Optional[float] = None
    chi: Optional[float] = None
    phi: Optional[float] = None


class Peak(NamedTuple):
    """
    Phsical/virtual monitor positions

    hkl: miller indice (h,k,l)
    angles: moter angles
    """

    hkl: tuple[float, float, float]
    angles: MotorAngles


class UBConf(NamedTuple):
    """Logs for UB matrix determination


    ub_peaks (tuple of Peaks): peaks used to determine the UB matrix
    u_mat (np.adarray): U matrix
    b_mat (np.adarray): B matrix
    ub_matrix (np.adarray): UB matrix
    plane_normal (np.adarray): normal vector in Qsample frame, goniometers at zero
    in_plane_ref (np.adarray): in plane vector in Qsample frame, goniometers at zero
    """

    ub_mat: np.ndarray
    plane_normal: Optional[np.ndarray] = None
    in_plane_ref: Optional[np.ndarray] = None
    u_mat: Optional[np.ndarray] = None
    b_mat: Optional[np.ndarray] = None
    ub_peaks: Optional[tuple[Peak]] = None


# --------------------------------------------------------------------------
# helper functions
# --------------------------------------------------------------------------
def en2q(en: float) -> float:
    """convert energy en in meV to momontum transfer q in inverse Angstrom"""
    if en < 0:
        raise ValueError(f"Cannot convert negative energy en={en} to momentum transfer q.")
    q = np.sqrt(en / ksq2eng)
    return q


def q2en(q: float) -> float:
    """convert momontum transfer q in inverse Angstrom to energy en in mev"""
    if q < 0:
        raise ValueError(f"Converting negative momentum transfer q={q} to energy.")
    en = ksq2eng * q**2
    return en


def get_angle_from_triangle(a: float, b: float, c: float) -> Optional[float]:
    """In a triangle with sides a,b and c, get angle between a and b in radian
    Note:
        return value in [0,pi]"""
    acos = (a**2 + b**2 - c**2) / (2 * a * b)
    if acos > 1 or acos < -1:
        angle = None
    else:
        angle = np.arccos(acos)
    return angle


def get_angle_vec(v1, v2):
    """Get the angle in degress between two vectors v1 and v2"""
    return np.arccos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)) / np.pi * 180


def get_angle_bragg(
    neutron_momentum: float,
    sample_d_spaceing: float,
):
    """return angle based on Bragg's law, in radian
    2d sin(theta) = lambda = 2 pi /q
    """
    theta = np.arcsin(np.pi / (sample_d_spaceing * neutron_momentum))
    return theta


def rotation_matrix_2d(phi):
    """rotate the coordination system by angle of phi about z-axis

    Args:
        phi (float): angle in radian

    Note:
        This is to rotate the coordination system, NOT the vector!


    """
    s = np.sin(phi)
    c = np.cos(phi)
    mat = np.array(
        [
            [c, s, 0],
            [-s, c, 0],
            [0, 0, 1],
        ]
    )
    return mat


def spice_to_mantid(vec):
    """suffle the order from spice convention to mantid convention"""
    return np.array([vec[0], vec[2], -vec[1]])


def mantid_to_spice(vec):
    """suffle the order from mantid convention to spice convention"""
    return np.array([vec[0], -vec[2], vec[1]])


def rot_x(nu):
    """rotation matrix about y-axis by angle nu

    Args:
        nu (float): angle in degrees

    Note:
        Using Mantid convention, beam along z, y is up, x in plane
    """

    angle = np.deg2rad(nu)
    c = np.cos(angle)
    s = np.sin(angle)
    mat = np.array(
        [
            [1, 0, 0],
            [0, c, -s],
            [0, s, c],
        ]
    )
    return mat


def rot_y(omega):
    """rotation matrix about y-axis by angle omega

    Args:
        omega (float): angle in degrees

    Note:
        Using Mantid convention, beam along z, y is up, x in plane
    """

    angle = np.deg2rad(omega)
    c = np.cos(angle)
    s = np.sin(angle)
    mat = np.array(
        [
            [c, 0, s],
            [0, 1, 0],
            [-s, 0, c],
        ]
    )
    return mat


def rot_z(mu):
    """rotation matrix about z-axis by angle mu

    Args:
        mu (float): angle in degrees

    Note:
        Using Mantid convention, beam along z, y is up, x in plane
    """

    angle = np.deg2rad(mu)
    c = np.cos(angle)
    s = np.sin(angle)
    mat = np.array(
        [
            [c, -s, 0],
            [s, c, 0],
            [0, 0, 1],
        ]
    )
    return mat
