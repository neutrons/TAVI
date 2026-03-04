"""General utilities for tas related functions and classes."""

import numpy as np
from scipy.constants import e, hbar, m_n

from tavi.library.component.goniometer import Goniometer
from tavi.library.experiment.peak import Peak


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


def q_norm_from_hkl(hkl: tuple[float, float, float], b_mat: np.ndarray) -> np.ndarray:
    """
    Return norm of q for given (h,k,l).

    Note:
        Either b_mat or ub_mat would work, since U^T.U=U^-1.U=1

    """
    q_norm = 2 * np.pi * np.linalg.norm(b_mat @ np.array(hkl))

    return q_norm


# -----------R matrix-----------
def r_matrix_with_minimal_tilt(
    peak: Peak,
    ei: float,
    ef: float,
    sense: str,
    ub_mat: np.ndarray,
    plane_normal: np.ndarray,
    in_plane_ref: np.ndarray,
) -> np.ndarray:
    """
    Calculate R matrix when the tilt from the scattering plane is minimal.

    Args:
        peak: Peak
        ei: incident energy, in meV
        ef: final energy, in meV
        sense: either + or -, indicate if phi =0 or pi
        ub_mat: UB matrix
        plane_normal: plane_normal
        in_plane_ref: in plane reflection.

    """
    ki = SE2K(ei)
    kf = SE2K(ef)
    q_norm = q_norm_from_hkl(peak.hkl, ub_mat)

    # Eq.112
    Q_squared = 4 * np.pi**2 * (np.array(peak.hkl).T @ ub_mat.T @ ub_mat @ peak.hkl)
    # Eq.113
    analyzer_position = 1 if sense == "+" else -1
    two_theta = analyzer_position * np.arccos((ki**2 + kf**2 - Q_squared) / (2 * ki * kf))
    # with minimal tilt, we are considering scenario described above Eq.114
    q_lab1 = np.array([-kf * np.sin(two_theta), 0, ki - kf * np.cos(two_theta)]) / q_norm
    q_lab2 = np.array([ki - kf * np.cos(two_theta), 0, kf * np.sin(two_theta)]) / q_norm
    q_lab3 = np.array([0, 1, 0])

    Q_lab = np.array([q_lab1, q_lab2, q_lab3]).T
    tol = 1e-5

    # Eq.106
    t1 = ub_mat @ np.array(peak.hkl)
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


# -----------Calculate UB Matrix-----------


def find_u_from_two_peaks(
    peaks: tuple[Peak, Peak], b_mat: np.ndarray, r_mat: Goniometer.r_mat, ei: float, ef: float
) -> np.ndarray:
    """
    Calculate U matrix from two peaks.

    r_mat can be removed later when goniometer is implemented as it can be calculated from
    peak.angles.

    Follow Eq.76-81 and Eq.83-88. We assume q_3 is perpendicular from the two peaks
    """
    peak1, peak2 = peaks
    B = b_mat

    t1_c = B @ peak1.hkl
    t3_c = np.cross(B @ peak1.hkl, B @ peak2.hkl)
    t2_c = np.cross(t3_c, t1_c)

    T_c = np.array(
        [
            t1_c / np.linalg.norm(t1_c),
            t2_c / np.linalg.norm(t2_c),
            t3_c / np.linalg.norm(t3_c),
        ]
    ).T

    # We need to create vectors t1_v, t2_v, t3_v
    q_lab_1 = q_lab(ei, ef, peak1.angles.two_theta)
    q_lab_2 = q_lab(ei, ef, peak2.angles.two_theta)

    # In identical fashion as described above Eq.79
    t1_v = np.linalg.inv(r_mat(peak1.angles)) @ q_lab_1
    t3_v = np.cross(np.linalg.inv(r_mat(peak1.angles)) @ q_lab_1, np.linalg.inv(r_mat(peak2.angles)) @ q_lab_2)
    t2_v = np.cross(t3_v, t1_v)
    T_v = np.array(
        [
            t1_v / np.linalg.norm(t1_v),
            t2_v / np.linalg.norm(t2_v),
            t3_v / np.linalg.norm(t3_v),
        ]
    ).T

    u_mat = T_v @ T_c.T
    return u_mat


def find_ub_from_three_peaks(peaks: tuple[Peak, Peak, Peak], r_mat: np.ndarray, ei: float, ef: float) -> np.ndarray:
    """
    Calculate U matrix from three peaks.

    r_mat can be removed later when goniometer is implemented as it can be calculated from
    peak.angles.

    Follow Eq.83-90.
    """
    peak1, peak2, peak3 = peaks

    # we directly use the three peaks as t1_c, t2_c and t3_c
    V = np.array([peak1.hkl, peak2.hkl, peak3.hkl]).T

    q_lab_1 = q_lab(ei, ef, peak1.angles.two_theta)
    q_lab_2 = q_lab(ei, ef, peak2.angles.two_theta)
    q_lab_3 = q_lab(ei, ef, peak3.angles.two_theta)

    q1_v = np.linalg.inv(r_mat(peak1.angles)) @ q_lab_1 / (2 * np.pi)
    q2_v = np.linalg.inv(r_mat(peak2.angles)) @ q_lab_2 / (2 * np.pi)
    q3_v = np.linalg.inv(r_mat(peak3.angles)) @ q_lab_3 / (2 * np.pi)
    Q_v = np.array([q1_v, q2_v, q3_v]).T

    ub_mat = Q_v @ np.linalg.inv(V)
    return ub_mat


def find_ub_from_multiple_peaks(peaks: tuple[Peak, ...], r_mat: np.ndarray, ei: float, ef: float) -> np.ndarray:
    """
    Calculate U matrix from three peaks.

    r_mat can be removed later when goniometer is implemented as it can be calculated from
    peak.angles.

    Follow Eq.89-98.
    """
    n = len(peaks)
    Q_v = np.zeros((3, 3))
    VV = np.zeros((3, 3))

    for i in range(n):
        hkl = peaks[i].hkl
        q_lab_i = q_lab(ei, ef, peaks[i].angles.two_theta)
        q_v_i = np.linalg.inv(r_mat(peaks[i].angles)) @ q_lab_i / (2 * np.pi)
        for j in range(3):
            for k in range(3):
                Q_v[j, k] += q_v_i[k] * hkl[j]
                VV[j, k] += hkl[k] * hkl[j]
    ub_mat = Q_v.T @ np.linalg.inv(VV).T
    return ub_mat


def plane_normal_from_two_peaks(
    u_mat: np.ndarray, b_mat: np.ndarray, peaks: tuple[Peak, Peak]
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate plane_normal and in_plane reflection.

    Both are vectors representing peaks in Q_lab.
    """
    peak1, peak2 = peaks
    t1_c = b_mat @ peak1.hkl
    t3_c = np.cross(b_mat @ peak1.hkl, b_mat @ peak2.hkl)

    # Eq. 79(3) calculate t3_v from U and t3_c
    plane_normal = u_mat @ t3_c / np.linalg.norm(t3_c)
    # if y is pointing down, set it to point up
    plane_normal = -plane_normal if plane_normal[1] < 0 else plane_normal
    in_plane_ref = u_mat @ t1_c / np.linalg.norm(t1_c)
    return (plane_normal, in_plane_ref)


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
