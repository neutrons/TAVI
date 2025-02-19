from typing import Optional

import numpy as np

from tavi.lattice_algorithm import b_mat_from_lattice, lattice_params_from_g_star_mat
from tavi.utilities import Peak, UBConf, en2q, get_angle_from_triangle

# -----------------------------------------------------
# Angle and Q math
# -----------------------------------------------------


def q_lab(ei: float, ef: float, theta: float, phi: float = 0) -> np.ndarray:
    """
    Reutrn momentum transfer vector q in lab frame

    Args:
        ei: incident energy, in meV
        ef: final energy, in meV
        theta: In degree. Same as two theta for TAS with a single detector
                in the scattering plane.
        phi: In degree. Always zero for TAS with a single detector in the scattering plane.

    Return:
        np.ndarry of size 3

    Notes:
        Using convention in Mantid/Internaltional Crystallography Table

    """

    ki = en2q(ei)
    kf = en2q(ef)
    theta_rad = np.deg2rad(theta)
    phi_rad = np.deg2rad(phi)
    q_lab = np.array(
        [
            -kf * np.sin(theta_rad) * np.cos(phi_rad),
            -kf * np.sin(theta_rad) * np.sin(phi_rad),
            ki - kf * np.cos(theta_rad),
        ],
    )
    return q_lab


def q_norm_from_hkl(hkl: tuple[float, float, float], b_mat: np.ndarray):
    """return norm of q for given (h,k,l)

    Note:
        Either b_mat or ub_mat would work, since U^T.U=U^-1.U=1"""

    b_hkl = np.matmul(b_mat, np.array(hkl))
    q_squared = np.matmul(b_hkl.T, b_hkl)
    q_norm = np.sqrt(q_squared) * 2 * np.pi
    return q_norm


def two_theta_from_hkle(
    hkl: tuple[float, float, float],
    ei: float,
    ef: float,
    b_mat: np.ndarray,
) -> Optional[float]:
    """Return None is (h,k,l) can't be reached.

    Note:
        Either b_mat or ub_mat would work, since U^T.U=U^-1.U=1"""

    ki = en2q(ei)
    kf = en2q(ef)
    q_norm = q_norm_from_hkl(hkl, b_mat)
    two_theta_radian = get_angle_from_triangle(ki, kf, q_norm)
    return two_theta_radian


# -----------------------------------------------------
# R matrix math
# -----------------------------------------------------


def r_matrix_with_minimal_tilt(
    hkl: tuple[float, float, float],
    ei: float,
    ef: float,
    two_theta: float,
    ub_conf: UBConf,
) -> np.ndarray:
    """Calculate R matrix when the tilt from the scattering plane is minimal
    Args:
        hkl: tuple of Miller indices
        ei: incident energy, in meV
        ef: finial energy, in meV
        two_theta: two_theta angle, in degrees
        ub_conf: need to have ub_mat, plane_normal and in_plane_ref
    """

    ub_mat = ub_conf.ub_mat
    plane_normal = ub_conf.plane_normal
    in_plane_ref = ub_conf.in_plane_ref

    ki = en2q(ei)
    kf = en2q(ef)

    q_norm = q_norm_from_hkl(hkl, ub_mat)
    tt = np.deg2rad(two_theta)
    q_lab1 = np.array([-kf * np.sin(tt), 0, ki - kf * np.cos(tt)]) / q_norm
    q_lab2 = np.array([ki - kf * np.cos(tt), 0, kf * np.sin(tt)]) / q_norm
    q_lab3 = np.array([0, 1, 0])
    q_lab_mat = np.array([q_lab1, q_lab2, q_lab3]).T

    ZERO = 1e-6
    t1 = ub_mat @ np.array(hkl)
    if np.abs(np.dot(t1, plane_normal)) < ZERO:
        # t1 in plane
        t3 = plane_normal
        t2 = np.cross(t3, t1)
    elif np.linalg.norm(np.cross(plane_normal, t1)) < ZERO:
        # oops, t1 along plane_normal
        t2 = in_plane_ref
        t3 = np.cross(t1, t2)
    else:
        # t1 not in plane, need to change tilts
        # t2p = np.cross(plane_normal, t1)
        # t3 = np.cross(t1, t2p)
        # t2 = np.cross(t3, t1)
        t2 = np.cross(plane_normal, t1)
        t3 = np.cross(t1, t2)

    t_mat = np.array(
        [
            t1 / np.linalg.norm(t1),
            t2 / np.linalg.norm(t2),
            t3 / np.linalg.norm(t3),
        ]
    ).T
    t_mat_inv = np.linalg.inv(t_mat)
    r_mat = q_lab_mat @ t_mat_inv
    return r_mat


# -----------------------------------------------------
# UB matrix convsersion to u and v vectors
# -----------------------------------------------------
def ub_matrix_to_uv(ub_matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Calculate u and v vector from UB matrix

    Note:
        u vector, in reciprocal lattice unit, along beam
        v vector, in reciprocal lattice unit,in the horizaontal scattering plane
    """

    inv_ub_matrix = np.linalg.inv(ub_matrix)
    u = np.matmul(inv_ub_matrix, np.array([0, 0, 1]))
    v = np.matmul(inv_ub_matrix, np.array([1, 0, 0]))
    return (u, v)


def ub_matrix_to_lattice_params(
    ub_matrix: np.ndarray,
) -> tuple[float, float, float, float, float, float]:
    """Calculate lattice parameters from UB matrix"""

    g_star_mat = np.matmul(np.transpose(ub_matrix), ub_matrix)
    return lattice_params_from_g_star_mat(g_star_mat)


def uv_to_ub_matrix(
    u: np.ndarray,
    v: np.ndarray,
    lattice_params: tuple[float, float, float, float, float, float],
) -> np.ndarray:
    """Calculate UB matrix from u and v vector, and lattice parameters"""

    b_mat = b_mat_from_lattice(lattice_params)
    t1 = np.matmul(b_mat, u)
    t2_prime = np.matmul(b_mat, v)
    t3 = np.cross(t1, t2_prime)
    t2 = np.cross(t3, t1)
    t_mat = np.array(
        [
            t1 / np.linalg.norm(t1),
            t2 / np.linalg.norm(t2),
            t3 / np.linalg.norm(t3),
        ]
    ).T
    q_mat = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]).T
    ub_matrix = np.matmul(
        q_mat,
        np.matmul(np.linalg.inv(t_mat), b_mat),
    )

    return ub_matrix


# -----------------------------------------------------
# UB matrix math
# -----------------------------------------------------


def find_u_from_two_peaks(
    peaks: tuple[Peak, Peak],
    b_mat: np.ndarray,
    r_mat_inv,
    ei: float,
    ef: float,
):
    """Calculate U matrix from two peaks, need to know B matrix"""

    peak1, peak2 = peaks
    q_hkl1 = np.dot(b_mat, np.array(peak1.hkl))
    q_hkl2p = np.dot(b_mat, np.array(peak2.hkl))
    q_hkl3 = np.cross(q_hkl1, q_hkl2p)
    q_hkl_2 = np.cross(q_hkl3, q_hkl1)

    q_hkl_mat = np.array(
        [
            q_hkl1 / np.linalg.norm(q_hkl1),
            q_hkl_2 / np.linalg.norm(q_hkl_2),
            q_hkl3 / np.linalg.norm(q_hkl3),
        ]
    ).T

    q_lab1 = q_lab(ei=ei, ef=ef, theta=peak1.angles.two_theta)
    q_lab2 = q_lab(ei=ei, ef=ef, theta=peak2.angles.two_theta)

    # Goniometer angles all zeros in q_sample frame
    q_sample1 = np.matmul(r_mat_inv(peak1.angles), q_lab1)
    q_sample2p = np.matmul(r_mat_inv(peak2.angles), q_lab2)
    q_sample3 = np.cross(q_sample1, q_sample2p)
    q_sample2 = np.cross(q_sample3, q_sample1)

    q_sample1 = q_sample1 / np.linalg.norm(q_sample1)
    q_sample2 = q_sample2 / np.linalg.norm(q_sample2)
    q_sample3 = q_sample3 / np.linalg.norm(q_sample3)

    q_sample_mat = np.array([q_sample1, q_sample2, q_sample3]).T

    u_mat = np.matmul(q_sample_mat, q_hkl_mat.T)

    # plane normal always up along +Y
    plane_normal = -q_sample3 if q_sample3[1] < 0 else q_sample3
    in_plane_ref = q_sample1
    return u_mat, plane_normal, in_plane_ref


# TODO
def find_ub_from_three_peaks(
    peaks: tuple[Peak, Peak, Peak],
    r_mat_inv,
    ei: float,
    ef: float,
):
    """Find UB matrix from three observed peaks for a given goniomete"""
    u_mat = None
    b_mat = None
    ub_mat = None
    plane_normal = None
    in_plane_ref = None
    return (u_mat, b_mat, ub_mat, plane_normal, in_plane_ref)


# TODO
def find_ub_from_multiple_peaks(
    peaks: tuple[Peak, ...],
    r_mat_inv,
    ei: float,
    ef: float,
):
    """Find UB matrix from more than three observed peaks for a given goniomete"""
    u_mat = None
    b_mat = None
    ub_mat = None
    plane_normal = None
    in_plane_ref = None
    return (u_mat, b_mat, ub_mat, plane_normal, in_plane_ref)
