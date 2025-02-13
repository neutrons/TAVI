from typing import Optional

import numpy as np

from tavi.lattice_algorithm import b_mat_from_lattice, lattice_params_from_g_star_mat
from tavi.utilities import Peak, en2q


def q_lab(
    theta: float,
    phi: float = 0,
    ei: Optional[float] = None,
    ef: Optional[float] = None,
) -> np.ndarray:
    """
    Reutrn momentum transfer vector q in lab frame

    Args:
        theta: In degree. Same as two theta for TAS with a single detector
                in the scattering plane.
        phi: In degree. Always zero for TAS with a single detector in the scattering plane.

    Return:
        np.ndarry of size 3

    Notes:
        Using convention in Mantid/Internaltional Crystallography Table

    """

    if (ei is None) and (ef is not None):
        ei = ef
    elif (ei is not None) and (ef is None):
        ef = ei
    else:
        raise ValueError("Ei and Ef needs to be provided to calculate q_lab.")

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


# -----------------------------------------------------
# U matrix math
# -----------------------------------------------------


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
    u = inv_ub_matrix @ np.array([0, 0, 1])
    v = inv_ub_matrix @ np.array([1, 0, 0])
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


def find_u_from_two_peaks(peaks: tuple[Peak, Peak], b_mat, r_mat_inv):
    """Calculate U matrix from two peaks, need to know B matrix"""

    peak1, peak2 = peaks
    q_hkl1 = np.matmul(b_mat, np.array(peak1.hkl))
    q_hkl2p = np.matmul(b_mat, np.array(peak2.hkl))
    q_hkl3 = np.cross(q_hkl1, q_hkl2p)
    q_hkl_2 = np.cross(q_hkl3, q_hkl1)

    q_hkl_mat = np.array(
        [
            q_hkl1 / np.linalg.norm(q_hkl1),
            q_hkl_2 / np.linalg.norm(q_hkl_2),
            q_hkl3 / np.linalg.norm(q_hkl3),
        ]
    ).T

    q_lab1 = q_lab(
        peak1.angles.two_theta,
        ei=peak1.ei,
        ef=peak1.ef if peak1.ef is not None else peak1.ei,
    )
    q_lab2 = q_lab(
        peak2.angles.two_theta,
        ei=peak2.ei,
        ef=peak2.ef if peak2.ef is not None else peak2.ei,
    )

    # Goniometer angles all zeros in q_sample frame
    q_sample1 = np.matmul(r_mat_inv(peak1.angles), q_lab1)
    q_sample2p = np.matmul(r_mat_inv(peak2.angles), q_lab2)
    q_sample3 = np.cross(q_sample1, q_sample2p)
    q_sample2 = np.cross(q_sample3, q_sample1)

    q_sample1 = q_sample1 / np.linalg.norm(q_sample1)
    q_sample2 = q_sample2 / np.linalg.norm(q_sample2)
    q_sample3 = q_sample3 / np.linalg.norm(q_sample3)

    q_sample_mat = np.array([q_sample1, q_sample2, q_sample3]).T

    u_mat = np.matmul(q_sample_mat, np.linalg.inv(q_hkl_mat))

    plane_normal = q_sample3
    if plane_normal[1] < 0:  # plane normal always up along +Y
        plane_normal = -plane_normal

    in_plane_ref = q_sample1

    return u_mat, plane_normal, in_plane_ref


# # TODO
# def _find_ub_from_three_peaks(
#     self,
#     peaks: tuple[Peak, Peak, Peak],
# ) -> UBConf:
#     """Find UB matrix from three observed peaks for a given goniomete"""
#     ubconf = UBConf()
#     return ubconf


# # TODO
# def _find_ub_from_multiple_peaks(
#     self,
#     peaks: tuple[Peak, ...],
# ) -> UBConf:
#     """Find UB matrix from more than three observed peaks for a given goniomete"""
#     ubconf = UBConf()
#     return ubconf


# @staticmethod
# def norm_mat(t1, t2, t3):
#     mat = np.array([t1 / np.linalg.norm(t1), t2 / np.linalg.norm(t2), t3 / np.linalg.norm(t3)]).T
#     return mat


# def _t_mat_minimal_tilt(self, hkl: np.ndarray):
#     """Build matrix T assuming minimal goniometer tilt angles"""

#     if not isinstance(self.sample, Xtal):
#         raise ValueError("Sample needs to be Xtal class for UB calculation.")
#     if self.sample.ub_mat is None:
#         raise ValueError("UB matrix is unknown.")
#     if self.sample.plane_normal is None:
#         raise ValueError("Plane normal vector is not known.")
#     if self.sample.in_plane_ref is None:
#         raise ValueError("In-plnae reference vector is not known.")

#     EPS = 1e-8  # zero

#     plane_normal = np.array(self.sample.plane_normal)
#     in_plane_ref = np.array(self.sample.in_plane_ref)
#     ub_mat = self.sample.ub_mat

#     if self.SPICE_CONVENTION:  # suffle the order following SPICE convention
#         plane_normal = spice_to_mantid(plane_normal)
#         in_plane_ref = spice_to_mantid(in_plane_ref)
#         ub_mat = spice_to_mantid(ub_mat)

#     q = ub_mat @ hkl
#     t1 = q / np.linalg.norm(q)

#     if np.dot(t1, plane_normal) < EPS:  # t1 in plane
#         t3 = plane_normal
#         t2 = np.cross(t3, t1)
#         return TAS.norm_mat(t1, t2, t3)

#     # t1 not in plane, need to change tilts
#     if np.linalg.norm(np.cross(plane_normal, t1)) < EPS:
#         # oops, t1 along plane_normal
#         t2 = in_plane_ref
#         t3 = np.cross(t1, t2)
#         return TAS.norm_mat(t1, t2, t3)

#     else:
#         t2p = np.cross(plane_normal, t1)
#         t3 = np.cross(t1, t2p)
#         t2 = np.cross(t3, t1)
#         return TAS.norm_mat(t1, t2, t3)
