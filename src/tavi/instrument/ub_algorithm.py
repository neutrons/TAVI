from typing import Optional

import numpy as np

from tavi.utilities import Peak, UBConf, en2q, spice_to_mantid


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
# real and reciprocal lattice math
# -----------------------------------------------------


def v_alpha_beta_gamma_calc(alpha: float, beta: float, gamma: float) -> float:
    """
    Calculate V_alpha_bet_gamma = Volume/(abc) where Volume = a * (b x c) is the volume of real space cell
    """
    cos_alpha = np.cos(np.deg2rad(alpha))
    cos_beta = np.cos(np.deg2rad(beta))
    cos_gamma = np.cos(np.deg2rad(gamma))
    v_alpha_beta_gamma = np.sqrt(1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 2 * cos_alpha * cos_beta * cos_gamma)
    return v_alpha_beta_gamma


def real_space_vectors(
    lattice_params: tuple[float, float, float, float, float, float]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the real space lattice vectors in Cartesian coordiantes

    Args:
        lattice_params: tuple of (a, b, c, alpha, beta, gamma) in Angstrom and degrees

    Return:
        tuple of vectors

    Note:
        a-vec is always along x-axis
    """
    a, b, c, alpha, beta, gamma = lattice_params
    cos_alpha = np.cos(np.deg2rad(alpha))
    cos_beta = np.cos(np.deg2rad(beta))
    cos_gamma = np.cos(np.deg2rad(gamma))
    sin_gamma = np.sin(np.deg2rad(gamma))

    ac = np.array([a, 0, 0])
    bc = np.array([b * cos_gamma, b * sin_gamma, 0])
    v_abg = v_alpha_beta_gamma_calc(alpha, beta, gamma)
    cc = np.array(
        [
            c * cos_beta,
            c * (cos_alpha - cos_gamma * cos_beta) / sin_gamma,
            c * v_abg / sin_gamma,
        ]
    )
    return (ac, bc, cc)


def reciprocal_latt_params(
    lattice_params: tuple[float, float, float, float, float, float]
) -> tuple[float, float, float, float, float, float]:
    """Calculate the reciprocal lattice parameter lengths and angles in inverse Angstrom and degrees"""

    a, b, c, alpha, beta, gamma = lattice_params
    alpha_rad = np.deg2rad(alpha)
    beta_rad = np.deg2rad(beta)
    gamma_rad = np.deg2rad(gamma)

    sin_alpha = np.sin(alpha_rad)
    cos_alpha = np.cos(alpha_rad)
    sin_beta = np.sin(beta_rad)
    cos_beta = np.cos(beta_rad)
    cos_gamma = np.cos(gamma_rad)
    sin_gamma = np.sin(gamma_rad)

    v_abg = v_alpha_beta_gamma_calc(alpha, beta, gamma)

    a_star = sin_alpha / a / v_abg * np.pi * 2
    b_star = sin_beta / b / v_abg * np.pi * 2
    c_star = sin_gamma / c / v_abg * np.pi * 2
    alpha_star = np.arccos((cos_beta * cos_gamma - cos_alpha) / sin_beta / sin_gamma)
    beta_star = np.arccos((cos_gamma * cos_alpha - cos_beta) / sin_alpha / sin_gamma)
    gamma_star = np.arccos((cos_alpha * cos_beta - cos_gamma) / sin_beta / sin_alpha)
    alpha_star = np.rad2deg(alpha_star)
    beta_star = np.rad2deg(beta_star)
    gamma_star = np.rad2deg(gamma_star)

    return (a_star, b_star, c_star, alpha_star, beta_star, gamma_star)


def reciprocal_space_vectors(
    lattice_params: tuple[float, float, float, float, float, float]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the reciprocal space lattice vectors in the Cartesian coordinates
    """
    a, b, c, alpha, beta, gamma = lattice_params
    v_abg = v_alpha_beta_gamma_calc(alpha, beta, gamma)
    v = v_abg * a * b * c
    prefactor = 2 * np.pi / v
    a_vec, b_vec, c_vec = real_space_vectors(lattice_params)
    a_star_vec = np.cross(b_vec, c_vec) * prefactor
    b_star_vec = np.cross(c_vec, a_vec) * prefactor
    c_star_vec = np.cross(a_vec, b_vec) * prefactor

    return (a_star_vec, b_star_vec, c_star_vec)


def reciprocal_basis(
    lattice_params: tuple[float, float, float, float, float, float]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate the reciprocal basis vectors i_star, j_star, k_star"""

    (a_star_vec, b_star_vec, c_star_vec) = reciprocal_space_vectors(lattice_params)
    i_star = a_star_vec / np.linalg.norm(a_star_vec)
    a_star_perp = np.cross(
        np.cross(a_star_vec, b_star_vec),
        a_star_vec,
    )
    j_star = a_star_perp / np.linalg.norm(a_star_perp)
    k_star = np.cross(i_star, j_star)
    return (i_star, j_star, k_star)


# -----------------------------------------------------
# B matrix math
# -----------------------------------------------------
def b_mat_from_lattice(lattice_params) -> np.ndarray:
    """Calculate the B matrix from lattice parameters"""

    _, _, c, alpha, _, _ = lattice_params
    (a_star, b_star, c_star, _, beta_star, gamma_star) = reciprocal_latt_params(lattice_params)

    alpha_rad = np.deg2rad(alpha)
    beta_star_rad = np.deg2rad(beta_star)
    gamma_star_rad = np.deg2rad(gamma_star)

    b_mat = np.array(
        [
            [a_star, b_star * np.cos(gamma_star_rad), c_star * np.cos(beta_star_rad)],
            [0, b_star * np.sin(gamma_star_rad), -c_star * np.sin(beta_star_rad) * np.cos(alpha_rad)],
            [0, 0, 2 * np.pi / c],
        ]
    )
    b_mat = b_mat / (2 * np.pi)

    return b_mat


def lattice_params_from_g_star_mat(g_star_mat: np.ndarray) -> tuple[float, float, float, float, float, float]:
    """Return lattice parameters from G* matrix

    Args:
        g_star_mat: G* = B^T * B, 3 by 3

    Return:
       lattice_params = (a, b, c, alpha, beta, gamma) in Angstrom and degrees
    """
    try:
        g_mat = np.linalg.inv(g_star_mat)
    except np.linalg.LinAlgError:
        print("G* matrix is singluar. Returnin lattice_param=(1, 1, 1, 90, 90, 90).")
        g_mat = np.identity(3)
    a = np.sqrt(g_mat[0, 0])
    b = np.sqrt(g_mat[1, 1])
    c = np.sqrt(g_mat[2, 2])
    alpha_rad = np.arccos((g_mat[1, 2] + g_mat[2, 1]) / 2)
    beta_rad = np.arccos((g_mat[0, 2] + g_mat[2, 0]) / 2)
    gamma_rad = np.arccos((g_mat[0, 1] + g_mat[1, 0]) / 2)

    lattice_params = (a, b, c, np.rad2deg(alpha_rad), np.rad2deg(beta_rad), np.rad2deg(gamma_rad))
    return lattice_params


def lattice_params_from_b_mat(b_mat: np.ndarray) -> tuple[float, float, float, float, float, float]:
    """Calculate lattice parameters from B matrix
    Args:
        b_mat: B matrix, 3 by 3

    Return:
       lattice_params = (a, b, c, alpha, beta, gamma) in Angstrom and degrees
    """
    g_star_mat = np.matmul(np.transpose(b_mat), b_mat)
    return lattice_params_from_g_star_mat(g_star_mat)


# -----------------------------------------------------
# U matrix math
# -----------------------------------------------------

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


# TODO
def _find_ub_from_three_peaks(
    self,
    peaks: tuple[Peak, Peak, Peak],
) -> UBConf:
    """Find UB matrix from three observed peaks for a given goniomete"""
    ubconf = UBConf()
    return ubconf


# TODO
def _find_ub_from_multiple_peaks(
    self,
    peaks: tuple[Peak, ...],
) -> UBConf:
    """Find UB matrix from more than three observed peaks for a given goniomete"""
    ubconf = UBConf()
    return ubconf


@staticmethod
def norm_mat(t1, t2, t3):
    mat = np.array([t1 / np.linalg.norm(t1), t2 / np.linalg.norm(t2), t3 / np.linalg.norm(t3)]).T
    return mat


def _t_mat_minimal_tilt(self, hkl: np.ndarray):
    """Build matrix T assuming minimal goniometer tilt angles"""

    if not isinstance(self.sample, Xtal):
        raise ValueError("Sample needs to be Xtal class for UB calculation.")
    if self.sample.ub_mat is None:
        raise ValueError("UB matrix is unknown.")
    if self.sample.plane_normal is None:
        raise ValueError("Plane normal vector is not known.")
    if self.sample.in_plane_ref is None:
        raise ValueError("In-plnae reference vector is not known.")

    EPS = 1e-8  # zero

    plane_normal = np.array(self.sample.plane_normal)
    in_plane_ref = np.array(self.sample.in_plane_ref)
    ub_mat = self.sample.ub_mat

    if self.SPICE_CONVENTION:  # suffle the order following SPICE convention
        plane_normal = spice_to_mantid(plane_normal)
        in_plane_ref = spice_to_mantid(in_plane_ref)
        ub_mat = spice_to_mantid(ub_mat)

    q = ub_mat @ hkl
    t1 = q / np.linalg.norm(q)

    if np.dot(t1, plane_normal) < EPS:  # t1 in plane
        t3 = plane_normal
        t2 = np.cross(t3, t1)
        return TAS.norm_mat(t1, t2, t3)

    # t1 not in plane, need to change tilts
    if np.linalg.norm(np.cross(plane_normal, t1)) < EPS:
        # oops, t1 along plane_normal
        t2 = in_plane_ref
        t3 = np.cross(t1, t2)
        return TAS.norm_mat(t1, t2, t3)

    else:
        t2p = np.cross(plane_normal, t1)
        t3 = np.cross(t1, t2p)
        t2 = np.cross(t3, t1)
        return TAS.norm_mat(t1, t2, t3)
