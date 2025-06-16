from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np

from tavi.lattice_algorithm import b_mat_from_lattice, lattice_params_from_g_star_mat
from tavi.utilities import MotorAngles, Peak, en2q, get_angle_from_triangle


def mantid_to_spice(v):
    t = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    return t.dot(v)


def spice_to_mantid(v):
    t = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
    return t.dot(v)


@dataclass
class UBConf:
    """Logs for UB matrix determination

    Attributes:
        ub_peaks (tuple of Peaks): peaks used to determine the UB matrix
        _u_mat (np.adarray): U matrix
        b_mat (np.adarray): B matrix
        _ub_matrix (np.adarray): UB matrix
        _plane_normal (np.adarray): normal vector in Qsample frame, goniometers at zero
        _in_plane_ref (np.adarray): in plane vector in Qsample frame, goniometers at zero

    Note:
        _u_mat, _ub_mat, _plnae_normal and _in_plane_ref uses Mandid/International
        Crystallography Table convention.
        When u_mat, ub_mat, plane_normal and in_plane_ref are requested, a conversion from
        Mantid is performed based on what convention is used
    """

    _ub_mat: np.ndarray
    convention: Literal["Mantid", "Spice"] = "Spice"
    _plane_normal: Optional[np.ndarray] = None
    _in_plane_ref: Optional[np.ndarray] = None
    _u_mat: Optional[np.ndarray] = None
    b_mat: Optional[np.ndarray] = None
    ub_peaks: Optional[tuple[Peak, ...]] = None

    def __init__(self, convention="Spice", **kwargs):
        self.convention = "Spice" if convention is None else convention
        for k, v in kwargs.items():
            if v is not None:
                self.__setattr__(k, v)

    def _from_mantid(self, v: np.ndarray):
        """Convert v from Mantid to other conventions"""
        match self.convention:
            case "Mantid":
                return v
            case "Spice":
                return mantid_to_spice(v)
            case _:
                raise ValueError("Unrecogonized convention: " + self.convention)

    def _to_mantid(self, v: np.ndarray):
        """Convert v from some convention to Mantid convention"""
        match self.convention:
            case "Mantid":
                return v
            case "Spice":
                return spice_to_mantid(v)
            case _:
                raise ValueError("Unrecogonized convention: " + self.convention)

    @property
    def ub_mat(self):
        return self._from_mantid(self._ub_mat)

    @ub_mat.setter
    def ub_mat(self, value):
        self._ub_mat = self._to_mantid(value)

    @property
    def u_mat(self):
        return self._from_mantid(self._u_mat)

    @u_mat.setter
    def u_mat(self, value):
        self._u_mat = self._to_mantid(value)

    @property
    def plane_normal(self):
        return self._from_mantid(self._plane_normal)

    @plane_normal.setter
    def plane_normal(self, value):
        self._plane_normal = self._to_mantid(value)

    @property
    def in_plane_ref(self):
        return self._from_mantid(self._in_plane_ref)

    @in_plane_ref.setter
    def in_plane_ref(self, value):
        self._in_plane_ref = self._to_mantid(value)

    def __repr__(self):
        ub_str = ""
        for name, value in self.__dict__.items():
            match name:
                case "convention":
                    conv_str = "UBconf using " + value + " convention.\n"
                    ub_str += conv_str
                case "ub_peaks":
                    ub_peaks_str = ""
                    sz = len(value)
                    for peak in value:
                        ub_peaks_str += str(peak.hkl) + ", "
                    ub_str += f"UB matrix determined from {sz} peaks: {ub_peaks_str}"
                case _:
                    if value is not None:
                        if name in ["_ub_mat", "_u_mat", "_plane_normal", "_in_plane_ref"]:
                            value = self._from_mantid(value)
                            name = name.strip("_")
                        value_str = np.array2string(
                            value,
                            precision=6,
                            suppress_small=True,
                            separator=",",
                        ).replace("\n", "")
                        ub_str += f"{name}=" + value_str + "\n"
        return ub_str


# -----------------------------------------------------
# Angle and Q math
# -----------------------------------------------------


def q_lab(ei: float, ef: float, theta: float, phi: float = 0) -> np.ndarray:
    """
    Return momentum transfer vector q in lab frame

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
    theta_rad = np.radians(theta)
    phi_rad = np.radians(phi)
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
) -> float:
    """Calculate the angle between ki and kf. .

    Note:
        Either b_mat or ub_mat would work, since U^T.U=U^-1.U=1"""

    ki = en2q(ei)
    kf = en2q(ef)
    q_norm = q_norm_from_hkl(hkl, b_mat)
    two_theta_radians = get_angle_from_triangle(ki, kf, q_norm)
    return two_theta_radians


def psi_from_hkle(
    hkl: tuple[float, float, float],
    ei: float,
    ef: float,
    b_mat: np.ndarray,
) -> float:
    """Calculate the angle between ki and Q=ki-kf.

    Note:
        Either b_mat or ub_mat would work, since U^T.U=U^-1.U=1"""

    ki = en2q(ei)
    kf = en2q(ef)
    q_norm = q_norm_from_hkl(hkl, b_mat)
    psi_radians = get_angle_from_triangle(ki, q_norm, kf)

    return psi_radians


def plane_normal_from_two_peaks(
    u_mat: np.ndarray,
    b_mat: np.ndarray,
    hkl1: tuple[float, float, float],
    hkl2: tuple[float, float, float],
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate plane_normal and in_plane_ref in Q_sample coordinate,
    giving Miller indices of two pekas in plane
    """
    t1 = b_mat.dot(hkl1)
    t2p = b_mat.dot(hkl2)
    t3 = np.cross(t1, t2p)
    t3 /= np.linalg.norm(t3)
    plane_normal = u_mat.dot(t3)
    plane_normal = -plane_normal if plane_normal[1] < 0 else plane_normal
    in_plane_ref = u_mat.dot(t1 / np.linalg.norm(t1))
    return (plane_normal, in_plane_ref)


def angle_between_two_hkls(
    hkl1: tuple[float, float, float],
    hkl2: tuple[float, float, float],
    b_mat: np.ndarray,
):
    q1 = b_mat.dot(hkl1)
    q2 = b_mat.dot(hkl2)
    c = np.dot(q1, q2) / np.linalg.norm(q1) / np.linalg.norm(q2)
    angle = np.arccos(np.clip(c, -1, 1))
    return np.degrees(angle)


def angle_between_two_motor_positions(
    angles1: MotorAngles,
    angles2: MotorAngles,
    r_mat_inv,
    ei: float,
    ef: float,
):
    q1 = r_mat_inv(angles1).dot(q_lab(ei, ef, angles1.two_theta))
    q2 = r_mat_inv(angles2).dot(q_lab(ei, ef, angles2.two_theta))
    c = np.dot(q1, q2) / np.linalg.norm(q1) / np.linalg.norm(q2)  # -> cosine of the angle
    angle = np.arccos(np.clip(c, -1, 1))
    return np.degrees(angle)


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

    ub_mat = ub_conf._ub_mat
    plane_normal = ub_conf._plane_normal
    in_plane_ref = ub_conf._in_plane_ref

    if (plane_normal is None) or (in_plane_ref is None):
        raise ValueError("plane_normal or in_plane_ref cannot be None.")

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


def plane_normal_from_one_peak(
    hkl: tuple[float, float, float],
    angles: MotorAngles,
    r_mat_inv,
    ub_mat: np.ndarray,
):
    plane_normal = r_mat_inv(angles).dot((0, 1, 0))
    q = ub_mat.dot(hkl)
    in_plane_ref = np.cross(plane_normal, q / np.linalg.norm(q))
    return plane_normal, in_plane_ref


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


def b_mat_from_ub_matrix(ub_matrix: np.ndarray) -> np.ndarray:
    """Calculate the B matrix from UB matrix"""
    g_star_mat = ub_matrix.T @ ub_matrix
    lattice_params = lattice_params_from_g_star_mat(g_star_mat)
    return b_mat_from_lattice(lattice_params)


def u_mat_from_ub_matrix(ub_matrix: np.ndarray) -> np.ndarray:
    """Calculate the U matrix from UB matrix"""
    b_mat = b_mat_from_ub_matrix(ub_matrix)
    inv_b_mat = np.linalg.inv(b_mat)
    return ub_matrix.dot(inv_b_mat)


# -----------------------------------------------------
# Find UB matrix from peaks
# -----------------------------------------------------


def find_u_from_one_peak_and_scattering_plane(
    peak: Peak,
    scattering_plane: tuple,
    b_mat: np.ndarray,
    r_mat_inv,
    ei: float,
    ef: float,
):
    """Calculate U matrix from one peak and a scattering plane defined by two vectors

    Note:
        v1, v2= scatttering plane needs to be right handed, making v3 = v1 x v3 pointing up.
    """
    q_hkl1 = b_mat.dot(peak.hkl)
    v1, v2 = scattering_plane
    coeff1, coeff2 = np.dot(v1, peak.hkl), np.dot(v2, peak.hkl)

    if np.abs(coeff1) > np.abs(coeff2):
        q_hkl2p = b_mat.dot(v2) * np.sign(coeff1)
    else:
        q_hkl2p = b_mat.dot(v1) * np.sign(coeff2)
    q_hkl3 = np.cross(q_hkl1, q_hkl2p)
    q_hkl2 = np.cross(q_hkl3, q_hkl1)

    q_hkl_mat = np.array(
        [
            q_hkl1 / np.linalg.norm(q_hkl1),
            q_hkl2 / np.linalg.norm(q_hkl2),
            q_hkl3 / np.linalg.norm(q_hkl3),
        ]
    ).T

    q_lab1 = q_lab(ei=ei, ef=ef, theta=peak.angles.two_theta)
    q_sample1 = r_mat_inv(peak.angles).dot(q_lab1)
    q_sample3 = r_mat_inv(peak.angles).dot((0, 1, 0))
    q_sample2 = np.cross(q_sample3, q_sample1)

    q_sample1 = q_sample1 / np.linalg.norm(q_sample1)
    q_sample2 = q_sample2 / np.linalg.norm(q_sample2)
    q_sample3 = q_sample3 / np.linalg.norm(q_sample3)

    q_sample_mat = np.array([q_sample1, q_sample2, q_sample3]).T

    u_mat = q_sample_mat.dot(np.linalg.inv(q_hkl_mat))

    return u_mat


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
    q_hkl2 = np.cross(q_hkl3, q_hkl1)

    q_hkl_mat = np.array(
        [
            q_hkl1 / np.linalg.norm(q_hkl1),
            q_hkl2 / np.linalg.norm(q_hkl2),
            q_hkl3 / np.linalg.norm(q_hkl3),
        ]
    ).T

    q_lab1 = q_lab(ei=ei, ef=ef, theta=peak1.angles.two_theta)
    q_lab2 = q_lab(ei=ei, ef=ef, theta=peak2.angles.two_theta)

    # Goniometer angles all zeros in q_sample frame
    q_sample1 = r_mat_inv(peak1.angles).dot(q_lab1)
    q_sample2p = r_mat_inv(peak2.angles).dot(q_lab2)
    q_sample3 = np.cross(q_sample1, q_sample2p)
    q_sample2 = np.cross(q_sample3, q_sample1)

    q_sample1 = q_sample1 / np.linalg.norm(q_sample1)
    q_sample2 = q_sample2 / np.linalg.norm(q_sample2)
    q_sample3 = q_sample3 / np.linalg.norm(q_sample3)

    q_sample_mat = np.array([q_sample1, q_sample2, q_sample3]).T

    u_mat = q_sample_mat.dot(np.linalg.inv(q_hkl_mat))

    return u_mat


def find_ub_from_three_peaks(
    peaks: tuple[Peak, Peak, Peak],
    r_mat_inv,
    ei: float,
    ef: float,
) -> np.ndarray:
    """Find UB matrix from three non-coplanar peaks for a given goniometer"""
    peak1, peak2, peak3 = peaks

    hkl_mat = np.array([peak1.hkl, peak2.hkl, peak3.hkl]).T

    q_lab1 = q_lab(ei=ei, ef=ef, theta=peak1.angles.two_theta)
    q_lab2 = q_lab(ei=ei, ef=ef, theta=peak2.angles.two_theta)
    q_lab3 = q_lab(ei=ei, ef=ef, theta=peak3.angles.two_theta)
    # Goniometer angles all zeros in q_sample frame
    q_sample1 = r_mat_inv(peak1.angles).dot(q_lab1) / (2 * np.pi)
    q_sample2 = r_mat_inv(peak2.angles).dot(q_lab2) / (2 * np.pi)
    q_sample3 = r_mat_inv(peak3.angles).dot(q_lab3) / (2 * np.pi)

    q_sample_mat = np.array([q_sample1, q_sample2, q_sample3]).T

    ub_mat = q_sample_mat.dot(np.linalg.inv(hkl_mat))

    return ub_mat


def find_ub_from_multiple_peaks(
    peaks: tuple[Peak, ...],
    r_mat_inv,
    ei: float,
    ef: float,
) -> Optional[np.ndarray]:
    """Find UB matrix from more than three observed peaks for a given goniomete"""

    num = len(peaks)
    qv_mat = np.zeros((3, 3))
    vv_mat = np.zeros((3, 3))
    for i in range(num):
        hkl = peaks[i].hkl
        angles = peaks[i].angles
        qi_lab = q_lab(ei=ei, ef=ef, theta=angles.two_theta)
        qi = r_mat_inv(angles).dot(qi_lab) / (2 * np.pi)
        for j in range(3):
            for k in range(3):
                qv_mat[j, k] += qi[k] * hkl[j]
                vv_mat[j, k] += hkl[k] * hkl[j]

    ub_mat = np.matmul(qv_mat.T, np.linalg.inv(vv_mat).T)

    return ub_mat
