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
    Return qlab matrix. Follows mantid convention. Eq.8.

    Args:
        ei: incident energy, in meV.
        ef: final energy, in meV.
        theta: angle between scattered beam and z axis, in degree.
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


def proj_with_constraint(mat: np.ndarray, b: np.ndarray, k: int) -> np.ndarray:
    """
    Reduce a quadratic form under the constraint x_k = b·x'' (R. A. Robinson, Analytical Techniques, Eq. 18).

        |a_00, a_01, a_02, a_03|
    A = |a_10, a_11, a_12, a_13|
        |a_20, a_21, a_22, a_23|
        |a_30, a_31, a_32, a_33|

    The indexing of A''' is the same as A but with k removed. Take k = 2 as an example:

           |a'''_00, a'''_01, a'''_03|
    A''' = |a'''_10, a'''_11, a'''_13| 
           |a'''_30, a'''_30, a'''_33|.

    Python does not allow leaving a gap in numbering, so wee need to account for that with row_k, col_k and reduced.

    Again, k = 2 as an example:

    row_k = |a_20, a_21, a_23|
    col_k = |a_02, a_12, a_32|
    a_kk = |a_kk|.

    mat: n x n resolution matrix.
    b: 1-D array of the n-1 coefficients of x'' (all components of x except k), per Eq. 16.
    k: index of the component eliminated by the constraint.
    """
    n = len(mat)
    # Row and column are taken separately so Eq. 18 holds even if mat is not symmetric.
    row_k = np.delete(mat[k, :], k)  # a_kj for j != k, re-indexed to x''
    col_k = np.delete(mat[:, k], k)  # a_ik for i != k, re-indexed to x''
    a_kk = mat[k, k]
    reduced = np.delete(np.delete(mat, k, axis=0), k, axis=1)  # a_ij for i, j != k

    new_mat = np.zeros((n - 1, n - 1))
    for i in range(n - 1):
        for j in range(n - 1):
            new_mat[i][j] = reduced[i][j] + b[i] * row_k[j] + b[j] * col_k[i] + b[i] * b[j] * a_kk
    return new_mat


def get_angle_vec(v1: np.ndarray, v2: np.ndarray) -> float:
    """Get the angle in degrees between two vectors v1 and v2."""
    return np.degrees(np.arccos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))


def mantid_to_spice(v: np.ndarray, version: str = "new") -> np.ndarray:
    """Convert a vector from the Mantid frame to the SPICE frame."""
    if version == "new":
        t = np.array([[-1, 0, 0], [0, 0, 1], [0, 1, 0]])
    else:
        t = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    return t.dot(v)


def spice_to_mantid(v: np.ndarray, version: str = "new") -> np.ndarray:
    """Convert a vector from the SPICE frame to the Mantid frame."""
    if version == "new":
        t = np.array([[-1, 0, 0], [0, 0, 1], [0, 1, 0]])
    else:
        t = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
    return t.dot(v)


def get_side_from_triangle(a: float, b: float, angle: float) -> float:
    """In a triangle with sides a,b and angle between a and b, get the third side c."""
    if (np.abs(a) < 1e-6) or (np.abs(b) < 1e-6):
        raise ValueError("Triangle cannot be closed.")
    return np.sqrt(a**2 + b**2 - 2 * a * b * np.cos(angle))


def translation(number: float, character: str) -> str:
    """Used in projection. Borrowed from SHIVER histogram_parameters.py. Thanks to Andrei."""
    if number == 0:
        return "0"
    if number == 1:
        return character
    if number == -1:
        return "-" + character
    return str(number) + character


def update_dimension_names(resolution_frame: tuple) -> list:
    """Update the combo box dimension selection items based on the projection values.Borrowed from SHIVER histogram_parameters.py. Thanks to Andrei."""
    chars = ["H", "K", "L"]
    combo_dimensions = []
    for i in range(3):
        index_max = np.argmax(np.abs(resolution_frame[i]))
        combo_dimensions.append("[" + ",".join([translation(x, chars[index_max]) for x in resolution_frame[i]]) + "]")
    combo_dimensions.append("DeltaE")
    return combo_dimensions
