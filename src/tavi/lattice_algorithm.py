import numpy as np

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
    lattice_params: tuple[float, float, float, float, float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the real space lattice vectors in Cartesian coordinates

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
    lattice_params: tuple[float, float, float, float, float, float],
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
    lattice_params: tuple[float, float, float, float, float, float],
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
    lattice_params: tuple[float, float, float, float, float, float],
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


def lattice_params_from_g_star_mat(
    g_star_mat: np.ndarray,
) -> tuple[float, float, float, float, float, float]:
    """Return lattice parameters from G* matrix

    Args:
        g_star_mat: G* = B^T * B, 3 by 3

    Return:
       lattice_params = (a, b, c, alpha, beta, gamma) in Angstrom and degrees
    """
    try:
        g_mat = np.linalg.inv(g_star_mat)
    except np.linalg.LinAlgError:
        print("G* matrix is singular. Returning lattice_param=(1, 1, 1, 90, 90, 90).")
        g_mat = np.identity(3)

    a = np.sqrt(g_mat[0, 0])
    b = np.sqrt(g_mat[1, 1])
    c = np.sqrt(g_mat[2, 2])
    alpha_rad = np.arccos((g_mat[1, 2] + g_mat[2, 1]) / (2 * b * c))
    beta_rad = np.arccos((g_mat[0, 2] + g_mat[2, 0]) / (2 * a * c))
    gamma_rad = np.arccos((g_mat[0, 1] + g_mat[1, 0]) / (2 * a * b))
    alpha = np.rad2deg(alpha_rad)
    beta = np.rad2deg(beta_rad)
    gamma = np.rad2deg(gamma_rad)
    lattice_params = (a, b, c, alpha, beta, gamma)
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
