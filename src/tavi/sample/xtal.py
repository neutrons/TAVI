import numpy as np
from tavi.utilities import *
from tavi.sample.sample import Sample


class Xtal(Sample):
    """Singel crystal sample

    Attibutes:
        type (str): "xtal"
        inv_ub_matrix
        in_plane_ref: in plane vector in Qsample frame, goniometers at zero
        plane_normal: normal vector in Qsample frame, goniometers at zero
        u (tuple)
        v (tuple)
    Methods:


    """

    def __init__(self, lattice_params=(1, 1, 1, 90, 90, 90)):
        super().__init__(lattice_params)
        self.type = "xtal"

        self.ub_peaks = None
        self.ub_angles = None
        self.ub_matrix = None
        self.inv_ub_matrix = None

        self.plane_normal = None
        self.in_plane_ref = None

        self.i_star, self.j_star, self.k_star = self.reciprocal_basis()

    @property
    def u(self):
        """u vector, in reciprocal lattice unit, along beam"""
        return self.ub_matrix_to_uv(self.ub_matrix)[0]

    @property
    def v(self):
        """
        v vector, in reciprocal lattice unit,
        in the horizaontal scattering plane
        """
        return self.ub_matrix_to_uv(self.ub_matrix)[1]

    @staticmethod
    def ub_matrix_to_uv(ub_matrix):
        """ "Calculate u and v vector from UB matrix"""
        inv_ub_matrix = np.linalg.inv(ub_matrix)
        u = inv_ub_matrix @ np.array([0, 0, 1])
        v = inv_ub_matrix @ np.array([1, 0, 0])
        return (u, v)

    @staticmethod
    def ub_matrix_to_lattice_params(ub_matrix):
        """Calculate lattice parameters from UB matrix"""
        g_mat = np.linalg.inv(ub_matrix.T @ ub_matrix)
        a = np.sqrt(g_mat[0, 0])
        b = np.sqrt(g_mat[1, 1])
        c = np.sqrt(g_mat[2, 2])
        alpha = np.arccos((g_mat[1, 2] + g_mat[2, 1]) / (2 * b * c)) * rad2deg
        beta = np.arccos((g_mat[0, 2] + g_mat[2, 0]) / (2 * a * c)) * rad2deg
        gamma = np.arccos((g_mat[0, 1] + g_mat[1, 0]) / (2 * a * b)) * rad2deg
        return (a, b, c, alpha, beta, gamma)

    def uv_to_ub_matrix(self, u, v):
        """Calculate UB matrix from u and v vector, and lattice parameters"""

        b_mat = self.b_mat()
        t1 = b_mat @ u
        t2_prime = b_mat @ v
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
        ub_matrix = q_mat @ np.linalg.inv(t_mat) @ b_mat

        return ub_matrix