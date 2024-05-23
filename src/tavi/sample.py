import numpy as np


class Sample(object):
    """
    Attributes:
        a, b, c                 lattice constants in Angstrom
        alpha, beta, gamma      angles in degrees
        a_vec, b_vec, c_veßc    real sapce lattice vector
        v_abg                   V_alpha_beta_gamma = unit_cell_volume/(abc)
        a_star, b_star, c_star  lattice constants in inverse Angstrom
        alpha_star, beta_star, gamma_star       reciprocal angles in degrees
        i_star, j_star, k_star  bases for the reciprocal space lattice vectors
        a_star_vec, b_star_vec, c_star_vec      reciprocal lattice vector


    Methods:
        real_vec_cart
        reciprocal_vec_cart
        reciprocal_latt_params
        reciprocal_basis
        b_mat

    Static Methods:
        v_alpha_beta_gamma_calc(alpha,m beta, gamma)


    """

    def __init__(self, lattice_params=(1, 1, 1, 90, 90, 90)):
        a, b, c, alpha, beta, gamma = lattice_params
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.v_abg = Sample.v_alpha_beta_gamma_calc(alpha, beta, gamma)
        self.a_vec, self.b_vec, self.c_vec = self.real_vec_cart()
        (
            self.a_star,
            self.b_star,
            self.c_star,
            self.alpha_star,
            self.beta_star,
            self.gamma_star,
        ) = self.reciprocal_latt_params()
        self.a_star_vec, self.b_star_vec, self.c_star_vec = self.reciprocal_vec_cart()
        self.i_star, self.j_star, self.k_star = self.reciprocal_basis()

    @staticmethod
    def v_alpha_beta_gamma_calc(alpha, beta, gamma):
        """
        Calculate V_alpha_bet_gamma = Volume/(abc)
        Volume = a * (b x c)
        """
        cos_alpha = np.cos(alpha / 180 * np.pi)
        cos_beta = np.cos(beta / 180 * np.pi)
        cos_gamma = np.cos(gamma / 180 * np.pi)
        v_alpha_beta_gamma = np.sqrt(
            1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 2 * cos_alpha * cos_beta * cos_gamma
        )
        return v_alpha_beta_gamma

    def real_vec_cart(self):
        """
        Calculate the real space lattice vectors in Cartesian coordiantes
        """
        cos_alpha = np.cos(self.alpha / 180 * np.pi)
        cos_beta = np.cos(self.beta / 180 * np.pi)
        cos_gamma = np.cos(self.gamma / 180 * np.pi)
        sin_gamma = np.sin(self.gamma / 180 * np.pi)

        ac = np.array([self.a, 0, 0])
        bc = np.array(
            [
                self.b * cos_gamma,
                self.b * sin_gamma,
                0,
            ]
        )

        cc = np.array(
            [
                self.c * cos_beta,
                self.c * (cos_alpha - cos_gamma * cos_beta) / sin_gamma,
                self.c * self.v_abg / sin_gamma,
            ]
        )
        ac = np.round(ac, 8)
        bc = np.round(bc, 8)
        cc = np.round(cc, 8)
        return (ac, bc, cc)

    def reciprocal_latt_params(self):
        """Calculate the reciprocal lattice parameters and angles"""
        sin_alpha = np.sin(self.alpha / 180 * np.pi)
        cos_alpha = np.cos(self.alpha / 180 * np.pi)
        sin_beta = np.sin(self.beta / 180 * np.pi)
        cos_beta = np.cos(self.beta / 180 * np.pi)
        cos_gamma = np.cos(self.gamma / 180 * np.pi)
        sin_gamma = np.sin(self.gamma / 180 * np.pi)

        a_star = sin_alpha / self.a / self.v_abg * np.pi * 2
        b_star = sin_beta / self.b / self.v_abg * np.pi * 2
        c_star = sin_gamma / self.c / self.v_abg * np.pi * 2
        alpha_star = np.arccos((cos_beta * cos_gamma - cos_alpha) / sin_beta / sin_gamma) / np.pi * 180
        beta_star = np.arccos((cos_gamma * cos_alpha - cos_beta) / sin_alpha / sin_gamma) / np.pi * 180
        gamma_star = np.arccos((cos_alpha * cos_beta - cos_gamma) / sin_beta / sin_alpha) / np.pi * 180
        a_star = np.round(a_star, 8)
        b_star = np.round(b_star, 8)
        c_star = np.round(c_star, 8)
        alpha_star = np.round(alpha_star, 8)
        beta_star = np.round(beta_star, 8)
        gamma_star = np.round(gamma_star, 8)

        return (a_star, b_star, c_star, alpha_star, beta_star, gamma_star)

    def reciprocal_vec_cart(self):
        """
        Calculate the reciprocal space lattice vectors in the Cartesian coordinates
        """
        v = self.v_abg * self.a * self.b * self.c
        prefactor = 2 * np.pi / v
        a_star_vec = np.cross(self.b_vec, self.c_vec) * prefactor
        b_star_vec = np.cross(self.c_vec, self.a_vec) * prefactor
        c_star_vec = np.cross(self.a_vec, self.b_vec) * prefactor

        return (a_star_vec, b_star_vec, c_star_vec)

    def b_mat(self):
        """
        Calculate the B matrix
        B * (h,k,l) gives Q in terms of i_star, j_star, k_star
        """
        b_mat = np.array(
            [
                [
                    self.a_star,
                    self.b_star * np.cos(self.gamma_star / 180 * np.pi),
                    self.c_star * np.cos(self.beta_star / 180 * np.pi),
                ],
                [
                    0,
                    self.b_star * np.sin(self.gamma_star / 180 * np.pi),
                    -self.c_star * np.sin(self.beta_star / 180 * np.pi) * np.cos(self.alpha / 180 * np.pi),
                ],
                [0, 0, 2 * np.pi / self.c],
            ]
        )
        b_mat = np.round(b_mat, 8)
        return b_mat

    def reciprocal_basis(self):
        """Calculate the reciprocal basis vectors i_star, j_star, k_star"""
        i_star = self.a_star_vec / np.linalg.norm(self.a_star_vec)
        a_star_perp = np.cross(np.cross(self.a_star_vec, self.b_star_vec), self.a_star_vec)
        j_star = a_star_perp / np.linalg.norm(a_star_perp)
        k_star = np.cross(i_star, j_star)
        return (i_star, j_star, k_star)

    def set_UB(self, u=[1, 0, 0], v=[0, 0, 1]):
        self.u = u
        self.v = v

    def set_mosaic(self, mosaic_params):
        pass

    def set_shape(self, shape_params):
        pass


if __name__ == "__main__":
    sample = Sample(lattice_params=(1, 1, 1, 90, 90, 120))
    # print(sample.a_star / 2 / np.pi)
    # nprint(sample.b_mat() / (2 * np.pi) @ np.array([1, 0, 0]))
    # print(sample.b_star_vec / np.pi / 2)
    # print(sample.reciprocal_basis())
    print(sample.b_mat() / 2 / np.pi)
    print(sample.b_mat() @ np.array([0, 1, 0]) / 2 / np.pi)
