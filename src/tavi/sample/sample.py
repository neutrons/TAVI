# -*- coding: utf-8 -*-
import json

import numpy as np


class Sample(object):
    """
    Sample class

    Attributes:
        shape (str): "cuboid" or "cylindrical"
        width (float): in units of cm
        height (float): in units of cm
        depth (float): in units of cm

        mosaic_h (fload): in units of minutes of arc
        mosaic_v (fload): verital mosaic if anisotropic, in units of minutes of arc

        a, b, c                 lattice constants in Angstrom
        alpha, beta, gamma      angles in degrees
        a_vec, b_vec, c_vec     real sapce lattice vector
        a_star, b_star, c_star  lattice constants in inverse Angstrom
        alpha_star, beta_star, gamma_star       reciprocal angles in degrees
        a_star_vec, b_star_vec, c_star_vec      reciprocal lattice vector
        i_star, j_star, k_star  bases for the reciprocal space lattice vectors

    Methods:
        real_vec_cart
        reciprocal_vec_cart
        reciprocal_latt_params
        reciprocal_basis
        b_mat

    Static Methods:
        v_alpha_beta_gamma_calc(alpha,m beta, gamma)

    Class Methods:
        from_json(path_to_json): construct sample from a json file


    """

    def __init__(
        self,
        lattice_params: tuple[float] = (1, 1, 1, 90, 90, 90),
    ) -> None:
        assert len(lattice_params) == 6, "Incomplete lattice parameters."
        for length in lattice_params[:3]:
            assert length > 0, "Lattice parameters smaller than zero."
        for angle in lattice_params[3:]:
            assert 0.0 < angle < 180.0, "Lattice angles out of range."

        self.update_lattice_parameters(lattice_params)
        # parameters for resolution calculation
        self.set_mosaic()  # with defalt values
        self.set_shape()  # with defalt values

    # TODO
    @classmethod
    def from_json(cls, path_to_json):
        """Alternate constructor from json"""

        with open(path_to_json, "r", encoding="utf-8") as file:
            sample_params = json.load(file)

        lattice_params = (
            sample_params["a"],
            sample_params["b"],
            sample_params["c"],
            sample_params["alpha"],
            sample_params["beta"],
            sample_params["gamma"],
        )
        sample = cls(lattice_params=lattice_params)

        ub_matrix = sample_params.get("ub_matrix")
        if ub_matrix is not None:
            sample.ub_matrix = np.array(ub_matrix).reshape(3, 3)

        plane_normal = sample_params.get("plane_normal")
        if plane_normal is not None:
            sample.plane_normal = np.array(plane_normal)

        shape = sample_params.get("shape")
        width = sample_params.get("width")
        height = sample_params.get("height")
        depth = sample_params.get("depth")
        if all([shape, width, height, depth]):
            sample.set_shape(shape, width, height, depth)

        mosaic_h = sample_params.get("mosaic_h")
        mosaic_v = sample_params.get("mosaic_v")
        if all([mosaic_h, mosaic_v]):
            sample.set_mosaic(mosaic_h, mosaic_v)

        return sample

    def set_shape(
        self,
        shape: str = "cuboid",
        width: float = 1.0,  # in cm
        height: float = 1.0,  # in cm
        depth: float = 1.0,  # in cm
    ) -> None:
        """set sample shape"""
        self.shape = shape
        self.width = width
        self.height = height
        self.depth = depth

    def set_mosaic(
        self,
        mosaic_h: float = 30,  # horizontal mosaic
        mosaic_v: float = 30,  # vertical mosaic
    ) -> None:
        """Set horizontal and vertical mosaic in units of minitues of arc"""
        self.mosaic_h = mosaic_h  # * min2rad
        self.mosaic_v = mosaic_v  # * min2rad

    def update_lattice_parameters(
        self,
        lattice_params: tuple[float] = (1, 1, 1, 90, 90, 90),
    ):
        """update real and reciprocal space lattice parameters and vectors"""

        a, b, c, alpha, beta, gamma = lattice_params
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        (
            self.a_vec,
            self.b_vec,
            self.c_vec,
        ) = self._real_space_vectors()
        (
            self.a_star,
            self.b_star,
            self.c_star,
            self.alpha_star,
            self.beta_star,
            self.gamma_star,
        ) = self.reciprocal_latt_params()

        (
            self.a_star_vec,
            self.b_star_vec,
            self.c_star_vec,
        ) = self._reciprocal_space_vectors()

        (
            self.i_star,
            self.j_star,
            self.k_star,
        ) = self._reciprocal_basis()

    @staticmethod
    def v_alpha_beta_gamma_calc(alpha, beta, gamma) -> float:
        """
        Calculate V_alpha_bet_gamma = Volume/(abc)
        Volume = a * (b x c)
        """
        cos_alpha = np.cos(np.deg2rad(alpha))
        cos_beta = np.cos(np.deg2rad(beta))
        cos_gamma = np.cos(np.deg2rad(gamma))
        v_alpha_beta_gamma = np.sqrt(
            1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 2 * cos_alpha * cos_beta * cos_gamma
        )
        return v_alpha_beta_gamma

    def _real_space_vectors(self) -> tuple[np.ndarray]:
        """
        Calculate the real space lattice vectors in Cartesian coordiantes
        """
        cos_alpha = np.cos(np.deg2rad(self.alpha))
        cos_beta = np.cos(np.deg2rad(self.beta))
        cos_gamma = np.cos(np.deg2rad(self.gamma))
        sin_gamma = np.sin(np.deg2rad(self.gamma))

        ac = np.array([self.a, 0, 0])
        bc = np.array(
            [
                self.b * cos_gamma,
                self.b * sin_gamma,
                0,
            ]
        )
        v_abg = Sample.v_alpha_beta_gamma_calc(self.alpha, self.beta, self.gamma)
        cc = np.array(
            [
                self.c * cos_beta,
                self.c * (cos_alpha - cos_gamma * cos_beta) / sin_gamma,
                self.c * v_abg / sin_gamma,
            ]
        )
        return (ac, bc, cc)

    def reciprocal_latt_params(self):
        """Calculate the reciprocal lattice parameter lengths and angles"""
        sin_alpha = np.sin(np.deg2rad(self.alpha))
        cos_alpha = np.cos(np.deg2rad(self.alpha))
        sin_beta = np.sin(np.deg2rad(self.beta))
        cos_beta = np.cos(np.deg2rad(self.beta))
        cos_gamma = np.cos(np.deg2rad(self.gamma))
        sin_gamma = np.sin(np.deg2rad(self.gamma))

        v_abg = Sample.v_alpha_beta_gamma_calc(self.alpha, self.beta, self.gamma)

        a_star = sin_alpha / self.a / v_abg * np.pi * 2
        b_star = sin_beta / self.b / v_abg * np.pi * 2
        c_star = sin_gamma / self.c / v_abg * np.pi * 2
        alpha_star = np.arccos((cos_beta * cos_gamma - cos_alpha) / sin_beta / sin_gamma)
        beta_star = np.arccos((cos_gamma * cos_alpha - cos_beta) / sin_alpha / sin_gamma)
        gamma_star = np.arccos((cos_alpha * cos_beta - cos_gamma) / sin_beta / sin_alpha)
        alpha_star = np.rad2deg(alpha_star)
        beta_star = np.rad2deg(beta_star)
        gamma_star = np.rad2deg(gamma_star)

        return (a_star, b_star, c_star, alpha_star, beta_star, gamma_star)

    def _reciprocal_space_vectors(self) -> tuple[np.ndarray]:
        """
        Calculate the reciprocal space lattice vectors in the Cartesian coordinates
        """
        v_abg = Sample.v_alpha_beta_gamma_calc(self.alpha, self.beta, self.gamma)
        v = v_abg * self.a * self.b * self.c
        prefactor = 2 * np.pi / v
        a_star_vec = np.cross(self.b_vec, self.c_vec) * prefactor
        b_star_vec = np.cross(self.c_vec, self.a_vec) * prefactor
        c_star_vec = np.cross(self.a_vec, self.b_vec) * prefactor

        return (a_star_vec, b_star_vec, c_star_vec)

    def hkl2q(self, hkl: tuple[float]) -> float:
        """Convert (h,k,l) to q, in units of inverse Angstrom"""
        assert len(hkl) == 3, "Length of (h,k,l) is not 3."
        (qh, qk, ql) = hkl
        q_vec = qh * self.a_star_vec + qk * self.b_star_vec + ql * self.c_star_vec
        q_norm = np.linalg.norm(q_vec)
        return q_norm

    def b_mat(self) -> np.ndarray:
        """
        Calculate the B matrix
        B * (h,k,l) gives Q in terms of i_star, j_star, k_star
        """
        alpha_star_deg = np.deg2rad(self.alpha_star)
        beta_star_deg = np.deg2rad(self.beta_star)
        gamma_star_deg = np.deg2rad(self.gamma_star)
        b_mat = np.array(
            [
                [
                    self.a_star,
                    self.b_star * np.cos(gamma_star_deg),
                    self.c_star * np.cos(beta_star_deg),
                ],
                [
                    0,
                    self.b_star * np.sin(gamma_star_deg),
                    -self.c_star * np.sin(beta_star_deg) * np.cos(alpha_star_deg),
                ],
                [0, 0, 2 * np.pi / self.c],
            ]
        )
        b_mat = b_mat / (2 * np.pi)
        # b_mat = np.round(b_mat, 8)
        return b_mat

    def _reciprocal_basis(self) -> tuple[np.ndarray]:
        """Calculate the reciprocal basis vectors i_star, j_star, k_star"""
        i_star = self.a_star_vec / np.linalg.norm(self.a_star_vec)
        a_star_perp = np.cross(
            np.cross(self.a_star_vec, self.b_star_vec),
            self.a_star_vec,
        )
        j_star = a_star_perp / np.linalg.norm(a_star_perp)
        k_star = np.cross(i_star, j_star)
        return (i_star, j_star, k_star)
