# -*- coding: utf-8 -*-
import json
from typing import Literal, Optional

import numpy as np

from tavi.lattice_algorithm import (
    b_mat_from_lattice,
    lattice_params_from_b_mat,
    lattice_params_from_g_star_mat,
    real_space_vectors,
    reciprocal_latt_params,
    reciprocal_space_vectors,
)
from tavi.ub_algorithm import UBConf


class Sample(object):
    """
    Sample class

    Attributes:
        type (str): "crystal" or "powder"
        json_dict (dict): dictionary from json file
        shape (str): "cuboid" or "cylindrical"
        width (float): in units of cm
        height (float): in units of cm
        depth (float): in units of cm

        mosaic_h (fload): in units of minutes of arc
        mosaic_v (fload): verital mosaic if anisotropic, in units of minutes of arc

        a, b, c                 lattice constants in Angstrom
        alpha, beta, gamma      angles in degrees

        ub_conf (UBConf): latest UB configuration information
        u (tuple): u vector, (h, k, l) along the beam direction when all goniometer angles are zero
        v (tuple): v vector, (h ,k, l) in the scattering plane

    Methods:
        update_lattice_parameters: recommended method to set lattice parameters.
                                   This method update all lattice related vectors
        reciprocal_latt_params: returns (a_star, b_star, c_star, alpha_star, beta_star, gamma_star)
        lattice_params: return (a, b, c, alpha, beta, gamma)
        real_space_vectors: return (a_vec, b_vec, c_vec)
        reciprocal_latt_params: return (a_star, b_star, c_star, alpha_star, beta_star, gamma_star)
            lattice constants in inverse Angstrom, reciprocal angles in degrees
        reciprocal_space_vectors: return (a_star_vec, b_star_vec, c_star_vec)
        b_mat: B matrix calculated from lattice parameters
        q_norm: length of q given (h, k ,l)


    Class Methods:
        from_json(path_to_json): construct sample from a json file


    """

    def __init__(
        self,
        lattice_params: tuple[float, float, float, float, float, float] = (1, 1, 1, 90, 90, 90),
    ) -> None:
        """Initialization from lattice parameters"""
        self.type: Literal["crystal", "powder"] = "crystal"
        self.json_dict: Optional[dict] = None
        # lattice parameters
        self.a: float
        self.b: float
        self.c: float
        self.alpha: float
        self.beta: float
        self.gamma: float

        # parameters for resolution calculation
        self.shape: str = "cuboid"
        self.width: float = 1.0  # in cm
        self.height: float = 1.0  # in cm
        self.depth: float = 1.0  # in cm
        self.mosaic_h: float = 30  # horizontal mosaic in minutes of arc
        self.mosaic_v: float = 30  # vertical mosaic in minutes of arc

        self.ub_conf: Optional[UBConf] = None
        self.update_lattice_parameters(lattice_params)

    def __repr__(self):
        cls = self.__class__.__name__
        cls_str = (
            f"{cls} class, a={self.a:.4f}, b={self.b:.4f}, c={self.c:.4f}, "
            + f"alpha={self.alpha:.2f}, beta={self.beta:.2f}, gamma={self.gamma:.2f},"
        )
        return cls_str

    @classmethod
    def from_json(cls, path_to_json):
        """Alternate constructor from json"""

        with open(path_to_json, "r", encoding="utf-8") as file:
            sample_params_dict = json.load(file)

        lattice_params = (
            sample_params_dict["a"],
            sample_params_dict["b"],
            sample_params_dict["c"],
            sample_params_dict["alpha"],
            sample_params_dict["beta"],
            sample_params_dict["gamma"],
        )
        sample = cls(lattice_params=lattice_params)
        # sample.json_dict = sample_params_dict

        # setting sample shape
        shape = sample_params_dict.get("shape")
        width = sample_params_dict.get("width")
        height = sample_params_dict.get("height")
        depth = sample_params_dict.get("depth")
        if None not in [shape, width, height, depth]:
            sample.set_shape(shape, width, height, depth)

        # setting sample mosaic
        mosaic_h = sample_params_dict.get("mosaic_h")
        mosaic_v = sample_params_dict.get("mosaic_v")
        if None not in [mosaic_h, mosaic_v]:
            sample.set_mosaic(mosaic_h, mosaic_v)

        # setting UB matrix
        if (ub_matrix := sample_params_dict.get("ub_matrix")) is not None:
            ub_matrix = np.array(ub_matrix).reshape(3, 3)
        if (plane_normal := sample_params_dict.get("plane_normal")) is not None:
            plane_normal = np.array(plane_normal)
        if (in_plane_ref := sample_params_dict.get("in_plane_ref")) is not None:
            in_plane_ref = np.array(in_plane_ref)

        convention = sample_params_dict.get("convention")
        sample.ub_conf = UBConf(
            convention=convention,
            ub_mat=ub_matrix,
            plane_normal=plane_normal,
            in_plane_ref=in_plane_ref,
        )

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
        mosaic_h: float = 30,  # horizontal mosaic in minutes of arc
        mosaic_v: float = 30,  # vertical mosaic in minutes of arc
    ) -> None:
        """Set horizontal and vertical mosaic in units of minitues of arc"""
        self.mosaic_h = mosaic_h
        self.mosaic_v = mosaic_v

    @property
    def _mosaic_h(self) -> float:
        """horizontal mosaic in radian, reserved for resolution calculation"""
        return np.deg2rad(self.mosaic_h / 60)

    @property
    def _mosaic_v(self) -> float:
        """vertical mosaic in radian, reserved for resolution calculation"""
        return np.deg2rad(self.mosaic_v / 60)

    def update_lattice_parameters(
        self,
        lattice_params: tuple[float, float, float, float, float, float] = (1, 1, 1, 90, 90, 90),
    ):
        """update real and reciprocal space lattice parameters and vectors"""

        try:
            a, b, c, alpha, beta, gamma = lattice_params
        except ValueError:
            print(f"Trying to set lattice parameters={lattice_params}.")

        for length in (a, b, c):
            if length <= 0:
                raise ValueError(f"Lattice parameters={lattice_params[0:3]} smaller than zero.")

        for angle in (alpha, beta, gamma):
            if 0.0 > angle or angle > 180.0:
                raise ValueError(f"Lattice angles= {lattice_params[3:6]}out of range.")

        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

    @property
    def lattice_params(self):
        return (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)

    @property
    def real_space_vectors(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Reutrun real space lattice vectors in Cartesian coordiantes"""
        return real_space_vectors(self.lattice_params)

    @property
    def reciprocal_latt_params(self) -> tuple[float, float, float, float, float, float]:
        """Return reciprocal lattice parameter lengths and angles in inverse Angstrom and degrees"""
        return reciprocal_latt_params(self.lattice_params)

    @property
    def reciprocal_space_vectors(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return reciprocal space lattice vectors in the Cartesian coordinates"""
        return reciprocal_space_vectors(self.lattice_params)

    @property
    def b_mat(self):
        return b_mat_from_lattice(self.lattice_params)

    def get_q_norm(self, hkl: tuple[float, float, float]) -> float:
        """Convert (h,k,l) to q, in units of inverse Angstrom"""
        try:
            qh, qk, ql = hkl
        except ValueError:
            print("hkl needs to be a tuple of length 3")

        (a_star_vec, b_star_vec, c_star_vec) = reciprocal_space_vectors(self.lattice_params)
        q_vec = qh * a_star_vec + qk * b_star_vec + ql * c_star_vec
        q_norm = float(np.linalg.norm(q_vec))
        return q_norm

    def update_lattice_parametres_from_b_mat(self, b_mat: np.ndarray):
        """Update lattice paramters based on B matrix"""
        lattice_params = lattice_params_from_b_mat(b_mat)
        self.update_lattice_parameters(lattice_params)

    def update_lattice_parametres_from_g_star_mat(self, g_star_mat: np.ndarray):
        lattice_params = lattice_params_from_g_star_mat(g_star_mat)
        self.update_lattice_parameters(lattice_params)
