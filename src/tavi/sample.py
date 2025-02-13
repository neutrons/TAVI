# -*- coding: utf-8 -*-
import json
from typing import Optional

import numpy as np

from tavi.instrument.ub_algorithm import (
    b_mat_from_lattice,
    real_space_vectors,
    reciprocal_latt_params,
    reciprocal_space_vectors,
)


class Sample(object):
    """
    Sample class

    Attributes:
        json_dict (dict): dictionary from json file
        shape (str): "cuboid" or "cylindrical"
        width (float): in units of cm
        height (float): in units of cm
        depth (float): in units of cm

        mosaic_h (fload): in units of minutes of arc
        mosaic_v (fload): verital mosaic if anisotropic, in units of minutes of arc

        a, b, c                 lattice constants in Angstrom
        alpha, beta, gamma      angles in degrees

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
        self.type = "generic"
        self.json_dict: Optional[dict] = None
        # lattice parameters
        self.a: float
        self.b: float
        self.c: float
        self.alpha: float
        self.beta: float
        self.gamma: float
        self.b_mat: np.ndarray
        # parameters for resolution calculation
        self.shape: str = "cuboid"
        self.width: float = 1.0  # in cm
        self.height: float = 1.0  # in cm
        self.depth: float = 1.0  # in cm
        self.mosaic_h: float = 30  # horizontal mosaic in minutes of arc
        self.mosaic_v: float = 30  # vertical mosaic in minutes of arc

        try:
            a, b, c, alpha, beta, gamma = lattice_params
        except ValueError:
            print("Incomplete lattice parameters.")

        for length in (a, b, c):
            if length <= 0:
                raise ValueError("Lattice parameters smaller than zero.")

        for angle in (alpha, beta, gamma):
            if 0.0 > angle or angle > 180.0:
                raise ValueError("Lattice angles out of range.")

        self.update_lattice_parameters(lattice_params)

    def _unpack_json_parameters(self):
        sample_params = self.json_dict

        shape = sample_params.get("shape")
        width = sample_params.get("width")
        height = sample_params.get("height")
        depth = sample_params.get("depth")
        if None not in [shape, width, height, depth]:
            self.set_shape(shape, width, height, depth)

        mosaic_h = sample_params.get("mosaic_h")
        mosaic_v = sample_params.get("mosaic_v")
        if None not in [mosaic_h, mosaic_v]:
            self.set_mosaic(mosaic_h, mosaic_v)

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
        sample.json_dict = sample_params
        sample._unpack_json_parameters()

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

    @property
    def _mosaic_h(self) -> float:
        """horizontal mosaic in radian, for resolution calculation"""
        return np.deg2rad(self.mosaic_h / 60)

    @property
    def _mosaic_v(self) -> float:
        """vertical mosaic in radian, for resolution calculation"""
        return np.deg2rad(self.mosaic_v / 60)

    def update_lattice_parameters(
        self,
        lattice_params: tuple[float, float, float, float, float, float] = (1, 1, 1, 90, 90, 90),
    ):
        """update real and reciprocal space lattice parameters and vectors"""

        a, b, c, alpha, beta, gamma = lattice_params
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

    def q_norm(self, hkl: tuple[float, float, float]) -> float:
        """Convert (h,k,l) to q, in units of inverse Angstrom"""
        try:
            qh, qk, ql = hkl
        except ValueError:
            print("hkl needs to be a tuple of length 3")

        (a_star_vec, b_star_vec, c_star_vec) = reciprocal_space_vectors(self.lattice_params)
        q_vec = qh * a_star_vec + qk * b_star_vec + ql * c_star_vec
        q_norm = float(np.linalg.norm(q_vec))
        return q_norm
