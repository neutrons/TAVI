# -*- coding: utf-8 -*-
from typing import Literal, Optional, Union

import numpy as np

from tavi.utilities import cm2angstrom


class TASComponent(object):
    """
    A component of the triple-axis spectrometer

    Attributes:
        type (str): identify which component
    Methods:
        update(param, value): update a parameter
    """

    def __init__(
        self,
        param_dict: Optional[dict] = None,
        component_name: str = "",
    ):
        """Load parameters if provided."""
        self.component_name = component_name
        if param_dict is not None:
            for key, val in param_dict.items():
                setattr(self, key, val)

    def update(
        self,
        param: str,
        value: Union[str, float, int],
    ) -> None:
        """update a paramter if exist"""
        if hasattr(self, param):
            setattr(self, param, value)
            print(f"Setting {self.component_name} parameter {param} to {value}")
        else:
            print(f"{self.component_name} has no attribute {param}")

    @staticmethod
    def _min2rad(
        angle_in_minutes: Optional[float],
        param_name: str = "Angle",
    ) -> Optional[float]:
        """Convert from minutes to radian is angle is not None"""
        if angle_in_minutes is None:
            print(f"{param_name} is None.")
            return None
        return np.deg2rad(angle_in_minutes / 60)

    @staticmethod
    def _cm2angstrom_given_shape(
        length_in_cm: Optional[float],
        shape: str,
        param_name: str = "Length",
    ) -> Optional[float]:
        """
        Convert length from centimeter to angstrom.

        Apply correction based on shape, for resolution calculation.
        Divide length by np.sqrt(12) if rectangular,
        Divive diameter D by 4 if circular or spherical.
        """
        if length_in_cm is None:
            print(f"{param_name} is None.")
            return None
        match shape:
            case "rectangular":
                return length_in_cm / np.sqrt(12) * cm2angstrom
            case "circular" | "spherical":
                return length_in_cm / 4 * cm2angstrom
            case _:
                print("Unrecognized shape. Needs to be rectangular or circular/spherical.")
                return None

    @staticmethod
    def _cm2angstrom(
        length_in_cm: Optional[float],
        param_name: str = "Length",
    ) -> Optional[float]:
        """
        Convert length from centimeter to angstrom.
        """
        if length_in_cm is None:
            print(f"{param_name} is None.")
            return None

        return length_in_cm * cm2angstrom


class Source(TASComponent):
    """Neutron source"""

    def __init__(
        self,
        param_dict: Optional[dict] = None,
        component_name="source",
    ):
        self.name: str = "HFIR"
        self.shape: Literal["rectangular", "circular"] = "rectangular"
        self.width: Optional[float] = None
        self.height: Optional[float] = None

        super().__init__(param_dict, component_name)

    @property
    def _width(self):
        """width in angstrom, with correction based on shape, for resolution calculation"""
        return TASComponent._cm2angstrom_given_shape(self.width, self.shape, "Source width")

    @property
    def _height(self):
        """Height in angstrom, with correction based on shape, for resolution calculation"""
        return TASComponent._cm2angstrom_given_shape(self.height, self.shape, "Source height")


class Guide(TASComponent):
    """Neutron guide in the monochromator drum"""

    def __init__(
        self,
        param_dict: Optional[dict] = None,
        component_name: str = "guide",
    ):
        self.in_use: bool = False
        self.div_h: Optional[float] = None
        self.div_v: Optional[float] = None

        super().__init__(param_dict, component_name)

    @property
    def _div_h(self):
        """Horizontal divergence in radian"""
        return TASComponent._min2rad(self.div_h, "Guide horizontal divergence")

    @property
    def _div_v(self):
        """Vertical divergence in radian"""
        return TASComponent._min2rad(self.div_v, "Guide vertical divergence")


class Collimators(TASComponent):
    """collimitor divergences, in mins of arc"""

    def __init__(
        self,
        param_dict: Optional[dict] = None,
        component_name: str = "collimitors",
    ):
        # defalut values
        self.h_pre_mono: float = 30.0  # mins of arc
        self.h_pre_sample: float = 30.0
        self.h_post_sample: float = 30.0
        self.h_post_ana: float = 30.0
        self.v_pre_mono: float = 30.0
        self.v_pre_sample: float = 30.0
        self.v_post_sample: float = 30.0
        self.v_post_ana: float = 30.0

        super().__init__(param_dict, component_name)

    @property
    def horizontal_divergence(self) -> list:
        """list of horizontal divergence in minitus of arc"""
        return [
            self.h_pre_mono,
            self.h_pre_sample,
            self.h_post_sample,
            self.h_post_ana,
        ]

    @property
    def _horizontal_divergence(self) -> list:
        """list of horizontal divergence in radian"""
        return [
            np.deg2rad(self.h_pre_mono / 60),
            np.deg2rad(self.h_pre_sample / 60),
            np.deg2rad(self.h_post_sample / 60),
            np.deg2rad(self.h_post_ana / 60),
        ]

    @property
    def vertical_divergence(self) -> list:
        """list of vertical divergence in minutes of arcs"""
        return [
            self.v_pre_mono,
            self.v_pre_sample,
            self.v_post_sample,
            self.v_post_ana,
        ]

    @property
    def _vertical_divergence(self) -> list:
        """list of vertical divergence in radian"""
        return [
            np.deg2rad(self.v_pre_mono / 60),
            np.deg2rad(self.v_pre_sample / 60),
            np.deg2rad(self.v_post_sample / 60),
            np.deg2rad(self.v_post_ana / 60),
        ]


# TODO not implemented in resolution calculation
class Monitor(TASComponent):
    """
    Neutron Monitor

    Note:
        Needs to think about monitor efficiency
    """

    def __init__(
        self,
        param_dict: Optional[dict] = None,
        component_name: str = "monitor",
    ):
        self.shape: Literal["rectangular", "spherical"] = "rectangular"
        self.width: Optional[float] = None
        self.height: Optional[float] = None

        super().__init__(param_dict, component_name)

    @property
    def _width(self):
        """width in angstrom, with correction based on shape, for resolution calculation"""
        return TASComponent._cm2angstrom_given_shape(self.width, self.shape, "Monitor width")

    @property
    def _height(self):
        """Height in angstrom, with correction based on shape, for resolution calculation"""
        return TASComponent._cm2angstrom_given_shape(self.height, self.shape, "Monitor height")


# TODO not implemented in resolution calculation
class Detector(TASComponent):
    """
    Neutron detector

    Note:
        Need to think about detector efficiency, saturation
    """

    def __init__(
        self,
        param_dict: Optional[dict] = None,
        component_name: str = "detector",
    ):
        self.shape: Literal["rectangular", "spherical"] = "rectangular"
        self.width: Optional[float] = None
        self.height: Optional[float] = None

        super().__init__(param_dict, component_name)

    @property
    def _width(self):
        """width in angstrom, with correction based on shape, for resolution calculation"""
        return TASComponent._cm2angstrom_given_shape(self.width, self.shape, "Detector width")

    @property
    def _height(self):
        """Height in angstrom, with correction based on shape, for resolution calculation"""
        return TASComponent._cm2angstrom_given_shape(self.height, self.shape, "Detector height")


class Distances(TASComponent):
    """Distance between components, in units of centimeters"""

    def __init__(
        self,
        param_dict: Optional[dict] = None,
        component_name: str = "distances",
    ):
        self.src_mono: Optional[float] = None
        self.mono_sample: Optional[float] = None
        self.sample_ana: Optional[float] = None
        self.ana_det: Optional[float] = None
        self.mono_monitor: Optional[float] = None

        super().__init__(param_dict, component_name)

    @property
    def _src_mono(self):
        return TASComponent._cm2angstrom(self.src_mono, "Source to monochromator distance")

    @property
    def _mono_sample(self):
        return TASComponent._cm2angstrom(self.src_mono, "Monochromator to sample distance")

    @property
    def _sample_ana(self):
        return TASComponent._cm2angstrom(self.src_mono, "Sample to analyzer distance")

    @property
    def _ana_det(self):
        return TASComponent._cm2angstrom(self.src_mono, "Analyzer to detector distance")

    @property
    def _mono_monitor(self):
        return TASComponent._cm2angstrom(self.src_mono, "Monochromator to monitor distance")
