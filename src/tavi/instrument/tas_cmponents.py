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
    ) -> float:
        """Convert from minutes to radian is angle is not None"""
        if angle_in_minutes is None:
            print(f"{param_name} is None.")
            return None
        return np.deg2rad(angle_in_minutes / 60)

    @staticmethod
    def _cm2angstrom(
        length_in_cm: Optional[float],
        shape: str,
        param_name: str = "Length",
    ) -> float:
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


class Source(TASComponent):
    """Neutron source"""

    def __init__(
        self,
        param_dict: Optional[dict] = None,
        component_name="source",
    ):
        self.name: str = "HFIR"
        self.shape: Literal["rectangular", "circular"] = "rectangular"
        self.width: Optional[str] = None
        self.height: Optional[str] = None

        super().__init__(param_dict, component_name)
        pass

    @property
    def _width(self):
        """width in angstrom, with correction based on shape, for resolution calculation"""
        return TASComponent._cm2angstrom(self.width, self.shape, "Source width")

    @property
    def _height(self):
        """Height in angstrom, with correction based on shape, for resolution calculation"""
        return TASComponent._cm2angstrom(self.height, self.shape, "Source height")


class Guide(TASComponent):
    """Neutron guide in the monochromator drum"""

    def __init__(
        self,
        param_dict: Optional[dict] = None,
        component_name: str = "guide",
    ):
        self.in_use: bool = False
        self.div_h: Optional[str] = None
        self.div_v: Optional[str] = None

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
# TODO =========================================================
class Monitor(object):
    def __init__(self, param_dict, component_name):
        self.shape = "rectangular"  # or spherical
        self.width = 5 * cm2angstrom
        self.height = 12 * cm2angstrom

        for key, val in param_dict.items():
            match key:
                case "width" | "height":
                    # divide by np.sqrt(12) if rectangular
                    # Diameter D/4 if spherical
                    if param_dict["shape"] == "rectangular":
                        setattr(self, key, val / np.sqrt(12) * cm2angstrom)
                    elif param_dict["shape"] == "spherical":
                        setattr(self, key, val / 4 * cm2angstrom)
                    else:
                        print("Monitor shape needs to be either rectangular or spherical.")
                case _:
                    setattr(self, key, val)


# TODO
class Detector(object):
    def __init__(self, param_dict, component_name):
        self.shape = "rectangular"
        self.width = 1.5  # in cm
        self.height = 5.0  # in cm

        for key, val in param_dict.items():
            match key:
                case "width" | "height":
                    # divide by np.sqrt(12) if rectangular
                    # Diameter D/4 if spherical
                    if param_dict["shape"] == "rectangular":
                        setattr(self, key, val / np.sqrt(12) * cm2angstrom)
                    elif param_dict["shape"] == "spherical":
                        setattr(self, key, val / 4 * cm2angstrom)
                    else:
                        print("Detector shape needs to be either rectangular or spherical.")
                case _:
                    setattr(self, key, val)


class Arms(object):
    """lengths of arms, in units of cm"""

    def __init__(self, param_dict, component_name):
        # in units of cm
        self.src_mono = 0
        self.mono_sample = 0
        self.sample_ana = 0
        self.ana_det = 0
        self.mono_monitor = 0

        for key, val in param_dict.items():
            setattr(self, key, val * cm2angstrom)
