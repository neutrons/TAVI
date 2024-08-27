from typing import Literal, Optional

import numpy as np

from tavi.instrument.tas_cmponents import TASComponent

# ---------------------------------------------------------------
# d_spacing table from Shirane Appendix 3, in units of Angstrom
# ---------------------------------------------------------------
mono_ana_xtal = {
    "PG002": 3.35416,
    "Pg002": 3.35416,
    "PG004": 1.67708,
    "Cu111": 2.08717,
    "Cu220": 1.27813,
    "Ge111": 3.26627,
    "Ge220": 2.00018,
    "Ge311": 1.70576,
    "Ge331": 1.29789,
    "Be002": 1.79160,
    "Be110": 1.14280,
    "Heusler": 3.435,  # Cu2MnAl(111)
}


class MonoAna(TASComponent):
    """
    Monochromator and analyzer class

    Note:
        protected attributes are reserved for resolution calculation
    """

    def __init__(
        self,
        param_dict: Optional[dict] = None,
        component_name: str = "",
    ):
        # defalut values
        self.type: str = "PG002"
        self.d_spacing: float = mono_ana_xtal["PG002"]
        self.mosaic_h: float = 45.0  # horizontal mosaic, min of arc
        self.mosaic_v: float = 45.0  # vertical mosaic, if anisotropic
        self.sense: Literal[-1, 1] = -1  # +1 for counter-clockwise, -1 for clockwise

        # divide by np.sqrt(12) if rectangular
        # Diameter D/4 if spherical
        self.shape: Literal["rectangular", "spherical"] = "rectangular"
        self.width: float = 12.0
        self.height: float = 8.0
        self.depth: float = 0.15
        # horizontal focusing
        self.curved_h: bool = False
        self.curvh: float = 0.0
        self.optimally_curved_h: bool = False
        # vertical focusing
        self.curved_v: bool = False
        self.curvv: float = 0.0
        self.optimally_curved_v: bool = False

        super().__init__(param_dict, component_name)
        self.d_spacing = mono_ana_xtal[self.type]

    @property
    def _mosaic_h(self):
        return np.deg2rad(self.mosaic_h / 60)

    @property
    def _mosaic_v(self):
        return np.deg2rad(self.mosaic_v / 60)

    @property
    def _width(self):
        """width in angstrom, with correction based on shape, for resolution calculation"""
        return TASComponent._cm2angstrom_given_shape(self.width, self.shape, f"{self.component_name} width")

    @property
    def _height(self):
        """height in angstrom, with correction based on shape, for resolution calculation"""
        return TASComponent._cm2angstrom_given_shape(self.height, self.shape, f"{self.component_name} height")

    @property
    def _depth(self):
        """depth in angstrom, with correction based on shape, for resolution calculation"""
        return TASComponent._cm2angstrom_given_shape(self.depth, self.shape, f"{self.component_name} depth")
