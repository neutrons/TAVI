# -*- coding: utf-8 -*-
from typing import Literal, Optional

from tavi.instrument.components.component_base import TASComponent


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
