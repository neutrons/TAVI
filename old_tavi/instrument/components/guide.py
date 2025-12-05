# -*- coding: utf-8 -*-
from typing import Optional

from tavi.instrument.components.component_base import TASComponent


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
