# -*- coding: utf-8 -*-
from typing import Optional

import numpy as np

from tavi.instrument.components.component_base import TASComponent


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

    def __repr__(self):
        return "Collimator"

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
