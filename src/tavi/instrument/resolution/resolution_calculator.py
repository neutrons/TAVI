from typing import Union

import numpy as np


class ResolutionCalculator:
    """
    Resoluatoin Calculator base class

    Attributes:
        instrument (TAS):

    Methods:
        generate_hkle_list
        validate_instrument_parameters

    """

    def __init__(self, instrument):
        self.instrument = instrument

    @staticmethod
    def require_existance(obj, attr_name):
        val = getattr(obj, attr_name, None)
        if val is None:
            raise AttributeError(attr_name + f" is missing in {repr(obj)}.")
        return val

    @staticmethod
    def require_positive(obj, attr_name):
        val = getattr(obj, attr_name, None)
        if val is None:
            raise AttributeError(attr_name + f" is missing in {repr(obj)}.")
        if isinstance(val, (float, int)):
            if val < 0:
                raise ValueError(attr_name + f" = {val} in {repr(obj)} cannot be negative.")
        elif isinstance(val, (list, tuple, np.ndarray)):
            if any(v < 0 for v in val):
                raise ValueError(attr_name + f" = {val} in {repr(obj)} cannot be negative.")

    @staticmethod
    def check_sense(obj):
        val = getattr(obj, "sense", None)
        if val is None:
            raise AttributeError(f"sense is missing in {repr(obj)}.")
        if val not in ("+", "-"):
            raise ValueError(f"sense of {repr(obj)} needs to be either '+' or '-'")

    @staticmethod
    def generate_hkle_from_projection(u, v, w, en, projection):
        variable_str = ("u", "v", "w", "en")
        variables = (u, v, w, en)

        # Generate value arrays
        variable_list = []
        for name, var in zip(variable_str, variables):
            if isinstance(var, tuple):
                if len(var) != 3:
                    raise ValueError(f"{name}={var} must be a tuple of size 3")
                variable_list.append(np.arange(*var))
            elif isinstance(var, (int, float)):
                variable_list.append(np.array([var]))
            else:
                raise ValueError(f"{name}={var} must be a tuple of size 3 or a numeric scalar")

        # Convert projection vectors to arrays once
        proj = [np.asarray(p) for p in projection]

        # Generate hkle list
        return [
            tuple(proj[0] * u0 + proj[1] * v0 + proj[2] * w0) + (en0,)
            for u0 in variable_list[0]
            for v0 in variable_list[1]
            for w0 in variable_list[2]
            for en0 in variable_list[3]
        ]

    def generate_hkleief_list(
        self,
        hkle: Union[tuple, list[tuple]],
    ) -> tuple[tuple[tuple[float, float, float], float, float], ...]:
        """Generate a list containing tuple ((h, k, l), ei, ef)

        Arguments:
            hkle (tuple | list(tuple)): (h,k,l, en)

        Return:
            hkleief_list (tuple): list of (((h, k, l), ei, ef), ...)
        """

        hkle_list = [hkle] if not isinstance(hkle, list | np.ndarray) else hkle

        hkleief_list = []
        for h, k, l, en in hkle_list:
            ei, ef = self.instrument._get_ei_ef(en=en)
            hkleief_list.append(((h, k, l), ei, ef))

        return tuple(hkleief_list)

    def validate_instrument_parameters(self):
        """Check if enough instrument parameters are provided for Cooper-Nathans method"""

        require_existance = ResolutionCalculator.require_existance
        require_positive = ResolutionCalculator.require_positive
        check_sense = ResolutionCalculator.check_sense

        instrument = self.instrument
        mono = require_existance(instrument, "monochromator")
        require_positive(mono, "mosaic_h")
        require_positive(mono, "mosaic_v")
        check_sense(mono)

        ana = require_existance(instrument, "analyzer")
        require_positive(ana, "mosaic_h")
        require_positive(ana, "mosaic_v")
        check_sense(ana)

        coll = require_existance(instrument, "collimators")
        require_positive(coll, "horizontal_divergence")
        require_positive(coll, "vertical_divergence")

        goni = require_existance(instrument, "goniometer")
        check_sense(goni)

        sample = require_existance(instrument, "sample")
