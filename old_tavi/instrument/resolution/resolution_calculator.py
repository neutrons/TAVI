from typing import Union

import numpy as np


class ResolutionCalculator:
    """
    Resoluatoin Calculator base class

    Attributes:
        instrument (TAS):

    Methods:
        generate_hkleief_list
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
    def generate_hkle(grid: tuple, axes: tuple = ((1, 0, 0), (0, 1, 0), (0, 0, 1), "en")):
        """
        Generate a grid of (h, k, l, en) points based on given axes and grid parameters.

        Args:
            grid (tuple): A tuple of length 4. Each element can be:
                        - a scalar -> number of points along the corresponding axis
                        - a tuple (min, max, step) defining range and resolution
            axes (tuple): Axes definition. Should be a tuple of 4 elements:
                        - 3 numeric vectors for reciprocal space directions
                        - "en" for energy dimension

        Returns:
            np.ndarray: Array of shape (N, 4) for (h, k, l, en) grid points
        """
        # validate length of grid and axes
        if len(grid) != 4:
            raise ValueError("Grid must have 4 elements.")

        if len(axes) != 4 or "en" not in axes:
            raise ValueError('Axes must contain 3 Q vectors and "en".')

        q_vectors = [np.array(ax) for ax in axes if ax != "en"]

        # Generate value arrays
        q_list = []
        for name, var in zip(axes, grid):
            if isinstance(var, tuple):
                if len(var) != 3:
                    raise ValueError(f"{name} = {var} must be (min, max, step)")
                vals = np.arange(*var)
            elif isinstance(var, (int, float)):
                vals = np.array([var])
            else:
                raise ValueError(f"{name} = {var} must be (min, max, step) or scalar")

            if name == "en":
                e_list = vals
            else:
                q_list.append(vals)

        # Generate hkle list
        return [
            tuple(q_vectors[0] * u + q_vectors[1] * v + q_vectors[2] * w) + (en,)
            for u in q_list[0]
            for v in q_list[1]
            for w in q_list[2]
            for en in e_list
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
