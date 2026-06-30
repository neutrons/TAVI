"""Enumeration of fixed-energy modes for triple-axis experiments."""

from enum import Enum


class FixedEnergyMode(Enum):
    """Which energy is held fixed in a triple-axis experiment."""

    FIX_Ei = "fix_ei"
    FIX_Ef = "fix_ef"
