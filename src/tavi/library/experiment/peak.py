"""Representing a Bragg peak."""

from dataclasses import field
from typing import Optional

from pydantic.dataclasses import dataclass


@dataclass(frozen=True)
class MotorAngles:
    """
    Motor anlges.

    angles_dict: currently only stores {name: value}. But in future it will
                stores {name: (direction, angle)} for specific motor.
                Direction follows mantid convention. i.e. {"sgl", ("+y", 30)}.
                This should be done in the same PR restructure goneiometer.

    Note:
        use angles = (two_theta, omega, sgl, sgu) for a Huber table,
        angles = (two_theta, omega, chi, phi) for a four-circle in the bisect mode.

    """

    angles_dict: Optional[dict] = field(default_factory=dict)


@dataclass(frozen=True)
class DataPoint:
    """
    Description of an experiment data point that contains hkl, ei, ef, motor angles.

    Phyical/virtual monitor positions

    hkl: miller indice (h,k,l)
    ei: incident energy in meV
    ef: scattered energy in meV
    angles:
    """

    hkl: tuple[float, float, float]
    ei: Optional[float] = None
    ef: Optional[float] = None
    angles: Optional[MotorAngles] = None
