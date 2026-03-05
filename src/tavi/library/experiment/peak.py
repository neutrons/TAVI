"""Representing a Bragg peak."""

from dataclasses import field
from typing import Optional

from pydantic.dataclasses import dataclass


@dataclass(frozen=True)
class Peak:
    """
    Description of a peak and motor angles.

    Phyical/virtual monitor positions

    hkl: miller indice (h,k,l)
    angles:
    """

    hkl: tuple[float, float, float]
    angles: Optional[MotorAngles] = None


@dataclass(frozen=True)
class MotorAngles:
    """
    Moter anlges.

    Using angles_dict to store motor angles as key,val pairs.

    Recogonized keyword:

    two_theta: s2 angle, in degree
    omega: s1 angle, in degree
    sgl: sample goniometer lower, in degree
    sgu: sample goniometer upper, in degree
    chi: chi angle for a four-circle goniometer, in degree
    phi: phi angle for a four-circle goniometer, in degree

    Note:
        use angles = (two_theta, omega, sgl, sgu) for a Huber table,
        angles = (two_theta, omega, chi, phi) for a four-circle in the bisect mode.

    """

    angles_dict: Optional[dict] = field(default_factory=dict)
