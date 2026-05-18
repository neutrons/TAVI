"""
Goniometer component.

TODO: add 4 circle mode, express rotation with quaternion.
"""

from typing import Literal, Optional

import numpy as np

from tavi.library.experiment.peak import MotorAngles
from tavi.library.experiment.utilities import rot_x, rot_y, rot_z


class Goniometer:
    """
    Goniometer.

    For Huber table, use type Y,-Z,X or Y,Z,-X
    For Four-Cricle in bisect mode, use type ?Y,Z,Y,bisect

    Attributes:
        type (str):
        sense (str): "+" if s2 is right-hand
        limits (dict): (min, max) of s2 angle

    Methods:
        _sense
        stacking_order
        mode (str): Mode for Four-Circle. "bisect" or "azimuthal"
        motor_senses: signs in the goniometer type string
        r_mat
        r_mat_inv
        set_limit
        validate_motor_positions

    """

    def __init__(self, type: str = "Y, -Z, X", s2_sense: Optional[Literal["-", "+"]] = None) -> None:
        """Init."""
        self.type = type
        self._sense = s2_sense

    @property
    def sense(self) -> int:
        """Get s2 sense of goniometer."""
        return 1 if self._sense == "+" else -1

    @sense.setter
    def sense(self, val: str) -> None:
        """Set sense."""
        self._sense = val

    def _get_motor_senses(self) -> tuple[int, int, int]:
        ax0, ax1, ax2, *__ = self.type.split(",")
        signs = []
        for ax in (ax0, ax1, ax2):
            signs.append(-1 if ax.strip().startswith("-") else 1)
        return tuple(signs)

    def r_mat(self, angles: MotorAngles) -> np.ndarray:
        """
        Goniometer rotation matrix R.

        Args:
            angles (MotorAngle): two_theta, omega, sgl, sgu, chi, phi.

        Note:
            type refers to the axis of rotation for each motor,
            + for CCW, - for CW.
            The Cartesian coordinate Z is along the incoming beam,
            X in plane, Y up, right handed. This is the convention used
            in Mantid/International Crystallography Table.

        """
        signs = self._get_motor_senses()
        # initial implementation of huber table YZX mode
        r_mat = (
            rot_y(angles.angles_dict["omega"] * signs[0])
            @ rot_z(angles.angles_dict["sgl"] * signs[1])
            @ rot_x(angles.angles_dict["sgu"] * signs[2])
        )
        return r_mat

    def angles_from_r_mat(self, r_mat: np.ndarray) -> tuple[float, float, float]:
        """
        Calculate goniometer angles from the R matrix. Eq. 105.

        Note:
            In bisect mode, omega angle is half of two_theta for diffrction. But for inelastic scattering,
            chi ring should be rotated where the axis of chi is perpendicular to Q.

            range of np.arcsin is -pi/2 to pi/2.
            range of np.arctan2 is -pi to pi.

        """
        signs = self._get_motor_senses()
        # for YZX mode, Eq.100-105. Notice they share a common denominator.

        denomitor = np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2)
        # omega rotate along y, omega ~ omega
        omega = np.arctan2(-r_mat[2, 0], r_mat[0, 0]) * signs[0]
        # sgl rotate along z, sgl ~ mu
        sgl = np.arctan2(r_mat[1, 0], denomitor) * signs[1]
        # sgu rotate along x, sgu ~ epsilon
        sgu = np.arctan2(-r_mat[1, 2], r_mat[1, 1]) * signs[2]

        angles_dict = dict(omega=np.degrees(omega), sgl=np.degrees(sgl), sgu=np.degrees(sgu))
        return MotorAngles(angles_dict)
