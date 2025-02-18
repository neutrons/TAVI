from typing import List, Literal, Optional

import numpy as np

from tavi.instrument.components.component_base import TASComponent
from tavi.utilities import MotorAngles, rot_x, rot_y, rot_z


class Goniometer(TASComponent):
    """
    Goniometer
    For Huber table, use type Y,-Z,X or Y,Z,-X
    For Four-Cricle in bisect mode, use type ?Y,Z,Y,bisect
    """

    def __init__(
        self,
        param_dict: Optional[dict] = None,
        component_name: str = "goniometer",
    ):
        self.type: str = "Y,-Z,X"  # Y-mZ-X for (s1, sgl, sgu), Huber table at HB1A and HB3,
        self.sense: Literal["-", "+"] = "-"  # determines the sign of s2
        self.s2_limit: Optional[tuple] = None
        self.omega_limit: Optional[tuple] = None
        self.sgl_limit: Optional[tuple] = (-10, 10)
        self.sgu_limit: Optional[tuple] = (-10, 10)
        self.chi_limit: Optional[tuple] = None
        self.phi_limit: Optional[tuple] = None

        super().__init__(param_dict)
        self.component_name = component_name

    @property
    def stacking_order(self) -> str:
        """Return motors staking order without the signs"""

        ax0, ax1, ax2, _ = self.type.split(",")
        stacking_order = ""
        for ax in [ax0, ax1, ax2]:
            if ax.startswith("-"):
                stacking_order += ax.strip("-")
            else:
                stacking_order += ax
        return stacking_order

    @property
    def mode(self) -> Optional[str]:
        """Return mode if running a Four-Circle"""
        _, _, _, *mode_list = self.type.split(",")
        mode = None if not mode_list else mode_list[0]
        return mode

    @property
    def motor_senses(self) -> List[Literal[+1, -1]]:
        ax0, ax1, ax2, _ = self.type.split(",")
        signs = []
        for ax in [ax0, ax1, ax2]:
            sign = -1 if ax.startswith("-") else +1
            signs.append(sign)
        return signs

    @property
    def _sense(self):
        """Sense of goniometer is the sign of s2"""
        match self.sense:
            case "+":
                return +1
            case "-":
                return -1
            case _:
                raise ValueError(f"sense {self.sense} needs to be either '+' or '-'.")

    def r_mat(self, angles: MotorAngles) -> np.ndarray:
        """Goniometer rotation matrix R

        Args:
            angles (MotorAngle): two_theta, omega, sgl, sgu, chi, phi

        Note:
            type refers to the axis of rotation for each motor,
            + for CCW, - for CW.
            The Cartesian coordinate Z is along the incoming beam,
            X in plane, Y up, right handed. This is the convention used
            in Mantid/International Crystallography Table.
        """
        signs = self.motor_senses
        operating_mode = self.stacking_order + self.mode if self.mode is not None else self.stacking_order
        match operating_mode:
            case "YZX":  # Huber table
                r_mat = np.matmul(
                    rot_y(angles.omega * signs[0]),
                    np.matmul(
                        rot_z(angles.sgl * signs[1]),
                        rot_x(angles.sgu * signs[2]),
                    ),
                )
            case "YZY_bisect":  # Four-Circle in s1 bisect mode
                r_mat = np.matmul(
                    rot_y(angles.omega * signs[0]),
                    np.matmul(
                        rot_z(angles.chi * signs[1]),
                        rot_y(angles.phi * signs[2]),
                    ),
                )
            case _:
                raise ValueError(f"Unknow goniometert={operating_mode}")

        return r_mat

    def r_mat_inv(self, angles: MotorAngles) -> np.ndarray:
        """inverse of rotation matrix, equivalent to transpose since R is unitary"""
        # return np.linalg.inv(self.r_mat(angles))
        return self.r_mat(angles).T

    def angles_from_r_mat(self, r_mat: np.ndarray, two_theta: Optional[float] = None) -> MotorAngles:
        """Calculate goniometer angles from the R matrix

        Note:
            range of np.arcsin is -pi/2 to pi/2
            range of np.arctan2 is -pi to pi
        """

        signs = self.motor_senses
        operating_mode = self.stacking_order + self.mode if self.mode is not None else self.stacking_order

        match operating_mode:
            case "YZX":  # Y-mZ-X (s1, sgl, sgu) for HB1A and HB3,
                denominator = np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2)
                sgl_rad = np.arctan2(r_mat[1, 0], denominator)
                sgu_rad = np.arctan2(-r_mat[1, 2] / denominator, r_mat[1, 1] / denominator)
                omega_rad = np.arctan2(-r_mat[2, 0] / denominator, r_mat[0, 0] / denominator)
                omega, sgl, sgu = np.rad2deg(
                    [
                        signs[0] * omega_rad,
                        signs[1] * sgl_rad,
                        signs[2] * sgu_rad,
                    ]
                )
                angles = MotorAngles(two_theta, omega, sgl, sgu, None, None)
                self.validate_motor_positions(angles)

            case "YZY_bisect":
                # omega, chi, phi = stacking_order_yzy(r_mat)
                angles = MotorAngles(two_theta, signs[0] * omega, None, None, signs[1] * chi, signs[2] * phi)
                self.validate_motor_positions(angles)

            case _:
                raise ValueError("Unknow goniometer type.  Curruntly support Y-ZX and YZ-X.")

        return angles

    def set_limit(
        self,
        motor_name: Literal["s2", "omega", "sgl", "sgu", "chi", "phi"],
        motor_range: tuple[float, float] = (-180, 180),
    ):
        "set goiometer motor limt"
        setattr(self, motor_name + "_limit", motor_range)

    # TODO
    def validate_motor_positions(self, angles: MotorAngles):
        "check if all goiometer motors are within the limits"
        pass
